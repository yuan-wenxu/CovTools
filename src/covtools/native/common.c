#include "common.h"
#include <ctype.h>
#include <fcntl.h>
#include <locale.h>
#include <pthread.h>
#include <stdarg.h>

static locale_t numeric_locale;
static pthread_once_t locale_once = PTHREAD_ONCE_INIT;

static void make_locale(void) {
    numeric_locale = newlocale(LC_NUMERIC_MASK, "C", (locale_t)0);
}

int value_error(Context *ctx, const char *format, ...) {
    if (!ctx->kind) {
        va_list args;
        ctx->kind = 1;
        va_start(args, format);
        vsnprintf(ctx->message, sizeof(ctx->message), format, args);
        va_end(args);
    }
    return 0;
}

int io_error(Context *ctx, const char *path) {
    if (!ctx->kind) {
        ctx->kind = 2;
        ctx->code = errno ? errno : EIO;
        snprintf(ctx->message, sizeof(ctx->message), "%s", path);
    }
    return 0;
}

int memory_error(Context *ctx) {
    if (!ctx->kind) ctx->kind = 3;
    return 0;
}

PyObject *raise_error(Context *ctx) {
    if (ctx->kind == 1) PyErr_SetString(PyExc_ValueError, ctx->message);
    else if (ctx->kind == 2) {
        errno = ctx->code;
        PyErr_SetFromErrnoWithFilename(PyExc_OSError, ctx->message);
    } else if (ctx->kind == 3) PyErr_NoMemory();
    else if (!PyErr_Occurred())
        PyErr_SetString(PyExc_RuntimeError, "Native processing failed");
    return NULL;
}

int tick(Context *ctx, uint64_t bytes, int force) {
    ctx->pending += bytes;
    if (!force && ctx->pending < 1024 * 1024) return 1;
    PyGILState_STATE state = PyGILState_Ensure();
    int ok = PyErr_CheckSignals() == 0;
    if (ok && ctx->callback && ctx->callback != Py_None && ctx->pending) {
        PyObject *result = PyObject_CallFunction(
            ctx->callback, "K", (unsigned long long)ctx->pending);
        ok = result != NULL;
        Py_XDECREF(result);
    }
    if (!ok) ctx->kind = 4;
    PyGILState_Release(state);
    ctx->pending = 0;
    return ok;
}

int fields(Context *ctx, Record *record, const char *line, size_t length) {
    while (length && (line[length - 1] == '\n' || line[length - 1] == '\r'))
        length--;
    record->length = length;
    record->count = 0;
    size_t start = 0;
    for (size_t i = 0; i <= length; i++) {
        if (i != length && line[i] != '\t') continue;
        if (record->count == record->capacity) {
            size_t capacity = record->capacity ? record->capacity * 2 : 16;
            if (capacity > SIZE_MAX / sizeof(Span)) return memory_error(ctx);
            Span *resized = realloc(record->fields, capacity * sizeof(Span));
            if (!resized) return memory_error(ctx);
            record->fields = resized;
            record->capacity = capacity;
        }
        record->fields[record->count++] = (Span){line + start, i - start};
        start = i + 1;
    }
    return 1;
}

static Span trim(Span span) {
    while (span.len && isspace((unsigned char)*span.data)) {
        span.data++; span.len--;
    }
    while (span.len && isspace((unsigned char)span.data[span.len - 1]))
        span.len--;
    return span;
}

int unsigned_value(Span span, uint64_t *result) {
    span = trim(span);
    int negative = 0;
    if (span.len && (*span.data == '+' || *span.data == '-')) {
        negative = *span.data == '-';
        span.data++; span.len--;
    }
    if (!span.len) return 0;
    uint64_t value = 0;
    int previous_digit = 0;
    for (size_t i = 0; i < span.len; i++) {
        unsigned char c = span.data[i];
        if (c == '_' && previous_digit && i + 1 < span.len &&
            span.data[i + 1] >= '0' && span.data[i + 1] <= '9') {
            previous_digit = 0; continue;
        }
        if (c < '0' || c > '9' || value > (UINT64_MAX - (c - '0')) / 10)
            return 0;
        value = value * 10 + c - '0';
        previous_digit = 1;
    }
    if (!previous_digit || (negative && value)) return 0;
    *result = value;
    return 1;
}

static int nonnegative_integer(Span span) {
    span = trim(span);
    int negative = 0, nonzero = 0, previous_digit = 0;
    if (span.len && (*span.data == '+' || *span.data == '-')) {
        negative = *span.data == '-'; span.data++; span.len--;
    }
    if (!span.len) return 0;
    for (size_t i = 0; i < span.len; i++) {
        unsigned char c = span.data[i];
        if (c == '_' && previous_digit && i + 1 < span.len &&
            span.data[i + 1] >= '0' && span.data[i + 1] <= '9') {
            previous_digit = 0; continue;
        }
        if (c < '0' || c > '9') return 0;
        if (c != '0') nonzero = 1;
        previous_digit = 1;
    }
    return previous_digit && !(negative && nonzero);
}

int validate(Context *ctx, Record *record, uint64_t offset) {
    if (record->count < 6 || record->count > 9)
        return value_error(ctx, "Invalid record at byte %" PRIu64
                           ": expected six to nine columns", offset);
    for (size_t i = 0; i < record->count; i++)
        if (!record->fields[i].len)
            return value_error(ctx, "Invalid record at byte %" PRIu64
                               ": empty field; use '.' for missing metadata",
                               offset);
    uint64_t end;
    if (!unsigned_value(record->fields[1], &record->position) ||
        !unsigned_value(record->fields[2], &end) ||
        !record->position || end != record->position)
        return value_error(ctx, "Invalid record at byte %" PRIu64
                           ": expected 1-based point coordinates: start = end"
                           " (unsigned 64-bit integers)", offset);
    Span percentage = trim(record->fields[3]);
    char local[256];
    char *number = percentage.len < sizeof(local) ? local :
                   malloc(percentage.len + 1);
    if (!number) return memory_error(ctx);
    size_t digits = 0;
    int valid_number = 1;
    for (size_t i = 0; i < percentage.len; i++) {
        char c = percentage.data[i];
        if (c == '_' && i && i + 1 < percentage.len &&
            percentage.data[i - 1] >= '0' && percentage.data[i - 1] <= '9' &&
            percentage.data[i + 1] >= '0' && percentage.data[i + 1] <= '9')
            continue;
        if (!c || !strchr("0123456789+-.eE", c)) valid_number = 0;
        number[digits++] = c;
    }
    number[digits] = 0;
    char *tail;
    pthread_once(&locale_once, make_locale);
    if (!numeric_locale) {
        if (number != local) free(number);
        return memory_error(ctx);
    }
    double value = strtod_l(number, &tail, numeric_locale);
    int valid = valid_number && tail != number && !*tail && isfinite(value) &&
                value >= 0 && value <= 100;
    if (number != local) free(number);
    if (!valid)
        return value_error(ctx, "Invalid record at byte %" PRIu64
                           ": methylation percentage must be in [0, 100]",
                           offset);
    if (!nonnegative_integer(record->fields[4]) ||
        !nonnegative_integer(record->fields[5]))
        return value_error(ctx, "Invalid record at byte %" PRIu64
                           ": counts must be nonnegative integers", offset);
    if (record->count >= 8 && !(record->fields[7].len == 1 &&
        (*record->fields[7].data == '+' || *record->fields[7].data == '-' ||
         *record->fields[7].data == '.')))
        return value_error(ctx, "Invalid record at byte %" PRIu64
                           ": strand in column 8 must be '+', '-' or '.'",
                           offset);
    return 1;
}

int bytes_compare(Span a, Span b) {
    size_t length = a.len < b.len ? a.len : b.len;
    int comparison = memcmp(a.data, b.data, length);
    return comparison ? comparison : (a.len > b.len) - (a.len < b.len);
}

int span_equal(Span span, const char *data, size_t length) {
    return span.len == length && (!length || !memcmp(span.data, data, length));
}

int utf8_valid(Span span) {
    size_t i = 0;
    while (i < span.len) {
        unsigned char c = span.data[i++];
        if (c < 128) continue;
        int n = c >= 0xc2 && c <= 0xdf ? 1 :
                c >= 0xe0 && c <= 0xef ? 2 :
                c >= 0xf0 && c <= 0xf4 ? 3 : -1;
        if (n < 0 || i + (size_t)n > span.len) return 0;
        unsigned char next = span.data[i];
        if ((c == 0xe0 && next < 0xa0) || (c == 0xed && next >= 0xa0) ||
            (c == 0xf0 && next < 0x90) || (c == 0xf4 && next >= 0x90))
            return 0;
        while (n--) if (((unsigned char)span.data[i++] & 0xc0) != 0x80)
            return 0;
    }
    return 1;
}

int64_t stat_mtime(const struct stat *status) {
    return (int64_t)status->st_mtim.tv_sec * 1000000000LL +
           status->st_mtim.tv_nsec;
}

int reader_open(Context *ctx, Reader *reader, const char *path,
                uint64_t start, uint64_t end, uint64_t size, int64_t mtime) {
    reader->path = path;
    reader->size = size; reader->mtime = mtime;
    reader->offset = start; reader->end = end;
    reader->file = fopen(path, "rb");
    if (!reader->file) return io_error(ctx, path);
    setvbuf(reader->file, NULL, _IOFBF, 1024 * 1024);
    struct stat status;
    if (fstat(fileno(reader->file), &status)) return io_error(ctx, path);
    if ((uint64_t)status.st_size != size || stat_mtime(&status) != mtime)
        return value_error(ctx, "Source changed before processing: %s", path);
    if (start > end || end > size || end > INT64_MAX)
        return value_error(ctx, "Invalid source segment: %s", path);
    if (!start) {
        unsigned char magic[2];
        size_t count = fread(magic, 1, 2, reader->file);
        if (count == 2 && magic[0] == 0x1f && magic[1] == 0x8b)
            return value_error(ctx, "Input must be uncompressed: %s", path);
    }
    if (fseeko(reader->file, (off_t)start, SEEK_SET))
        return io_error(ctx, path);
    return 1;
}

int reader_next(Context *ctx, Reader *reader) {
    if (reader->offset >= reader->end) return 0;
    ssize_t length = getline(&reader->line, &reader->capacity, reader->file);
    if (length < 0) {
        if (ferror(reader->file)) io_error(ctx, reader->path);
        else value_error(ctx, "Source changed while reading: %s", reader->path);
        return -1;
    }
    if ((uint64_t)length > reader->end - reader->offset) {
        value_error(ctx, "Invalid source segment boundary: %s", reader->path);
        return -1;
    }
    reader->length = (size_t)length;
    return 1;
}

int reader_check(Context *ctx, Reader *reader) {
    struct stat fd_status, path_status;
    if (fstat(fileno(reader->file), &fd_status) ||
        stat(reader->path, &path_status)) return io_error(ctx, reader->path);
    if ((uint64_t)fd_status.st_size != reader->size ||
        (uint64_t)path_status.st_size != reader->size ||
        stat_mtime(&fd_status) != reader->mtime ||
        stat_mtime(&path_status) != reader->mtime)
        return value_error(ctx, "Source changed during processing: %s",
                           reader->path);
    return 1;
}

void reader_close(Reader *reader) {
    if (reader->file) fclose(reader->file);
    free(reader->line); free(reader->record.fields);
    reader->file = NULL;
    reader->line = NULL; reader->record.fields = NULL;
    reader->capacity = reader->record.capacity = 0;
}

FILE *output_fd(Context *ctx, int fd) {
    int copy = dup(fd);
    if (copy < 0) { io_error(ctx, "output"); return NULL; }
    FILE *file = fdopen(copy, "wb");
    if (!file) { close(copy); io_error(ctx, "output"); return NULL; }
    setvbuf(file, NULL, _IOFBF, 1024 * 1024);
    return file;
}

FILE *output_new(Context *ctx, const char *path) {
    int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd < 0) { io_error(ctx, path); return NULL; }
    FILE *file = fdopen(fd, "wb");
    if (!file) {
        int saved = errno; close(fd); unlink(path); errno = saved;
        io_error(ctx, path); return NULL;
    }
    setvbuf(file, NULL, _IOFBF, 1024 * 1024);
    return file;
}

int write_bytes(Context *ctx, FILE *file, const void *data, size_t length) {
    if (length && fwrite_unlocked(data, 1, length, file) != length)
        return io_error(ctx, "output");
    return 1;
}

int finish_output(Context *ctx, FILE *file) {
    if (file && fclose(file)) return io_error(ctx, "output");
    return !ctx->kind;
}

int write_block(Context *ctx, FILE *file, Span chromosome, uint64_t window,
                uint64_t start, uint64_t end) {
    if (!write_bytes(ctx, file, chromosome.data, chromosome.len)) return 0;
    if (fprintf(file, "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                window, start, end) < 0) return io_error(ctx, "index");
    return 1;
}

PyObject *cov_copy(PyObject *self, PyObject *args) {
    (void)self;
    const char *path; int fd;
    if (!PyArg_ParseTuple(args, "yi", &path, &fd)) return NULL;
    Context ctx = {0};
    uint64_t total = 0;
    Py_BEGIN_ALLOW_THREADS
    char *buffer = malloc(1024 * 1024);
    FILE *input = NULL, *output = NULL;
    if (!buffer) memory_error(&ctx);
    if (!ctx.kind) {
        input = fopen(path, "rb");
        if (!input) io_error(&ctx, path);
    }
    if (!ctx.kind) output = output_fd(&ctx, fd);
    size_t length;
    while (!ctx.kind && (length = fread(buffer, 1, 1024 * 1024, input))) {
        if (!write_bytes(&ctx, output, buffer, length) ||
            !tick(&ctx, length, 0)) break;
        total += length;
    }
    if (input) {
        if (ferror(input)) io_error(&ctx, path);
        fclose(input);
    }
    if (output) finish_output(&ctx, output);
    if (!ctx.kind) tick(&ctx, 0, 1);
    free(buffer);
    Py_END_ALLOW_THREADS
    if (ctx.kind) return raise_error(&ctx);
    return PyLong_FromUnsignedLongLong(total);
}

PyObject *cov_parse(PyObject *self, PyObject *args) {
    (void)self;
    const char *line;
    Py_ssize_t length;
    unsigned long long offset;
    if (!PyArg_ParseTuple(args, "y#K", &line, &length, &offset)) return NULL;
    Context ctx = {0}; Record record = {0};
    if (!fields(&ctx, &record, line, (size_t)length) ||
        !validate(&ctx, &record, offset)) {
        free(record.fields); return raise_error(&ctx);
    }
    PyObject *values = PyList_New((Py_ssize_t)record.count);
    if (!values) { free(record.fields); return NULL; }
    for (size_t i = 0; i < record.count; i++) {
        PyObject *field = PyBytes_FromStringAndSize(
            record.fields[i].data, (Py_ssize_t)record.fields[i].len);
        if (!field) {
            free(record.fields); Py_DECREF(values); return NULL;
        }
        PyList_SET_ITEM(values, (Py_ssize_t)i, field);
    }
    PyObject *result = Py_BuildValue("NK", values,
                                   (unsigned long long)record.position);
    free(record.fields);
    return result;
}
