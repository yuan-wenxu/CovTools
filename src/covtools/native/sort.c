#include "common.h"

typedef struct {
    char *data;
    size_t offset, length, chromosome_length, barcode_offset, barcode_length;
    uint64_t position, order;
} SortRow;

static int compare_rows(const void *left, const void *right) {
    const SortRow *a = left, *b = right;
    int comparison = bytes_compare(
        (Span){a->data, a->chromosome_length},
        (Span){b->data, b->chromosome_length});
    if (comparison) return comparison;
    if (a->position != b->position)
        return a->position > b->position ? 1 : -1;
    comparison = bytes_compare(
        (Span){a->data + a->barcode_offset, a->barcode_length},
        (Span){b->data + b->barcode_offset, b->barcode_length});
    if (comparison) return comparison;
    return (a->order > b->order) - (a->order < b->order);
}

PyObject *cov_sort_chunk(PyObject *self, PyObject *args) {
    (void)self;
    const char *path, *destination;
    unsigned long long start, end, size;
    long long mtime;
    if (!PyArg_ParseTuple(args, "yyKKKL", &path, &destination, &start,
                          &end, &size, &mtime)) return NULL;
    Context ctx = {0}; Reader reader = {0};
    char *arena = NULL;
    SortRow *rows = NULL;
    size_t count = 0, capacity = 0, used = 0;
    FILE *output = NULL;
    int created = 0;
    Py_BEGIN_ALLOW_THREADS
    if (reader_open(&ctx, &reader, path, start, end, size, mtime)) {
        if (end - start >= SIZE_MAX) memory_error(&ctx);
        else {
            arena = malloc((size_t)(end - start) + 1);
            if (!arena) memory_error(&ctx);
        }
    }
    while (!ctx.kind && reader_next(&ctx, &reader) > 0) {
        Record *record = &reader.record;
        if (!fields(&ctx, record, reader.line, reader.length) ||
            !validate(&ctx, record, reader.offset)) break;
        if (count == capacity) {
            size_t next = capacity ? capacity * 2 : 1024;
            if (next > SIZE_MAX / sizeof(SortRow)) {
                memory_error(&ctx); break;
            }
            SortRow *resized = realloc(rows, next * sizeof(SortRow));
            if (!resized) { memory_error(&ctx); break; }
            rows = resized; capacity = next;
        }
        SortRow *row = &rows[count];
        *row = (SortRow){.offset = used, .length = record->length + 1,
                         .chromosome_length = record->fields[0].len,
                         .position = record->position, .order = count};
        if (record->count >= 7) {
            row->barcode_offset =
                (size_t)(record->fields[6].data - reader.line);
            row->barcode_length = record->fields[6].len;
        }
        memcpy(arena + used, reader.line, record->length);
        arena[used + record->length] = '\n';
        used += row->length; count++;
        reader.offset += reader.length;
        if (!tick(&ctx, reader.length, 0)) break;
    }
    if (!ctx.kind && !count) value_error(&ctx, "The coverage file is empty");
    if (!ctx.kind) reader_check(&ctx, &reader);
    if (!ctx.kind) {
        for (size_t i = 0; i < count; i++)
            rows[i].data = arena + rows[i].offset;
        qsort(rows, count, sizeof(SortRow), compare_rows);
        if (tick(&ctx, 0, 1)) {
            output = output_new(&ctx, destination);
            created = output != NULL;
        }
    }
    for (size_t i = 0; output && !ctx.kind && i < count; i++) {
        write_bytes(&ctx, output, rows[i].data, rows[i].length);
        if (!tick(&ctx, rows[i].length, 0)) break;
    }
    if (output) finish_output(&ctx, output);
    if (!ctx.kind) reader_check(&ctx, &reader);
    if (!ctx.kind) tick(&ctx, 0, 1);
    if (ctx.kind && created) unlink(destination);
    free(rows); free(arena); reader_close(&reader);
    Py_END_ALLOW_THREADS
    if (ctx.kind) return raise_error(&ctx);
    return PyLong_FromUnsignedLongLong(end - start);
}

typedef struct {
    FILE *file;
    const char *path;
    char *line;
    size_t capacity, length;
    Span chromosome, barcode;
    uint64_t position;
    Record record;
} Cursor;

static int cursor_next(Context *ctx, Cursor *cursor) {
    ssize_t length = getline(&cursor->line, &cursor->capacity, cursor->file);
    if (length < 0) {
        if (ferror(cursor->file)) io_error(ctx, cursor->path);
        return 0;
    }
    cursor->length = (size_t)length;
    Record *record = &cursor->record;
    if (!fields(ctx, record, cursor->line, cursor->length)) return 0;
    if (record->count < 6 || record->count > 9 ||
        !unsigned_value(record->fields[1], &cursor->position)) {
        value_error(ctx, "Invalid sorted run: %s", cursor->path); return 0;
    }
    cursor->chromosome = record->fields[0];
    cursor->barcode = record->count >= 7 ? record->fields[6] : (Span){"", 0};
    return 1;
}

static int cursor_less(Cursor *cursors, size_t a, size_t b) {
    int comparison = bytes_compare(cursors[a].chromosome,
                                   cursors[b].chromosome);
    if (comparison) return comparison < 0;
    if (cursors[a].position != cursors[b].position)
        return cursors[a].position < cursors[b].position;
    comparison = bytes_compare(cursors[a].barcode, cursors[b].barcode);
    return comparison ? comparison < 0 : a < b;
}

static void sift(Cursor *cursors, size_t *heap, size_t count, size_t parent) {
    for (;;) {
        size_t child = parent * 2 + 1;
        if (child >= count) return;
        if (child + 1 < count &&
            cursor_less(cursors, heap[child + 1], heap[child])) child++;
        if (!cursor_less(cursors, heap[child], heap[parent])) return;
        size_t swap = heap[parent]; heap[parent] = heap[child];
        heap[child] = swap; parent = child;
    }
}

PyObject *cov_merge_runs(PyObject *self, PyObject *args) {
    (void)self;
    PyObject *paths, *callback = Py_None;
    int fd;
    if (!PyArg_ParseTuple(args, "Oi|O", &paths, &fd, &callback)) return NULL;
    PyObject *sequence = PySequence_Fast(paths, "Expected a list of run paths");
    if (!sequence) return NULL;
    size_t count = (size_t)PySequence_Fast_GET_SIZE(sequence);
    Cursor *cursors = calloc(count ? count : 1, sizeof(Cursor));
    size_t *heap = malloc((count ? count : 1) * sizeof(size_t));
    Context ctx = {.callback = callback};
    if (!cursors || !heap) memory_error(&ctx);
    for (size_t i = 0; !ctx.kind && i < count; i++) {
        char *path = PyBytes_AsString(PySequence_Fast_GET_ITEM(sequence, i));
        if (!path) { ctx.kind = 4; break; }
        cursors[i].path = strdup(path);
        if (!cursors[i].path) memory_error(&ctx);
    }
    FILE *output = NULL;
    size_t heap_count = 0;
    Py_BEGIN_ALLOW_THREADS
    if (!ctx.kind) output = output_fd(&ctx, fd);
    for (size_t i = 0; !ctx.kind && i < count; i++) {
        cursors[i].file = fopen(cursors[i].path, "rb");
        if (!cursors[i].file) { io_error(&ctx, cursors[i].path); break; }
        setvbuf(cursors[i].file, NULL, _IOFBF, 256 * 1024);
        if (cursor_next(&ctx, &cursors[i])) heap[heap_count++] = i;
    }
    for (size_t i = heap_count / 2; i > 0; i--)
        sift(cursors, heap, heap_count, i - 1);
    while (!ctx.kind && heap_count) {
        size_t index = heap[0];
        Cursor *cursor = &cursors[index];
        if (!write_bytes(&ctx, output, cursor->line, cursor->length) ||
            !tick(&ctx, cursor->length, 0)) break;
        if (!cursor_next(&ctx, cursor)) heap[0] = heap[--heap_count];
        if (heap_count) sift(cursors, heap, heap_count, 0);
    }
    if (output) finish_output(&ctx, output);
    if (!ctx.kind) tick(&ctx, 0, 1);
    for (size_t i = 0; cursors && i < count; i++) {
        if (cursors[i].file) fclose(cursors[i].file);
        free(cursors[i].line); free(cursors[i].record.fields);
        free((char *)cursors[i].path);
    }
    free(cursors); free(heap);
    Py_END_ALLOW_THREADS
    Py_DECREF(sequence);
    if (ctx.kind) return raise_error(&ctx);
    Py_RETURN_NONE;
}
