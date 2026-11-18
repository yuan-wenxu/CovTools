#include "common.h"
#include <fcntl.h>

static int metadata_options(PyObject *columns, PyObject *values,
                            int selected[3], Span literals[3]) {
    const char *names[] = {"barcode", "strand", "context"};
    for (int i = 0; i < 3; i++) {
        selected[i] = -1;
        literals[i] = (Span){NULL, 0};
    }
    for (int i = 0; i < 3; i++) {
        PyObject *column = PyDict_GetItemString(columns, names[i]);
        if (column) {
            long number = PyLong_AsLong(column);
            if (PyErr_Occurred()) goto failed;
            if (number < 1 || number > INT_MAX) {
                PyErr_SetString(PyExc_ValueError, "Invalid metadata column");
                goto failed;
            }
            selected[i] = (int)number - 1;
        }
        PyObject *value = PyDict_GetItemString(values, names[i]);
        if (value) {
            char *data; Py_ssize_t length;
            if (PyBytes_AsStringAndSize(value, &data, &length)) goto failed;
            char *copy = malloc(length ? (size_t)length : 1);
            if (!copy) { PyErr_NoMemory(); goto failed; }
            memcpy(copy, data, (size_t)length);
            literals[i] = (Span){copy, (size_t)length};
        }
    }
    return 1;
failed:
    for (int i = 0; i < 3; i++) free((char *)literals[i].data);
    return 0;
}

PyObject *cov_convert(PyObject *self, PyObject *args) {
    (void)self;
    const char *path;
    int fd;
    unsigned long long start, end, size;
    long long mtime;
    PyObject *columns, *values, *callback = Py_None;
    if (!PyArg_ParseTuple(args, "yiKKKLOO|O", &path, &fd, &start, &end,
                          &size, &mtime, &columns, &values, &callback))
        return NULL;
    int selected[3]; Span literals[3];
    if (!metadata_options(columns, values, selected, literals)) return NULL;
    Context ctx = {.callback = callback}; Reader reader = {0};
    FILE *output = NULL;
    uint64_t rows = 0;
    const char *names[] = {"barcode", "strand", "context"};
    Py_BEGIN_ALLOW_THREADS
    if (reader_open(&ctx, &reader, path, start, end, size, mtime))
        output = output_fd(&ctx, fd);
    while (output && !ctx.kind && reader_next(&ctx, &reader) > 0) {
        if (!fields(&ctx, &reader.record, reader.line, reader.length)) break;
        Record *input = &reader.record;
        if (input->count < 6) {
            value_error(&ctx, "Expected at least six columns at byte %"
                        PRIu64, reader.offset); break;
        }
        Span normalized[9];
        memcpy(normalized, input->fields, 6 * sizeof(Span));
        for (int i = 0; i < 3; i++) {
            Span value = {".", 1};
            if (selected[i] >= 0) {
                if ((size_t)selected[i] >= input->count) {
                    value_error(&ctx, "Expected at least %d columns at byte "
                                "%" PRIu64, selected[i] + 1, reader.offset);
                    break;
                }
                value = input->fields[selected[i]];
            }
            if (literals[i].data) {
                if (!span_equal(value, ".", 1) &&
                    bytes_compare(value, literals[i])) {
                    value_error(&ctx, "Conflicting %s at byte %" PRIu64,
                                names[i], reader.offset); break;
                }
                value = literals[i];
            }
            normalized[i + 6] = value;
        }
        if (ctx.kind) break;
        Record record = {.fields = normalized, .count = 9};
        if (!validate(&ctx, &record, reader.offset)) break;
        size_t prefix = (size_t)(normalized[5].data - reader.line) +
                        normalized[5].len;
        if (!write_bytes(&ctx, output, reader.line, prefix)) break;
        for (int i = 6; i < 9 && !ctx.kind; i++) {
            write_bytes(&ctx, output, "\t", 1);
            write_bytes(&ctx, output, normalized[i].data, normalized[i].len);
        }
        write_bytes(&ctx, output, "\n", 1);
        reader.offset += reader.length; rows++;
        if (!tick(&ctx, reader.length, 0)) break;
    }
    if (!ctx.kind && !rows) value_error(&ctx, "The coverage file is empty");
    if (!ctx.kind) reader_check(&ctx, &reader);
    if (output) finish_output(&ctx, output);
    if (!ctx.kind) tick(&ctx, 0, 1);
    reader_close(&reader);
    for (int i = 0; i < 3; i++) free((char *)literals[i].data);
    Py_END_ALLOW_THREADS
    if (ctx.kind) return raise_error(&ctx);
    return PyLong_FromUnsignedLongLong(rows);
}

PyObject *cov_index(PyObject *self, PyObject *args) {
    (void)self;
    const char *path;
    int fd;
    unsigned long long start, end, size, bin_size;
    long long mtime;
    PyObject *callback = Py_None;
    if (!PyArg_ParseTuple(args, "yiKKKLK|O", &path, &fd, &start, &end,
                          &size, &mtime, &bin_size, &callback)) return NULL;
    if (!bin_size) {
        PyErr_SetString(PyExc_ValueError, "bin_size must be positive");
        return NULL;
    }
    Context ctx = {.callback = callback}; Reader reader = {0};
    FILE *output = NULL;
    char *chromosome = NULL, *first_chromosome = NULL;
    size_t chromosome_length = 0, first_chromosome_length = 0;
    uint64_t window = 0, block_start = start, rows = 0, blocks = 0;
    uint64_t first_position = 0, previous_position = 0;
    Py_BEGIN_ALLOW_THREADS
    if (reader_open(&ctx, &reader, path, start, end, size, mtime))
        output = output_fd(&ctx, fd);
    while (output && !ctx.kind && reader_next(&ctx, &reader) > 0) {
        Record *record = &reader.record;
        if (!fields(&ctx, record, reader.line, reader.length) ||
            !validate(&ctx, record, reader.offset)) break;
        Span current = record->fields[0];
        if (!utf8_valid(current)) {
            value_error(&ctx, "Invalid UTF-8 chromosome at byte %" PRIu64,
                        reader.offset); break;
        }
        int order = rows ? bytes_compare(
            (Span){chromosome, chromosome_length}, current) : 0;
        if (order > 0 || (!order && rows &&
                          record->position < previous_position)) {
            value_error(&ctx, "Coverage must be sorted by chromosome and "
                        "position; run covtools sort (byte %" PRIu64 ")",
                        reader.offset); break;
        }
        if (!rows) {
            first_chromosome = malloc(current.len);
            if (!first_chromosome) { memory_error(&ctx); break; }
            memcpy(first_chromosome, current.data, current.len);
            first_chromosome_length = current.len;
            first_position = record->position;
        }
        uint64_t current_window = (record->position - 1) / bin_size;
        if (!rows || window != current_window ||
            !span_equal(current, chromosome, chromosome_length)) {
            if (rows) {
                write_block(&ctx, output,
                            (Span){chromosome, chromosome_length}, window,
                            block_start, reader.offset);
                blocks++;
            }
            char *copy = realloc(chromosome, current.len ? current.len : 1);
            if (!copy) { memory_error(&ctx); break; }
            chromosome = copy; chromosome_length = current.len;
            memcpy(chromosome, current.data, current.len);
            window = current_window; block_start = reader.offset;
        }
        previous_position = record->position;
        reader.offset += reader.length; rows++;
        if (!tick(&ctx, reader.length, 0)) break;
    }
    if (!ctx.kind && !rows) value_error(&ctx, "The coverage file is empty");
    if (!ctx.kind) {
        write_block(&ctx, output, (Span){chromosome, chromosome_length},
                    window, block_start, reader.offset);
        blocks++;
        reader_check(&ctx, &reader);
    }
    if (output) finish_output(&ctx, output);
    if (!ctx.kind) tick(&ctx, 0, 1);
    reader_close(&reader);
    Py_END_ALLOW_THREADS
    PyObject *result = ctx.kind ? raise_error(&ctx) : Py_BuildValue(
        "KK(y#K)(y#K)", (unsigned long long)rows,
        (unsigned long long)blocks,
        first_chromosome, (Py_ssize_t)first_chromosome_length,
        (unsigned long long)first_position,
        chromosome, (Py_ssize_t)chromosome_length,
        (unsigned long long)previous_position);
    free(first_chromosome); free(chromosome);
    return result;
}

typedef struct {
    PyObject_HEAD
    FILE *output;
    char *chromosome;
    size_t chromosome_length;
    uint64_t window, start, end, blocks;
    int pending;
} IndexJoin;

static PyObject *index_add(IndexJoin *join, PyObject *args) {
    const char *path;
    if (!PyArg_ParseTuple(args, "y", &path)) return NULL;
    if (!join->output) {
        PyErr_SetString(PyExc_ValueError, "Index writer is closed");
        return NULL;
    }
    Context ctx = {0}; Record record = {0};
    char *line = NULL; size_t capacity = 0;
    Py_BEGIN_ALLOW_THREADS
    FILE *input = fopen(path, "rb");
    if (!input) io_error(&ctx, path);
    ssize_t length;
    while (input && !ctx.kind &&
           (length = getline(&line, &capacity, input)) >= 0) {
        uint64_t window, start, end;
        if (!fields(&ctx, &record, line, (size_t)length)) break;
        if (record.count != 4 ||
            !unsigned_value(record.fields[1], &window) ||
            !unsigned_value(record.fields[2], &start) ||
            !unsigned_value(record.fields[3], &end)) {
            value_error(&ctx, "Invalid partial index"); break;
        }
        Span chromosome = record.fields[0];
        int matches = join->pending && window == join->window &&
            span_equal(chromosome, join->chromosome, join->chromosome_length);
        if (!matches) {
            if (join->pending) {
                write_block(&ctx, join->output,
                            (Span){join->chromosome, join->chromosome_length},
                            join->window, join->start, join->end);
                join->blocks++;
            }
            char *copy = realloc(join->chromosome, chromosome.len);
            if (!copy) { memory_error(&ctx); break; }
            join->chromosome = copy;
            memcpy(copy, chromosome.data, chromosome.len);
            join->chromosome_length = chromosome.len;
            join->window = window; join->start = start; join->pending = 1;
        } else if (join->end != start) {
            value_error(&ctx, "Noncontiguous partial index"); break;
        }
        join->end = end;
        if (!tick(&ctx, (uint64_t)length, 0)) break;
    }
    if (input) {
        if (ferror(input)) io_error(&ctx, path);
        fclose(input);
    }
    free(line); free(record.fields);
    if (!ctx.kind) tick(&ctx, 0, 1);
    Py_END_ALLOW_THREADS
    if (ctx.kind) return raise_error(&ctx);
    Py_RETURN_NONE;
}

static PyObject *index_finish(IndexJoin *join, PyObject *unused) {
    (void)unused;
    if (!join->output) {
        PyErr_SetString(PyExc_ValueError, "Index writer is closed");
        return NULL;
    }
    Context ctx = {0};
    Py_BEGIN_ALLOW_THREADS
    if (join->pending) {
        write_block(&ctx, join->output,
                    (Span){join->chromosome, join->chromosome_length},
                    join->window, join->start, join->end);
        join->blocks++; join->pending = 0;
    }
    finish_output(&ctx, join->output); join->output = NULL;
    Py_END_ALLOW_THREADS
    if (ctx.kind) return raise_error(&ctx);
    return PyLong_FromUnsignedLongLong(join->blocks);
}

static void index_destroy(IndexJoin *join) {
    if (join->output) fclose(join->output);
    free(join->chromosome);
    Py_TYPE(join)->tp_free((PyObject *)join);
}

static PyObject *index_close(IndexJoin *join, PyObject *unused) {
    (void)unused;
    if (join->output) fclose(join->output);
    join->output = NULL; join->pending = 0;
    Py_RETURN_NONE;
}

static PyMethodDef index_methods[] = {
    {"add", (PyCFunction)index_add, METH_VARARGS, "Append partial index."},
    {"finish", (PyCFunction)index_finish, METH_NOARGS, "Flush pending window."},
    {"close", (PyCFunction)index_close, METH_NOARGS, "Close partial writer."},
    {NULL, NULL, 0, NULL}
};

PyTypeObject CovIndexJoinType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "covtools._native.IndexJoin",
    .tp_basicsize = sizeof(IndexJoin),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)index_destroy,
    .tp_methods = index_methods,
};

PyObject *cov_index_join(PyObject *self, PyObject *args) {
    (void)self;
    int fd;
    if (!PyArg_ParseTuple(args, "i", &fd)) return NULL;
    IndexJoin *join = PyObject_New(IndexJoin, &CovIndexJoinType);
    if (!join) return NULL;
    join->output = NULL; join->chromosome = NULL;
    join->chromosome_length = 0; join->blocks = 0; join->pending = 0;
    Context ctx = {0};
    join->output = output_fd(&ctx, fd);
    if (!join->output) { Py_DECREF(join); return raise_error(&ctx); }
    return (PyObject *)join;
}

typedef struct Cell {
    struct Cell *hash_next, *all_next, *previous, *next;
    char *barcode, *name;
    size_t length;
    uint64_t hash;
    FILE *file;
    int created;
} Cell;

typedef struct {
    Cell **buckets, *all, *last, *head, *tail;
    size_t capacity, count, opened, max_open;
    int directory;
} Cells;

static uint64_t hash_span(Span key) {
    uint64_t hash = 14695981039346656037ULL;
    for (size_t i = 0; i < key.len; i++)
        hash = (hash ^ (unsigned char)key.data[i]) * 1099511628211ULL;
    return hash;
}

static char *encoded_name(Context *ctx, Span barcode) {
    if (barcode.len > (SIZE_MAX - 5) / 3) {
        memory_error(ctx); return NULL;
    }
    char *name = malloc(barcode.len * 3 + 5);
    if (!name) { memory_error(ctx); return NULL; }
    const char *hex = "0123456789ABCDEF";
    size_t length = 0;
    for (size_t i = 0; i < barcode.len; i++) {
        unsigned char c = barcode.data[i];
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || (c && strchr("_.-~", c)))
            name[length++] = (char)c;
        else {
            name[length++] = '%'; name[length++] = hex[c >> 4];
            name[length++] = hex[c & 15];
        }
    }
    memcpy(name + length, ".cov", 5);
    return name;
}

static void unlink_lru(Cells *cells, Cell *cell) {
    if (cell->previous) cell->previous->next = cell->next;
    else cells->head = cell->next;
    if (cell->next) cell->next->previous = cell->previous;
    else cells->tail = cell->previous;
    cell->previous = cell->next = NULL;
}

static Cell *cell_get(Context *ctx, Cells *cells, Span barcode) {
    uint64_t hash = hash_span(barcode);
    size_t bucket = hash % cells->capacity;
    Cell *cell = cells->buckets[bucket];
    while (cell && (cell->hash != hash || cell->length != barcode.len ||
           memcmp(cell->barcode, barcode.data, barcode.len)))
        cell = cell->hash_next;
    if (!cell) {
        if (!utf8_valid(barcode)) {
            value_error(ctx, "Invalid UTF-8 barcode"); return NULL;
        }
        cell = calloc(1, sizeof(Cell));
        if (!cell) { memory_error(ctx); return NULL; }
        cell->barcode = malloc(barcode.len);
        cell->name = encoded_name(ctx, barcode);
        if (!cell->barcode || !cell->name) {
            free(cell->barcode); free(cell->name); free(cell);
            memory_error(ctx); return NULL;
        }
        memcpy(cell->barcode, barcode.data, barcode.len);
        cell->length = barcode.len; cell->hash = hash;
        cell->hash_next = cells->buckets[bucket];
        cells->buckets[bucket] = cell;
        if (cells->last) cells->last->all_next = cell;
        else cells->all = cell;
        cells->last = cell; cells->count++;
        if (cells->count > cells->capacity * 3 / 4) {
            size_t capacity = cells->capacity * 2;
            Cell **buckets = calloc(capacity, sizeof(Cell *));
            if (!buckets) { memory_error(ctx); return NULL; }
            for (Cell *entry = cells->all; entry; entry = entry->all_next) {
                size_t index = entry->hash % capacity;
                entry->hash_next = buckets[index]; buckets[index] = entry;
            }
            free(cells->buckets); cells->buckets = buckets;
            cells->capacity = capacity;
        }
    }
    if (cell->file) unlink_lru(cells, cell);
    else {
        if (cells->opened == cells->max_open) {
            Cell *oldest = cells->tail;
            unlink_lru(cells, oldest);
            finish_output(ctx, oldest->file); oldest->file = NULL;
            cells->opened--;
            if (ctx->kind) return NULL;
        }
        int flags = O_WRONLY | (cell->created ? O_APPEND : O_CREAT | O_EXCL);
        int fd = openat(cells->directory, cell->name, flags, 0666);
        if (fd < 0) { io_error(ctx, cell->name); return NULL; }
        cell->created = 1;
        cell->file = fdopen(fd, "wb");
        if (!cell->file) {
            close(fd); io_error(ctx, cell->name); return NULL;
        }
        setvbuf(cell->file, NULL, _IOFBF, 64 * 1024);
        cells->opened++;
    }
    cell->next = cells->head;
    if (cells->head) cells->head->previous = cell;
    else cells->tail = cell;
    cells->head = cell;
    return cell;
}

static void cells_free(Cells *cells, int remove_files) {
    Cell *cell = cells->all;
    while (cell) {
        Cell *next = cell->all_next;
        if (cell->file) fclose(cell->file);
        if (remove_files && cell->created)
            unlinkat(cells->directory, cell->name, 0);
        free(cell->barcode); free(cell->name); free(cell);
        cell = next;
    }
    free(cells->buckets);
    if (cells->directory >= 0) close(cells->directory);
}

PyObject *cov_split(PyObject *self, PyObject *args) {
    (void)self;
    const char *path, *directory;
    unsigned long long start, end, size;
    long long mtime;
    int max_open, names = 1;
    PyObject *callback = Py_None;
    if (!PyArg_ParseTuple(args, "yyKKKLi|pO", &path, &directory, &start,
                          &end, &size, &mtime, &max_open, &names, &callback))
        return NULL;
    if (max_open < 1) {
        PyErr_SetString(PyExc_ValueError, "max_open_files must be positive");
        return NULL;
    }
    Context ctx = {.callback = callback}; Reader reader = {0};
    Cells cells = {.capacity = 1024, .max_open = (size_t)max_open,
                   .directory = -1};
    uint64_t rows = 0;
    Py_BEGIN_ALLOW_THREADS
    cells.buckets = calloc(cells.capacity, sizeof(Cell *));
    if (!cells.buckets) memory_error(&ctx);
    if (!ctx.kind) {
        cells.directory = open(directory, O_RDONLY | O_DIRECTORY);
        if (cells.directory < 0) io_error(&ctx, directory);
    }
    if (!ctx.kind)
        reader_open(&ctx, &reader, path, start, end, size, mtime);
    while (!ctx.kind && reader_next(&ctx, &reader) > 0) {
        Record *record = &reader.record;
        if (!fields(&ctx, record, reader.line, reader.length) ||
            !validate(&ctx, record, reader.offset)) break;
        if (record->count < 7 || span_equal(record->fields[6], ".", 1)) {
            value_error(&ctx, "Missing barcode at byte %" PRIu64,
                        reader.offset); break;
        }
        Span barcode = record->fields[6];
        Cell *cell = cell_get(&ctx, &cells, barcode);
        if (!cell) break;
        size_t prefix = (size_t)(barcode.data - reader.line);
        size_t suffix_start = prefix + barcode.len;
        write_bytes(&ctx, cell->file, reader.line, prefix);
        write_bytes(&ctx, cell->file, ".", 1);
        write_bytes(&ctx, cell->file, reader.line + suffix_start,
                    record->length - suffix_start);
        write_bytes(&ctx, cell->file, "\n", 1);
        reader.offset += reader.length; rows++;
        if (!tick(&ctx, reader.length, 0)) break;
    }
    if (!ctx.kind && !rows) value_error(&ctx, "The coverage file is empty");
    if (!ctx.kind) reader_check(&ctx, &reader);
    for (Cell *cell = cells.all; cell; cell = cell->all_next) {
        if (cell->file) {
            finish_output(&ctx, cell->file); cell->file = NULL;
        }
    }
    if (!ctx.kind) tick(&ctx, 0, 1);
    reader_close(&reader);
    Py_END_ALLOW_THREADS
    PyObject *result = NULL;
    if (!ctx.kind && names) {
        result = PyDict_New();
        if (!result) ctx.kind = 4;
        for (Cell *cell = cells.all; !ctx.kind && cell;
             cell = cell->all_next) {
            PyObject *barcode = PyUnicode_DecodeUTF8(
                cell->barcode, (Py_ssize_t)cell->length, "strict");
            PyObject *name = PyUnicode_FromString(cell->name);
            if (!barcode || !name || PyDict_SetItem(result, barcode, name))
                ctx.kind = 4;
            Py_XDECREF(barcode); Py_XDECREF(name);
        }
    } else if (!ctx.kind) { result = Py_None; Py_INCREF(result); }
    cells_free(&cells, ctx.kind != 0);
    if (ctx.kind) { Py_XDECREF(result); return raise_error(&ctx); }
    return result;
}
