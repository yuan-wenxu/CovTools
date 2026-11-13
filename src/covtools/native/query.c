#include "common.h"

typedef struct {
    PyObject_HEAD
    FILE *index;
    Reader source;
    char *path, *chromosome, *barcode, *index_line;
    size_t chromosome_length, barcode_length, index_capacity;
    uint64_t first_window, last_window, start, end;
    uint64_t previous_end, block_count, expected_blocks;
    Record index_record;
    int finished, busy;
} Query;

static void query_close(Query *query) {
    if (query->index) fclose(query->index);
    query->index = NULL;
    reader_close(&query->source);
    query->finished = 1;
}

static void query_destroy(Query *query) {
    query_close(query);
    free(query->path); free(query->chromosome); free(query->barcode);
    free(query->index_line); free(query->index_record.fields);
    Py_TYPE(query)->tp_free((PyObject *)query);
}

static int next_block(Context *ctx, Query *query) {
    ssize_t length;
    while ((length = getline(&query->index_line, &query->index_capacity,
                              query->index)) >= 0) {
        if (length >= 10 && !memcmp(query->index_line, "#complete\t", 10)) {
            if (query->block_count != query->expected_blocks ||
                query->previous_end != query->source.size ||
                fgetc(query->index) != EOF) {
                value_error(ctx, "Invalid index block count or source extent");
                return -1;
            }
            if (!reader_check(ctx, &query->source)) return -1;
            return 0;
        }
        Record *record = &query->index_record;
        uint64_t window, start, end;
        if (!fields(ctx, record, query->index_line, (size_t)length)) return -1;
        if (record->count != 4 || !record->fields[0].len ||
            !utf8_valid(record->fields[0]) ||
            !unsigned_value(record->fields[1], &window) ||
            !unsigned_value(record->fields[2], &start) ||
            !unsigned_value(record->fields[3], &end) ||
            start != query->previous_end || start >= end ||
            end > query->source.size) {
            value_error(ctx, "Invalid index window or byte offsets"); return -1;
        }
        query->previous_end = end; query->block_count++;
        if (!tick(ctx, (uint64_t)length, 0)) return -1;
        if (span_equal(record->fields[0], query->chromosome,
                       query->chromosome_length) &&
            window >= query->first_window && window <= query->last_window) {
            query->source.offset = start; query->source.end = end;
            if (fseeko(query->source.file, (off_t)start, SEEK_SET)) {
                io_error(ctx, query->path); return -1;
            }
            return 1;
        }
    }
    if (ferror(query->index)) io_error(ctx, "index");
    else value_error(ctx, "Index build is incomplete");
    return -1;
}

static PyObject *query_next(Query *query) {
    if (query->finished) return NULL;
    if (query->busy) {
        PyErr_SetString(PyExc_ValueError,
                        "Query iterator is already executing");
        return NULL;
    }
    query->busy = 1;
    Context ctx = {0};
    int found = 0;
    Py_BEGIN_ALLOW_THREADS
    while (!ctx.kind && !found) {
        Reader *source = &query->source;
        if (source->offset >= source->end) {
            int available = next_block(&ctx, query);
            if (available <= 0) { query_close(query); break; }
        }
        if (reader_next(&ctx, source) < 1) {
            if (!ctx.kind) value_error(&ctx, "Source has invalid offsets");
            break;
        }
        Record *record = &source->record;
        if (!fields(&ctx, record, source->line, source->length) ||
            !validate(&ctx, record, source->offset)) break;
        source->offset += source->length;
        if (!tick(&ctx, source->length, 0)) break;
        found = span_equal(record->fields[0], query->chromosome,
                           query->chromosome_length) &&
                record->position >= query->start &&
                record->position <= query->end &&
                (!query->barcode || (record->count >= 7 &&
                 span_equal(record->fields[6], query->barcode,
                            query->barcode_length)));
    }
    if (ctx.kind) query_close(query);
    Py_END_ALLOW_THREADS
    query->busy = 0;
    if (ctx.kind) return raise_error(&ctx);
    if (!found) return NULL;
    PyObject *result = PyUnicode_DecodeUTF8(query->source.line,
        (Py_ssize_t)query->source.length, "strict");
    if (!result) query_close(query);
    return result;
}

PyTypeObject CovQueryType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "covtools._native.Query",
    .tp_basicsize = sizeof(Query),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)query_destroy,
    .tp_iter = PyObject_SelfIter,
    .tp_iternext = (iternextfunc)query_next,
};

PyObject *cov_query(PyObject *self, PyObject *args) {
    (void)self;
    const char *path, *index_path, *chromosome;
    unsigned long long start, end, size, bin_size, blocks;
    long long expected_mtime;
    PyObject *barcode;
    if (!PyArg_ParseTuple(args, "yyyKKOKLKK", &path, &index_path, &chromosome,
                          &start, &end, &barcode, &size, &expected_mtime,
                          &bin_size, &blocks)) return NULL;
    if (!start || end < start || !bin_size || !*chromosome) {
        PyErr_SetString(PyExc_ValueError, "Invalid query region or bin size");
        return NULL;
    }
    Query *query = (Query *)CovQueryType.tp_alloc(&CovQueryType, 0);
    if (!query) return NULL;
    query->path = strdup(path); query->chromosome = strdup(chromosome);
    query->chromosome_length = strlen(chromosome);
    query->start = start; query->end = end;
    query->first_window = (start - 1) / bin_size;
    query->last_window = (end - 1) / bin_size;
    query->expected_blocks = blocks;
    if (!query->path || !query->chromosome) {
        PyErr_NoMemory(); Py_DECREF(query); return NULL;
    }
    if (barcode != Py_None) {
        char *data; Py_ssize_t length;
        if (PyBytes_AsStringAndSize(barcode, &data, &length)) {
            Py_DECREF(query); return NULL;
        }
        query->barcode = malloc((size_t)length + 1);
        if (!query->barcode) {
            PyErr_NoMemory(); Py_DECREF(query); return NULL;
        }
        memcpy(query->barcode, data, (size_t)length);
        query->barcode[length] = 0; query->barcode_length = (size_t)length;
    }
    struct stat status;
    if (stat(path, &status)) {
        PyErr_SetFromErrnoWithFilename(PyExc_OSError, path);
        Py_DECREF(query); return NULL;
    }
    if ((uint64_t)status.st_size != size) {
        PyErr_SetString(PyExc_ValueError,
                        "Source size changed; build a new index");
        Py_DECREF(query); return NULL;
    }
    int64_t mtime = stat_mtime(&status);
    if (mtime != expected_mtime && PyErr_WarnEx(PyExc_UserWarning,
        "Source modification time differs from the index; continuing because "
        "the file size matches. Ensure the contents are unchanged.", 1) < 0) {
        Py_DECREF(query); return NULL;
    }
    Context ctx = {0};
    if (!reader_open(&ctx, &query->source, query->path, 0, 0, size, mtime)) {
        Py_DECREF(query); return raise_error(&ctx);
    }
    query->index = fopen(index_path, "rb");
    if (!query->index) {
        PyErr_SetFromErrnoWithFilename(PyExc_OSError, index_path);
        Py_DECREF(query); return NULL;
    }
    for (int i = 0; i < 2; i++) {
        if (getline(&query->index_line, &query->index_capacity,
                    query->index) < 0) {
            PyErr_SetString(PyExc_ValueError, "Invalid index header");
            Py_DECREF(query); return NULL;
        }
    }
    return (PyObject *)query;
}
