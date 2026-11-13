#ifndef COVTOOLS_COMMON_H
#define COVTOOLS_COMMON_H

#define _GNU_SOURCE
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

typedef struct { const char *data; size_t len; } Span;
typedef struct {
    Span *fields;
    size_t count, capacity, length;
    uint64_t position;
} Record;
typedef struct {
    int kind, code;
    char message[512];
    PyObject *callback;
    uint64_t pending;
} Context;
typedef struct {
    FILE *file;
    const char *path;
    uint64_t offset, end, size;
    int64_t mtime;
    char *line;
    size_t capacity, length;
    Record record;
} Reader;

int value_error(Context *, const char *, ...);
int io_error(Context *, const char *);
int memory_error(Context *);
PyObject *raise_error(Context *);
int tick(Context *, uint64_t, int);
int fields(Context *, Record *, const char *, size_t);
int validate(Context *, Record *, uint64_t);
int unsigned_value(Span, uint64_t *);
int bytes_compare(Span, Span);
int span_equal(Span, const char *, size_t);
int utf8_valid(Span);
int64_t stat_mtime(const struct stat *);
int reader_open(Context *, Reader *, const char *, uint64_t, uint64_t,
                uint64_t, int64_t);
int reader_next(Context *, Reader *);
int reader_check(Context *, Reader *);
void reader_close(Reader *);
FILE *output_fd(Context *, int);
FILE *output_new(Context *, const char *);
int write_bytes(Context *, FILE *, const void *, size_t);
int finish_output(Context *, FILE *);
int write_block(Context *, FILE *, Span, uint64_t, uint64_t, uint64_t);

PyObject *cov_parse(PyObject *, PyObject *);
PyObject *cov_copy(PyObject *, PyObject *);
PyObject *cov_convert(PyObject *, PyObject *);
PyObject *cov_index(PyObject *, PyObject *);
PyObject *cov_index_join(PyObject *, PyObject *);
PyObject *cov_split(PyObject *, PyObject *);
PyObject *cov_sort_chunk(PyObject *, PyObject *);
PyObject *cov_merge_runs(PyObject *, PyObject *);
PyObject *cov_query(PyObject *, PyObject *);
extern PyTypeObject CovQueryType;
extern PyTypeObject CovIndexJoinType;

#endif
