#include "common.h"

static PyMethodDef methods[] = {
    {"parse_record", cov_parse, METH_VARARGS, "Validate one coverage record."},
    {"copy", cov_copy, METH_VARARGS, "Copy an owned temporary file in C."},
    {"convert", cov_convert, METH_VARARGS, "Convert a source byte range."},
    {"index", cov_index, METH_VARARGS, "Index a source byte range."},
    {"index_join", cov_index_join, METH_VARARGS, "Join partial indexes."},
    {"split", cov_split, METH_VARARGS, "Split a source byte range."},
    {"sort_chunk", cov_sort_chunk, METH_VARARGS, "Sort a source byte range."},
    {"merge_runs", cov_merge_runs, METH_VARARGS, "Merge sorted native runs."},
    {"query", cov_query, METH_VARARGS, "Open a lazy indexed query."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_native",
    .m_doc = "C coverage kernels; long-running operations release the GIL.",
    .m_size = -1,
    .m_methods = methods,
};

PyMODINIT_FUNC PyInit__native(void) {
    if (PyType_Ready(&CovQueryType) < 0 ||
        PyType_Ready(&CovIndexJoinType) < 0) return NULL;
    return PyModule_Create(&module);
}
