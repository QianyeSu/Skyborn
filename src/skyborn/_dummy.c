// Dummy C extension to force platform wheel generation
#include <Python.h>

static PyObject* dummy_function(PyObject* self, PyObject* args) {
    Py_RETURN_NONE;
}

static PyMethodDef dummy_methods[] = {
    {"dummy_function", dummy_function, METH_VARARGS, "Dummy function"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef dummy_module = {
    PyModuleDef_HEAD_INIT,
    "_dummy",
    "Dummy module for platform wheel generation",
    -1,
    dummy_methods
};

PyMODINIT_FUNC PyInit__dummy(void) {
    return PyModule_Create(&dummy_module);
}
