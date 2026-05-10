#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void dinterp_pressure_1d(
    void *ppin,
    void *xxin,
    void *ppout,
    void *xxout,
    int linlog,
    double xmsg,
    void *ier,
    int npin,
    int npout
);

static PyArrayObject *to_double_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
}

static PyObject *py_dinterp_pressure_1d(PyObject *self, PyObject *args) {
    PyObject *ppin_obj = NULL;
    PyObject *xxin_obj = NULL;
    PyObject *ppout_obj = NULL;
    PyArrayObject *ppin_arr = NULL;
    PyArrayObject *xxin_arr = NULL;
    PyArrayObject *ppout_arr = NULL;
    PyArrayObject *xxout_arr = NULL;
    int linlog = 0;
    int npin = 0;
    int npout = 0;
    int ier = 0;
    double xmsg = 0.0;
    npy_intp out_dims[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOid", &ppin_obj, &xxin_obj, &ppout_obj, &linlog, &xmsg)) {
        return NULL;
    }

    ppin_arr = to_double_1d(ppin_obj);
    xxin_arr = to_double_1d(xxin_obj);
    ppout_arr = to_double_1d(ppout_obj);
    if (ppin_arr == NULL || xxin_arr == NULL || ppout_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(ppin_arr) != 1 || PyArray_NDIM(xxin_arr) != 1 || PyArray_NDIM(ppout_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "dinterp_pressure_1d expects three 1D arrays");
        goto fail;
    }

    npin = (int) PyArray_DIM(ppin_arr, 0);
    npout = (int) PyArray_DIM(ppout_arr, 0);
    if ((int) PyArray_DIM(xxin_arr, 0) != npin) {
        PyErr_SetString(PyExc_ValueError, "ppin and xxin must have the same length");
        goto fail;
    }

    out_dims[0] = npout;
    xxout_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_FLOAT64, 1);
    if (xxout_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dinterp_pressure_1d(
        PyArray_DATA(ppin_arr),
        PyArray_DATA(xxin_arr),
        PyArray_DATA(ppout_arr),
        PyArray_DATA(xxout_arr),
        linlog,
        xmsg,
        &ier,
        npin,
        npout
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(ppin_arr);
    Py_DECREF(xxin_arr);
    Py_DECREF(ppout_arr);
    {
        PyObject *result = Py_BuildValue("Oi", xxout_arr, ier);
        Py_DECREF(xxout_arr);
        return result;
    }

fail:
    Py_XDECREF(ppin_arr);
    Py_XDECREF(xxin_arr);
    Py_XDECREF(ppout_arr);
    Py_XDECREF(xxout_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"dinterp_pressure_1d", py_dinterp_pressure_1d, METH_VARARGS, "Interpolate one pressure profile and return (xxout, ier)."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "int2p_kernels",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_int2p_kernels(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
