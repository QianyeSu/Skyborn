#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void dinterp_to_isentropic_corder_into(
    void *data_flat,
    void *temp_flat,
    void *pressure_flat,
    void *theta_levels,
    void *output_flat,
    double p0,
    double kappa,
    double spvl,
    int kxtrp,
    int nouter,
    int nlev,
    int ninner,
    int ntheta
);

static PyArrayObject *to_double_1d(PyObject *obj, int writable) {
    int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS;
    if (writable) {
        flags |= NPY_ARRAY_WRITEABLE;
    }
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, flags);
}

static PyObject *py_dinterp_to_isentropic_corder_into(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *theta_obj = NULL;
    PyObject *output_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *theta_arr = NULL;
    PyArrayObject *output_arr = NULL;
    double p0 = 0.0;
    double kappa = 0.0;
    double spvl = 0.0;
    int kxtrp = 0;
    int nouter = 0;
    int nlev = 0;
    int ninner = 0;
    int ntheta = 0;

    (void) self;

    if (!PyArg_ParseTuple(
            args,
            "OOOOOdddiiii",
            &data_obj,
            &temp_obj,
            &pressure_obj,
            &theta_obj,
            &output_obj,
            &p0,
            &kappa,
            &spvl,
            &kxtrp,
            &nouter,
            &nlev,
            &ninner
        )) {
        return NULL;
    }

    data_arr = to_double_1d(data_obj, 0);
    temp_arr = to_double_1d(temp_obj, 0);
    pressure_arr = to_double_1d(pressure_obj, 0);
    theta_arr = to_double_1d(theta_obj, 0);
    output_arr = to_double_1d(output_obj, 1);
    if (data_arr == NULL || temp_arr == NULL || pressure_arr == NULL || theta_arr == NULL || output_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(data_arr) != 1 || PyArray_NDIM(temp_arr) != 1 || PyArray_NDIM(pressure_arr) != 1 ||
        PyArray_NDIM(theta_arr) != 1 || PyArray_NDIM(output_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "isentropic kernel expects flat 1D buffers");
        goto fail;
    }

    ntheta = (int) PyArray_DIM(theta_arr, 0);
    if ((int) PyArray_DIM(data_arr, 0) != nouter * nlev * ninner ||
        (int) PyArray_DIM(temp_arr, 0) != nouter * nlev * ninner ||
        (int) PyArray_DIM(pressure_arr, 0) != nouter * nlev * ninner) {
        PyErr_SetString(PyExc_ValueError, "input flat buffers must have length nouter*nlev*ninner");
        goto fail;
    }
    if ((int) PyArray_DIM(output_arr, 0) != nouter * ntheta * ninner) {
        PyErr_SetString(PyExc_ValueError, "output_flat must have length nouter*ntheta*ninner");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dinterp_to_isentropic_corder_into(
        PyArray_DATA(data_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(theta_arr),
        PyArray_DATA(output_arr),
        p0,
        kappa,
        spvl,
        kxtrp,
        nouter,
        nlev,
        ninner,
        ntheta
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(data_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(theta_arr);
    Py_DECREF(output_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(theta_arr);
    Py_XDECREF(output_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"dinterp_to_isentropic_corder_into", py_dinterp_to_isentropic_corder_into, METH_VARARGS, "In-place C-order interpolation to isentropic surfaces."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "isentropic_kernels",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_isentropic_kernels(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
