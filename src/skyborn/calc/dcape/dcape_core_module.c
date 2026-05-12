#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void dcape_profile_c(
    int nlev,
    void *pressure_hpa,
    void *temperature_c,
    void *dewpoint_c,
    void *dcape
);
void dcape_grid_c(
    int nlev,
    int nlat,
    int nlon,
    void *pressure_3d,
    void *t_3d,
    void *td_3d,
    void *out
);

static PyObject *py_dcape_profile(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    double dcape = 0.0;
    int nlev;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOO",
            kwlist,
            &pressure_obj,
            &temperature_obj,
            &dewpoint_obj)) {
        return NULL;
    }

    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );

    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL) {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }

    if (PyArray_NDIM(pressure_arr) != 1 || PyArray_NDIM(temperature_arr) != 1 ||
        PyArray_NDIM(dewpoint_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must be 1D arrays");
        goto fail;
    }

    nlev = (int) PyArray_DIM(pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a length");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dcape_profile_c(
        nlev,
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temperature_arr),
        PyArray_DATA(dewpoint_arr),
        &dcape
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return PyFloat_FromDouble(dcape);

fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    return NULL;
}

static PyObject *py_dcape_grid(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *out_arr = NULL;
    npy_intp dims[2];
    int nlev, nlat, nlon;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOO",
            kwlist,
            &pressure_obj,
            &temperature_obj,
            &dewpoint_obj)) {
        return NULL;
    }

    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );

    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL) {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }

    if (PyArray_NDIM(pressure_arr) != 3 || PyArray_NDIM(temperature_arr) != 3 ||
        PyArray_NDIM(dewpoint_arr) != 3) {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must be 3D arrays");
        goto fail;
    }

    nlev = (int) PyArray_DIM(pressure_arr, 0);
    nlat = (int) PyArray_DIM(pressure_arr, 1);
    nlon = (int) PyArray_DIM(pressure_arr, 2);

    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(temperature_arr, 1) != nlat ||
        (int) PyArray_DIM(temperature_arr, 2) != nlon ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 1) != nlat ||
        (int) PyArray_DIM(dewpoint_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a shape");
        goto fail;
    }

    dims[0] = nlat;
    dims[1] = nlon;
    out_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (out_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dcape_grid_c(
        nlev,
        nlat,
        nlon,
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temperature_arr),
        PyArray_DATA(dewpoint_arr),
        PyArray_DATA(out_arr)
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return (PyObject *) out_arr;

fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    Py_XDECREF(out_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {
        "dcape_profile",
        (PyCFunction) (void (*)(void)) py_dcape_profile,
        METH_VARARGS | METH_KEYWORDS,
        "Compute DCAPE for a single profile."
    },
    {
        "dcape_grid",
        (PyCFunction) (void (*)(void)) py_dcape_grid,
        METH_VARARGS | METH_KEYWORDS,
        "Compute DCAPE for a 3D grid of profiles."
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "dcape_core",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_dcape_core(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
