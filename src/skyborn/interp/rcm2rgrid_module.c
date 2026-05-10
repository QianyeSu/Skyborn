#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void drcm2rgrid(
    int ngrd,
    int nyi,
    int nxi,
    void *yi,
    void *xi,
    void *fi,
    int nyo,
    void *yo,
    int nxo,
    void *xo,
    void *fo,
    double xmsg,
    int ncrit,
    int opt,
    void *ier
);

void drgrid2rcm(
    int ngrd,
    int nyi,
    int nxi,
    void *yi,
    void *xi,
    void *fi,
    int nyo,
    int nxo,
    void *yo,
    void *xo,
    void *fo,
    double xmsg,
    int ncrit,
    int opt,
    void *ier
);

static PyArrayObject *to_double_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
}

static PyArrayObject *to_double_2d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
}

static PyObject *py_drcm2rgrid(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"lat2d", "lon2d", "fi", "lat1d", "lon1d", "xmsg", "ncrit", "opt", NULL};
    PyObject *lat2d_obj = NULL;
    PyObject *lon2d_obj = NULL;
    PyObject *fi_obj = NULL;
    PyObject *lat1d_obj = NULL;
    PyObject *lon1d_obj = NULL;
    PyArrayObject *lat2d_arr = NULL;
    PyArrayObject *lon2d_arr = NULL;
    PyArrayObject *fi_arr = NULL;
    PyArrayObject *lat1d_arr = NULL;
    PyArrayObject *lon1d_arr = NULL;
    PyArrayObject *fo_arr = NULL;
    npy_intp fo_dims[3];
    int ngrd = 0;
    int nyi = 0;
    int nxi = 0;
    int nyo = 0;
    int nxo = 0;
    int ncrit = 1;
    int opt = 1;
    int ier = 0;
    double xmsg = -99.0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|dii",
            kwlist,
            &lat2d_obj,
            &lon2d_obj,
            &fi_obj,
            &lat1d_obj,
            &lon1d_obj,
            &xmsg,
            &ncrit,
            &opt
        )) {
        return NULL;
    }

    lat2d_arr = to_double_2d_fortran(lat2d_obj);
    lon2d_arr = to_double_2d_fortran(lon2d_obj);
    fi_arr = (PyArrayObject *) PyArray_FROM_OTF(fi_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    lat1d_arr = to_double_1d(lat1d_obj);
    lon1d_arr = to_double_1d(lon1d_obj);
    if (lat2d_arr == NULL || lon2d_arr == NULL || fi_arr == NULL || lat1d_arr == NULL || lon1d_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(lat2d_arr) != 2 || PyArray_NDIM(lon2d_arr) != 2 || PyArray_NDIM(fi_arr) < 2 ||
        PyArray_NDIM(lat1d_arr) != 1 || PyArray_NDIM(lon1d_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "rcm2rgrid expects 2D lat/lon grids, 1D target coords, and fi with >=2 dims");
        goto fail;
    }

    nxi = (int) PyArray_DIM(lat2d_arr, 0);
    nyi = (int) PyArray_DIM(lat2d_arr, 1);
    if ((int) PyArray_DIM(lon2d_arr, 0) != nxi || (int) PyArray_DIM(lon2d_arr, 1) != nyi) {
        PyErr_SetString(PyExc_ValueError, "lat2d and lon2d must have the same shape");
        goto fail;
    }
    if ((int) PyArray_DIM(fi_arr, 0) != nxi ||
        (int) PyArray_DIM(fi_arr, 1) != nyi) {
        PyErr_SetString(PyExc_ValueError, "fi leading dimensions must match lat2d/lon2d after transpose");
        goto fail;
    }

    nyo = (int) PyArray_DIM(lat1d_arr, 0);
    nxo = (int) PyArray_DIM(lon1d_arr, 0);
    if (nyo <= 0 || nxo <= 0) {
        PyErr_SetString(PyExc_ValueError, "target coordinate arrays must be non-empty");
        goto fail;
    }

    fo_dims[0] = nxo;
    fo_dims[1] = nyo;
    fo_dims[2] = PyArray_NDIM(fi_arr) > 2 ? (int) PyArray_DIM(fi_arr, 2) : 1;
    fo_arr = (PyArrayObject *) PyArray_EMPTY(3, fo_dims, NPY_FLOAT64, 1);
    if (fo_arr == NULL) {
        goto fail;
    }

    ngrd = PyArray_NDIM(fi_arr) > 2 ? (int) PyArray_DIM(fi_arr, 2) : 1;

    Py_BEGIN_ALLOW_THREADS
    drcm2rgrid(
        ngrd,
        nyi,
        nxi,
        PyArray_DATA(lat2d_arr),
        PyArray_DATA(lon2d_arr),
        PyArray_DATA(fi_arr),
        nyo,
        PyArray_DATA(lat1d_arr),
        nxo,
        PyArray_DATA(lon1d_arr),
        PyArray_DATA(fo_arr),
        xmsg,
        ncrit,
        opt,
        &ier
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(lat2d_arr);
    Py_DECREF(lon2d_arr);
    Py_DECREF(fi_arr);
    Py_DECREF(lat1d_arr);
    Py_DECREF(lon1d_arr);

    if (ier != 0) {
        Py_DECREF(fo_arr);
        PyErr_Format(PyExc_RuntimeError, "drcm2rgrid backend returned ier=%d", ier);
        return NULL;
    }

    return PyArray_Return(fo_arr);

fail:
    Py_XDECREF(lat2d_arr);
    Py_XDECREF(lon2d_arr);
    Py_XDECREF(fi_arr);
    Py_XDECREF(lat1d_arr);
    Py_XDECREF(lon1d_arr);
    Py_XDECREF(fo_arr);
    return NULL;
}

static PyObject *py_drgrid2rcm(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"lat1d", "lon1d", "fi", "lat2d", "lon2d", "xmsg", "ncrit", "opt", NULL};
    PyObject *lat1d_obj = NULL;
    PyObject *lon1d_obj = NULL;
    PyObject *fi_obj = NULL;
    PyObject *lat2d_obj = NULL;
    PyObject *lon2d_obj = NULL;
    PyArrayObject *lat1d_arr = NULL;
    PyArrayObject *lon1d_arr = NULL;
    PyArrayObject *fi_arr = NULL;
    PyArrayObject *lat2d_arr = NULL;
    PyArrayObject *lon2d_arr = NULL;
    PyArrayObject *fo_arr = NULL;
    npy_intp fo_dims[3];
    int ngrd = 0;
    int nyi = 0;
    int nxi = 0;
    int nyo = 0;
    int nxo = 0;
    int ncrit = 1;
    int opt = 1;
    int ier = 0;
    double xmsg = -99.0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|dii",
            kwlist,
            &lat1d_obj,
            &lon1d_obj,
            &fi_obj,
            &lat2d_obj,
            &lon2d_obj,
            &xmsg,
            &ncrit,
            &opt
        )) {
        return NULL;
    }

    lat1d_arr = to_double_1d(lat1d_obj);
    lon1d_arr = to_double_1d(lon1d_obj);
    fi_arr = (PyArrayObject *) PyArray_FROM_OTF(fi_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    lat2d_arr = to_double_2d_fortran(lat2d_obj);
    lon2d_arr = to_double_2d_fortran(lon2d_obj);
    if (lat1d_arr == NULL || lon1d_arr == NULL || fi_arr == NULL || lat2d_arr == NULL || lon2d_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(lat1d_arr) != 1 || PyArray_NDIM(lon1d_arr) != 1 || PyArray_NDIM(fi_arr) < 2 ||
        PyArray_NDIM(lat2d_arr) != 2 || PyArray_NDIM(lon2d_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "drgrid2rcm expects 1D source coords, 2D target coords, and fi with >=2 dims");
        goto fail;
    }

    nyi = (int) PyArray_DIM(lat1d_arr, 0);
    nxi = (int) PyArray_DIM(lon1d_arr, 0);
    nxo = (int) PyArray_DIM(lat2d_arr, 0);
    nyo = (int) PyArray_DIM(lat2d_arr, 1);
    if ((int) PyArray_DIM(lon2d_arr, 0) != nxo || (int) PyArray_DIM(lon2d_arr, 1) != nyo) {
        PyErr_SetString(PyExc_ValueError, "lat2d and lon2d must have the same shape");
        goto fail;
    }
    if ((int) PyArray_DIM(fi_arr, 0) != nxi ||
        (int) PyArray_DIM(fi_arr, 1) != nyi) {
        PyErr_SetString(PyExc_ValueError, "fi leading dimensions must match lon1d/lat1d after transpose");
        goto fail;
    }

    fo_dims[0] = nxo;
    fo_dims[1] = nyo;
    fo_dims[2] = PyArray_NDIM(fi_arr) > 2 ? (int) PyArray_DIM(fi_arr, 2) : 1;
    fo_arr = (PyArrayObject *) PyArray_EMPTY(3, fo_dims, NPY_FLOAT64, 1);
    if (fo_arr == NULL) {
        goto fail;
    }

    ngrd = PyArray_NDIM(fi_arr) > 2 ? (int) PyArray_DIM(fi_arr, 2) : 1;

    Py_BEGIN_ALLOW_THREADS
    drgrid2rcm(
        ngrd,
        nyi,
        nxi,
        PyArray_DATA(lat1d_arr),
        PyArray_DATA(lon1d_arr),
        PyArray_DATA(fi_arr),
        nyo,
        nxo,
        PyArray_DATA(lat2d_arr),
        PyArray_DATA(lon2d_arr),
        PyArray_DATA(fo_arr),
        xmsg,
        ncrit,
        opt,
        &ier
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(lat1d_arr);
    Py_DECREF(lon1d_arr);
    Py_DECREF(fi_arr);
    Py_DECREF(lat2d_arr);
    Py_DECREF(lon2d_arr);

    if (ier != 0) {
        Py_DECREF(fo_arr);
        PyErr_Format(PyExc_RuntimeError, "drgrid2rcm backend returned ier=%d", ier);
        return NULL;
    }

    return PyArray_Return(fo_arr);

fail:
    Py_XDECREF(lat1d_arr);
    Py_XDECREF(lon1d_arr);
    Py_XDECREF(fi_arr);
    Py_XDECREF(lat2d_arr);
    Py_XDECREF(lon2d_arr);
    Py_XDECREF(fo_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"drcm2rgrid", (PyCFunction) (void (*)(void)) py_drcm2rgrid, METH_VARARGS | METH_KEYWORDS, "Interpolate curvilinear data to a rectilinear grid."},
    {"drgrid2rcm", (PyCFunction) (void (*)(void)) py_drgrid2rcm, METH_VARARGS | METH_KEYWORDS, "Interpolate rectilinear data to a curvilinear grid."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "rcm2rgrid",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_rcm2rgrid(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
