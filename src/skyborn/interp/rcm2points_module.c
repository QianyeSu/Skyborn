#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void drcm2points(
    int ngrd,
    int nyi,
    int nxi,
    void *yi,
    void *xi,
    void *fi,
    int nxyo,
    void *yo,
    void *xo,
    void *fo,
    double xmsg,
    int opt,
    int ncrit,
    int kval,
    void *ier
);

static PyArrayObject *to_double_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
}

static PyArrayObject *to_double_2d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
}

static PyObject *py_drcm2points(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"lat2d", "lon2d", "fi", "lat1d", "lon1d", "xmsg", "opt", NULL};
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
    npy_intp fo_dims[2];
    int ngrd = 0;
    int nyi = 0;
    int nxi = 0;
    int nxyo = 0;
    int opt = 0;
    int ncrit = 1;
    int kval = 1;
    int ier = 0;
    double xmsg = -99.0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|di",
            kwlist,
            &lat2d_obj,
            &lon2d_obj,
            &fi_obj,
            &lat1d_obj,
            &lon1d_obj,
            &xmsg,
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
        PyErr_SetString(PyExc_ValueError, "rcm2points expects 2D lat/lon grids, 1D target coords, and fi with >=2 dims");
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

    nxyo = (int) PyArray_DIM(lat1d_arr, 0);
    if ((int) PyArray_DIM(lon1d_arr, 0) != nxyo) {
        PyErr_SetString(PyExc_ValueError, "lat1d and lon1d must have the same length");
        goto fail;
    }
    ngrd = PyArray_NDIM(fi_arr) > 2 ? (int) PyArray_DIM(fi_arr, 2) : 1;

    fo_dims[0] = nxyo;
    fo_dims[1] = ngrd;
    fo_arr = (PyArrayObject *) PyArray_EMPTY(2, fo_dims, NPY_FLOAT64, 1);
    if (fo_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    drcm2points(
        ngrd,
        nyi,
        nxi,
        PyArray_DATA(lat2d_arr),
        PyArray_DATA(lon2d_arr),
        PyArray_DATA(fi_arr),
        nxyo,
        PyArray_DATA(lat1d_arr),
        PyArray_DATA(lon1d_arr),
        PyArray_DATA(fo_arr),
        xmsg,
        opt,
        ncrit,
        kval,
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
        PyErr_Format(PyExc_RuntimeError, "drcm2points backend returned ier=%d", ier);
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

static PyMethodDef module_methods[] = {
    {"drcm2points", (PyCFunction) (void (*)(void)) py_drcm2points, METH_VARARGS | METH_KEYWORDS, "Interpolate curvilinear-grid data to point targets."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "rcm2points",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_rcm2points(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
