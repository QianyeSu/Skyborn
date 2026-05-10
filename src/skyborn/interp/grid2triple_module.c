#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void grid2triple(
    void *x,
    void *y,
    void *z,
    void *d,
    void *ld,
    double zmsg,
    void *ier,
    int mx,
    int ny,
    int ldmax
);

static PyArrayObject *to_double_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
}

static PyArrayObject *to_double_2d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
}

static PyObject *py_grid2triple(PyObject *self, PyObject *args) {
    PyObject *x_obj = NULL;
    PyObject *y_obj = NULL;
    PyObject *z_obj = NULL;
    PyArrayObject *x_arr = NULL;
    PyArrayObject *y_arr = NULL;
    PyArrayObject *z_arr = NULL;
    PyArrayObject *d_arr = NULL;
    npy_intp d_dims[2];
    int mx = 0;
    int ny = 0;
    int ldmax = 0;
    int ld = 0;
    int ier = 0;
    double zmsg = 0.0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOd", &x_obj, &y_obj, &z_obj, &zmsg)) {
        return NULL;
    }

    x_arr = to_double_1d(x_obj);
    y_arr = to_double_1d(y_obj);
    z_arr = to_double_2d_fortran(z_obj);
    if (x_arr == NULL || y_arr == NULL || z_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(x_arr) != 1 || PyArray_NDIM(y_arr) != 1 || PyArray_NDIM(z_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "grid2triple expects x(1D), y(1D), z(2D)");
        goto fail;
    }

    mx = (int) PyArray_DIM(x_arr, 0);
    ny = (int) PyArray_DIM(y_arr, 0);
    if ((int) PyArray_DIM(z_arr, 0) != mx || (int) PyArray_DIM(z_arr, 1) != ny) {
        PyErr_SetString(PyExc_ValueError, "z must have shape (mx, ny)");
        goto fail;
    }

    ldmax = mx * ny;
    d_dims[0] = ldmax;
    d_dims[1] = 3;
    d_arr = (PyArrayObject *) PyArray_EMPTY(2, d_dims, NPY_FLOAT64, 1);
    if (d_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    grid2triple(
        PyArray_DATA(x_arr),
        PyArray_DATA(y_arr),
        PyArray_DATA(z_arr),
        PyArray_DATA(d_arr),
        &ld,
        zmsg,
        &ier,
        mx,
        ny,
        ldmax
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(x_arr);
    Py_DECREF(y_arr);
    Py_DECREF(z_arr);
    if (ier != 0 && ier != -10) {
        Py_DECREF(d_arr);
        PyErr_Format(PyExc_RuntimeError, "grid2triple backend returned ier=%d", ier);
        return NULL;
    }
    {
        PyObject *result = Py_BuildValue("Oi", d_arr, ld);
        Py_DECREF(d_arr);
        return result;
    }

fail:
    Py_XDECREF(x_arr);
    Py_XDECREF(y_arr);
    Py_XDECREF(z_arr);
    Py_XDECREF(d_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"grid2triple", py_grid2triple, METH_VARARGS, "Pack a 2D grid into a legacy triple array and return (d, ld)."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "grid2triple",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_grid2triple(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
