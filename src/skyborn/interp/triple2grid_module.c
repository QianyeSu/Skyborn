#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void triple2grid1(
    int kz,
    void *xi,
    void *yi,
    void *zi,
    double zmsg,
    int mx,
    int ny,
    void *gx,
    void *gy,
    void *grid,
    double domain,
    int loop,
    int method,
    double distmx,
    int mx2,
    int ny2,
    void *x,
    void *y,
    void *z,
    void *gbigx,
    void *gbigy,
    void *gbigxy,
    void *ier
);

static PyArrayObject *to_double_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
}

static PyObject *py_triple2grid1(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"xi", "yi", "zi", "gx", "gy", "zmsg", "domain", "method", "distmx", NULL};
    PyObject *xi_obj = NULL;
    PyObject *yi_obj = NULL;
    PyObject *zi_obj = NULL;
    PyObject *gx_obj = NULL;
    PyObject *gy_obj = NULL;
    PyArrayObject *xi_arr = NULL;
    PyArrayObject *yi_arr = NULL;
    PyArrayObject *zi_arr = NULL;
    PyArrayObject *gx_arr = NULL;
    PyArrayObject *gy_arr = NULL;
    PyArrayObject *grid_arr = NULL;
    PyArrayObject *x_work = NULL;
    PyArrayObject *y_work = NULL;
    PyArrayObject *z_work = NULL;
    PyArrayObject *gbigx_arr = NULL;
    PyArrayObject *gbigy_arr = NULL;
    PyArrayObject *gbigxy_arr = NULL;
    npy_intp grid_dims[2];
    npy_intp bigx_dims[1];
    npy_intp bigy_dims[1];
    npy_intp bigxy_dims[2];
    int kz = 0;
    int mx = 0;
    int ny = 0;
    int mx2 = 0;
    int ny2 = 0;
    int method = 1;
    int loop = 0;
    int ier = 0;
    double zmsg = -99.0;
    double domain = 1.0;
    double distmx = 1.0e20;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|ddid",
            kwlist,
            &xi_obj,
            &yi_obj,
            &zi_obj,
            &gx_obj,
            &gy_obj,
            &zmsg,
            &domain,
            &method,
            &distmx
        )) {
        return NULL;
    }

    xi_arr = to_double_1d(xi_obj);
    yi_arr = to_double_1d(yi_obj);
    zi_arr = to_double_1d(zi_obj);
    gx_arr = to_double_1d(gx_obj);
    gy_arr = to_double_1d(gy_obj);
    if (xi_arr == NULL || yi_arr == NULL || zi_arr == NULL || gx_arr == NULL || gy_arr == NULL) {
        goto fail;
    }

    if (PyArray_NDIM(xi_arr) != 1 || PyArray_NDIM(yi_arr) != 1 || PyArray_NDIM(zi_arr) != 1 ||
        PyArray_NDIM(gx_arr) != 1 || PyArray_NDIM(gy_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "triple2grid1 expects 1D xi, yi, zi, gx, and gy arrays");
        goto fail;
    }

    kz = (int) PyArray_DIM(xi_arr, 0);
    if ((int) PyArray_DIM(yi_arr, 0) != kz || (int) PyArray_DIM(zi_arr, 0) != kz) {
        PyErr_SetString(PyExc_ValueError, "xi, yi, and zi must have the same length");
        goto fail;
    }

    mx = (int) PyArray_DIM(gx_arr, 0);
    ny = (int) PyArray_DIM(gy_arr, 0);
    mx2 = mx + 2;
    ny2 = ny + 2;

    grid_dims[0] = mx;
    grid_dims[1] = ny;
    bigx_dims[0] = mx2;
    bigy_dims[0] = ny2;
    bigxy_dims[0] = mx2;
    bigxy_dims[1] = ny2;

    grid_arr = (PyArrayObject *) PyArray_EMPTY(2, grid_dims, NPY_FLOAT64, 1);
    x_work = (PyArrayObject *) PyArray_EMPTY(1, (npy_intp[]) {kz}, NPY_FLOAT64, 1);
    y_work = (PyArrayObject *) PyArray_EMPTY(1, (npy_intp[]) {kz}, NPY_FLOAT64, 1);
    z_work = (PyArrayObject *) PyArray_EMPTY(1, (npy_intp[]) {kz}, NPY_FLOAT64, 1);
    gbigx_arr = (PyArrayObject *) PyArray_EMPTY(1, bigx_dims, NPY_FLOAT64, 1);
    gbigy_arr = (PyArrayObject *) PyArray_EMPTY(1, bigy_dims, NPY_FLOAT64, 1);
    gbigxy_arr = (PyArrayObject *) PyArray_EMPTY(2, bigxy_dims, NPY_FLOAT64, 1);
    if (grid_arr == NULL || x_work == NULL || y_work == NULL || z_work == NULL ||
        gbigx_arr == NULL || gbigy_arr == NULL || gbigxy_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    triple2grid1(
        kz,
        PyArray_DATA(xi_arr),
        PyArray_DATA(yi_arr),
        PyArray_DATA(zi_arr),
        zmsg,
        mx,
        ny,
        PyArray_DATA(gx_arr),
        PyArray_DATA(gy_arr),
        PyArray_DATA(grid_arr),
        domain,
        loop,
        method,
        distmx,
        mx2,
        ny2,
        PyArray_DATA(x_work),
        PyArray_DATA(y_work),
        PyArray_DATA(z_work),
        PyArray_DATA(gbigx_arr),
        PyArray_DATA(gbigy_arr),
        PyArray_DATA(gbigxy_arr),
        &ier
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(xi_arr);
    Py_DECREF(yi_arr);
    Py_DECREF(zi_arr);
    Py_DECREF(gx_arr);
    Py_DECREF(gy_arr);
    Py_DECREF(x_work);
    Py_DECREF(y_work);
    Py_DECREF(z_work);
    Py_DECREF(gbigx_arr);
    Py_DECREF(gbigy_arr);
    Py_DECREF(gbigxy_arr);

    if (ier != 0 && ier != -1 && ier != -11) {
        Py_DECREF(grid_arr);
        PyErr_Format(PyExc_RuntimeError, "triple2grid1 backend returned ier=%d", ier);
        return NULL;
    }
    return (PyObject *) grid_arr;

fail:
    Py_XDECREF(xi_arr);
    Py_XDECREF(yi_arr);
    Py_XDECREF(zi_arr);
    Py_XDECREF(gx_arr);
    Py_XDECREF(gy_arr);
    Py_XDECREF(grid_arr);
    Py_XDECREF(x_work);
    Py_XDECREF(y_work);
    Py_XDECREF(z_work);
    Py_XDECREF(gbigx_arr);
    Py_XDECREF(gbigy_arr);
    Py_XDECREF(gbigxy_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"triple2grid1", (PyCFunction) py_triple2grid1, METH_VARARGS | METH_KEYWORDS, "Map a triple list onto a target rectilinear grid."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "triple2grid",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_triple2grid(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
