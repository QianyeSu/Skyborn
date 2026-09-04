#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void liang_single_c(
    void *y1, void *y2, double *t21, double *tau21, double *z,
    double *dh1_star, double *dh1_noise, int nm, int npt, int *ierr
);
void liang_batch_c(
    void *y1, void *y2, void *t21, void *tau21,
    int nm, int nsim, int npt, int *ierr
);
void ar1_filter_batch_c(
    void *innovations, void *output, double g, int burnin,
    int nnoise, int nsim, int nout, int *ierr
);

static void set_liang_error(const char *context, int ierr)
{
    if (ierr == 1) {
        PyErr_SetString(
            PyExc_ValueError,
            "npt must leave at least two differenced samples"
        );
    } else if (ierr == 2) {
        PyErr_SetString(
            PyExc_ValueError,
            "Liang inputs must contain only finite values"
        );
    } else if (ierr == 3) {
        PyErr_SetString(
            PyExc_ValueError,
            "Liang information flow is undefined for singular covariance"
        );
    } else if (ierr == 4) {
        PyErr_SetString(
            PyExc_ValueError,
            "Liang information flow normalization is undefined"
        );
    } else if (ierr == 5) {
        PyErr_SetString(
            PyExc_ValueError,
            "AR(1) filter dimensions are invalid"
        );
    } else {
        PyErr_Format(
            PyExc_ValueError,
            "%s failed with error code %d",
            context,
            ierr
        );
    }
}

static PyArrayObject *as_double_1d(PyObject *obj)
{
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
}

static PyArrayObject *as_double_2d_fortran(PyObject *obj)
{
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
}

static PyObject *py_liang_single(PyObject *self, PyObject *args)
{
    PyObject *y1_obj = NULL;
    PyObject *y2_obj = NULL;
    PyArrayObject *y1 = NULL;
    PyArrayObject *y2 = NULL;
    int npt = 0;
    int nm = 0;
    int ierr = 0;
    double t21, tau21, z, dh1_star, dh1_noise;

    (void) self;
    if (!PyArg_ParseTuple(args, "OOi", &y1_obj, &y2_obj, &npt)) {
        return NULL;
    }

    y1 = as_double_1d(y1_obj);
    y2 = as_double_1d(y2_obj);
    if (y1 == NULL || y2 == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(y1) != 1 || PyArray_NDIM(y2) != 1) {
        PyErr_SetString(PyExc_ValueError, "Liang inputs must be 1D arrays");
        goto fail;
    }
    nm = (int) PyArray_DIM(y1, 0);
    if ((int) PyArray_DIM(y2, 0) != nm) {
        PyErr_SetString(PyExc_ValueError, "Liang inputs must share a length");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    liang_single_c(
        PyArray_DATA(y1), PyArray_DATA(y2), &t21, &tau21, &z,
        &dh1_star, &dh1_noise, nm, npt, &ierr
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(y1);
    Py_DECREF(y2);
    if (ierr != 0) {
        set_liang_error("Liang information-flow kernel", ierr);
        return NULL;
    }
    return Py_BuildValue("ddddd", t21, tau21, z, dh1_star, dh1_noise);

fail:
    Py_XDECREF(y1);
    Py_XDECREF(y2);
    return NULL;
}

static PyObject *py_liang_batch(PyObject *self, PyObject *args)
{
    PyObject *y1_obj = NULL;
    PyObject *y2_obj = NULL;
    PyArrayObject *y1 = NULL;
    PyArrayObject *y2 = NULL;
    PyArrayObject *t21 = NULL;
    PyArrayObject *tau21 = NULL;
    PyObject *result = NULL;
    npy_intp shape[1];
    int npt = 0;
    int nm = 0;
    int nsim = 0;
    int ierr = 0;

    (void) self;
    if (!PyArg_ParseTuple(args, "OOi", &y1_obj, &y2_obj, &npt)) {
        return NULL;
    }

    y1 = as_double_2d_fortran(y1_obj);
    y2 = as_double_2d_fortran(y2_obj);
    if (y1 == NULL || y2 == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(y1) != 2 || PyArray_NDIM(y2) != 2) {
        PyErr_SetString(PyExc_ValueError, "Liang batch inputs must be 2D arrays");
        goto fail;
    }
    if (PyArray_DIM(y1, 0) != PyArray_DIM(y2, 0) ||
        PyArray_DIM(y1, 1) != PyArray_DIM(y2, 1)) {
        PyErr_SetString(PyExc_ValueError, "Liang batch inputs must share a shape");
        goto fail;
    }

    nm = (int) PyArray_DIM(y1, 0);
    nsim = (int) PyArray_DIM(y1, 1);
    shape[0] = nsim;
    t21 = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 0);
    tau21 = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 0);
    if (t21 == NULL || tau21 == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    liang_batch_c(
        PyArray_DATA(y1), PyArray_DATA(y2),
        PyArray_DATA(t21), PyArray_DATA(tau21),
        nm, nsim, npt, &ierr
    );
    Py_END_ALLOW_THREADS

    if (ierr != 0) {
        set_liang_error("Liang batch kernel", ierr);
        goto fail;
    }

    result = Py_BuildValue("NN", (PyObject *) t21, (PyObject *) tau21);
    t21 = NULL;
    tau21 = NULL;

fail:
    Py_XDECREF(y1);
    Py_XDECREF(y2);
    Py_XDECREF(t21);
    Py_XDECREF(tau21);
    return result;
}

static PyObject *py_ar1_filter_batch(PyObject *self, PyObject *args)
{
    PyObject *innovations_obj = NULL;
    PyArrayObject *innovations = NULL;
    PyArrayObject *output = NULL;
    PyObject *result = NULL;
    npy_intp shape[2];
    double g = 0.0;
    int nout = 0;
    int burnin = 50;
    int nnoise = 0;
    int nsim = 0;
    int ierr = 0;

    (void) self;
    if (!PyArg_ParseTuple(args, "Odi", &innovations_obj, &g, &nout)) {
        return NULL;
    }

    innovations = as_double_2d_fortran(innovations_obj);
    if (innovations == NULL) {
        return NULL;
    }
    if (PyArray_NDIM(innovations) != 2) {
        PyErr_SetString(
            PyExc_ValueError,
            "AR(1) filter innovations must be a 2D array"
        );
        goto fail;
    }

    nnoise = (int) PyArray_DIM(innovations, 0);
    nsim = (int) PyArray_DIM(innovations, 1);
    if (nout < 1 || nnoise < burnin + nout || nsim < 1) {
        PyErr_SetString(
            PyExc_ValueError,
            "AR(1) filter dimensions are invalid"
        );
        goto fail;
    }

    shape[0] = nout;
    shape[1] = nsim;
    output = (PyArrayObject *) PyArray_EMPTY(2, shape, NPY_FLOAT64, 1);
    if (output == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    ar1_filter_batch_c(
        PyArray_DATA(innovations), PyArray_DATA(output), g, burnin,
        nnoise, nsim, nout, &ierr
    );
    Py_END_ALLOW_THREADS

    if (ierr != 0) {
        set_liang_error("AR(1) filter kernel", ierr + 4);
        goto fail;
    }

    result = (PyObject *) output;
    output = NULL;

fail:
    Py_XDECREF(innovations);
    Py_XDECREF(output);
    return result;
}

static PyMethodDef module_methods[] = {
    {"liang_single", py_liang_single, METH_VARARGS, "Calculate one Liang flow."},
    {"liang_batch", py_liang_batch, METH_VARARGS, "Calculate Liang flows by column."},
    {"ar1_filter_batch", py_ar1_filter_batch, METH_VARARGS, "Filter AR(1) innovations by column."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_causality_core",
    "Compiled Liang information-flow kernels.",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__causality_core(void)
{
    PyObject *module = PyModule_Create(&module_def);
    if (module == NULL) {
        return NULL;
    }
    import_array();
    return module;
}
