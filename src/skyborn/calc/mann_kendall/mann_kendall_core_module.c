#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void mk_score_var_batch(
    void *data,
    void *s_values,
    void *var_values,
    int modified,
    int ntime,
    int nseries
);
void sen_slope_batch(
    void *data,
    void *slopes,
    int ntime,
    int nseries
);
void grouped_sen_slope_batch(
    void *data,
    void *slopes,
    int period,
    int ntime,
    int nseries
);
void grouped_correlated_stats_batch(
    void *data,
    void *s_values,
    void *var_values,
    void *denom,
    int period,
    int ntime,
    int nseries
);
void partial_stats_batch(
    void *response,
    void *covariate,
    void *s_values,
    void *var_values,
    void *tau_values,
    int ntime,
    int nseries
);
void partial_stats_sen_batch(
    void *response,
    void *covariate,
    void *s_values,
    void *var_values,
    void *tau_values,
    void *slopes,
    int ntime,
    int nseries
);
void mk_score_var_sen_batch(
    void *data,
    void *s_values,
    void *var_values,
    void *slopes,
    int modified,
    int ntime,
    int nseries
);
void mk_yue_wang_score_var_sen_batch(
    void *data,
    void *s_values,
    void *var_values,
    void *slopes,
    int lag,
    int ntime,
    int nseries
);
void mk_hamed_rao_score_var_sen_batch(
    void *data,
    void *s_values,
    void *var_values,
    void *slopes,
    double interval,
    int lag,
    int ntime,
    int nseries
);

static PyArrayObject *to_double_2d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
}

static int require_2d(PyArrayObject *arr, const char *name) {
    if (PyArray_NDIM(arr) != 2) {
        PyErr_Format(PyExc_ValueError, "%s must be a 2D array", name);
        return -1;
    }
    return 0;
}

static int same_shape(PyArrayObject *a, PyArrayObject *b, const char *a_name, const char *b_name) {
    if (PyArray_DIM(a, 0) != PyArray_DIM(b, 0) || PyArray_DIM(a, 1) != PyArray_DIM(b, 1)) {
        PyErr_Format(PyExc_ValueError, "%s and %s must share the same 2D shape", a_name, b_name);
        return -1;
    }
    return 0;
}

static PyObject *mk_score_var_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    int modified = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Oi", &data_obj, &modified)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL) {
        goto fail;
    }

    mk_score_var_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        modified,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    {
        PyObject *result = PyTuple_Pack(2, (PyObject *) s_arr, (PyObject *) var_arr);
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        return result;
    }

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    return NULL;
}

static PyObject *sen_slope_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "O", &data_obj)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (slopes_arr == NULL) {
        goto fail;
    }

    sen_slope_batch(PyArray_DATA(data_arr), PyArray_DATA(slopes_arr), ntime, nseries);

    Py_DECREF(data_arr);
    return (PyObject *) slopes_arr;

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyObject *grouped_sen_slope_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    int period = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Oi", &data_obj, &period)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (slopes_arr == NULL) {
        goto fail;
    }

    grouped_sen_slope_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(slopes_arr),
        period,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    return (PyObject *) slopes_arr;

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyObject *grouped_correlated_stats_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    double denom = 0.0;
    int period = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Oi", &data_obj, &period)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL) {
        goto fail;
    }

    grouped_correlated_stats_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        &denom,
        period,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    {
        PyObject *result = Py_BuildValue("OOd", s_arr, var_arr, denom);
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        return result;
    }

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    return NULL;
}

static PyObject *partial_stats_batch_py(PyObject *self, PyObject *args) {
    PyObject *response_obj = NULL;
    PyObject *covariate_obj = NULL;
    PyArrayObject *response_arr = NULL;
    PyArrayObject *covariate_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    PyArrayObject *tau_arr = NULL;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "OO", &response_obj, &covariate_obj)) {
        return NULL;
    }

    response_arr = to_double_2d_fortran(response_obj);
    covariate_arr = to_double_2d_fortran(covariate_obj);
    if (response_arr == NULL || covariate_arr == NULL) {
        goto fail;
    }
    if (require_2d(response_arr, "response") != 0 ||
        require_2d(covariate_arr, "covariate") != 0 ||
        same_shape(response_arr, covariate_arr, "response", "covariate") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(response_arr, 0);
    nseries = (int) PyArray_DIM(response_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    tau_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL || tau_arr == NULL) {
        goto fail;
    }

    partial_stats_batch(
        PyArray_DATA(response_arr),
        PyArray_DATA(covariate_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        PyArray_DATA(tau_arr),
        ntime,
        nseries
    );

    Py_DECREF(response_arr);
    Py_DECREF(covariate_arr);
    {
        PyObject *result = PyTuple_Pack(
            3,
            (PyObject *) s_arr,
            (PyObject *) var_arr,
            (PyObject *) tau_arr
        );
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        Py_DECREF(tau_arr);
        return result;
    }

fail:
    Py_XDECREF(response_arr);
    Py_XDECREF(covariate_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    Py_XDECREF(tau_arr);
    return NULL;
}

static PyObject *partial_stats_sen_batch_py(PyObject *self, PyObject *args) {
    PyObject *response_obj = NULL;
    PyObject *covariate_obj = NULL;
    PyArrayObject *response_arr = NULL;
    PyArrayObject *covariate_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    PyArrayObject *tau_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "OO", &response_obj, &covariate_obj)) {
        return NULL;
    }

    response_arr = to_double_2d_fortran(response_obj);
    covariate_arr = to_double_2d_fortran(covariate_obj);
    if (response_arr == NULL || covariate_arr == NULL) {
        goto fail;
    }
    if (require_2d(response_arr, "response") != 0 ||
        require_2d(covariate_arr, "covariate") != 0 ||
        same_shape(response_arr, covariate_arr, "response", "covariate") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(response_arr, 0);
    nseries = (int) PyArray_DIM(response_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    tau_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL || tau_arr == NULL || slopes_arr == NULL) {
        goto fail;
    }

    partial_stats_sen_batch(
        PyArray_DATA(response_arr),
        PyArray_DATA(covariate_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        PyArray_DATA(tau_arr),
        PyArray_DATA(slopes_arr),
        ntime,
        nseries
    );

    Py_DECREF(response_arr);
    Py_DECREF(covariate_arr);
    {
        PyObject *result = PyTuple_Pack(
            4,
            (PyObject *) s_arr,
            (PyObject *) var_arr,
            (PyObject *) tau_arr,
            (PyObject *) slopes_arr
        );
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        Py_DECREF(tau_arr);
        Py_DECREF(slopes_arr);
        return result;
    }

fail:
    Py_XDECREF(response_arr);
    Py_XDECREF(covariate_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    Py_XDECREF(tau_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyObject *mk_score_var_sen_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    int modified = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Oi", &data_obj, &modified)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL || slopes_arr == NULL) {
        goto fail;
    }

    mk_score_var_sen_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        PyArray_DATA(slopes_arr),
        modified,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    {
        PyObject *result = PyTuple_Pack(
            3,
            (PyObject *) s_arr,
            (PyObject *) var_arr,
            (PyObject *) slopes_arr
        );
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        Py_DECREF(slopes_arr);
        return result;
    }

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyObject *mk_yue_wang_score_var_sen_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    int lag = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Oi", &data_obj, &lag)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL || slopes_arr == NULL) {
        goto fail;
    }

    mk_yue_wang_score_var_sen_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        PyArray_DATA(slopes_arr),
        lag,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    {
        PyObject *result = PyTuple_Pack(
            3,
            (PyObject *) s_arr,
            (PyObject *) var_arr,
            (PyObject *) slopes_arr
        );
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        Py_DECREF(slopes_arr);
        return result;
    }

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyObject *mk_hamed_rao_score_var_sen_batch_py(PyObject *self, PyObject *args) {
    PyObject *data_obj = NULL;
    PyArrayObject *data_arr = NULL;
    PyArrayObject *s_arr = NULL;
    PyArrayObject *var_arr = NULL;
    PyArrayObject *slopes_arr = NULL;
    double interval = 0.0;
    int lag = 0;
    int ntime = 0;
    int nseries = 0;
    npy_intp shape[1];

    (void) self;

    if (!PyArg_ParseTuple(args, "Odi", &data_obj, &interval, &lag)) {
        return NULL;
    }

    data_arr = to_double_2d_fortran(data_obj);
    if (data_arr == NULL) {
        return NULL;
    }
    if (require_2d(data_arr, "data") != 0) {
        goto fail;
    }

    ntime = (int) PyArray_DIM(data_arr, 0);
    nseries = (int) PyArray_DIM(data_arr, 1);
    shape[0] = nseries;
    s_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    var_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    slopes_arr = (PyArrayObject *) PyArray_EMPTY(1, shape, NPY_FLOAT64, 1);
    if (s_arr == NULL || var_arr == NULL || slopes_arr == NULL) {
        goto fail;
    }

    mk_hamed_rao_score_var_sen_batch(
        PyArray_DATA(data_arr),
        PyArray_DATA(s_arr),
        PyArray_DATA(var_arr),
        PyArray_DATA(slopes_arr),
        interval,
        lag,
        ntime,
        nseries
    );

    Py_DECREF(data_arr);
    {
        PyObject *result = PyTuple_Pack(
            3,
            (PyObject *) s_arr,
            (PyObject *) var_arr,
            (PyObject *) slopes_arr
        );
        Py_DECREF(s_arr);
        Py_DECREF(var_arr);
        Py_DECREF(slopes_arr);
        return result;
    }

fail:
    Py_XDECREF(data_arr);
    Py_XDECREF(s_arr);
    Py_XDECREF(var_arr);
    Py_XDECREF(slopes_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"mk_score_var_batch", mk_score_var_batch_py, METH_VARARGS, "Compute MK score and variance for a clean 2D batch."},
    {"sen_slope_batch", sen_slope_batch_py, METH_VARARGS, "Compute Sen slopes for a clean 2D batch."},
    {"grouped_sen_slope_batch", grouped_sen_slope_batch_py, METH_VARARGS, "Compute grouped Sen slopes for a clean 2D batch."},
    {"grouped_correlated_stats_batch", grouped_correlated_stats_batch_py, METH_VARARGS, "Compute correlated grouped MK statistics for a clean 2D batch."},
    {"partial_stats_batch", partial_stats_batch_py, METH_VARARGS, "Compute partial MK statistics for paired clean 2D batches."},
    {"partial_stats_sen_batch", partial_stats_sen_batch_py, METH_VARARGS, "Compute partial MK statistics and slopes for paired clean 2D batches."},
    {"mk_score_var_sen_batch", mk_score_var_sen_batch_py, METH_VARARGS, "Compute MK score, variance, and Sen slope for a clean 2D batch."},
    {"mk_yue_wang_score_var_sen_batch", mk_yue_wang_score_var_sen_batch_py, METH_VARARGS, "Compute Yue-Wang MK score, variance, and slope for a clean 2D batch."},
    {"mk_hamed_rao_score_var_sen_batch", mk_hamed_rao_score_var_sen_batch_py, METH_VARARGS, "Compute Hamed-Rao MK score, variance, and slope for a clean 2D batch."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mann_kendall_core",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_mann_kendall_core(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
