#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void tropopause_grid_2d(
    int nspatial,
    int nlev,
    int nlevm,
    void *pfull,
    void *tfull,
    double tmsg,
    double lapsec,
    int punit,
    void *ptrop_hpa,
    void *htrop_m,
    void *itrop,
    void *lapse_rate,
    void *success
);
void tropopause_grid_3d(
    int nlat,
    int nlon,
    int nlev,
    int nlevm,
    void *pfull,
    void *tfull,
    double tmsg,
    double lapsec,
    int punit,
    void *ptrop_hpa,
    void *htrop_m,
    void *itrop,
    void *lapse_rate,
    void *success
);
void tropopause_grid_4d(
    int nlat,
    int nlon,
    int nlev,
    int ntime,
    int nlevm,
    void *pfull,
    void *tfull,
    double tmsg,
    double lapsec,
    int punit,
    void *ptrop_hpa,
    void *htrop_m,
    void *itrop,
    void *lapse_rate,
    void *success
);
void tropopause_profile_1d(
    int nlev,
    int nlevm,
    void *pfull,
    void *tfull,
    double tmsg,
    double lapsec,
    int punit,
    void *ptrop_hpa,
    void *htrop_m,
    void *itrop,
    void *lapse_rate,
    void *success
);

static int parse_common_call(
    PyObject *args,
    PyObject *kwargs,
    int *nlevm,
    PyObject **pfull_obj,
    PyObject **tfull_obj,
    double *tmsg,
    double *lapsec,
    int *punit
) {
    static char *kwlist[] = {"nlevm", "pfull", "tfull", "tmsg", "lapsec", "punit", NULL};

    return PyArg_ParseTupleAndKeywords(
        args,
        kwargs,
        "iOOddi",
        kwlist,
        nlevm,
        pfull_obj,
        tfull_obj,
        tmsg,
        lapsec,
        punit
    );
}

static PyObject *build_grid_result(
    PyArrayObject *ptrop_hpa,
    PyArrayObject *htrop_m,
    PyArrayObject *itrop,
    PyArrayObject *lapse_rate,
    PyArrayObject *success
) {
    PyObject *result = PyTuple_Pack(
        5,
        (PyObject *) ptrop_hpa,
        (PyObject *) htrop_m,
        (PyObject *) itrop,
        (PyObject *) lapse_rate,
        (PyObject *) success
    );

    Py_DECREF(ptrop_hpa);
    Py_DECREF(htrop_m);
    Py_DECREF(itrop);
    Py_DECREF(lapse_rate);
    Py_DECREF(success);
    return result;
}

static int convert_success_to_bool(PyArrayObject *src, PyArrayObject *dst) {
    PyArrayObject *iter_src = NULL;
    PyArrayObject *iter_dst = NULL;
    NpyIter *src_iter = NULL;
    NpyIter *dst_iter = NULL;
    NpyIter_IterNextFunc *src_next = NULL;
    NpyIter_IterNextFunc *dst_next = NULL;
    char **src_data = NULL;
    char **dst_data = NULL;

    iter_src = (PyArrayObject *) PyArray_FromAny(
        (PyObject *) src,
        PyArray_DescrFromType(NPY_INT32),
        0,
        0,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS,
        NULL
    );
    if (iter_src == NULL) {
        return -1;
    }

    iter_dst = (PyArrayObject *) PyArray_FromAny(
        (PyObject *) dst,
        PyArray_DescrFromType(NPY_BOOL),
        0,
        0,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITEABLE,
        NULL
    );
    if (iter_dst == NULL) {
        Py_DECREF(iter_src);
        return -1;
    }

    src_iter = NpyIter_New(iter_src, NPY_ITER_READONLY, NPY_KEEPORDER, NPY_NO_CASTING, NULL);
    if (src_iter == NULL) {
        Py_DECREF(iter_src);
        Py_DECREF(iter_dst);
        return -1;
    }
    dst_iter = NpyIter_New(iter_dst, NPY_ITER_READWRITE, NPY_KEEPORDER, NPY_NO_CASTING, NULL);
    if (dst_iter == NULL) {
        NpyIter_Deallocate(src_iter);
        Py_DECREF(iter_src);
        Py_DECREF(iter_dst);
        return -1;
    }

    src_next = NpyIter_GetIterNext(src_iter, NULL);
    dst_next = NpyIter_GetIterNext(dst_iter, NULL);
    src_data = NpyIter_GetDataPtrArray(src_iter);
    dst_data = NpyIter_GetDataPtrArray(dst_iter);

    if (src_next == NULL || dst_next == NULL || src_data == NULL || dst_data == NULL) {
        NpyIter_Deallocate(src_iter);
        NpyIter_Deallocate(dst_iter);
        Py_DECREF(iter_src);
        Py_DECREF(iter_dst);
        return -1;
    }

    do {
        *((npy_bool *) *dst_data) = (*((int *) *src_data) != 0) ? NPY_TRUE : NPY_FALSE;
    } while (src_next(src_iter) && dst_next(dst_iter));

    NpyIter_Deallocate(src_iter);
    NpyIter_Deallocate(dst_iter);
    Py_DECREF(iter_src);
    Py_DECREF(iter_dst);
    return 0;
}

static PyObject *py_tropopause_grid_2d(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pfull_obj = NULL;
    PyObject *tfull_obj = NULL;
    PyArrayObject *pfull_arr = NULL;
    PyArrayObject *tfull_arr = NULL;
    PyArrayObject *ptrop_hpa_arr = NULL;
    PyArrayObject *htrop_m_arr = NULL;
    PyArrayObject *itrop_arr = NULL;
    PyArrayObject *lapse_rate_arr = NULL;
    PyArrayObject *success_tmp_arr = NULL;
    PyArrayObject *success_arr = NULL;
    npy_intp out_dims[1];
    int nspatial;
    int nlev;
    int nlevm;
    int punit;
    double tmsg;
    double lapsec;

    (void) self;

    if (!parse_common_call(args, kwargs, &nlevm, &pfull_obj, &tfull_obj, &tmsg, &lapsec, &punit)) {
        return NULL;
    }

    pfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    tfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        tfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    if (pfull_arr == NULL || tfull_arr == NULL) {
        Py_XDECREF(pfull_arr);
        Py_XDECREF(tfull_arr);
        return NULL;
    }

    if (PyArray_NDIM(pfull_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pfull must be a 1D array");
        goto fail;
    }
    if (PyArray_NDIM(tfull_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "tfull must be a 2D array");
        goto fail;
    }

    nspatial = (int) PyArray_DIM(tfull_arr, 0);
    nlev = (int) PyArray_DIM(tfull_arr, 1);
    if ((int) PyArray_DIM(pfull_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pfull length must match tfull.shape[1]");
        goto fail;
    }
    if (nlevm != nlev + 1) {
        PyErr_SetString(PyExc_ValueError, "nlevm must equal len(pfull) + 1");
        goto fail;
    }

    out_dims[0] = nspatial;
    ptrop_hpa_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_FLOAT64, 1);
    htrop_m_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_FLOAT64, 1);
    itrop_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_INT32, 1);
    lapse_rate_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_FLOAT64, 1);
    success_tmp_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_INT32, 1);
    success_arr = (PyArrayObject *) PyArray_EMPTY(1, out_dims, NPY_BOOL, 0);
    if (ptrop_hpa_arr == NULL || htrop_m_arr == NULL || itrop_arr == NULL ||
        lapse_rate_arr == NULL || success_tmp_arr == NULL || success_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    tropopause_grid_2d(
        nspatial,
        nlev,
        nlevm,
        PyArray_DATA(pfull_arr),
        PyArray_DATA(tfull_arr),
        tmsg,
        lapsec,
        punit,
        PyArray_DATA(ptrop_hpa_arr),
        PyArray_DATA(htrop_m_arr),
        PyArray_DATA(itrop_arr),
        PyArray_DATA(lapse_rate_arr),
        PyArray_DATA(success_tmp_arr)
    );
    Py_END_ALLOW_THREADS

    if (convert_success_to_bool(success_tmp_arr, success_arr) < 0) {
        goto fail;
    }

    Py_DECREF(pfull_arr);
    Py_DECREF(tfull_arr);
    Py_DECREF(success_tmp_arr);
    return build_grid_result(ptrop_hpa_arr, htrop_m_arr, itrop_arr, lapse_rate_arr, success_arr);

fail:
    Py_XDECREF(pfull_arr);
    Py_XDECREF(tfull_arr);
    Py_XDECREF(ptrop_hpa_arr);
    Py_XDECREF(htrop_m_arr);
    Py_XDECREF(itrop_arr);
    Py_XDECREF(lapse_rate_arr);
    Py_XDECREF(success_tmp_arr);
    Py_XDECREF(success_arr);
    return NULL;
}

static PyObject *py_tropopause_grid_3d(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pfull_obj = NULL;
    PyObject *tfull_obj = NULL;
    PyArrayObject *pfull_arr = NULL;
    PyArrayObject *tfull_arr = NULL;
    PyArrayObject *ptrop_hpa_arr = NULL;
    PyArrayObject *htrop_m_arr = NULL;
    PyArrayObject *itrop_arr = NULL;
    PyArrayObject *lapse_rate_arr = NULL;
    PyArrayObject *success_tmp_arr = NULL;
    PyArrayObject *success_arr = NULL;
    npy_intp out_dims[2];
    int nlat;
    int nlon;
    int nlev;
    int nlevm;
    int punit;
    double tmsg;
    double lapsec;

    (void) self;

    if (!parse_common_call(args, kwargs, &nlevm, &pfull_obj, &tfull_obj, &tmsg, &lapsec, &punit)) {
        return NULL;
    }

    pfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    tfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        tfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    if (pfull_arr == NULL || tfull_arr == NULL) {
        Py_XDECREF(pfull_arr);
        Py_XDECREF(tfull_arr);
        return NULL;
    }

    if (PyArray_NDIM(pfull_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pfull must be a 1D array");
        goto fail;
    }
    if (PyArray_NDIM(tfull_arr) != 3) {
        PyErr_SetString(PyExc_ValueError, "tfull must be a 3D array");
        goto fail;
    }

    nlat = (int) PyArray_DIM(tfull_arr, 0);
    nlon = (int) PyArray_DIM(tfull_arr, 1);
    nlev = (int) PyArray_DIM(tfull_arr, 2);
    if ((int) PyArray_DIM(pfull_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pfull length must match tfull.shape[2]");
        goto fail;
    }
    if (nlevm != nlev + 1) {
        PyErr_SetString(PyExc_ValueError, "nlevm must equal len(pfull) + 1");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    ptrop_hpa_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_FLOAT64, 1);
    htrop_m_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_FLOAT64, 1);
    itrop_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_INT32, 1);
    lapse_rate_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_FLOAT64, 1);
    success_tmp_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_INT32, 1);
    success_arr = (PyArrayObject *) PyArray_EMPTY(2, out_dims, NPY_BOOL, 0);
    if (ptrop_hpa_arr == NULL || htrop_m_arr == NULL || itrop_arr == NULL ||
        lapse_rate_arr == NULL || success_tmp_arr == NULL || success_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    tropopause_grid_3d(
        nlat,
        nlon,
        nlev,
        nlevm,
        PyArray_DATA(pfull_arr),
        PyArray_DATA(tfull_arr),
        tmsg,
        lapsec,
        punit,
        PyArray_DATA(ptrop_hpa_arr),
        PyArray_DATA(htrop_m_arr),
        PyArray_DATA(itrop_arr),
        PyArray_DATA(lapse_rate_arr),
        PyArray_DATA(success_tmp_arr)
    );
    Py_END_ALLOW_THREADS

    if (convert_success_to_bool(success_tmp_arr, success_arr) < 0) {
        goto fail;
    }

    Py_DECREF(pfull_arr);
    Py_DECREF(tfull_arr);
    Py_DECREF(success_tmp_arr);
    return build_grid_result(ptrop_hpa_arr, htrop_m_arr, itrop_arr, lapse_rate_arr, success_arr);

fail:
    Py_XDECREF(pfull_arr);
    Py_XDECREF(tfull_arr);
    Py_XDECREF(ptrop_hpa_arr);
    Py_XDECREF(htrop_m_arr);
    Py_XDECREF(itrop_arr);
    Py_XDECREF(lapse_rate_arr);
    Py_XDECREF(success_tmp_arr);
    Py_XDECREF(success_arr);
    return NULL;
}

static PyObject *py_tropopause_grid_4d(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pfull_obj = NULL;
    PyObject *tfull_obj = NULL;
    PyArrayObject *pfull_arr = NULL;
    PyArrayObject *tfull_arr = NULL;
    PyArrayObject *ptrop_hpa_arr = NULL;
    PyArrayObject *htrop_m_arr = NULL;
    PyArrayObject *itrop_arr = NULL;
    PyArrayObject *lapse_rate_arr = NULL;
    PyArrayObject *success_tmp_arr = NULL;
    PyArrayObject *success_arr = NULL;
    npy_intp out_dims[3];
    int nlat;
    int nlon;
    int nlev;
    int ntime;
    int nlevm;
    int punit;
    double tmsg;
    double lapsec;

    (void) self;

    if (!parse_common_call(args, kwargs, &nlevm, &pfull_obj, &tfull_obj, &tmsg, &lapsec, &punit)) {
        return NULL;
    }

    pfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    tfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        tfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    if (pfull_arr == NULL || tfull_arr == NULL) {
        Py_XDECREF(pfull_arr);
        Py_XDECREF(tfull_arr);
        return NULL;
    }

    if (PyArray_NDIM(pfull_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pfull must be a 1D array");
        goto fail;
    }
    if (PyArray_NDIM(tfull_arr) != 4) {
        PyErr_SetString(PyExc_ValueError, "tfull must be a 4D array");
        goto fail;
    }

    nlat = (int) PyArray_DIM(tfull_arr, 0);
    nlon = (int) PyArray_DIM(tfull_arr, 1);
    nlev = (int) PyArray_DIM(tfull_arr, 2);
    ntime = (int) PyArray_DIM(tfull_arr, 3);
    if ((int) PyArray_DIM(pfull_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pfull length must match tfull.shape[2]");
        goto fail;
    }
    if (nlevm != nlev + 1) {
        PyErr_SetString(PyExc_ValueError, "nlevm must equal len(pfull) + 1");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    out_dims[2] = ntime;
    ptrop_hpa_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_FLOAT64, 1);
    htrop_m_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_FLOAT64, 1);
    itrop_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_INT32, 1);
    lapse_rate_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_FLOAT64, 1);
    success_tmp_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_INT32, 1);
    success_arr = (PyArrayObject *) PyArray_EMPTY(3, out_dims, NPY_BOOL, 0);
    if (ptrop_hpa_arr == NULL || htrop_m_arr == NULL || itrop_arr == NULL ||
        lapse_rate_arr == NULL || success_tmp_arr == NULL || success_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    tropopause_grid_4d(
        nlat,
        nlon,
        nlev,
        ntime,
        nlevm,
        PyArray_DATA(pfull_arr),
        PyArray_DATA(tfull_arr),
        tmsg,
        lapsec,
        punit,
        PyArray_DATA(ptrop_hpa_arr),
        PyArray_DATA(htrop_m_arr),
        PyArray_DATA(itrop_arr),
        PyArray_DATA(lapse_rate_arr),
        PyArray_DATA(success_tmp_arr)
    );
    Py_END_ALLOW_THREADS

    if (convert_success_to_bool(success_tmp_arr, success_arr) < 0) {
        goto fail;
    }

    Py_DECREF(pfull_arr);
    Py_DECREF(tfull_arr);
    Py_DECREF(success_tmp_arr);
    return build_grid_result(ptrop_hpa_arr, htrop_m_arr, itrop_arr, lapse_rate_arr, success_arr);

fail:
    Py_XDECREF(pfull_arr);
    Py_XDECREF(tfull_arr);
    Py_XDECREF(ptrop_hpa_arr);
    Py_XDECREF(htrop_m_arr);
    Py_XDECREF(itrop_arr);
    Py_XDECREF(lapse_rate_arr);
    Py_XDECREF(success_tmp_arr);
    Py_XDECREF(success_arr);
    return NULL;
}

static PyObject *py_tropopause_profile_1d(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *pfull_obj = NULL;
    PyObject *tfull_obj = NULL;
    PyArrayObject *pfull_arr = NULL;
    PyArrayObject *tfull_arr = NULL;
    double ptrop_hpa;
    double htrop_m;
    double lapse_rate;
    int itrop;
    int success;
    int nlev;
    int nlevm;
    int punit;
    double tmsg;
    double lapsec;

    (void) self;

    if (!parse_common_call(args, kwargs, &nlevm, &pfull_obj, &tfull_obj, &tmsg, &lapsec, &punit)) {
        return NULL;
    }

    pfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    tfull_arr = (PyArrayObject *) PyArray_FROM_OTF(
        tfull_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    if (pfull_arr == NULL || tfull_arr == NULL) {
        Py_XDECREF(pfull_arr);
        Py_XDECREF(tfull_arr);
        return NULL;
    }

    if (PyArray_NDIM(pfull_arr) != 1 || PyArray_NDIM(tfull_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pfull and tfull must be 1D arrays");
        goto fail;
    }

    nlev = (int) PyArray_DIM(pfull_arr, 0);
    if ((int) PyArray_DIM(tfull_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pfull and tfull must have the same length");
        goto fail;
    }
    if (nlevm != nlev + 1) {
        PyErr_SetString(PyExc_ValueError, "nlevm must equal len(pfull) + 1");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    tropopause_profile_1d(
        nlev,
        nlevm,
        PyArray_DATA(pfull_arr),
        PyArray_DATA(tfull_arr),
        tmsg,
        lapsec,
        punit,
        &ptrop_hpa,
        &htrop_m,
        &itrop,
        &lapse_rate,
        &success
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(pfull_arr);
    Py_DECREF(tfull_arr);
    return PyTuple_Pack(
        5,
        PyFloat_FromDouble(ptrop_hpa),
        PyFloat_FromDouble(htrop_m),
        PyLong_FromLong((long) itrop),
        PyFloat_FromDouble(lapse_rate),
        PyBool_FromLong(success != 0)
    );

fail:
    Py_XDECREF(pfull_arr);
    Py_XDECREF(tfull_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {
        "tropopause_grid_2d",
        (PyCFunction) (void (*)(void)) py_tropopause_grid_2d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute WMO tropopause properties for 2D cross-sections."
    },
    {
        "tropopause_grid_3d",
        (PyCFunction) (void (*)(void)) py_tropopause_grid_3d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute WMO tropopause properties for 3D grids."
    },
    {
        "tropopause_grid_4d",
        (PyCFunction) (void (*)(void)) py_tropopause_grid_4d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute WMO tropopause properties for 4D time series."
    },
    {
        "tropopause_profile_1d",
        (PyCFunction) (void (*)(void)) py_tropopause_profile_1d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute WMO tropopause properties for a single profile."
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tropopause_height",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_tropopause_height(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
