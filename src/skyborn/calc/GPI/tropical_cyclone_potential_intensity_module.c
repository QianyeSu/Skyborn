#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void calculate_pi_gridded_data(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    int nlat,
    int nlon,
    int num_levels
);
void calculate_pi_gridded_with_missing(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    int nlat,
    int nlon,
    int num_levels
);
void calculate_pi_gridded_diagnostics(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    void *outflow_temp,
    void *outflow_level,
    void *lnpi,
    void *lneff,
    void *lndiseq,
    void *lnckcd,
    int outflow_source_flag,
    float ckcd_in,
    int nlat,
    int nlon,
    int num_levels
);
void calculate_pi_gridded_diagnostics_with_missing(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    void *outflow_temp,
    void *outflow_level,
    void *lnpi,
    void *lneff,
    void *lndiseq,
    void *lnckcd,
    int outflow_source_flag,
    float ckcd_in,
    int nlat,
    int nlon,
    int num_levels
);
void calculate_pi_single_profile(
    float sst_in,
    float psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    int actual_levels,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    int num_levels
);
void calculate_pi_profile_diagnostics(
    float sst_in,
    float psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    int actual_levels,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    void *outflow_temp,
    void *outflow_level,
    void *lnpi,
    void *lneff,
    void *lndiseq,
    void *lnckcd,
    int outflow_source_flag,
    float ckcd_in,
    int num_levels
);
void calculate_pi_4d_data(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    int nlat,
    int nlon,
    int num_levels,
    int num_times
);
void calculate_pi_4d_with_missing(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    int nlat,
    int nlon,
    int num_levels,
    int num_times
);
void calculate_pi_4d_diagnostics(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    void *outflow_temp,
    void *outflow_level,
    void *lnpi,
    void *lneff,
    void *lndiseq,
    void *lnckcd,
    int outflow_source_flag,
    float ckcd_in,
    int nlat,
    int nlon,
    int num_levels,
    int num_times
);
void calculate_pi_4d_diagnostics_with_missing(
    void *sst_in,
    void *psl_in,
    void *pressure_levels,
    void *temp_in,
    void *mixing_ratio_in,
    void *min_pressure,
    void *max_wind,
    void *error_flag,
    void *outflow_temp,
    void *outflow_level,
    void *lnpi,
    void *lneff,
    void *lndiseq,
    void *lnckcd,
    int outflow_source_flag,
    float ckcd_in,
    int nlat,
    int nlon,
    int num_levels,
    int num_times
);
void cape(
    float parcel_temp,
    float parcel_mixing_ratio,
    float parcel_pressure,
    void *temp_profile,
    void *mixing_ratio_profile,
    void *pressure_profile,
    int num_points,
    float buoyancy_param,
    void *cape_value,
    void *outflow_temp,
    void *outflow_level,
    void *error_flag,
    int array_size
);

static PyArrayObject *to_float32_1d(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT32,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
}

static PyArrayObject *to_float32_2d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT32,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
}

static PyArrayObject *to_float32_3d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT32,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
}

static PyArrayObject *to_float32_4d_fortran(PyObject *obj) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT32,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
}

static int require_ndim(PyArrayObject *arr, int ndim, const char *name) {
    if (PyArray_NDIM(arr) != ndim) {
        PyErr_Format(PyExc_ValueError, "%s must be a %dD array", name, ndim);
        return -1;
    }
    return 0;
}

static PyArrayObject *new_float32_array(int ndim, const npy_intp *dims) {
    return (PyArrayObject *) PyArray_EMPTY(ndim, (npy_intp *) dims, NPY_FLOAT32, 1);
}

static PyObject *build_gridded_result(PyArrayObject *min_pressure, PyArrayObject *max_wind, int error_flag) {
    PyObject *result = Py_BuildValue("OOi", min_pressure, max_wind, error_flag);
    Py_DECREF(min_pressure);
    Py_DECREF(max_wind);
    return result;
}

static PyObject *build_gridded_diag_result(
    PyArrayObject *min_pressure,
    PyArrayObject *max_wind,
    int error_flag,
    PyArrayObject *outflow_temp,
    PyArrayObject *outflow_level,
    PyArrayObject *lnpi,
    PyArrayObject *lneff,
    PyArrayObject *lndiseq,
    float lnckcd
) {
    PyObject *result = Py_BuildValue(
        "OOiOOOOOf",
        min_pressure,
        max_wind,
        error_flag,
        outflow_temp,
        outflow_level,
        lnpi,
        lneff,
        lndiseq,
        lnckcd
    );
    Py_DECREF(min_pressure);
    Py_DECREF(max_wind);
    Py_DECREF(outflow_temp);
    Py_DECREF(outflow_level);
    Py_DECREF(lnpi);
    Py_DECREF(lneff);
    Py_DECREF(lndiseq);
    return result;
}

static PyObject *py_calculate_pi_gridded_data(PyObject *self, PyObject *args) {
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    npy_intp out_dims[2];
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOO", &sst_obj, &psl_obj, &pressure_obj, &temp_obj, &mixing_obj)) {
        return NULL;
    }

    sst_arr = to_float32_2d_fortran(sst_obj);
    psl_arr = to_float32_2d_fortran(psl_obj);
    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_3d_fortran(temp_obj);
    mixing_arr = to_float32_3d_fortran(mixing_obj);
    if (sst_arr == NULL || psl_arr == NULL || pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(sst_arr, 2, "sst_in") != 0 ||
        require_ndim(psl_arr, 2, "psl_in") != 0 ||
        require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 3, "temp_in") != 0 ||
        require_ndim(mixing_arr, 3, "mixing_ratio_in") != 0) {
        goto fail;
    }

    nlat = (int) PyArray_DIM(sst_arr, 0);
    nlon = (int) PyArray_DIM(sst_arr, 1);
    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(psl_arr, 0) != nlat || PyArray_DIM(psl_arr, 1) != nlon) {
        PyErr_SetString(PyExc_ValueError, "psl_in must match sst_in shape");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(temp_arr, 1) != nlat || PyArray_DIM(temp_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "temp_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }
    if (PyArray_DIM(mixing_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 1) != nlat || PyArray_DIM(mixing_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "mixing_ratio_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    min_pressure_arr = new_float32_array(2, out_dims);
    max_wind_arr = new_float32_array(2, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_gridded_data(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        nlat,
        nlon,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_result(min_pressure_arr, max_wind_arr, error_flag);

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    return NULL;
}

static PyObject *py_calculate_pi_gridded_with_missing(PyObject *self, PyObject *args) {
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    npy_intp out_dims[2];
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOO", &sst_obj, &psl_obj, &pressure_obj, &temp_obj, &mixing_obj)) {
        return NULL;
    }

    sst_arr = to_float32_2d_fortran(sst_obj);
    psl_arr = to_float32_2d_fortran(psl_obj);
    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_3d_fortran(temp_obj);
    mixing_arr = to_float32_3d_fortran(mixing_obj);
    if (sst_arr == NULL || psl_arr == NULL || pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(sst_arr, 2, "sst_in") != 0 ||
        require_ndim(psl_arr, 2, "psl_in") != 0 ||
        require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 3, "temp_in") != 0 ||
        require_ndim(mixing_arr, 3, "mixing_ratio_in") != 0) {
        goto fail;
    }

    nlat = (int) PyArray_DIM(sst_arr, 0);
    nlon = (int) PyArray_DIM(sst_arr, 1);
    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(psl_arr, 0) != nlat || PyArray_DIM(psl_arr, 1) != nlon) {
        PyErr_SetString(PyExc_ValueError, "psl_in must match sst_in shape");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(temp_arr, 1) != nlat || PyArray_DIM(temp_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "temp_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }
    if (PyArray_DIM(mixing_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 1) != nlat || PyArray_DIM(mixing_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "mixing_ratio_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    min_pressure_arr = new_float32_array(2, out_dims);
    max_wind_arr = new_float32_array(2, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_gridded_with_missing(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        nlat,
        nlon,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_result(min_pressure_arr, max_wind_arr, error_flag);

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    return NULL;
}

static PyObject *py_calculate_pi_gridded_diagnostics(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "sst_in",
        "psl_in",
        "pressure_levels",
        "temp_in",
        "mixing_ratio_in",
        "outflow_source_flag",
        "ckcd_in",
        NULL,
    };
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    PyArrayObject *outflow_temp_arr = NULL;
    PyArrayObject *outflow_level_arr = NULL;
    PyArrayObject *lnpi_arr = NULL;
    PyArrayObject *lneff_arr = NULL;
    PyArrayObject *lndiseq_arr = NULL;
    npy_intp out_dims[2];
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;
    int outflow_source_flag = 0;
    float ckcd_in = 0.9f;
    float lnckcd = 0.0f;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|if",
            kwlist,
            &sst_obj,
            &psl_obj,
            &pressure_obj,
            &temp_obj,
            &mixing_obj,
            &outflow_source_flag,
            &ckcd_in
        )) {
        return NULL;
    }

    sst_arr = to_float32_2d_fortran(sst_obj);
    psl_arr = to_float32_2d_fortran(psl_obj);
    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_3d_fortran(temp_obj);
    mixing_arr = to_float32_3d_fortran(mixing_obj);
    if (sst_arr == NULL || psl_arr == NULL || pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(sst_arr, 2, "sst_in") != 0 ||
        require_ndim(psl_arr, 2, "psl_in") != 0 ||
        require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 3, "temp_in") != 0 ||
        require_ndim(mixing_arr, 3, "mixing_ratio_in") != 0) {
        goto fail;
    }

    nlat = (int) PyArray_DIM(sst_arr, 0);
    nlon = (int) PyArray_DIM(sst_arr, 1);
    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(psl_arr, 0) != nlat || PyArray_DIM(psl_arr, 1) != nlon) {
        PyErr_SetString(PyExc_ValueError, "psl_in must match sst_in shape");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(temp_arr, 1) != nlat || PyArray_DIM(temp_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "temp_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }
    if (PyArray_DIM(mixing_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 1) != nlat || PyArray_DIM(mixing_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "mixing_ratio_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    min_pressure_arr = new_float32_array(2, out_dims);
    max_wind_arr = new_float32_array(2, out_dims);
    outflow_temp_arr = new_float32_array(2, out_dims);
    outflow_level_arr = new_float32_array(2, out_dims);
    lnpi_arr = new_float32_array(2, out_dims);
    lneff_arr = new_float32_array(2, out_dims);
    lndiseq_arr = new_float32_array(2, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL || outflow_temp_arr == NULL ||
        outflow_level_arr == NULL || lnpi_arr == NULL || lneff_arr == NULL || lndiseq_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_gridded_diagnostics(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        PyArray_DATA(outflow_temp_arr),
        PyArray_DATA(outflow_level_arr),
        PyArray_DATA(lnpi_arr),
        PyArray_DATA(lneff_arr),
        PyArray_DATA(lndiseq_arr),
        &lnckcd,
        outflow_source_flag,
        ckcd_in,
        nlat,
        nlon,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_diag_result(
        min_pressure_arr,
        max_wind_arr,
        error_flag,
        outflow_temp_arr,
        outflow_level_arr,
        lnpi_arr,
        lneff_arr,
        lndiseq_arr,
        lnckcd
    );

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    Py_XDECREF(outflow_temp_arr);
    Py_XDECREF(outflow_level_arr);
    Py_XDECREF(lnpi_arr);
    Py_XDECREF(lneff_arr);
    Py_XDECREF(lndiseq_arr);
    return NULL;
}

static PyObject *py_calculate_pi_gridded_diagnostics_with_missing(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "sst_in",
        "psl_in",
        "pressure_levels",
        "temp_in",
        "mixing_ratio_in",
        "outflow_source_flag",
        "ckcd_in",
        NULL,
    };
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    PyArrayObject *outflow_temp_arr = NULL;
    PyArrayObject *outflow_level_arr = NULL;
    PyArrayObject *lnpi_arr = NULL;
    PyArrayObject *lneff_arr = NULL;
    PyArrayObject *lndiseq_arr = NULL;
    npy_intp out_dims[2];
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;
    int outflow_source_flag = 0;
    float ckcd_in = 0.9f;
    float lnckcd = 0.0f;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|if",
            kwlist,
            &sst_obj,
            &psl_obj,
            &pressure_obj,
            &temp_obj,
            &mixing_obj,
            &outflow_source_flag,
            &ckcd_in
        )) {
        return NULL;
    }

    sst_arr = to_float32_2d_fortran(sst_obj);
    psl_arr = to_float32_2d_fortran(psl_obj);
    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_3d_fortran(temp_obj);
    mixing_arr = to_float32_3d_fortran(mixing_obj);
    if (sst_arr == NULL || psl_arr == NULL || pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(sst_arr, 2, "sst_in") != 0 ||
        require_ndim(psl_arr, 2, "psl_in") != 0 ||
        require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 3, "temp_in") != 0 ||
        require_ndim(mixing_arr, 3, "mixing_ratio_in") != 0) {
        goto fail;
    }

    nlat = (int) PyArray_DIM(sst_arr, 0);
    nlon = (int) PyArray_DIM(sst_arr, 1);
    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(psl_arr, 0) != nlat || PyArray_DIM(psl_arr, 1) != nlon) {
        PyErr_SetString(PyExc_ValueError, "psl_in must match sst_in shape");
        goto fail;
    }
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(temp_arr, 1) != nlat || PyArray_DIM(temp_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "temp_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }
    if (PyArray_DIM(mixing_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 1) != nlat || PyArray_DIM(mixing_arr, 2) != nlon) {
        PyErr_SetString(PyExc_ValueError, "mixing_ratio_in must have shape (num_levels, nlat, nlon)");
        goto fail;
    }

    out_dims[0] = nlat;
    out_dims[1] = nlon;
    min_pressure_arr = new_float32_array(2, out_dims);
    max_wind_arr = new_float32_array(2, out_dims);
    outflow_temp_arr = new_float32_array(2, out_dims);
    outflow_level_arr = new_float32_array(2, out_dims);
    lnpi_arr = new_float32_array(2, out_dims);
    lneff_arr = new_float32_array(2, out_dims);
    lndiseq_arr = new_float32_array(2, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL || outflow_temp_arr == NULL ||
        outflow_level_arr == NULL || lnpi_arr == NULL || lneff_arr == NULL || lndiseq_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_gridded_diagnostics_with_missing(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        PyArray_DATA(outflow_temp_arr),
        PyArray_DATA(outflow_level_arr),
        PyArray_DATA(lnpi_arr),
        PyArray_DATA(lneff_arr),
        PyArray_DATA(lndiseq_arr),
        &lnckcd,
        outflow_source_flag,
        ckcd_in,
        nlat,
        nlon,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_diag_result(
        min_pressure_arr,
        max_wind_arr,
        error_flag,
        outflow_temp_arr,
        outflow_level_arr,
        lnpi_arr,
        lneff_arr,
        lndiseq_arr,
        lnckcd
    );

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    Py_XDECREF(outflow_temp_arr);
    Py_XDECREF(outflow_level_arr);
    Py_XDECREF(lnpi_arr);
    Py_XDECREF(lneff_arr);
    Py_XDECREF(lndiseq_arr);
    return NULL;
}

static PyObject *py_calculate_pi_single_profile(PyObject *self, PyObject *args) {
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    float sst_in = 0.0f;
    float psl_in = 0.0f;
    float min_pressure = 0.0f;
    float max_wind = 0.0f;
    int error_flag = 0;
    int actual_levels = 0;
    int num_levels = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "ffOOOi", &sst_in, &psl_in, &pressure_obj, &temp_obj, &mixing_obj, &actual_levels)) {
        return NULL;
    }

    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_1d(temp_obj);
    mixing_arr = to_float32_1d(mixing_obj);
    if (pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 1, "temp_in") != 0 ||
        require_ndim(mixing_arr, 1, "mixing_ratio_in") != 0) {
        goto fail;
    }

    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 0) != num_levels) {
        PyErr_SetString(PyExc_ValueError, "profile inputs must share the same length");
        goto fail;
    }
    if (actual_levels < 1 || actual_levels > num_levels) {
        PyErr_SetString(PyExc_ValueError, "actual_levels must be between 1 and the profile length");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_single_profile(
        sst_in,
        psl_in,
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        actual_levels,
        &min_pressure,
        &max_wind,
        &error_flag,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return Py_BuildValue("ffi", min_pressure, max_wind, error_flag);

fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    return NULL;
}

static PyObject *py_calculate_pi_profile_diagnostics(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "sst_in",
        "psl_in",
        "pressure_levels",
        "temp_in",
        "mixing_ratio_in",
        "actual_levels",
        "outflow_source_flag",
        "ckcd_in",
        NULL,
    };
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    float sst_in = 0.0f;
    float psl_in = 0.0f;
    float min_pressure = 0.0f;
    float max_wind = 0.0f;
    float outflow_temp = 0.0f;
    float outflow_level = 0.0f;
    float lnpi = 0.0f;
    float lneff = 0.0f;
    float lndiseq = 0.0f;
    float lnckcd = 0.0f;
    int error_flag = 0;
    int actual_levels = 0;
    int outflow_source_flag = 0;
    int num_levels = 0;
    float ckcd_in = 0.9f;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "ffOOOi|if",
            kwlist,
            &sst_in,
            &psl_in,
            &pressure_obj,
            &temp_obj,
            &mixing_obj,
            &actual_levels,
            &outflow_source_flag,
            &ckcd_in
        )) {
        return NULL;
    }

    pressure_arr = to_float32_1d(pressure_obj);
    temp_arr = to_float32_1d(temp_obj);
    mixing_arr = to_float32_1d(mixing_obj);
    if (pressure_arr == NULL || temp_arr == NULL || mixing_arr == NULL) {
        goto fail;
    }

    if (require_ndim(pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(temp_arr, 1, "temp_in") != 0 ||
        require_ndim(mixing_arr, 1, "mixing_ratio_in") != 0) {
        goto fail;
    }

    num_levels = (int) PyArray_DIM(pressure_arr, 0);
    if (PyArray_DIM(temp_arr, 0) != num_levels || PyArray_DIM(mixing_arr, 0) != num_levels) {
        PyErr_SetString(PyExc_ValueError, "profile inputs must share the same length");
        goto fail;
    }
    if (actual_levels < 1 || actual_levels > num_levels) {
        PyErr_SetString(PyExc_ValueError, "actual_levels must be between 1 and the profile length");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_profile_diagnostics(
        sst_in,
        psl_in,
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        actual_levels,
        &min_pressure,
        &max_wind,
        &error_flag,
        &outflow_temp,
        &outflow_level,
        &lnpi,
        &lneff,
        &lndiseq,
        &lnckcd,
        outflow_source_flag,
        ckcd_in,
        num_levels
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return Py_BuildValue(
        "ffiffffff",
        min_pressure,
        max_wind,
        error_flag,
        outflow_temp,
        outflow_level,
        lnpi,
        lneff,
        lndiseq,
        lnckcd
    );

fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    return NULL;
}

static int prepare_4d_inputs(
    PyObject *sst_obj,
    PyObject *psl_obj,
    PyObject *pressure_obj,
    PyObject *temp_obj,
    PyObject *mixing_obj,
    PyArrayObject **sst_arr,
    PyArrayObject **psl_arr,
    PyArrayObject **pressure_arr,
    PyArrayObject **temp_arr,
    PyArrayObject **mixing_arr,
    int *num_times,
    int *nlat,
    int *nlon,
    int *num_levels
) {
    *sst_arr = to_float32_3d_fortran(sst_obj);
    *psl_arr = to_float32_3d_fortran(psl_obj);
    *pressure_arr = to_float32_1d(pressure_obj);
    *temp_arr = to_float32_4d_fortran(temp_obj);
    *mixing_arr = to_float32_4d_fortran(mixing_obj);
    if (*sst_arr == NULL || *psl_arr == NULL || *pressure_arr == NULL || *temp_arr == NULL || *mixing_arr == NULL) {
        return -1;
    }

    if (require_ndim(*sst_arr, 3, "sst_in") != 0 ||
        require_ndim(*psl_arr, 3, "psl_in") != 0 ||
        require_ndim(*pressure_arr, 1, "pressure_levels") != 0 ||
        require_ndim(*temp_arr, 4, "temp_in") != 0 ||
        require_ndim(*mixing_arr, 4, "mixing_ratio_in") != 0) {
        return -1;
    }

    *num_times = (int) PyArray_DIM(*sst_arr, 0);
    *nlat = (int) PyArray_DIM(*sst_arr, 1);
    *nlon = (int) PyArray_DIM(*sst_arr, 2);
    *num_levels = (int) PyArray_DIM(*pressure_arr, 0);
    if (PyArray_DIM(*psl_arr, 0) != *num_times || PyArray_DIM(*psl_arr, 1) != *nlat || PyArray_DIM(*psl_arr, 2) != *nlon) {
        PyErr_SetString(PyExc_ValueError, "psl_in must match sst_in shape");
        return -1;
    }
    if (PyArray_DIM(*temp_arr, 0) != *num_times || PyArray_DIM(*temp_arr, 1) != *num_levels ||
        PyArray_DIM(*temp_arr, 2) != *nlat || PyArray_DIM(*temp_arr, 3) != *nlon) {
        PyErr_SetString(PyExc_ValueError, "temp_in must have shape (num_times, num_levels, nlat, nlon)");
        return -1;
    }
    if (PyArray_DIM(*mixing_arr, 0) != *num_times || PyArray_DIM(*mixing_arr, 1) != *num_levels ||
        PyArray_DIM(*mixing_arr, 2) != *nlat || PyArray_DIM(*mixing_arr, 3) != *nlon) {
        PyErr_SetString(PyExc_ValueError, "mixing_ratio_in must have shape (num_times, num_levels, nlat, nlon)");
        return -1;
    }

    return 0;
}

static PyObject *py_calculate_pi_4d_data(PyObject *self, PyObject *args) {
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    npy_intp out_dims[3];
    int num_times = 0;
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOO", &sst_obj, &psl_obj, &pressure_obj, &temp_obj, &mixing_obj)) {
        return NULL;
    }

    if (prepare_4d_inputs(
            sst_obj,
            psl_obj,
            pressure_obj,
            temp_obj,
            mixing_obj,
            &sst_arr,
            &psl_arr,
            &pressure_arr,
            &temp_arr,
            &mixing_arr,
            &num_times,
            &nlat,
            &nlon,
            &num_levels
        ) != 0) {
        goto fail;
    }

    out_dims[0] = num_times;
    out_dims[1] = nlat;
    out_dims[2] = nlon;
    min_pressure_arr = new_float32_array(3, out_dims);
    max_wind_arr = new_float32_array(3, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_4d_data(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        nlat,
        nlon,
        num_levels,
        num_times
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_result(min_pressure_arr, max_wind_arr, error_flag);

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    return NULL;
}

static PyObject *py_calculate_pi_4d_with_missing(PyObject *self, PyObject *args) {
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    npy_intp out_dims[3];
    int num_times = 0;
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOO", &sst_obj, &psl_obj, &pressure_obj, &temp_obj, &mixing_obj)) {
        return NULL;
    }

    if (prepare_4d_inputs(
            sst_obj,
            psl_obj,
            pressure_obj,
            temp_obj,
            mixing_obj,
            &sst_arr,
            &psl_arr,
            &pressure_arr,
            &temp_arr,
            &mixing_arr,
            &num_times,
            &nlat,
            &nlon,
            &num_levels
        ) != 0) {
        goto fail;
    }

    out_dims[0] = num_times;
    out_dims[1] = nlat;
    out_dims[2] = nlon;
    min_pressure_arr = new_float32_array(3, out_dims);
    max_wind_arr = new_float32_array(3, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_4d_with_missing(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        nlat,
        nlon,
        num_levels,
        num_times
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_result(min_pressure_arr, max_wind_arr, error_flag);

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    return NULL;
}

static PyObject *py_calculate_pi_4d_diagnostics(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "sst_in",
        "psl_in",
        "pressure_levels",
        "temp_in",
        "mixing_ratio_in",
        "outflow_source_flag",
        "ckcd_in",
        NULL,
    };
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    PyArrayObject *outflow_temp_arr = NULL;
    PyArrayObject *outflow_level_arr = NULL;
    PyArrayObject *lnpi_arr = NULL;
    PyArrayObject *lneff_arr = NULL;
    PyArrayObject *lndiseq_arr = NULL;
    npy_intp out_dims[3];
    int num_times = 0;
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;
    int outflow_source_flag = 0;
    float ckcd_in = 0.9f;
    float lnckcd = 0.0f;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|if",
            kwlist,
            &sst_obj,
            &psl_obj,
            &pressure_obj,
            &temp_obj,
            &mixing_obj,
            &outflow_source_flag,
            &ckcd_in
        )) {
        return NULL;
    }

    if (prepare_4d_inputs(
            sst_obj,
            psl_obj,
            pressure_obj,
            temp_obj,
            mixing_obj,
            &sst_arr,
            &psl_arr,
            &pressure_arr,
            &temp_arr,
            &mixing_arr,
            &num_times,
            &nlat,
            &nlon,
            &num_levels
        ) != 0) {
        goto fail;
    }

    out_dims[0] = num_times;
    out_dims[1] = nlat;
    out_dims[2] = nlon;
    min_pressure_arr = new_float32_array(3, out_dims);
    max_wind_arr = new_float32_array(3, out_dims);
    outflow_temp_arr = new_float32_array(3, out_dims);
    outflow_level_arr = new_float32_array(3, out_dims);
    lnpi_arr = new_float32_array(3, out_dims);
    lneff_arr = new_float32_array(3, out_dims);
    lndiseq_arr = new_float32_array(3, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL || outflow_temp_arr == NULL ||
        outflow_level_arr == NULL || lnpi_arr == NULL || lneff_arr == NULL || lndiseq_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_4d_diagnostics(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        PyArray_DATA(outflow_temp_arr),
        PyArray_DATA(outflow_level_arr),
        PyArray_DATA(lnpi_arr),
        PyArray_DATA(lneff_arr),
        PyArray_DATA(lndiseq_arr),
        &lnckcd,
        outflow_source_flag,
        ckcd_in,
        nlat,
        nlon,
        num_levels,
        num_times
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_diag_result(
        min_pressure_arr,
        max_wind_arr,
        error_flag,
        outflow_temp_arr,
        outflow_level_arr,
        lnpi_arr,
        lneff_arr,
        lndiseq_arr,
        lnckcd
    );

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    Py_XDECREF(outflow_temp_arr);
    Py_XDECREF(outflow_level_arr);
    Py_XDECREF(lnpi_arr);
    Py_XDECREF(lneff_arr);
    Py_XDECREF(lndiseq_arr);
    return NULL;
}

static PyObject *py_calculate_pi_4d_diagnostics_with_missing(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "sst_in",
        "psl_in",
        "pressure_levels",
        "temp_in",
        "mixing_ratio_in",
        "outflow_source_flag",
        "ckcd_in",
        NULL,
    };
    PyObject *sst_obj = NULL;
    PyObject *psl_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyArrayObject *sst_arr = NULL;
    PyArrayObject *psl_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *min_pressure_arr = NULL;
    PyArrayObject *max_wind_arr = NULL;
    PyArrayObject *outflow_temp_arr = NULL;
    PyArrayObject *outflow_level_arr = NULL;
    PyArrayObject *lnpi_arr = NULL;
    PyArrayObject *lneff_arr = NULL;
    PyArrayObject *lndiseq_arr = NULL;
    npy_intp out_dims[3];
    int num_times = 0;
    int nlat = 0;
    int nlon = 0;
    int num_levels = 0;
    int error_flag = 0;
    int outflow_source_flag = 0;
    float ckcd_in = 0.9f;
    float lnckcd = 0.0f;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOO|if",
            kwlist,
            &sst_obj,
            &psl_obj,
            &pressure_obj,
            &temp_obj,
            &mixing_obj,
            &outflow_source_flag,
            &ckcd_in
        )) {
        return NULL;
    }

    if (prepare_4d_inputs(
            sst_obj,
            psl_obj,
            pressure_obj,
            temp_obj,
            mixing_obj,
            &sst_arr,
            &psl_arr,
            &pressure_arr,
            &temp_arr,
            &mixing_arr,
            &num_times,
            &nlat,
            &nlon,
            &num_levels
        ) != 0) {
        goto fail;
    }

    out_dims[0] = num_times;
    out_dims[1] = nlat;
    out_dims[2] = nlon;
    min_pressure_arr = new_float32_array(3, out_dims);
    max_wind_arr = new_float32_array(3, out_dims);
    outflow_temp_arr = new_float32_array(3, out_dims);
    outflow_level_arr = new_float32_array(3, out_dims);
    lnpi_arr = new_float32_array(3, out_dims);
    lneff_arr = new_float32_array(3, out_dims);
    lndiseq_arr = new_float32_array(3, out_dims);
    if (min_pressure_arr == NULL || max_wind_arr == NULL || outflow_temp_arr == NULL ||
        outflow_level_arr == NULL || lnpi_arr == NULL || lneff_arr == NULL || lndiseq_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    calculate_pi_4d_diagnostics_with_missing(
        PyArray_DATA(sst_arr),
        PyArray_DATA(psl_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(min_pressure_arr),
        PyArray_DATA(max_wind_arr),
        &error_flag,
        PyArray_DATA(outflow_temp_arr),
        PyArray_DATA(outflow_level_arr),
        PyArray_DATA(lnpi_arr),
        PyArray_DATA(lneff_arr),
        PyArray_DATA(lndiseq_arr),
        &lnckcd,
        outflow_source_flag,
        ckcd_in,
        nlat,
        nlon,
        num_levels,
        num_times
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(sst_arr);
    Py_DECREF(psl_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    return build_gridded_diag_result(
        min_pressure_arr,
        max_wind_arr,
        error_flag,
        outflow_temp_arr,
        outflow_level_arr,
        lnpi_arr,
        lneff_arr,
        lndiseq_arr,
        lnckcd
    );

fail:
    Py_XDECREF(sst_arr);
    Py_XDECREF(psl_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(min_pressure_arr);
    Py_XDECREF(max_wind_arr);
    Py_XDECREF(outflow_temp_arr);
    Py_XDECREF(outflow_level_arr);
    Py_XDECREF(lnpi_arr);
    Py_XDECREF(lneff_arr);
    Py_XDECREF(lndiseq_arr);
    return NULL;
}

static PyObject *py_cape(PyObject *self, PyObject *args) {
    PyObject *temp_obj = NULL;
    PyObject *mixing_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *mixing_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    float parcel_temp = 0.0f;
    float parcel_mixing_ratio = 0.0f;
    float parcel_pressure = 0.0f;
    float buoyancy_param = 0.0f;
    float cape_value_out = 0.0f;
    float outflow_temp_out = 0.0f;
    float outflow_level_out = 0.0f;
    int error_flag = 0;
    int num_points = 0;
    int array_size = 0;

    (void) self;

    if (!PyArg_ParseTuple(
            args,
            "fffOOOif",
            &parcel_temp,
            &parcel_mixing_ratio,
            &parcel_pressure,
            &temp_obj,
            &mixing_obj,
            &pressure_obj,
            &num_points,
            &buoyancy_param
        )) {
        return NULL;
    }

    temp_arr = to_float32_1d(temp_obj);
    mixing_arr = to_float32_1d(mixing_obj);
    pressure_arr = to_float32_1d(pressure_obj);
    if (temp_arr == NULL || mixing_arr == NULL || pressure_arr == NULL) {
        goto fail;
    }

    if (require_ndim(temp_arr, 1, "temp_profile") != 0 ||
        require_ndim(mixing_arr, 1, "mixing_ratio_profile") != 0 ||
        require_ndim(pressure_arr, 1, "pressure_profile") != 0) {
        goto fail;
    }

    array_size = (int) PyArray_DIM(temp_arr, 0);
    if (PyArray_DIM(mixing_arr, 0) != array_size || PyArray_DIM(pressure_arr, 0) != array_size) {
        PyErr_SetString(PyExc_ValueError, "CAPE profile inputs must share the same length");
        goto fail;
    }
    if (num_points < 1 || num_points > array_size) {
        PyErr_SetString(PyExc_ValueError, "num_points must be between 1 and the profile length");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    cape(
        parcel_temp,
        parcel_mixing_ratio,
        parcel_pressure,
        PyArray_DATA(temp_arr),
        PyArray_DATA(mixing_arr),
        PyArray_DATA(pressure_arr),
        num_points,
        buoyancy_param,
        &cape_value_out,
        &outflow_temp_out,
        &outflow_level_out,
        &error_flag,
        array_size
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(temp_arr);
    Py_DECREF(mixing_arr);
    Py_DECREF(pressure_arr);
    return Py_BuildValue("fffi", cape_value_out, outflow_temp_out, outflow_level_out, error_flag);

fail:
    Py_XDECREF(temp_arr);
    Py_XDECREF(mixing_arr);
    Py_XDECREF(pressure_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"calculate_pi_gridded_data", py_calculate_pi_gridded_data, METH_VARARGS, "Compute PI on a 3D gridded field."},
    {"calculate_pi_gridded_with_missing", py_calculate_pi_gridded_with_missing, METH_VARARGS, "Compute PI on a 3D gridded field with missing values."},
    {"calculate_pi_gridded_diagnostics", (PyCFunction) py_calculate_pi_gridded_diagnostics, METH_VARARGS | METH_KEYWORDS, "Compute PI diagnostics on a 3D gridded field."},
    {"calculate_pi_gridded_diagnostics_with_missing", (PyCFunction) py_calculate_pi_gridded_diagnostics_with_missing, METH_VARARGS | METH_KEYWORDS, "Compute PI diagnostics on a 3D gridded field with missing values."},
    {"calculate_pi_single_profile", py_calculate_pi_single_profile, METH_VARARGS, "Compute PI for one vertical profile."},
    {"calculate_pi_profile_diagnostics", (PyCFunction) py_calculate_pi_profile_diagnostics, METH_VARARGS | METH_KEYWORDS, "Compute PI diagnostics for one vertical profile."},
    {"calculate_pi_4d_data", py_calculate_pi_4d_data, METH_VARARGS, "Compute PI on a 4D time-varying field."},
    {"calculate_pi_4d_with_missing", py_calculate_pi_4d_with_missing, METH_VARARGS, "Compute PI on a 4D time-varying field with missing values."},
    {"calculate_pi_4d_diagnostics", (PyCFunction) py_calculate_pi_4d_diagnostics, METH_VARARGS | METH_KEYWORDS, "Compute PI diagnostics on a 4D time-varying field."},
    {"calculate_pi_4d_diagnostics_with_missing", (PyCFunction) py_calculate_pi_4d_diagnostics_with_missing, METH_VARARGS | METH_KEYWORDS, "Compute PI diagnostics on a 4D time-varying field with missing values."},
    {"cape", py_cape, METH_VARARGS, "Compute CAPE for one parcel/environment profile pair."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tropical_cyclone_potential_intensity",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_tropical_cyclone_potential_intensity(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
