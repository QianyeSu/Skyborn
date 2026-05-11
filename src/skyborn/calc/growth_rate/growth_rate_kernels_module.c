#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void dbarot_growth_rate_1d_c(
    void *lat,
    void *u,
    void *max_growth,
    void *ier,
    int nlat
);
void dbaroc_growth_rate_1d_c(
    void *u,
    void *theta,
    void *pressure,
    void *temperature,
    double f_cor,
    double beta,
    int smooth_window,
    int wavenumber_mode,
    int wavenumber_count,
    double zonal_length,
    void *max_growth,
    void *ier,
    int nlev
);
void dbaroc_growth_rate_profiles_c(
    void *u_input,
    void *temperature_input,
    void *source_pressure,
    void *target_pressure,
    void *f_cor,
    void *beta,
    int interp_kind,
    int smooth_window,
    int wavenumber_mode,
    int wavenumber_count,
    void *zonal_length,
    void *growth,
    void *ier,
    int nlev_in,
    int nprofile,
    int nlev_out
);

static PyArrayObject *to_double_1d(PyObject *obj, int contig_flag) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | contig_flag
    );
}

static PyObject *py_dbarot_growth_rate_1d(PyObject *self, PyObject *args) {
    PyObject *lat_obj = NULL;
    PyObject *u_obj = NULL;
    PyArrayObject *lat_arr = NULL;
    PyArrayObject *u_arr = NULL;
    double max_growth = 0.0;
    int ier = 0;
    int nlat = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OO", &lat_obj, &u_obj)) {
        return NULL;
    }

    lat_arr = to_double_1d(lat_obj, NPY_ARRAY_C_CONTIGUOUS);
    u_arr = to_double_1d(u_obj, NPY_ARRAY_C_CONTIGUOUS);
    if (lat_arr == NULL || u_arr == NULL) {
        Py_XDECREF(lat_arr);
        Py_XDECREF(u_arr);
        return NULL;
    }

    if (PyArray_NDIM(lat_arr) != 1 || PyArray_NDIM(u_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "lat and u must be 1D arrays");
        goto fail;
    }
    nlat = (int) PyArray_DIM(lat_arr, 0);
    if ((int) PyArray_DIM(u_arr, 0) != nlat) {
        PyErr_SetString(PyExc_ValueError, "lat and u must have the same length");
        goto fail;
    }

    dbarot_growth_rate_1d_c(
        PyArray_DATA(lat_arr),
        PyArray_DATA(u_arr),
        &max_growth,
        &ier,
        nlat
    );

    Py_DECREF(lat_arr);
    Py_DECREF(u_arr);
    return Py_BuildValue("di", max_growth, ier);

fail:
    Py_XDECREF(lat_arr);
    Py_XDECREF(u_arr);
    return NULL;
}

static PyObject *py_dbaroc_growth_rate_1d(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *theta_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *theta_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    double f_cor = 0.0;
    double beta = 0.0;
    double zonal_length = 0.0;
    double max_growth = 0.0;
    int smooth_window = 0;
    int wavenumber_mode = 0;
    int wavenumber_count = 0;
    int ier = 0;
    int nlev = 0;

    (void) self;

    if (!PyArg_ParseTuple(
            args,
            "OOOOddiiid",
            &u_obj,
            &theta_obj,
            &pressure_obj,
            &temperature_obj,
            &f_cor,
            &beta,
            &smooth_window,
            &wavenumber_mode,
            &wavenumber_count,
            &zonal_length
        )) {
        return NULL;
    }

    u_arr = to_double_1d(u_obj, NPY_ARRAY_C_CONTIGUOUS);
    theta_arr = to_double_1d(theta_obj, NPY_ARRAY_C_CONTIGUOUS);
    pressure_arr = to_double_1d(pressure_obj, NPY_ARRAY_C_CONTIGUOUS);
    temperature_arr = to_double_1d(temperature_obj, NPY_ARRAY_C_CONTIGUOUS);
    if (u_arr == NULL || theta_arr == NULL || pressure_arr == NULL || temperature_arr == NULL) {
        Py_XDECREF(u_arr);
        Py_XDECREF(theta_arr);
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        return NULL;
    }

    nlev = (int) PyArray_DIM(u_arr, 0);
    if ((int) PyArray_DIM(theta_arr, 0) != nlev ||
        (int) PyArray_DIM(pressure_arr, 0) != nlev ||
        (int) PyArray_DIM(temperature_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "all 1D inputs must have the same length");
        goto fail;
    }

    dbaroc_growth_rate_1d_c(
        PyArray_DATA(u_arr),
        PyArray_DATA(theta_arr),
        PyArray_DATA(pressure_arr),
        PyArray_DATA(temperature_arr),
        f_cor,
        beta,
        smooth_window,
        wavenumber_mode,
        wavenumber_count,
        zonal_length,
        &max_growth,
        &ier,
        nlev
    );

    Py_DECREF(u_arr);
    Py_DECREF(theta_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    return Py_BuildValue("di", max_growth, ier);

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(theta_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    return NULL;
}

static PyObject *py_dbaroc_growth_rate_profiles(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *source_pressure_obj = NULL;
    PyObject *target_pressure_obj = NULL;
    PyObject *f_cor_obj = NULL;
    PyObject *beta_obj = NULL;
    PyObject *zonal_length_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *source_pressure_arr = NULL;
    PyArrayObject *target_pressure_arr = NULL;
    PyArrayObject *f_cor_arr = NULL;
    PyArrayObject *beta_arr = NULL;
    PyArrayObject *zonal_length_arr = NULL;
    PyArrayObject *growth_arr = NULL;
    PyArrayObject *ier_arr = NULL;
    int interp_kind = 0;
    int smooth_window = 0;
    int wavenumber_mode = 0;
    int wavenumber_count = 0;
    int nlev_in = 0;
    int nprofile = 0;
    int nlev_out = 0;

    (void) self;

    if (!PyArg_ParseTuple(
            args,
            "OOOOOOiiiiO",
            &u_obj,
            &temperature_obj,
            &source_pressure_obj,
            &target_pressure_obj,
            &f_cor_obj,
            &beta_obj,
            &interp_kind,
            &smooth_window,
            &wavenumber_mode,
            &wavenumber_count,
            &zonal_length_obj
        )) {
        return NULL;
    }

    u_arr = (PyArrayObject *) PyArray_FROM_OTF(u_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    source_pressure_arr = to_double_1d(source_pressure_obj, NPY_ARRAY_C_CONTIGUOUS);
    target_pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        target_pressure_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    f_cor_arr = to_double_1d(f_cor_obj, NPY_ARRAY_C_CONTIGUOUS);
    beta_arr = to_double_1d(beta_obj, NPY_ARRAY_C_CONTIGUOUS);
    zonal_length_arr = to_double_1d(zonal_length_obj, NPY_ARRAY_C_CONTIGUOUS);
    if (u_arr == NULL || temperature_arr == NULL || source_pressure_arr == NULL ||
        target_pressure_arr == NULL || f_cor_arr == NULL || beta_arr == NULL ||
        zonal_length_arr == NULL) {
        Py_XDECREF(u_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(source_pressure_arr);
        Py_XDECREF(target_pressure_arr);
        Py_XDECREF(f_cor_arr);
        Py_XDECREF(beta_arr);
        Py_XDECREF(zonal_length_arr);
        return NULL;
    }

    if (PyArray_NDIM(u_arr) != 2 || PyArray_NDIM(temperature_arr) != 2 ||
        PyArray_NDIM(target_pressure_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "profile inputs must be 2D arrays");
        goto fail;
    }

    nlev_in = (int) PyArray_DIM(u_arr, 0);
    nprofile = (int) PyArray_DIM(u_arr, 1);
    nlev_out = (int) PyArray_DIM(target_pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev_in ||
        (int) PyArray_DIM(target_pressure_arr, 1) != nprofile ||
        (int) PyArray_DIM(temperature_arr, 1) != nprofile) {
        PyErr_SetString(PyExc_ValueError, "all batched inputs must share the profile dimension");
        goto fail;
    }

    growth_arr = (PyArrayObject *) PyArray_EMPTY(1, (npy_intp[]) {nprofile}, NPY_FLOAT64, 1);
    ier_arr = (PyArrayObject *) PyArray_EMPTY(1, (npy_intp[]) {nprofile}, NPY_INT32, 1);
    if (growth_arr == NULL || ier_arr == NULL) {
        Py_XDECREF(growth_arr);
        Py_XDECREF(ier_arr);
        goto fail;
    }

    dbaroc_growth_rate_profiles_c(
        PyArray_DATA(u_arr),
        PyArray_DATA(temperature_arr),
        PyArray_DATA(source_pressure_arr),
        PyArray_DATA(target_pressure_arr),
        PyArray_DATA(f_cor_arr),
        PyArray_DATA(beta_arr),
        interp_kind,
        smooth_window,
        wavenumber_mode,
        wavenumber_count,
        PyArray_DATA(zonal_length_arr),
        PyArray_DATA(growth_arr),
        PyArray_DATA(ier_arr),
        nlev_in,
        nprofile,
        nlev_out
    );

    Py_DECREF(u_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(source_pressure_arr);
    Py_DECREF(target_pressure_arr);
    Py_DECREF(f_cor_arr);
    Py_DECREF(beta_arr);
    Py_DECREF(zonal_length_arr);
    {
        PyObject *result = PyTuple_Pack(2, (PyObject *) growth_arr, (PyObject *) ier_arr);
        Py_DECREF(growth_arr);
        Py_DECREF(ier_arr);
        return result;
    }

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(source_pressure_arr);
    Py_XDECREF(target_pressure_arr);
    Py_XDECREF(f_cor_arr);
    Py_XDECREF(beta_arr);
    Py_XDECREF(zonal_length_arr);
    Py_XDECREF(growth_arr);
    Py_XDECREF(ier_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {
        "dbarot_growth_rate_1d",
        py_dbarot_growth_rate_1d,
        METH_VARARGS,
        "Compute the barotropic growth rate and return (max_growth, ier)."
    },
    {
        "dbaroc_growth_rate_1d",
        py_dbaroc_growth_rate_1d,
        METH_VARARGS,
        "Compute the baroclinic growth rate and return (max_growth, ier)."
    },
    {
        "dbaroc_growth_rate_profiles",
        py_dbaroc_growth_rate_profiles,
        METH_VARARGS,
        "Compute batched baroclinic growth rates."
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "growth_rate_kernels",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_growth_rate_kernels(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
