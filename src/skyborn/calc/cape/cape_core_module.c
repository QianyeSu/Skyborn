#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

/* Fortran kernels (bind(C)) */
void cape_profile_c(
    int nlev, void *pressure_hpa, void *temperature_c, void *dewpoint_c,
    void *cape, void *cin, void *lfc_p, void *el_p);
void cape_grid_c(
    int nlev, int nlat, int nlon,
    void *pressure_3d, void *t_3d, void *td_3d,
    void *cape2, void *cin2, void *lfc2, void *el2);
void parcel_profile_c(
    int nlev, void *pressure_hpa, void *temperature_c, void *dewpoint_c, void *t_par_c);
void parcel_profile_grid_c(
    int nlev, int nlat, int nlon,
    void *pressure_3d, void *t_3d, void *td_3d, void *out_3d);
void most_unstable_parcel_c(
    int nlev, void *pressure_hpa, void *temperature_c, void *dewpoint_c,
    double depth_hpa, void *mup_p, void *mup_t, void *mup_td, int *mup_idx);
void most_unstable_parcel_grid_c(
    int nlev, int nlat, int nlon,
    void *pressure_3d, void *t_3d, void *td_3d, double depth_hpa,
    void *out_p3, void *out_t3, void *out_td3, int *out_idx3);
void mucape_c(
    int nlev, void *pressure_hpa, void *temperature_c, void *dewpoint_c,
    double depth_hpa, void *cape, void *cin, void *lfc_p, void *el_p);
void mucape_grid_c(
    int nlev, int nlat, int nlon,
    void *pressure_3d, void *t_3d, void *td_3d, double depth_hpa,
    void *cape2, void *cin2, void *lfc2, void *el2);

static int check_3d_shapes(PyArrayObject *p, PyArrayObject *t, PyArrayObject *td)
{
    npy_intp nlev, nlat, nlon;
    nlev = PyArray_DIM(p, 0);
    nlat = PyArray_DIM(p, 1);
    nlon = PyArray_DIM(p, 2);
    if (PyArray_DIM(t, 0) != nlev || PyArray_DIM(t, 1) != nlat || PyArray_DIM(t, 2) != nlon ||
        PyArray_DIM(td, 0) != nlev || PyArray_DIM(td, 1) != nlat || PyArray_DIM(td, 2) != nlon)
    {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a shape");
        return 0;
    }
    return 1;
}

/* ---------------- cape_profile: single vertical profile ---------------- */
static PyObject *py_cape_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    double cape = 0.0, cin = 0.0, lfc_p = 0.0, el_p = 0.0;
    int nlev;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 1 || PyArray_NDIM(temperature_arr) != 1 ||
        PyArray_NDIM(dewpoint_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 1D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev)
    {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a length");
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    cape_profile_c(nlev, PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
                   PyArray_DATA(dewpoint_arr), &cape, &cin, &lfc_p, &el_p);
    Py_END_ALLOW_THREADS
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return Py_BuildValue("dddd", cape, cin, lfc_p, el_p);
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    return NULL;
}

/* ---------------- cape_grid: 3D grid of profiles ---------------- */
static PyObject *py_cape_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *cape_arr = NULL;
    PyArrayObject *cin_arr = NULL;
    PyArrayObject *lfc_arr = NULL;
    PyArrayObject *el_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[2];
    int nlev, nlat, nlon;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 3 || PyArray_NDIM(temperature_arr) != 3 ||
        PyArray_NDIM(dewpoint_arr) != 3)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 3D arrays");
        goto fail;
    }
    if (!check_3d_shapes(pressure_arr, temperature_arr, dewpoint_arr))
    {
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    nlat = (int) PyArray_DIM(pressure_arr, 1);
    nlon = (int) PyArray_DIM(pressure_arr, 2);
    dims[0] = nlat;
    dims[1] = nlon;
    cape_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    cin_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    lfc_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    el_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (cape_arr == NULL || cin_arr == NULL || lfc_arr == NULL || el_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    cape_grid_c(nlev, nlat, nlon,
                PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
                PyArray_DATA(dewpoint_arr), PyArray_DATA(cape_arr),
                PyArray_DATA(cin_arr), PyArray_DATA(lfc_arr), PyArray_DATA(el_arr));
    Py_END_ALLOW_THREADS
    result = Py_BuildValue("NNNN", (PyObject *) cape_arr, (PyObject *) cin_arr,
                           (PyObject *) lfc_arr, (PyObject *) el_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return result;
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    Py_XDECREF(cape_arr);
    Py_XDECREF(cin_arr);
    Py_XDECREF(lfc_arr);
    Py_XDECREF(el_arr);
    return NULL;
}

/* ---------------- parcel_profile: single profile ---------------- */
static PyObject *py_parcel_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *out_arr = NULL;
    npy_intp dims[1];
    int nlev;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 1 || PyArray_NDIM(temperature_arr) != 1 ||
        PyArray_NDIM(dewpoint_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 1D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev)
    {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a length");
        goto fail;
    }
    dims[0] = nlev;
    out_arr = (PyArrayObject *) PyArray_EMPTY(1, dims, NPY_FLOAT64, 0);
    if (out_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    parcel_profile_c(nlev, PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
                     PyArray_DATA(dewpoint_arr), PyArray_DATA(out_arr));
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

/* ---------------- parcel_profile_grid: 3D grid ---------------- */
static PyObject *py_parcel_profile_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *out_arr = NULL;
    npy_intp dims[3];
    int nlev, nlat, nlon;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 3 || PyArray_NDIM(temperature_arr) != 3 ||
        PyArray_NDIM(dewpoint_arr) != 3)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 3D arrays");
        goto fail;
    }
    if (!check_3d_shapes(pressure_arr, temperature_arr, dewpoint_arr))
    {
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    nlat = (int) PyArray_DIM(pressure_arr, 1);
    nlon = (int) PyArray_DIM(pressure_arr, 2);
    dims[0] = nlev;
    dims[1] = nlat;
    dims[2] = nlon;
    out_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    if (out_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    parcel_profile_grid_c(nlev, nlat, nlon, PyArray_DATA(pressure_arr),
                          PyArray_DATA(temperature_arr), PyArray_DATA(dewpoint_arr),
                          PyArray_DATA(out_arr));
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

/* ---------------- most_unstable_parcel: single profile ---------------- */
static PyObject *py_most_unstable_parcel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    double depth = 300.0, mup_p = 0.0, mup_t = 0.0, mup_td = 0.0;
    int mup_idx = -1, nlev;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", "depth", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj, &depth))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 1 || PyArray_NDIM(temperature_arr) != 1 ||
        PyArray_NDIM(dewpoint_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 1D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev)
    {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a length");
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    most_unstable_parcel_c(nlev, PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
                           PyArray_DATA(dewpoint_arr), depth, &mup_p, &mup_t, &mup_td, &mup_idx);
    Py_END_ALLOW_THREADS
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return Py_BuildValue("dddi", mup_p, mup_t, mup_td, mup_idx);
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    return NULL;
}

/* ---------------- most_unstable_parcel: 3D grid ---------------- */
static PyObject *py_most_unstable_parcel_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *p_arr = NULL, *t_arr = NULL, *td_arr = NULL, *idx_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[3];
    double depth = 300.0;
    int nlev, nlat, nlon;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", "depth", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj, &depth))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 3 || PyArray_NDIM(temperature_arr) != 3 ||
        PyArray_NDIM(dewpoint_arr) != 3)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 3D arrays");
        goto fail;
    }
    if (!check_3d_shapes(pressure_arr, temperature_arr, dewpoint_arr))
    {
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    nlat = (int) PyArray_DIM(pressure_arr, 1);
    nlon = (int) PyArray_DIM(pressure_arr, 2);
    dims[0] = nlev;
    dims[1] = nlat;
    dims[2] = nlon;
    p_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    t_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    td_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    idx_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_INT32, 1);
    if (p_arr == NULL || t_arr == NULL || td_arr == NULL || idx_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    most_unstable_parcel_grid_c(nlev, nlat, nlon, PyArray_DATA(pressure_arr),
                                PyArray_DATA(temperature_arr), PyArray_DATA(dewpoint_arr),
                                depth, PyArray_DATA(p_arr), PyArray_DATA(t_arr),
                                PyArray_DATA(td_arr), (int *) PyArray_DATA(idx_arr));
    Py_END_ALLOW_THREADS
    result = Py_BuildValue("NNNN", (PyObject *) p_arr, (PyObject *) t_arr,
                           (PyObject *) td_arr, (PyObject *) idx_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return result;
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    Py_XDECREF(p_arr);
    Py_XDECREF(t_arr);
    Py_XDECREF(td_arr);
    Py_XDECREF(idx_arr);
    return NULL;
}

/* ---------------- most unstable CAPE/CIN: single profile ---------------- */
static PyObject *py_mucape_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    double depth = 300.0, cape = 0.0, cin = 0.0, lfc_p = 0.0, el_p = 0.0;
    int nlev;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", "depth", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj, &depth))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 1 || PyArray_NDIM(temperature_arr) != 1 ||
        PyArray_NDIM(dewpoint_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 1D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    if ((int) PyArray_DIM(temperature_arr, 0) != nlev ||
        (int) PyArray_DIM(dewpoint_arr, 0) != nlev)
    {
        PyErr_SetString(PyExc_ValueError, "pressure, temperature, and dewpoint must share a length");
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    mucape_c(nlev, PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
             PyArray_DATA(dewpoint_arr), depth, &cape, &cin, &lfc_p, &el_p);
    Py_END_ALLOW_THREADS
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return Py_BuildValue("dddd", cape, cin, lfc_p, el_p);
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    return NULL;
}

/* ---------------- most unstable CAPE/CIN: 3D grid ---------------- */
static PyObject *py_mucape_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *pressure_obj = NULL;
    PyObject *temperature_obj = NULL;
    PyObject *dewpoint_obj = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *temperature_arr = NULL;
    PyArrayObject *dewpoint_arr = NULL;
    PyArrayObject *cape_arr = NULL, *cin_arr = NULL, *lfc_arr = NULL, *el_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[2];
    double depth = 300.0;
    int nlev, nlat, nlon;
    static char *kwlist[] = {"pressure", "temperature", "dewpoint", "depth", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|d", kwlist,
                                     &pressure_obj, &temperature_obj, &dewpoint_obj, &depth))
    {
        return NULL;
    }
    pressure_arr = (PyArrayObject *) PyArray_FROM_OTF(
        pressure_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    temperature_arr = (PyArrayObject *) PyArray_FROM_OTF(
        temperature_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    dewpoint_arr = (PyArrayObject *) PyArray_FROM_OTF(
        dewpoint_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (pressure_arr == NULL || temperature_arr == NULL || dewpoint_arr == NULL)
    {
        Py_XDECREF(pressure_arr);
        Py_XDECREF(temperature_arr);
        Py_XDECREF(dewpoint_arr);
        return NULL;
    }
    if (PyArray_NDIM(pressure_arr) != 3 || PyArray_NDIM(temperature_arr) != 3 ||
        PyArray_NDIM(dewpoint_arr) != 3)
    {
        PyErr_SetString(PyExc_ValueError,
                        "pressure, temperature, and dewpoint must be 3D arrays");
        goto fail;
    }
    if (!check_3d_shapes(pressure_arr, temperature_arr, dewpoint_arr))
    {
        goto fail;
    }
    nlev = (int) PyArray_DIM(pressure_arr, 0);
    nlat = (int) PyArray_DIM(pressure_arr, 1);
    nlon = (int) PyArray_DIM(pressure_arr, 2);
    dims[0] = nlat;
    dims[1] = nlon;
    cape_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    cin_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    lfc_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    el_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (cape_arr == NULL || cin_arr == NULL || lfc_arr == NULL || el_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    mucape_grid_c(nlev, nlat, nlon, PyArray_DATA(pressure_arr), PyArray_DATA(temperature_arr),
                  PyArray_DATA(dewpoint_arr), depth, PyArray_DATA(cape_arr),
                  PyArray_DATA(cin_arr), PyArray_DATA(lfc_arr), PyArray_DATA(el_arr));
    Py_END_ALLOW_THREADS
    result = Py_BuildValue("NNNN", (PyObject *) cape_arr, (PyObject *) cin_arr,
                           (PyObject *) lfc_arr, (PyObject *) el_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(temperature_arr);
    Py_DECREF(dewpoint_arr);
    return result;
fail:
    Py_XDECREF(pressure_arr);
    Py_XDECREF(temperature_arr);
    Py_XDECREF(dewpoint_arr);
    Py_XDECREF(cape_arr);
    Py_XDECREF(cin_arr);
    Py_XDECREF(lfc_arr);
    Py_XDECREF(el_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"cape_profile", (PyCFunction) (void (*)(void)) py_cape_profile,
     METH_VARARGS | METH_KEYWORDS,
     "Compute CAPE, CIN, LFC pressure, and EL pressure for a single profile."},
    {"cape_grid", (PyCFunction) (void (*)(void)) py_cape_grid,
     METH_VARARGS | METH_KEYWORDS,
     "Compute CAPE, CIN, LFC pressure, and EL pressure for a 3D grid of profiles."},
    {"parcel_profile", (PyCFunction) (void (*)(void)) py_parcel_profile,
     METH_VARARGS | METH_KEYWORDS,
     "Compute the surface-based parcel temperature profile."},
    {"parcel_profile_grid", (PyCFunction) (void (*)(void)) py_parcel_profile_grid,
     METH_VARARGS | METH_KEYWORDS,
     "Compute parcel temperature profiles for a 3D grid."},
    {"most_unstable_parcel", (PyCFunction) (void (*)(void)) py_most_unstable_parcel,
     METH_VARARGS | METH_KEYWORDS,
     "Find the most unstable parcel (max theta-e) within the bottom-depth layer."},
    {"most_unstable_parcel_grid", (PyCFunction) (void (*)(void)) py_most_unstable_parcel_grid,
     METH_VARARGS | METH_KEYWORDS,
     "Find the most unstable parcel for every column of a 3D grid."},
    {"mucape_profile", (PyCFunction) (void (*)(void)) py_mucape_profile,
     METH_VARARGS | METH_KEYWORDS,
     "Compute most unstable CAPE, CIN, LFC pressure, and EL pressure."},
    {"mucape_grid", (PyCFunction) (void (*)(void)) py_mucape_grid,
     METH_VARARGS | METH_KEYWORDS,
     "Compute most unstable CAPE/CIN for a 3D grid of profiles."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cape_core",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_cape_core(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
