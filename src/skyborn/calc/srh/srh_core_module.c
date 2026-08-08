#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

void srh_profile_c(
    int nlev, void *height_m, void *u_ms, void *v_ms,
    double depth_m, double bottom_m, double storm_u_ms, double storm_v_ms,
    void *srh_pos, void *srh_neg, void *srh_total);
void srh_grid_c(
    int nlev, int nlat, int nlon,
    void *height_3d, void *u_3d, void *v_3d,
    double depth_m, double bottom_m, double storm_u_ms, double storm_v_ms,
    void *pos2, void *neg2, void *tot2);

/* ---------------- srh_profile: single vertical profile ---------------- */
static PyObject *py_srh_profile(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *height_obj = NULL;
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *height_arr = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    double depth = 0.0, bottom = 0.0, storm_u = 0.0, storm_v = 0.0;
    double srh_pos = 0.0, srh_neg = 0.0, srh_total = 0.0;
    int nlev;
    static char *kwlist[] = {"height", "u", "v", "depth", "bottom",
                             "storm_u", "storm_v", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOd|ddd", kwlist,
                                     &height_obj, &u_obj, &v_obj, &depth,
                                     &bottom, &storm_u, &storm_v))
    {
        return NULL;
    }
    height_arr = (PyArrayObject *) PyArray_FROM_OTF(
        height_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    u_arr = (PyArrayObject *) PyArray_FROM_OTF(
        u_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    v_arr = (PyArrayObject *) PyArray_FROM_OTF(
        v_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (height_arr == NULL || u_arr == NULL || v_arr == NULL)
    {
        Py_XDECREF(height_arr);
        Py_XDECREF(u_arr);
        Py_XDECREF(v_arr);
        return NULL;
    }
    if (PyArray_NDIM(height_arr) != 1 || PyArray_NDIM(u_arr) != 1 ||
        PyArray_NDIM(v_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError,
                        "height, u, and v must be 1D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(height_arr, 0);
    if ((int) PyArray_DIM(u_arr, 0) != nlev || (int) PyArray_DIM(v_arr, 0) != nlev)
    {
        PyErr_SetString(PyExc_ValueError, "height, u, and v must share a length");
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    srh_profile_c(nlev, PyArray_DATA(height_arr), PyArray_DATA(u_arr),
                  PyArray_DATA(v_arr), depth, bottom, storm_u, storm_v,
                  &srh_pos, &srh_neg, &srh_total);
    Py_END_ALLOW_THREADS
    Py_DECREF(height_arr);
    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    return Py_BuildValue("ddd", srh_pos, srh_neg, srh_total);
fail:
    Py_XDECREF(height_arr);
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    return NULL;
}

/* ---------------- srh_grid: 3D grid of profiles ---------------- */
static PyObject *py_srh_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *height_obj = NULL;
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *height_arr = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *pos_arr = NULL, *neg_arr = NULL, *tot_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[2];
    double depth = 0.0, bottom = 0.0, storm_u = 0.0, storm_v = 0.0;
    int nlev, nlat, nlon;
    static char *kwlist[] = {"height", "u", "v", "depth", "bottom",
                             "storm_u", "storm_v", NULL};
    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOd|ddd", kwlist,
                                     &height_obj, &u_obj, &v_obj, &depth,
                                     &bottom, &storm_u, &storm_v))
    {
        return NULL;
    }
    height_arr = (PyArrayObject *) PyArray_FROM_OTF(
        height_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    u_arr = (PyArrayObject *) PyArray_FROM_OTF(
        u_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    v_arr = (PyArrayObject *) PyArray_FROM_OTF(
        v_obj, NPY_FLOAT64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (height_arr == NULL || u_arr == NULL || v_arr == NULL)
    {
        Py_XDECREF(height_arr);
        Py_XDECREF(u_arr);
        Py_XDECREF(v_arr);
        return NULL;
    }
    if (PyArray_NDIM(height_arr) != 3 || PyArray_NDIM(u_arr) != 3 ||
        PyArray_NDIM(v_arr) != 3)
    {
        PyErr_SetString(PyExc_ValueError,
                        "height, u, and v must be 3D arrays");
        goto fail;
    }
    nlev = (int) PyArray_DIM(height_arr, 0);
    nlat = (int) PyArray_DIM(height_arr, 1);
    nlon = (int) PyArray_DIM(height_arr, 2);
    if (PyArray_DIM(u_arr, 0) != nlev || PyArray_DIM(u_arr, 1) != nlat ||
        PyArray_DIM(u_arr, 2) != nlon || PyArray_DIM(v_arr, 0) != nlev ||
        PyArray_DIM(v_arr, 1) != nlat || PyArray_DIM(v_arr, 2) != nlon)
    {
        PyErr_SetString(PyExc_ValueError, "height, u, and v must share a shape");
        goto fail;
    }
    dims[0] = nlat;
    dims[1] = nlon;
    pos_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    neg_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    tot_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (pos_arr == NULL || neg_arr == NULL || tot_arr == NULL)
    {
        goto fail;
    }
    Py_BEGIN_ALLOW_THREADS
    srh_grid_c(nlev, nlat, nlon, PyArray_DATA(height_arr), PyArray_DATA(u_arr),
               PyArray_DATA(v_arr), depth, bottom, storm_u, storm_v,
               PyArray_DATA(pos_arr), PyArray_DATA(neg_arr), PyArray_DATA(tot_arr));
    Py_END_ALLOW_THREADS
    result = Py_BuildValue("NNN", (PyObject *) pos_arr, (PyObject *) neg_arr,
                           (PyObject *) tot_arr);
    Py_DECREF(height_arr);
    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    return result;
fail:
    Py_XDECREF(height_arr);
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(pos_arr);
    Py_XDECREF(neg_arr);
    Py_XDECREF(tot_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"srh_profile", (PyCFunction) (void (*)(void)) py_srh_profile,
     METH_VARARGS | METH_KEYWORDS,
     "Compute (positive, negative, total) storm-relative helicity for one profile."},
    {"srh_grid", (PyCFunction) (void (*)(void)) py_srh_grid,
     METH_VARARGS | METH_KEYWORDS,
     "Compute storm-relative helicity for a 3D grid of profiles."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "srh_core",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_srh_core(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
