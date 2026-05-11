#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void z2geouv_c(
    void *z,
    int nlat,
    int mlon,
    double zmsg,
    void *glon,
    void *glat,
    void *ug,
    void *vg,
    int iopt
);
void z2geouv_3d_c(
    void *z,
    int nlat,
    int mlon,
    int n3rd,
    double zmsg,
    void *glon,
    void *glat,
    void *ug,
    void *vg,
    int iopt
);
void z2geouv_4d_c(
    void *z,
    int nlat,
    int mlon,
    int n3rd,
    int n4th,
    double zmsg,
    void *glon,
    void *glat,
    void *ug,
    void *vg,
    int iopt
);
void zuvnew_c(
    void *z,
    int nlat,
    int mlon,
    double zmsg,
    void *glon,
    void *glat,
    void *ug,
    void *vg,
    int iopt
);
void z2guv_c(
    void *z,
    int nlat,
    int mlon,
    double zmsg,
    void *glon,
    void *glat,
    void *ug,
    void *vg,
    int iopt
);

typedef struct {
    PyArrayObject *z_arr;
    PyArrayObject *glon_arr;
    PyArrayObject *glat_arr;
} GeostrophicInputs;

typedef void (*geostrophic_2d_fn)(
    void *,
    int,
    int,
    double,
    void *,
    void *,
    void *,
    void *,
    int
);

typedef void (*geostrophic_3d_fn)(
    void *,
    int,
    int,
    int,
    double,
    void *,
    void *,
    void *,
    void *,
    int
);

typedef void (*geostrophic_4d_fn)(
    void *,
    int,
    int,
    int,
    int,
    double,
    void *,
    void *,
    void *,
    void *,
    int
);

static int get_optional_int(PyObject *kwargs, const char *name, int *value) {
    PyObject *item = NULL;

    if (kwargs == NULL) {
        return 0;
    }

    item = PyDict_GetItemString(kwargs, name);
    if (item == NULL) {
        return 0;
    }

    *value = (int) PyLong_AsLong(item);
    if (PyErr_Occurred()) {
        return -1;
    }

    return 1;
}

static int extract_common_keywords(
    PyObject *kwargs,
    double *zmsg,
    PyObject **glon_obj,
    PyObject **glat_obj,
    int *iopt
) {
    PyObject *item = NULL;

    if (kwargs == NULL) {
        PyErr_SetString(
            PyExc_TypeError,
            "keyword arguments zmsg, glon, glat, and iopt are required"
        );
        return -1;
    }

    item = PyDict_GetItemString(kwargs, "zmsg");
    if (item == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing required keyword argument: 'zmsg'");
        return -1;
    }
    *zmsg = PyFloat_AsDouble(item);
    if (PyErr_Occurred()) {
        return -1;
    }

    *glon_obj = PyDict_GetItemString(kwargs, "glon");
    if (*glon_obj == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing required keyword argument: 'glon'");
        return -1;
    }

    *glat_obj = PyDict_GetItemString(kwargs, "glat");
    if (*glat_obj == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing required keyword argument: 'glat'");
        return -1;
    }

    item = PyDict_GetItemString(kwargs, "iopt");
    if (item == NULL) {
        PyErr_SetString(PyExc_TypeError, "missing required keyword argument: 'iopt'");
        return -1;
    }
    *iopt = (int) PyLong_AsLong(item);
    if (PyErr_Occurred()) {
        return -1;
    }

    return 0;
}

static int parse_2d_call(
    PyObject *args,
    PyObject *kwargs,
    PyObject **z_obj,
    int *nlat,
    int *mlon,
    double *zmsg,
    PyObject **glon_obj,
    PyObject **glat_obj,
    int *iopt
) {
    Py_ssize_t argc = PyTuple_GET_SIZE(args);

    if (argc == 5) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 1));
        *glon_obj = PyTuple_GET_ITEM(args, 2);
        *glat_obj = PyTuple_GET_ITEM(args, 3);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 4));
        if (PyErr_Occurred()) {
            return -1;
        }
        *nlat = -1;
        *mlon = -1;
        return 0;
    }

    if (argc == 7) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 1));
        *mlon = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 2));
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 3));
        *glon_obj = PyTuple_GET_ITEM(args, 4);
        *glat_obj = PyTuple_GET_ITEM(args, 5);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 6));
        if (PyErr_Occurred()) {
            return -1;
        }
        return 0;
    }

    if (argc == 1) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = -1;
        *mlon = -1;
        return extract_common_keywords(kwargs, zmsg, glon_obj, glat_obj, iopt);
    }

    PyErr_SetString(
        PyExc_TypeError,
        "expected arguments (z, zmsg, glon, glat, iopt) or "
        "(z, nlat, mlon, zmsg, glon, glat, iopt)"
    );
    return -1;
}

static int parse_3d_call(
    PyObject *args,
    PyObject *kwargs,
    PyObject **z_obj,
    int *nlat,
    int *mlon,
    int *n3rd,
    double *zmsg,
    PyObject **glon_obj,
    PyObject **glat_obj,
    int *iopt
) {
    Py_ssize_t argc = PyTuple_GET_SIZE(args);

    if (argc == 5) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 1));
        *glon_obj = PyTuple_GET_ITEM(args, 2);
        *glat_obj = PyTuple_GET_ITEM(args, 3);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 4));
        if (PyErr_Occurred()) {
            return -1;
        }
        *nlat = -1;
        *mlon = -1;
        *n3rd = -1;
        return 0;
    }

    if (argc == 8) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 1));
        *mlon = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 2));
        *n3rd = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 3));
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 4));
        *glon_obj = PyTuple_GET_ITEM(args, 5);
        *glat_obj = PyTuple_GET_ITEM(args, 6);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 7));
        if (PyErr_Occurred()) {
            return -1;
        }
        return 0;
    }

    if (argc == 1) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = -1;
        *mlon = -1;
        *n3rd = -1;
        return extract_common_keywords(kwargs, zmsg, glon_obj, glat_obj, iopt);
    }

    PyErr_SetString(
        PyExc_TypeError,
        "expected arguments (z, zmsg, glon, glat, iopt) or "
        "(z, nlat, mlon, n3rd, zmsg, glon, glat, iopt)"
    );
    return -1;
}

static int parse_4d_call(
    PyObject *args,
    PyObject *kwargs,
    PyObject **z_obj,
    int *nlat,
    int *mlon,
    int *n3rd,
    int *n4th,
    double *zmsg,
    PyObject **glon_obj,
    PyObject **glat_obj,
    int *iopt
) {
    Py_ssize_t argc = PyTuple_GET_SIZE(args);

    if (argc == 5) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 1));
        *glon_obj = PyTuple_GET_ITEM(args, 2);
        *glat_obj = PyTuple_GET_ITEM(args, 3);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 4));
        if (PyErr_Occurred()) {
            return -1;
        }
        *nlat = -1;
        *mlon = -1;
        *n3rd = -1;
        *n4th = -1;
        return 0;
    }

    if (argc == 9) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 1));
        *mlon = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 2));
        *n3rd = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 3));
        *n4th = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 4));
        *zmsg = PyFloat_AsDouble(PyTuple_GET_ITEM(args, 5));
        *glon_obj = PyTuple_GET_ITEM(args, 6);
        *glat_obj = PyTuple_GET_ITEM(args, 7);
        *iopt = (int) PyLong_AsLong(PyTuple_GET_ITEM(args, 8));
        if (PyErr_Occurred()) {
            return -1;
        }
        return 0;
    }

    if (argc == 1) {
        *z_obj = PyTuple_GET_ITEM(args, 0);
        *nlat = -1;
        *mlon = -1;
        *n3rd = -1;
        *n4th = -1;
        return extract_common_keywords(kwargs, zmsg, glon_obj, glat_obj, iopt);
    }

    PyErr_SetString(
        PyExc_TypeError,
        "expected arguments (z, zmsg, glon, glat, iopt) or "
        "(z, nlat, mlon, n3rd, n4th, zmsg, glon, glat, iopt)"
    );
    return -1;
}

static int prepare_common_inputs(
    PyObject *z_obj,
    PyObject *glon_obj,
    PyObject *glat_obj,
    int expected_ndim,
    GeostrophicInputs *inputs
) {
    inputs->z_arr = NULL;
    inputs->glon_arr = NULL;
    inputs->glat_arr = NULL;

    inputs->z_arr = (PyArrayObject *) PyArray_FROM_OTF(
        z_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS
    );
    if (inputs->z_arr == NULL) {
        return -1;
    }

    inputs->glon_arr = (PyArrayObject *) PyArray_FROM_OTF(
        glon_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    if (inputs->glon_arr == NULL) {
        return -1;
    }

    inputs->glat_arr = (PyArrayObject *) PyArray_FROM_OTF(
        glat_obj,
        NPY_FLOAT64,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
    if (inputs->glat_arr == NULL) {
        return -1;
    }

    if (PyArray_NDIM(inputs->z_arr) != expected_ndim) {
        PyErr_Format(PyExc_ValueError, "z must be a %dD array", expected_ndim);
        return -1;
    }
    if (PyArray_NDIM(inputs->glon_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "glon must be a 1D array");
        return -1;
    }
    if (PyArray_NDIM(inputs->glat_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "glat must be a 1D array");
        return -1;
    }

    return 0;
}

static void release_inputs(GeostrophicInputs *inputs) {
    Py_XDECREF(inputs->z_arr);
    Py_XDECREF(inputs->glon_arr);
    Py_XDECREF(inputs->glat_arr);
}

static int validate_dims_2d(
    GeostrophicInputs *inputs,
    int *nlat,
    int *mlon,
    PyObject *kwargs
) {
    int kw_nlat = -1;
    int kw_mlon = -1;
    int status = 0;

    if (*nlat < 0) {
        *nlat = (int) PyArray_DIM(inputs->z_arr, 0);
    }
    if (*mlon < 0) {
        *mlon = (int) PyArray_DIM(inputs->z_arr, 1);
    }

    status = get_optional_int(kwargs, "nlat", &kw_nlat);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_nlat != *nlat) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }

    status = get_optional_int(kwargs, "mlon", &kw_mlon);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_mlon != *mlon) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }

    if (*nlat != (int) PyArray_DIM(inputs->z_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->z_arr, 1)) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->glon_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glon length does not match z.shape[1]");
        return -1;
    }
    if (*nlat != (int) PyArray_DIM(inputs->glat_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glat length does not match z.shape[0]");
        return -1;
    }

    return 0;
}

static int validate_dims_3d(
    GeostrophicInputs *inputs,
    int *nlat,
    int *mlon,
    int *n3rd,
    PyObject *kwargs
) {
    int kw_nlat = -1;
    int kw_mlon = -1;
    int kw_n3rd = -1;
    int status = 0;

    if (*nlat < 0) {
        *nlat = (int) PyArray_DIM(inputs->z_arr, 0);
    }
    if (*mlon < 0) {
        *mlon = (int) PyArray_DIM(inputs->z_arr, 1);
    }
    if (*n3rd < 0) {
        *n3rd = (int) PyArray_DIM(inputs->z_arr, 2);
    }

    status = get_optional_int(kwargs, "nlat", &kw_nlat);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_nlat != *nlat) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }

    status = get_optional_int(kwargs, "mlon", &kw_mlon);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_mlon != *mlon) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }

    status = get_optional_int(kwargs, "n3rd", &kw_n3rd);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_n3rd != *n3rd) {
        PyErr_SetString(PyExc_ValueError, "n3rd does not match z.shape[2]");
        return -1;
    }

    if (*nlat != (int) PyArray_DIM(inputs->z_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->z_arr, 1)) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }
    if (*n3rd != (int) PyArray_DIM(inputs->z_arr, 2)) {
        PyErr_SetString(PyExc_ValueError, "n3rd does not match z.shape[2]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->glon_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glon length does not match z.shape[1]");
        return -1;
    }
    if (*nlat != (int) PyArray_DIM(inputs->glat_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glat length does not match z.shape[0]");
        return -1;
    }

    return 0;
}

static int validate_dims_4d(
    GeostrophicInputs *inputs,
    int *nlat,
    int *mlon,
    int *n3rd,
    int *n4th,
    PyObject *kwargs
) {
    int kw_nlat = -1;
    int kw_mlon = -1;
    int kw_n3rd = -1;
    int kw_n4th = -1;
    int status = 0;

    if (*nlat < 0) {
        *nlat = (int) PyArray_DIM(inputs->z_arr, 0);
    }
    if (*mlon < 0) {
        *mlon = (int) PyArray_DIM(inputs->z_arr, 1);
    }
    if (*n3rd < 0) {
        *n3rd = (int) PyArray_DIM(inputs->z_arr, 2);
    }
    if (*n4th < 0) {
        *n4th = (int) PyArray_DIM(inputs->z_arr, 3);
    }

    status = get_optional_int(kwargs, "nlat", &kw_nlat);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_nlat != *nlat) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }

    status = get_optional_int(kwargs, "mlon", &kw_mlon);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_mlon != *mlon) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }

    status = get_optional_int(kwargs, "n3rd", &kw_n3rd);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_n3rd != *n3rd) {
        PyErr_SetString(PyExc_ValueError, "n3rd does not match z.shape[2]");
        return -1;
    }

    status = get_optional_int(kwargs, "n4th", &kw_n4th);
    if (status < 0) {
        return -1;
    }
    if (status > 0 && kw_n4th != *n4th) {
        PyErr_SetString(PyExc_ValueError, "n4th does not match z.shape[3]");
        return -1;
    }

    if (*nlat != (int) PyArray_DIM(inputs->z_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "nlat does not match z.shape[0]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->z_arr, 1)) {
        PyErr_SetString(PyExc_ValueError, "mlon does not match z.shape[1]");
        return -1;
    }
    if (*n3rd != (int) PyArray_DIM(inputs->z_arr, 2)) {
        PyErr_SetString(PyExc_ValueError, "n3rd does not match z.shape[2]");
        return -1;
    }
    if (*n4th != (int) PyArray_DIM(inputs->z_arr, 3)) {
        PyErr_SetString(PyExc_ValueError, "n4th does not match z.shape[3]");
        return -1;
    }
    if (*mlon != (int) PyArray_DIM(inputs->glon_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glon length does not match z.shape[1]");
        return -1;
    }
    if (*nlat != (int) PyArray_DIM(inputs->glat_arr, 0)) {
        PyErr_SetString(PyExc_ValueError, "glat length does not match z.shape[0]");
        return -1;
    }

    return 0;
}

static PyObject *call_2d_backend(
    PyObject *args,
    PyObject *kwargs,
    geostrophic_2d_fn backend
) {
    GeostrophicInputs inputs;
    PyObject *z_obj = NULL;
    PyObject *glon_obj = NULL;
    PyObject *glat_obj = NULL;
    PyArrayObject *ug_arr = NULL;
    PyArrayObject *vg_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[2];
    int nlat = -1;
    int mlon = -1;
    int iopt = 0;
    double zmsg = 0.0;

    if (parse_2d_call(
            args, kwargs, &z_obj, &nlat, &mlon, &zmsg, &glon_obj, &glat_obj, &iopt
        ) < 0) {
        return NULL;
    }

    if (prepare_common_inputs(z_obj, glon_obj, glat_obj, 2, &inputs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    if (validate_dims_2d(&inputs, &nlat, &mlon, kwargs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    dims[0] = nlat;
    dims[1] = mlon;
    ug_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    vg_arr = (PyArrayObject *) PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (ug_arr == NULL || vg_arr == NULL) {
        release_inputs(&inputs);
        Py_XDECREF(ug_arr);
        Py_XDECREF(vg_arr);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    backend(
        PyArray_DATA(inputs.z_arr),
        nlat,
        mlon,
        zmsg,
        PyArray_DATA(inputs.glon_arr),
        PyArray_DATA(inputs.glat_arr),
        PyArray_DATA(ug_arr),
        PyArray_DATA(vg_arr),
        iopt
    );
    Py_END_ALLOW_THREADS

    release_inputs(&inputs);
    result = PyTuple_Pack(2, (PyObject *) ug_arr, (PyObject *) vg_arr);
    Py_DECREF(ug_arr);
    Py_DECREF(vg_arr);
    return result;
}

static PyObject *call_3d_backend(
    PyObject *args,
    PyObject *kwargs,
    geostrophic_3d_fn backend
) {
    GeostrophicInputs inputs;
    PyObject *z_obj = NULL;
    PyObject *glon_obj = NULL;
    PyObject *glat_obj = NULL;
    PyArrayObject *ug_arr = NULL;
    PyArrayObject *vg_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[3];
    int nlat = -1;
    int mlon = -1;
    int n3rd = -1;
    int iopt = 0;
    double zmsg = 0.0;

    if (parse_3d_call(
            args,
            kwargs,
            &z_obj,
            &nlat,
            &mlon,
            &n3rd,
            &zmsg,
            &glon_obj,
            &glat_obj,
            &iopt
        ) < 0) {
        return NULL;
    }

    if (prepare_common_inputs(z_obj, glon_obj, glat_obj, 3, &inputs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    if (validate_dims_3d(&inputs, &nlat, &mlon, &n3rd, kwargs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    dims[0] = nlat;
    dims[1] = mlon;
    dims[2] = n3rd;
    ug_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    vg_arr = (PyArrayObject *) PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    if (ug_arr == NULL || vg_arr == NULL) {
        release_inputs(&inputs);
        Py_XDECREF(ug_arr);
        Py_XDECREF(vg_arr);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    backend(
        PyArray_DATA(inputs.z_arr),
        nlat,
        mlon,
        n3rd,
        zmsg,
        PyArray_DATA(inputs.glon_arr),
        PyArray_DATA(inputs.glat_arr),
        PyArray_DATA(ug_arr),
        PyArray_DATA(vg_arr),
        iopt
    );
    Py_END_ALLOW_THREADS

    release_inputs(&inputs);
    result = PyTuple_Pack(2, (PyObject *) ug_arr, (PyObject *) vg_arr);
    Py_DECREF(ug_arr);
    Py_DECREF(vg_arr);
    return result;
}

static PyObject *call_4d_backend(
    PyObject *args,
    PyObject *kwargs,
    geostrophic_4d_fn backend
) {
    GeostrophicInputs inputs;
    PyObject *z_obj = NULL;
    PyObject *glon_obj = NULL;
    PyObject *glat_obj = NULL;
    PyArrayObject *ug_arr = NULL;
    PyArrayObject *vg_arr = NULL;
    PyObject *result = NULL;
    npy_intp dims[4];
    int nlat = -1;
    int mlon = -1;
    int n3rd = -1;
    int n4th = -1;
    int iopt = 0;
    double zmsg = 0.0;

    if (parse_4d_call(
            args,
            kwargs,
            &z_obj,
            &nlat,
            &mlon,
            &n3rd,
            &n4th,
            &zmsg,
            &glon_obj,
            &glat_obj,
            &iopt
        ) < 0) {
        return NULL;
    }

    if (prepare_common_inputs(z_obj, glon_obj, glat_obj, 4, &inputs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    if (validate_dims_4d(&inputs, &nlat, &mlon, &n3rd, &n4th, kwargs) < 0) {
        release_inputs(&inputs);
        return NULL;
    }

    dims[0] = nlat;
    dims[1] = mlon;
    dims[2] = n3rd;
    dims[3] = n4th;
    ug_arr = (PyArrayObject *) PyArray_EMPTY(4, dims, NPY_FLOAT64, 1);
    vg_arr = (PyArrayObject *) PyArray_EMPTY(4, dims, NPY_FLOAT64, 1);
    if (ug_arr == NULL || vg_arr == NULL) {
        release_inputs(&inputs);
        Py_XDECREF(ug_arr);
        Py_XDECREF(vg_arr);
        return NULL;
    }

    Py_BEGIN_ALLOW_THREADS
    backend(
        PyArray_DATA(inputs.z_arr),
        nlat,
        mlon,
        n3rd,
        n4th,
        zmsg,
        PyArray_DATA(inputs.glon_arr),
        PyArray_DATA(inputs.glat_arr),
        PyArray_DATA(ug_arr),
        PyArray_DATA(vg_arr),
        iopt
    );
    Py_END_ALLOW_THREADS

    release_inputs(&inputs);
    result = PyTuple_Pack(2, (PyObject *) ug_arr, (PyObject *) vg_arr);
    Py_DECREF(ug_arr);
    Py_DECREF(vg_arr);
    return result;
}

static PyObject *py_z2geouv(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void) self;
    return call_2d_backend(args, kwargs, z2geouv_c);
}

static PyObject *py_z2geouv_3d(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void) self;
    return call_3d_backend(args, kwargs, z2geouv_3d_c);
}

static PyObject *py_z2geouv_4d(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void) self;
    return call_4d_backend(args, kwargs, z2geouv_4d_c);
}

static PyObject *py_zuvnew(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void) self;
    return call_2d_backend(args, kwargs, zuvnew_c);
}

static PyObject *py_z2guv(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void) self;
    return call_2d_backend(args, kwargs, z2guv_c);
}

static PyMethodDef module_methods[] = {
    {
        "z2geouv",
        (PyCFunction) (void (*)(void)) py_z2geouv,
        METH_VARARGS | METH_KEYWORDS,
        "Compute 2D geostrophic wind components from geopotential height."
    },
    {
        "z2geouv_3d",
        (PyCFunction) (void (*)(void)) py_z2geouv_3d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute 3D geostrophic wind components from geopotential height."
    },
    {
        "z2geouv_4d",
        (PyCFunction) (void (*)(void)) py_z2geouv_4d,
        METH_VARARGS | METH_KEYWORDS,
        "Compute 4D geostrophic wind components from geopotential height."
    },
    {
        "zuvnew",
        (PyCFunction) (void (*)(void)) py_zuvnew,
        METH_VARARGS | METH_KEYWORDS,
        "Reorder north-to-south latitude input before geostrophic computation."
    },
    {
        "z2guv",
        (PyCFunction) (void (*)(void)) py_z2guv,
        METH_VARARGS | METH_KEYWORDS,
        "Compute the low-level 2D geostrophic wind kernel directly."
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "geostrophicwind",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_geostrophicwind(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
