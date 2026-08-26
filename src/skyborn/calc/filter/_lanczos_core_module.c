#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <numpy/arrayobject.h>
#include <string.h>

/* External Fortran subroutines */
extern void compute_lanczos_weights_(double *cutoff_freq, int *nwt, int *pass_type,
                                     double *weights, int *ierr);
extern void apply_filter_1d_(double *data, int *n, double *weights, int *nwt,
                             double *filtered);
extern void apply_filter_1d_periodic_(double *data, int *n, double *weights, int *nwt,
                                      double *filtered);
extern void apply_filter_2d_(double *data, int *ny, int *nx,
                             double *weights_y, int *nwt_y,
                             double *weights_x, int *nwt_x,
                             double *filtered);
extern void apply_filter_3d_(double *data, int *nz, int *ny, int *nx,
                             double *weights_z, int *nwt_z,
                             double *weights_y, int *nwt_y,
                             double *weights_x, int *nwt_x,
                             double *filtered);

/* Python wrapper for compute_lanczos_weights */
static PyObject *py_compute_lanczos_weights(PyObject *self, PyObject *args) {
    double cutoff_freq;
    int nwt, pass_type;
    int ierr = 0;
    npy_intp dims[1];
    PyArrayObject *weights_arr = NULL;
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "dii", &cutoff_freq, &nwt, &pass_type)) {
        return NULL;
    }

    /* Allocate output array */
    dims[0] = nwt;
    weights_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 0);
    if (weights_arr == NULL) {
        return NULL;
    }

    /* Call Fortran subroutine */
    compute_lanczos_weights_(&cutoff_freq, &nwt, &pass_type,
                            (double *)PyArray_DATA(weights_arr), &ierr);

    if (ierr != 0) {
        Py_DECREF(weights_arr);
        if (ierr == 1) {
            PyErr_SetString(PyExc_ValueError, "window size must be >= 3");
        } else if (ierr == 2) {
            PyErr_SetString(PyExc_ValueError, "window size must be odd");
        } else if (ierr == 3) {
            PyErr_SetString(PyExc_ValueError, "cutoff_freq must be in (0, 0.5)");
        } else if (ierr == 4) {
            PyErr_SetString(PyExc_RuntimeError, "weight normalization failed");
        } else if (ierr == 5) {
            PyErr_SetString(PyExc_ValueError, "pass_type must be 1 (low) or 2 (high)");
        } else {
            PyErr_Format(PyExc_RuntimeError, "Fortran routine failed with error code %d", ierr);
        }
        return NULL;
    }

    /* Return tuple (weights, ierr) */
    result = Py_BuildValue("Oi", weights_arr, ierr);
    Py_DECREF(weights_arr);
    return result;
}

/* Python wrapper for apply_filter_1d */
static PyObject *py_apply_filter_1d(PyObject *self, PyObject *args) {
    PyObject *data_obj, *weights_obj;
    PyArrayObject *data_arr = NULL, *weights_arr = NULL, *filtered_arr = NULL;
    int n, nwt;
    npy_intp dims[1];
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "OO", &data_obj, &weights_obj)) {
        return NULL;
    }

    /* Convert input arrays to Fortran-contiguous float64 */
    data_arr = (PyArrayObject *)PyArray_FROM_OTF(data_obj, NPY_FLOAT64,
                                                  NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_obj, NPY_FLOAT64,
                                                     NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);

    if (data_arr == NULL || weights_arr == NULL) {
        Py_XDECREF(data_arr);
        Py_XDECREF(weights_arr);
        return NULL;
    }

    /* Check dimensions */
    if (PyArray_NDIM(data_arr) != 1 || PyArray_NDIM(weights_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "data and weights must be 1D arrays");
        Py_DECREF(data_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }

    n = (int)PyArray_DIM(data_arr, 0);
    nwt = (int)PyArray_DIM(weights_arr, 0);

    /* Allocate output array */
    dims[0] = n;
    filtered_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 1);
    if (filtered_arr == NULL) {
        Py_DECREF(data_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }

    /* Call Fortran subroutine */
    apply_filter_1d_((double *)PyArray_DATA(data_arr), &n,
                     (double *)PyArray_DATA(weights_arr), &nwt,
                     (double *)PyArray_DATA(filtered_arr));

    result = (PyObject *)filtered_arr;
    Py_DECREF(data_arr);
    Py_DECREF(weights_arr);
    return result;
}

/* Python wrapper for apply_filter_1d_periodic */
static PyObject *py_apply_filter_1d_periodic(PyObject *self, PyObject *args) {
    PyObject *data_obj, *weights_obj;
    PyArrayObject *data_arr = NULL, *weights_arr = NULL, *filtered_arr = NULL;
    int n, nwt;
    npy_intp dims[1];
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "OO", &data_obj, &weights_obj)) {
        return NULL;
    }

    data_arr = (PyArrayObject *)PyArray_FROM_OTF(data_obj, NPY_FLOAT64,
                                                  NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_obj, NPY_FLOAT64,
                                                     NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);

    if (data_arr == NULL || weights_arr == NULL) {
        Py_XDECREF(data_arr);
        Py_XDECREF(weights_arr);
        return NULL;
    }

    if (PyArray_NDIM(data_arr) != 1 || PyArray_NDIM(weights_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "data and weights must be 1D arrays");
        Py_DECREF(data_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }

    n = (int)PyArray_DIM(data_arr, 0);
    nwt = (int)PyArray_DIM(weights_arr, 0);

    dims[0] = n;
    filtered_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 1);
    if (filtered_arr == NULL) {
        Py_DECREF(data_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }

    apply_filter_1d_periodic_((double *)PyArray_DATA(data_arr), &n,
                              (double *)PyArray_DATA(weights_arr), &nwt,
                              (double *)PyArray_DATA(filtered_arr));

    result = (PyObject *)filtered_arr;
    Py_DECREF(data_arr);
    Py_DECREF(weights_arr);
    return result;
}

/* Python wrapper for apply_filter_2d */
static PyObject *py_apply_filter_2d(PyObject *self, PyObject *args) {
    PyObject *data_obj, *weights_y_obj, *weights_x_obj;
    PyArrayObject *data_arr = NULL, *weights_y_arr = NULL, *weights_x_arr = NULL;
    PyArrayObject *filtered_arr = NULL;
    int ny, nx, nwt_y, nwt_x;
    npy_intp dims[2];
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "OOO", &data_obj, &weights_y_obj, &weights_x_obj)) {
        return NULL;
    }

    data_arr = (PyArrayObject *)PyArray_FROM_OTF(data_obj, NPY_FLOAT64,
                                                  NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_y_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_y_obj, NPY_FLOAT64,
                                                       NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_x_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_x_obj, NPY_FLOAT64,
                                                       NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);

    if (data_arr == NULL || weights_y_arr == NULL || weights_x_arr == NULL) {
        Py_XDECREF(data_arr);
        Py_XDECREF(weights_y_arr);
        Py_XDECREF(weights_x_arr);
        return NULL;
    }

    if (PyArray_NDIM(data_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "data must be 2D array");
        Py_DECREF(data_arr);
        Py_DECREF(weights_y_arr);
        Py_DECREF(weights_x_arr);
        return NULL;
    }

    ny = (int)PyArray_DIM(data_arr, 0);
    nx = (int)PyArray_DIM(data_arr, 1);
    nwt_y = (int)PyArray_DIM(weights_y_arr, 0);
    nwt_x = (int)PyArray_DIM(weights_x_arr, 0);

    dims[0] = ny;
    dims[1] = nx;
    filtered_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT64, 1);
    if (filtered_arr == NULL) {
        Py_DECREF(data_arr);
        Py_DECREF(weights_y_arr);
        Py_DECREF(weights_x_arr);
        return NULL;
    }

    apply_filter_2d_((double *)PyArray_DATA(data_arr), &ny, &nx,
                     (double *)PyArray_DATA(weights_y_arr), &nwt_y,
                     (double *)PyArray_DATA(weights_x_arr), &nwt_x,
                     (double *)PyArray_DATA(filtered_arr));

    result = (PyObject *)filtered_arr;
    Py_DECREF(data_arr);
    Py_DECREF(weights_y_arr);
    Py_DECREF(weights_x_arr);
    return result;
}

/* Python wrapper for apply_filter_3d */
static PyObject *py_apply_filter_3d(PyObject *self, PyObject *args) {
    PyObject *data_obj, *weights_z_obj, *weights_y_obj, *weights_x_obj;
    PyArrayObject *data_arr = NULL, *weights_z_arr = NULL;
    PyArrayObject *weights_y_arr = NULL, *weights_x_arr = NULL;
    PyArrayObject *filtered_arr = NULL;
    int nz, ny, nx, nwt_z, nwt_y, nwt_x;
    npy_intp dims[3];
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOO", &data_obj, &weights_z_obj,
                          &weights_y_obj, &weights_x_obj)) {
        return NULL;
    }

    data_arr = (PyArrayObject *)PyArray_FROM_OTF(data_obj, NPY_FLOAT64,
                                                  NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_z_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_z_obj, NPY_FLOAT64,
                                                       NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_y_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_y_obj, NPY_FLOAT64,
                                                       NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    weights_x_arr = (PyArrayObject *)PyArray_FROM_OTF(weights_x_obj, NPY_FLOAT64,
                                                       NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);

    if (data_arr == NULL || weights_z_arr == NULL ||
        weights_y_arr == NULL || weights_x_arr == NULL) {
        Py_XDECREF(data_arr);
        Py_XDECREF(weights_z_arr);
        Py_XDECREF(weights_y_arr);
        Py_XDECREF(weights_x_arr);
        return NULL;
    }

    if (PyArray_NDIM(data_arr) != 3) {
        PyErr_SetString(PyExc_ValueError, "data must be 3D array");
        Py_DECREF(data_arr);
        Py_DECREF(weights_z_arr);
        Py_DECREF(weights_y_arr);
        Py_DECREF(weights_x_arr);
        return NULL;
    }

    nz = (int)PyArray_DIM(data_arr, 0);
    ny = (int)PyArray_DIM(data_arr, 1);
    nx = (int)PyArray_DIM(data_arr, 2);
    nwt_z = (int)PyArray_DIM(weights_z_arr, 0);
    nwt_y = (int)PyArray_DIM(weights_y_arr, 0);
    nwt_x = (int)PyArray_DIM(weights_x_arr, 0);

    dims[0] = nz;
    dims[1] = ny;
    dims[2] = nx;
    filtered_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT64, 1);
    if (filtered_arr == NULL) {
        Py_DECREF(data_arr);
        Py_DECREF(weights_z_arr);
        Py_DECREF(weights_y_arr);
        Py_DECREF(weights_x_arr);
        return NULL;
    }

    apply_filter_3d_((double *)PyArray_DATA(data_arr), &nz, &ny, &nx,
                     (double *)PyArray_DATA(weights_z_arr), &nwt_z,
                     (double *)PyArray_DATA(weights_y_arr), &nwt_y,
                     (double *)PyArray_DATA(weights_x_arr), &nwt_x,
                     (double *)PyArray_DATA(filtered_arr));

    result = (PyObject *)filtered_arr;
    Py_DECREF(data_arr);
    Py_DECREF(weights_z_arr);
    Py_DECREF(weights_y_arr);
    Py_DECREF(weights_x_arr);
    return result;
}

/* Method definitions */
static PyMethodDef lanczos_methods[] = {
    {"compute_lanczos_weights", py_compute_lanczos_weights, METH_VARARGS,
     "Compute Lanczos filter weights\n\n"
     "Parameters:\n"
     "  cutoff_freq : float - cutoff frequency (0 to 0.5)\n"
     "  nwt : int - number of weights (must be odd)\n"
     "  pass_type : int - 1 for low-pass, 2 for high-pass\n\n"
     "Returns:\n"
     "  (weights, ierr) : tuple of (ndarray, int)"},

    {"apply_filter_1d", py_apply_filter_1d, METH_VARARGS,
     "Apply 1D Lanczos filter with reflect boundary\n\n"
     "Parameters:\n"
     "  data : ndarray - 1D input data\n"
     "  weights : ndarray - filter weights\n\n"
     "Returns:\n"
     "  filtered : ndarray - filtered data"},

    {"apply_filter_1d_periodic", py_apply_filter_1d_periodic, METH_VARARGS,
     "Apply 1D Lanczos filter with periodic boundary\n\n"
     "Parameters:\n"
     "  data : ndarray - 1D input data\n"
     "  weights : ndarray - filter weights\n\n"
     "Returns:\n"
     "  filtered : ndarray - filtered data"},

    {"apply_filter_2d", py_apply_filter_2d, METH_VARARGS,
     "Apply 2D separable Lanczos filter\n\n"
     "Parameters:\n"
     "  data : ndarray - 2D input data (ny, nx)\n"
     "  weights_y : ndarray - filter weights for y-direction\n"
     "  weights_x : ndarray - filter weights for x-direction\n\n"
     "Returns:\n"
     "  filtered : ndarray - filtered data"},

    {"apply_filter_3d", py_apply_filter_3d, METH_VARARGS,
     "Apply 3D separable Lanczos filter\n\n"
     "Parameters:\n"
     "  data : ndarray - 3D input data (nz, ny, nx)\n"
     "  weights_z : ndarray - filter weights for z-direction\n"
     "  weights_y : ndarray - filter weights for y-direction\n"
     "  weights_x : ndarray - filter weights for x-direction\n\n"
     "Returns:\n"
     "  filtered : ndarray - filtered data"},

    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef lanczos_module = {
    PyModuleDef_HEAD_INIT,
    "_lanczos_core",
    "Fortran-accelerated Lanczos filtering routines",
    -1,
    lanczos_methods
};

/* Module initialization */
PyMODINIT_FUNC PyInit__lanczos_core(void) {
    PyObject *module;

    import_array();
    if (PyErr_Occurred()) {
        return NULL;
    }

    module = PyModule_Create(&lanczos_module);
    if (module == NULL) {
        return NULL;
    }

    return module;
}
