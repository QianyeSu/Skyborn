#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <numpy/arrayobject.h>

int reduced_import_array(void) {
    import_array1(-1);
    return 0;
}

void reduced_gaqd_c(int nlat, void *theta, void *wts, int *ierror);
void reduced_lap_c(void *dataspec, void *dataspec_lap, int nmdim, int nt, float rsphere);
void reduced_invlap_c(void *dataspec, void *dataspec_ilap, int nmdim, int nt, float rsphere);
void reduced_multsmoothfact_c(void *dataspec, void *dataspec_smooth, void *smooth, int nlat, int nmdim, int nt);
void reduced_gaussian_legendre_basis_c(void *basis, int nlat, int ntrunc, int *ierror);
void reduced_gaussian_legendre_derivative_from_basis_c(void *basis, void *dbasis, int nlat, int ntrunc, int *ierror);
void reduced_gaussian_legendre_derivative_basis_c(void *dbasis, int nlat, int ntrunc, int *ierror);
void reduced_gaussian_grdtospec_c(void *datagrid, void *pl, void *weights, void *basis, void *dataspec, int ngptot, int nlat, int ntrunc, int nt, int *ierror);
void reduced_gaussian_spectogrd_c(void *dataspec, void *pl, void *basis, void *datagrid, int nmdim, int nlat, int ntrunc, int nt, int ngptot, int *ierror);
void reduced_gaussian_spectogrd_pair_c(void *speca, void *specb, void *pl, void *basis, void *grida, void *gridb, int nmdim, int nlat, int ntrunc, int nt, int ngptot, int *ierror);
void reduced_gaussian_getgrad_c(void *dataspec, void *pl, void *basis, void *dbasis, void *sin_theta, void *ugrad, void *vgrad, int nmdim, int nlat, int ntrunc, int nt, int ngptot, float rsphere, int *ierror);
void reduced_gaussian_getgrad_pair_c(void *speca, void *specb, void *pl, void *basis, void *dbasis, void *sin_theta, void *a_ugrad, void *a_vgrad, void *b_ugrad, void *b_vgrad, int nmdim, int nlat, int ntrunc, int nt, int ngptot, float rsphere, int *ierror);
void reduced_gaussian_getvrtdivspec_c(void *ugrid, void *vgrid, void *pl, void *weights, void *basis, void *dbasis, void *sin_theta, void *vrtspec, void *divspec, int ngptot, int nlat, int ntrunc, int nt, float rsphere, int *ierror);
void reduced_gaussian_getvrtspec_c(void *ugrid, void *vgrid, void *pl, void *weights, void *basis, void *dbasis, void *sin_theta, void *vrtspec, int ngptot, int nlat, int ntrunc, int nt, float rsphere, int *ierror);
void reduced_gaussian_getdivspec_c(void *ugrid, void *vgrid, void *pl, void *weights, void *basis, void *dbasis, void *sin_theta, void *divspec, int ngptot, int nlat, int ntrunc, int nt, float rsphere, int *ierror);

static int as_complex64_2d_fortran(PyObject *obj, PyArrayObject **arr, int *nmdim, int *nt) {
    *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_COMPLEX64, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (*arr == NULL) {
        return -1;
    }
    if (PyArray_NDIM(*arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "expected a 2D complex64 Fortran-contiguous array");
        Py_DECREF(*arr);
        *arr = NULL;
        return -1;
    }
    *nmdim = (int)PyArray_DIM(*arr, 0);
    *nt = (int)PyArray_DIM(*arr, 1);
    return 0;
}

static int as_real32_2d_fortran(PyObject *obj, PyArrayObject **arr, int *d0, int *d1) {
    *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (*arr == NULL) {
        return -1;
    }
    if (PyArray_NDIM(*arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "expected a 2D float32 Fortran-contiguous array");
        Py_DECREF(*arr);
        *arr = NULL;
        return -1;
    }
    *d0 = (int)PyArray_DIM(*arr, 0);
    *d1 = (int)PyArray_DIM(*arr, 1);
    return 0;
}

static int as_real32_1d_contig(PyObject *obj, PyArrayObject **arr, int *n) {
    *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (*arr == NULL) {
        return -1;
    }
    if (PyArray_NDIM(*arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "expected a 1D float32 array");
        Py_DECREF(*arr);
        *arr = NULL;
        return -1;
    }
    *n = (int)PyArray_DIM(*arr, 0);
    return 0;
}

static int as_int32_1d_contig(PyObject *obj, PyArrayObject **arr, int *n) {
    *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_INT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (*arr == NULL) {
        return -1;
    }
    if (PyArray_NDIM(*arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "expected a 1D int32 array");
        Py_DECREF(*arr);
        *arr = NULL;
        return -1;
    }
    *n = (int)PyArray_DIM(*arr, 0);
    return 0;
}

static PyObject *py_reduced_gaqd(PyObject *self, PyObject *args) {
    int nlat;
    int ldwork;
    int ierror = 0;
    npy_intp dims[1];
    PyArrayObject *theta_arr = NULL;
    PyArrayObject *wts_arr = NULL;
    PyArrayObject *dwork_arr = NULL;
    PyObject *result = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "i", &nlat)) {
        return NULL;
    }
    ldwork = nlat * (nlat + 2);
    dims[0] = nlat;
    theta_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 1);
    wts_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 1);
    dims[0] = ldwork;
    dwork_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT64, 1);
    if (theta_arr == NULL || wts_arr == NULL || dwork_arr == NULL) {
        goto fail;
    }

    reduced_gaqd_c(nlat, PyArray_DATA(theta_arr), PyArray_DATA(wts_arr), &ierror);
    result = Py_BuildValue("NNi", theta_arr, wts_arr, ierror);
    Py_DECREF(dwork_arr);
    return result;

fail:
    Py_XDECREF(theta_arr);
    Py_XDECREF(wts_arr);
    Py_XDECREF(dwork_arr);
    return NULL;
}

static PyObject *py_reduced_lap(PyObject *self, PyObject *args) {
    PyObject *spec_obj = NULL;
    PyArrayObject *spec_arr = NULL;
    PyArrayObject *out_arr = NULL;
    float rsphere;
    int nmdim, nt;
    (void)self;

    if (!PyArg_ParseTuple(args, "Of", &spec_obj, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(spec_obj, &spec_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    out_arr = (PyArrayObject *)PyArray_EMPTY(PyArray_NDIM(spec_arr), PyArray_DIMS(spec_arr), NPY_COMPLEX64, 1);
    if (out_arr == NULL) {
        Py_DECREF(spec_arr);
        return NULL;
    }
    reduced_lap_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), nmdim, nt, rsphere);
    Py_DECREF(spec_arr);
    return (PyObject *)out_arr;
}

static PyObject *py_reduced_invlap(PyObject *self, PyObject *args) {
    PyObject *spec_obj = NULL;
    PyArrayObject *spec_arr = NULL;
    PyArrayObject *out_arr = NULL;
    float rsphere;
    int nmdim, nt;
    (void)self;

    if (!PyArg_ParseTuple(args, "Of", &spec_obj, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(spec_obj, &spec_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    out_arr = (PyArrayObject *)PyArray_EMPTY(PyArray_NDIM(spec_arr), PyArray_DIMS(spec_arr), NPY_COMPLEX64, 1);
    if (out_arr == NULL) {
        Py_DECREF(spec_arr);
        return NULL;
    }
    reduced_invlap_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), nmdim, nt, rsphere);
    Py_DECREF(spec_arr);
    return (PyObject *)out_arr;
}

static PyObject *py_reduced_multsmoothfact(PyObject *self, PyObject *args) {
    PyObject *spec_obj = NULL;
    PyObject *smooth_obj = NULL;
    PyArrayObject *spec_arr = NULL;
    PyArrayObject *smooth_arr = NULL;
    PyArrayObject *out_arr = NULL;
    int nmdim, nt, nlat;
    (void)self;

    if (!PyArg_ParseTuple(args, "OO", &spec_obj, &smooth_obj)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(spec_obj, &spec_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    if (as_real32_1d_contig(smooth_obj, &smooth_arr, &nlat) != 0) {
        Py_DECREF(spec_arr);
        return NULL;
    }
    out_arr = (PyArrayObject *)PyArray_EMPTY(PyArray_NDIM(spec_arr), PyArray_DIMS(spec_arr), NPY_COMPLEX64, 1);
    if (out_arr == NULL) {
        Py_DECREF(spec_arr);
        Py_DECREF(smooth_arr);
        return NULL;
    }
    reduced_multsmoothfact_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), PyArray_DATA(smooth_arr), nlat, nmdim, nt);
    Py_DECREF(spec_arr);
    Py_DECREF(smooth_arr);
    return (PyObject *)out_arr;
}

static PyObject *py_reduced_gaussian_legendre_basis(PyObject *self, PyObject *args) {
    int nlat, ntrunc, nmdim, ierror = 0;
    npy_intp dims[2];
    PyArrayObject *basis_arr = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "ii", &nlat, &ntrunc)) {
        return NULL;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nlat;
    dims[1] = nmdim;
    basis_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (basis_arr == NULL) {
        return NULL;
    }
    reduced_gaussian_legendre_basis_c(PyArray_DATA(basis_arr), nlat, ntrunc, &ierror);
    return Py_BuildValue("Ni", basis_arr, ierror);
}

static PyObject *py_reduced_gaussian_legendre_derivative_from_basis(PyObject *self, PyObject *args) {
    PyObject *basis_obj = NULL;
    PyArrayObject *basis_arr = NULL;
    PyArrayObject *dbasis_arr = NULL;
    int nlat, nmdim, ntrunc, expected_nmdim, ierror = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "Oi", &basis_obj, &ntrunc)) {
        return NULL;
    }
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &nlat, &nmdim) != 0) {
        return NULL;
    }
    expected_nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    if (nmdim != expected_nmdim) {
        PyErr_SetString(PyExc_ValueError, "basis shape does not match ntrunc");
        Py_DECREF(basis_arr);
        return NULL;
    }
    dims[0] = nlat;
    dims[1] = nmdim;
    dbasis_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (dbasis_arr == NULL) {
        Py_DECREF(basis_arr);
        return NULL;
    }
    reduced_gaussian_legendre_derivative_from_basis_c(
        PyArray_DATA(basis_arr), PyArray_DATA(dbasis_arr), nlat, ntrunc, &ierror);
    Py_DECREF(basis_arr);
    return Py_BuildValue("Ni", dbasis_arr, ierror);
}

static PyObject *py_reduced_gaussian_legendre_derivative_basis(PyObject *self, PyObject *args) {
    int nlat, ntrunc, nmdim, ierror = 0;
    npy_intp dims[2];
    PyArrayObject *dbasis_arr = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "ii", &nlat, &ntrunc)) {
        return NULL;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nlat;
    dims[1] = nmdim;
    dbasis_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (dbasis_arr == NULL) {
        return NULL;
    }
    reduced_gaussian_legendre_derivative_basis_c(PyArray_DATA(dbasis_arr), nlat, ntrunc, &ierror);
    return Py_BuildValue("Ni", dbasis_arr, ierror);
}

static PyObject *py_reduced_gaussian_grdtospec(PyObject *self, PyObject *args) {
    PyObject *datagrid_obj = NULL, *pl_obj = NULL, *weights_obj = NULL, *basis_obj = NULL;
    PyArrayObject *datagrid_arr = NULL, *pl_arr = NULL, *weights_arr = NULL, *basis_arr = NULL, *dataspec_arr = NULL;
    int ngptot, nt, nlat, nlat_w, basis_nlat, basis_nmdim, ntrunc, ierror = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOi", &datagrid_obj, &pl_obj, &weights_obj, &basis_obj, &ntrunc)) {
        return NULL;
    }
    if (as_real32_2d_fortran(datagrid_obj, &datagrid_arr, &ngptot, &nt) != 0) return NULL;
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_1d_contig(weights_obj, &weights_arr, &nlat_w) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nlat, &basis_nmdim) != 0) goto fail;
    if (nlat_w != nlat || basis_nlat != nlat) {
        PyErr_SetString(PyExc_ValueError, "pl, weights, and basis dimensions must match");
        goto fail;
    }
    dims[0] = basis_nmdim;
    dims[1] = nt;
    dataspec_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (dataspec_arr == NULL) goto fail;

    reduced_gaussian_grdtospec_c(
        PyArray_DATA(datagrid_arr), PyArray_DATA(pl_arr), PyArray_DATA(weights_arr),
        PyArray_DATA(basis_arr), PyArray_DATA(dataspec_arr), ngptot, nlat, ntrunc, nt, &ierror);

    Py_DECREF(datagrid_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(weights_arr);
    Py_DECREF(basis_arr);
    return Py_BuildValue("Ni", dataspec_arr, ierror);

fail:
    Py_XDECREF(datagrid_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(weights_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dataspec_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_spectogrd(PyObject *self, PyObject *args) {
    PyObject *dataspec_obj = NULL, *pl_obj = NULL, *basis_obj = NULL;
    PyArrayObject *dataspec_arr = NULL, *pl_arr = NULL, *basis_arr = NULL, *datagrid_arr = NULL;
    int nmdim, nt, nlat, basis_nmdim, basis_nlat, ntrunc, ngptot, ierror = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOii", &dataspec_obj, &pl_obj, &basis_obj, &ntrunc, &ngptot)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(dataspec_obj, &dataspec_arr, &nmdim, &nt) != 0) return NULL;
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (basis_nmdim != nmdim || basis_nlat != nlat) {
        PyErr_SetString(PyExc_ValueError, "basis shape must match spectrum and pl");
        goto fail;
    }
    dims[0] = ngptot;
    dims[1] = nt;
    datagrid_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (datagrid_arr == NULL) goto fail;

    reduced_gaussian_spectogrd_c(
        PyArray_DATA(dataspec_arr), PyArray_DATA(pl_arr), PyArray_DATA(basis_arr),
        PyArray_DATA(datagrid_arr), nmdim, nlat, ntrunc, nt, ngptot, &ierror);

    Py_DECREF(dataspec_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(basis_arr);
    return Py_BuildValue("Ni", datagrid_arr, ierror);

fail:
    Py_XDECREF(dataspec_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(datagrid_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_spectogrd_pair(PyObject *self, PyObject *args) {
    PyObject *speca_obj = NULL, *specb_obj = NULL, *pl_obj = NULL, *basis_obj = NULL;
    PyArrayObject *speca_arr = NULL, *specb_arr = NULL, *pl_arr = NULL, *basis_arr = NULL;
    PyArrayObject *grida_arr = NULL, *gridb_arr = NULL;
    int nmdim_a, nt_a, nmdim_b, nt_b, nlat, basis_nmdim, basis_nlat, ntrunc, ngptot, ierror = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOii", &speca_obj, &specb_obj, &pl_obj, &basis_obj, &ntrunc, &ngptot)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(speca_obj, &speca_arr, &nmdim_a, &nt_a) != 0) return NULL;
    if (as_complex64_2d_fortran(specb_obj, &specb_arr, &nmdim_b, &nt_b) != 0) goto fail;
    if (nmdim_a != nmdim_b || nt_a != nt_b) {
        PyErr_SetString(PyExc_ValueError, "paired spectra must have the same shape");
        goto fail;
    }
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (basis_nmdim != nmdim_a || basis_nlat != nlat) {
        PyErr_SetString(PyExc_ValueError, "basis shape must match spectra and pl");
        goto fail;
    }
    dims[0] = ngptot;
    dims[1] = nt_a;
    grida_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    gridb_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (grida_arr == NULL || gridb_arr == NULL) goto fail;

    reduced_gaussian_spectogrd_pair_c(
        PyArray_DATA(speca_arr), PyArray_DATA(specb_arr), PyArray_DATA(pl_arr), PyArray_DATA(basis_arr),
        PyArray_DATA(grida_arr), PyArray_DATA(gridb_arr), nmdim_a, nlat, ntrunc, nt_a, ngptot, &ierror);

    Py_DECREF(speca_arr);
    Py_DECREF(specb_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(basis_arr);
    return Py_BuildValue("NNi", grida_arr, gridb_arr, ierror);

fail:
    Py_XDECREF(speca_arr);
    Py_XDECREF(specb_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(grida_arr);
    Py_XDECREF(gridb_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_getgrad(PyObject *self, PyObject *args) {
    PyObject *dataspec_obj = NULL, *pl_obj = NULL, *basis_obj = NULL, *dbasis_obj = NULL, *sin_theta_obj = NULL;
    PyArrayObject *dataspec_arr = NULL, *pl_arr = NULL, *basis_arr = NULL, *dbasis_arr = NULL, *sin_theta_arr = NULL;
    PyArrayObject *ugrad_arr = NULL, *vgrad_arr = NULL;
    int nmdim, nt, nlat, basis_nmdim, basis_nlat, dbasis_nmdim, dbasis_nlat, sin_nlat, ntrunc, ngptot, ierror = 0;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOiif", &dataspec_obj, &pl_obj, &basis_obj, &dbasis_obj, &sin_theta_obj, &ntrunc, &ngptot, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(dataspec_obj, &dataspec_arr, &nmdim, &nt) != 0) return NULL;
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (as_real32_2d_fortran(dbasis_obj, &dbasis_arr, &dbasis_nmdim, &dbasis_nlat) != 0) goto fail;
    if (as_real32_1d_contig(sin_theta_obj, &sin_theta_arr, &sin_nlat) != 0) goto fail;
    if (basis_nmdim != nmdim || dbasis_nmdim != nmdim || basis_nlat != nlat || dbasis_nlat != nlat || sin_nlat != nlat) {
        PyErr_SetString(PyExc_ValueError, "basis, dbasis, sin_theta, spectrum, and pl dimensions must match");
        goto fail;
    }
    dims[0] = ngptot;
    dims[1] = nt;
    ugrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    vgrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (ugrad_arr == NULL || vgrad_arr == NULL) goto fail;

    reduced_gaussian_getgrad_c(
        PyArray_DATA(dataspec_arr), PyArray_DATA(pl_arr), PyArray_DATA(basis_arr), PyArray_DATA(dbasis_arr),
        PyArray_DATA(sin_theta_arr), PyArray_DATA(ugrad_arr), PyArray_DATA(vgrad_arr),
        nmdim, nlat, ntrunc, nt, ngptot, rsphere, &ierror);

    Py_DECREF(dataspec_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(dbasis_arr);
    Py_DECREF(sin_theta_arr);
    return Py_BuildValue("NNi", ugrad_arr, vgrad_arr, ierror);

fail:
    Py_XDECREF(dataspec_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dbasis_arr);
    Py_XDECREF(sin_theta_arr);
    Py_XDECREF(ugrad_arr);
    Py_XDECREF(vgrad_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_getgrad_pair(PyObject *self, PyObject *args) {
    PyObject *speca_obj = NULL, *specb_obj = NULL, *pl_obj = NULL, *basis_obj = NULL, *dbasis_obj = NULL, *sin_theta_obj = NULL;
    PyArrayObject *speca_arr = NULL, *specb_arr = NULL, *pl_arr = NULL, *basis_arr = NULL, *dbasis_arr = NULL, *sin_theta_arr = NULL;
    PyArrayObject *a_ugrad_arr = NULL, *a_vgrad_arr = NULL, *b_ugrad_arr = NULL, *b_vgrad_arr = NULL;
    int nmdim_a, nt_a, nmdim_b, nt_b, nlat, basis_nmdim, basis_nlat, dbasis_nmdim, dbasis_nlat, sin_nlat, ntrunc, ngptot, ierror = 0;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOiif", &speca_obj, &specb_obj, &pl_obj, &basis_obj, &dbasis_obj, &sin_theta_obj, &ntrunc, &ngptot, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(speca_obj, &speca_arr, &nmdim_a, &nt_a) != 0) return NULL;
    if (as_complex64_2d_fortran(specb_obj, &specb_arr, &nmdim_b, &nt_b) != 0) goto fail;
    if (nmdim_a != nmdim_b || nt_a != nt_b) {
        PyErr_SetString(PyExc_ValueError, "paired spectra must have the same shape");
        goto fail;
    }
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (as_real32_2d_fortran(dbasis_obj, &dbasis_arr, &dbasis_nmdim, &dbasis_nlat) != 0) goto fail;
    if (as_real32_1d_contig(sin_theta_obj, &sin_theta_arr, &sin_nlat) != 0) goto fail;
    if (basis_nmdim != nmdim_a || dbasis_nmdim != nmdim_a || basis_nlat != nlat || dbasis_nlat != nlat || sin_nlat != nlat) {
        PyErr_SetString(PyExc_ValueError, "basis, dbasis, sin_theta, spectra, and pl dimensions must match");
        goto fail;
    }
    dims[0] = ngptot;
    dims[1] = nt_a;
    a_ugrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    a_vgrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    b_ugrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    b_vgrad_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_FLOAT32, 1);
    if (a_ugrad_arr == NULL || a_vgrad_arr == NULL || b_ugrad_arr == NULL || b_vgrad_arr == NULL) goto fail;

    reduced_gaussian_getgrad_pair_c(
        PyArray_DATA(speca_arr), PyArray_DATA(specb_arr), PyArray_DATA(pl_arr), PyArray_DATA(basis_arr),
        PyArray_DATA(dbasis_arr), PyArray_DATA(sin_theta_arr), PyArray_DATA(a_ugrad_arr), PyArray_DATA(a_vgrad_arr),
        PyArray_DATA(b_ugrad_arr), PyArray_DATA(b_vgrad_arr), nmdim_a, nlat, ntrunc, nt_a, ngptot, rsphere, &ierror);

    Py_DECREF(speca_arr);
    Py_DECREF(specb_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(dbasis_arr);
    Py_DECREF(sin_theta_arr);
    return Py_BuildValue("NNNNi", a_ugrad_arr, a_vgrad_arr, b_ugrad_arr, b_vgrad_arr, ierror);

fail:
    Py_XDECREF(speca_arr);
    Py_XDECREF(specb_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dbasis_arr);
    Py_XDECREF(sin_theta_arr);
    Py_XDECREF(a_ugrad_arr);
    Py_XDECREF(a_vgrad_arr);
    Py_XDECREF(b_ugrad_arr);
    Py_XDECREF(b_vgrad_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_getvrtdivspec(PyObject *self, PyObject *args) {
    PyObject *ugrid_obj = NULL, *vgrid_obj = NULL, *pl_obj = NULL, *weights_obj = NULL, *basis_obj = NULL, *dbasis_obj = NULL, *sin_theta_obj = NULL;
    PyArrayObject *ugrid_arr = NULL, *vgrid_arr = NULL, *pl_arr = NULL, *weights_arr = NULL, *basis_arr = NULL, *dbasis_arr = NULL, *sin_theta_arr = NULL;
    PyArrayObject *vrtspec_arr = NULL, *divspec_arr = NULL;
    int ngptot_u, nt_u, ngptot_v, nt_v, nlat, nlat_w, basis_nmdim, basis_nlat, dbasis_nmdim, dbasis_nlat, sin_nlat, ntrunc, ierror = 0;
    int nmdim;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOOif", &ugrid_obj, &vgrid_obj, &pl_obj, &weights_obj, &basis_obj, &dbasis_obj, &sin_theta_obj, &ntrunc, &rsphere)) {
        return NULL;
    }
    if (as_real32_2d_fortran(ugrid_obj, &ugrid_arr, &ngptot_u, &nt_u) != 0) return NULL;
    if (as_real32_2d_fortran(vgrid_obj, &vgrid_arr, &ngptot_v, &nt_v) != 0) goto fail;
    if (ngptot_u != ngptot_v || nt_u != nt_v) {
        PyErr_SetString(PyExc_ValueError, "ugrid and vgrid must have the same shape");
        goto fail;
    }
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_1d_contig(weights_obj, &weights_arr, &nlat_w) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (as_real32_2d_fortran(dbasis_obj, &dbasis_arr, &dbasis_nmdim, &dbasis_nlat) != 0) goto fail;
    if (as_real32_1d_contig(sin_theta_obj, &sin_theta_arr, &sin_nlat) != 0) goto fail;
    if (nlat_w != nlat || basis_nlat != nlat || dbasis_nlat != nlat || sin_nlat != nlat || basis_nmdim != dbasis_nmdim) {
        PyErr_SetString(PyExc_ValueError, "pl, weights, basis, dbasis, and sin_theta dimensions must match");
        goto fail;
    }
    nmdim = basis_nmdim;
    dims[0] = nmdim;
    dims[1] = nt_u;
    vrtspec_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    divspec_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (vrtspec_arr == NULL || divspec_arr == NULL) goto fail;

    reduced_gaussian_getvrtdivspec_c(
        PyArray_DATA(ugrid_arr), PyArray_DATA(vgrid_arr), PyArray_DATA(pl_arr), PyArray_DATA(weights_arr),
        PyArray_DATA(basis_arr), PyArray_DATA(dbasis_arr), PyArray_DATA(sin_theta_arr),
        PyArray_DATA(vrtspec_arr), PyArray_DATA(divspec_arr), ngptot_u, nlat, ntrunc, nt_u, rsphere, &ierror);

    Py_DECREF(ugrid_arr);
    Py_DECREF(vgrid_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(weights_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(dbasis_arr);
    Py_DECREF(sin_theta_arr);
    return Py_BuildValue("NNi", vrtspec_arr, divspec_arr, ierror);

fail:
    Py_XDECREF(ugrid_arr);
    Py_XDECREF(vgrid_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(weights_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dbasis_arr);
    Py_XDECREF(sin_theta_arr);
    Py_XDECREF(vrtspec_arr);
    Py_XDECREF(divspec_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_getvrtspec(PyObject *self, PyObject *args) {
    PyObject *ugrid_obj = NULL, *vgrid_obj = NULL, *pl_obj = NULL, *weights_obj = NULL, *basis_obj = NULL, *dbasis_obj = NULL, *sin_theta_obj = NULL;
    PyArrayObject *ugrid_arr = NULL, *vgrid_arr = NULL, *pl_arr = NULL, *weights_arr = NULL, *basis_arr = NULL, *dbasis_arr = NULL, *sin_theta_arr = NULL, *vrtspec_arr = NULL;
    int ngptot_u, nt_u, ngptot_v, nt_v, nlat, nlat_w, basis_nmdim, basis_nlat, dbasis_nmdim, dbasis_nlat, sin_nlat, ntrunc, ierror = 0;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOOif", &ugrid_obj, &vgrid_obj, &pl_obj, &weights_obj, &basis_obj, &dbasis_obj, &sin_theta_obj, &ntrunc, &rsphere)) {
        return NULL;
    }
    if (as_real32_2d_fortran(ugrid_obj, &ugrid_arr, &ngptot_u, &nt_u) != 0) return NULL;
    if (as_real32_2d_fortran(vgrid_obj, &vgrid_arr, &ngptot_v, &nt_v) != 0) goto fail;
    if (ngptot_u != ngptot_v || nt_u != nt_v) {
        PyErr_SetString(PyExc_ValueError, "ugrid and vgrid must have the same shape");
        goto fail;
    }
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_1d_contig(weights_obj, &weights_arr, &nlat_w) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (as_real32_2d_fortran(dbasis_obj, &dbasis_arr, &dbasis_nmdim, &dbasis_nlat) != 0) goto fail;
    if (as_real32_1d_contig(sin_theta_obj, &sin_theta_arr, &sin_nlat) != 0) goto fail;
    if (nlat_w != nlat || basis_nlat != nlat || dbasis_nlat != nlat || sin_nlat != nlat || basis_nmdim != dbasis_nmdim) {
        PyErr_SetString(PyExc_ValueError, "pl, weights, basis, dbasis, and sin_theta dimensions must match");
        goto fail;
    }
    dims[0] = basis_nmdim;
    dims[1] = nt_u;
    vrtspec_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (vrtspec_arr == NULL) goto fail;

    reduced_gaussian_getvrtspec_c(
        PyArray_DATA(ugrid_arr), PyArray_DATA(vgrid_arr), PyArray_DATA(pl_arr), PyArray_DATA(weights_arr),
        PyArray_DATA(basis_arr), PyArray_DATA(dbasis_arr), PyArray_DATA(sin_theta_arr),
        PyArray_DATA(vrtspec_arr), ngptot_u, nlat, ntrunc, nt_u, rsphere, &ierror);

    Py_DECREF(ugrid_arr);
    Py_DECREF(vgrid_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(weights_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(dbasis_arr);
    Py_DECREF(sin_theta_arr);
    return Py_BuildValue("Ni", vrtspec_arr, ierror);

fail:
    Py_XDECREF(ugrid_arr);
    Py_XDECREF(vgrid_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(weights_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dbasis_arr);
    Py_XDECREF(sin_theta_arr);
    Py_XDECREF(vrtspec_arr);
    return NULL;
}

static PyObject *py_reduced_gaussian_getdivspec(PyObject *self, PyObject *args) {
    PyObject *ugrid_obj = NULL, *vgrid_obj = NULL, *pl_obj = NULL, *weights_obj = NULL, *basis_obj = NULL, *dbasis_obj = NULL, *sin_theta_obj = NULL;
    PyArrayObject *ugrid_arr = NULL, *vgrid_arr = NULL, *pl_arr = NULL, *weights_arr = NULL, *basis_arr = NULL, *dbasis_arr = NULL, *sin_theta_arr = NULL, *divspec_arr = NULL;
    int ngptot_u, nt_u, ngptot_v, nt_v, nlat, nlat_w, basis_nmdim, basis_nlat, dbasis_nmdim, dbasis_nlat, sin_nlat, ntrunc, ierror = 0;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOOif", &ugrid_obj, &vgrid_obj, &pl_obj, &weights_obj, &basis_obj, &dbasis_obj, &sin_theta_obj, &ntrunc, &rsphere)) {
        return NULL;
    }
    if (as_real32_2d_fortran(ugrid_obj, &ugrid_arr, &ngptot_u, &nt_u) != 0) return NULL;
    if (as_real32_2d_fortran(vgrid_obj, &vgrid_arr, &ngptot_v, &nt_v) != 0) goto fail;
    if (ngptot_u != ngptot_v || nt_u != nt_v) {
        PyErr_SetString(PyExc_ValueError, "ugrid and vgrid must have the same shape");
        goto fail;
    }
    if (as_int32_1d_contig(pl_obj, &pl_arr, &nlat) != 0) goto fail;
    if (as_real32_1d_contig(weights_obj, &weights_arr, &nlat_w) != 0) goto fail;
    if (as_real32_2d_fortran(basis_obj, &basis_arr, &basis_nmdim, &basis_nlat) != 0) goto fail;
    if (as_real32_2d_fortran(dbasis_obj, &dbasis_arr, &dbasis_nmdim, &dbasis_nlat) != 0) goto fail;
    if (as_real32_1d_contig(sin_theta_obj, &sin_theta_arr, &sin_nlat) != 0) goto fail;
    if (nlat_w != nlat || basis_nlat != nlat || dbasis_nlat != nlat || sin_nlat != nlat || basis_nmdim != dbasis_nmdim) {
        PyErr_SetString(PyExc_ValueError, "pl, weights, basis, dbasis, and sin_theta dimensions must match");
        goto fail;
    }
    dims[0] = basis_nmdim;
    dims[1] = nt_u;
    divspec_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (divspec_arr == NULL) goto fail;

    reduced_gaussian_getdivspec_c(
        PyArray_DATA(ugrid_arr), PyArray_DATA(vgrid_arr), PyArray_DATA(pl_arr), PyArray_DATA(weights_arr),
        PyArray_DATA(basis_arr), PyArray_DATA(dbasis_arr), PyArray_DATA(sin_theta_arr),
        PyArray_DATA(divspec_arr), ngptot_u, nlat, ntrunc, nt_u, rsphere, &ierror);

    Py_DECREF(ugrid_arr);
    Py_DECREF(vgrid_arr);
    Py_DECREF(pl_arr);
    Py_DECREF(weights_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(dbasis_arr);
    Py_DECREF(sin_theta_arr);
    return Py_BuildValue("Ni", divspec_arr, ierror);

fail:
    Py_XDECREF(ugrid_arr);
    Py_XDECREF(vgrid_arr);
    Py_XDECREF(pl_arr);
    Py_XDECREF(weights_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(dbasis_arr);
    Py_XDECREF(sin_theta_arr);
    Py_XDECREF(divspec_arr);
    return NULL;
}

PyMethodDef reduced_module_methods[] = {
    {"gaqd", py_reduced_gaqd, METH_VARARGS, "Compute Gaussian latitudes and weights."},
    {"lap", py_reduced_lap, METH_VARARGS, "Apply the spherical Laplacian to reduced spectra."},
    {"invlap", py_reduced_invlap, METH_VARARGS, "Apply the inverse spherical Laplacian to reduced spectra."},
    {"multsmoothfact", py_reduced_multsmoothfact, METH_VARARGS, "Apply spectral smoothing factors."},
    {"reduced_gaussian_legendre_basis", py_reduced_gaussian_legendre_basis, METH_VARARGS, "Build reduced Gaussian Legendre basis."},
    {"reduced_gaussian_legendre_derivative_from_basis", py_reduced_gaussian_legendre_derivative_from_basis, METH_VARARGS, "Build reduced Gaussian Legendre derivative basis from basis."},
    {"reduced_gaussian_legendre_derivative_basis", py_reduced_gaussian_legendre_derivative_basis, METH_VARARGS, "Build reduced Gaussian Legendre derivative basis."},
    {"reduced_gaussian_grdtospec", py_reduced_gaussian_grdtospec, METH_VARARGS, "Analyze packed reduced Gaussian grid to spectra."},
    {"reduced_gaussian_spectogrd", py_reduced_gaussian_spectogrd, METH_VARARGS, "Synthesize spectra to packed reduced Gaussian grid."},
    {"reduced_gaussian_spectogrd_pair", py_reduced_gaussian_spectogrd_pair, METH_VARARGS, "Synthesize paired spectra to packed reduced Gaussian grids."},
    {"reduced_gaussian_getgrad", py_reduced_gaussian_getgrad, METH_VARARGS, "Synthesize reduced Gaussian scalar gradient."},
    {"reduced_gaussian_getgrad_pair", py_reduced_gaussian_getgrad_pair, METH_VARARGS, "Synthesize paired reduced Gaussian scalar gradients."},
    {"reduced_gaussian_getvrtdivspec", py_reduced_gaussian_getvrtdivspec, METH_VARARGS, "Analyze packed reduced Gaussian winds to vorticity/divergence spectra."},
    {"reduced_gaussian_getvrtspec", py_reduced_gaussian_getvrtspec, METH_VARARGS, "Analyze packed reduced Gaussian winds to vorticity spectra."},
    {"reduced_gaussian_getdivspec", py_reduced_gaussian_getdivspec, METH_VARARGS, "Analyze packed reduced Gaussian winds to divergence spectra."},
    {NULL, NULL, 0, NULL},
};
