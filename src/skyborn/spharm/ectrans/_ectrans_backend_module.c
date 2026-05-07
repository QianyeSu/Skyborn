#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>
#include <math.h>

#define SKYBORN_ECTRANS_SETUP_CAPSULE "skyborn.spharm._ectrans_setup"

typedef struct {
    int ndgl;
    int ngptot;
    int max_nlon;
    int *nloen;
    int *lat_offsets;
    double *mu;
    double *weights;
} EctransSetup;

void validate_nloen(
    int ndgl,
    const int *nloen,
    int *ngptot,
    int *ierror
);

void build_grid_layout(
    int ndgl,
    const int *nloen,
    int ngptot,
    int *lat_offsets,
    int *max_nlon,
    int *ierror
);

void scalar_analysis_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int ntrunc,
    int nt,
    const double *datagrid,
    double *dataspec_packed,
    int *ierror
);

void scalar_synthesis_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int ntrunc,
    int nt,
    const double *dataspec_packed,
    double *datagrid,
    int *ierror
);

void vrtdiv_analysis_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int ntrunc,
    int nt,
    const double *ugrid,
    const double *vgrid,
    double *vrtspec_r,
    double *vrtspec_i,
    double *divspec_r,
    double *divspec_i,
    int *ierror
);

void uv_synthesis_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int ntrunc,
    int nt,
    const double *vrtspec_r,
    const double *vrtspec_i,
    const double *divspec_r,
    const double *divspec_i,
    double *ugrid,
    double *vgrid,
    int *ierror
);

void gradient_synthesis_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int ntrunc,
    int nt,
    const double *chispec_r,
    const double *chispec_i,
    double *ugrad,
    double *vgrad,
    int *ierror
);

static int infer_ntrunc_from_ncoeff(npy_intp ncoeff) {
    long long value = (long long) ncoeff;
    int ntrunc = (int) (-1.5 + 0.5 * sqrt(1.0 + 8.0 * (double) value));
    long long check = ((long long) (ntrunc + 1) * (long long) (ntrunc + 2)) / 2;
    return check == value ? ntrunc : -1;
}

static PyArrayObject *require_array(PyObject *obj, int typenum) {
    return (PyArrayObject *) PyArray_FROM_OTF(
        obj,
        typenum,
        NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS
    );
}

static void free_setup_contents(EctransSetup *setup) {
    if (setup == NULL) {
        return;
    }
    if (setup->nloen != NULL) {
        PyMem_Free(setup->nloen);
        setup->nloen = NULL;
    }
    if (setup->lat_offsets != NULL) {
        PyMem_Free(setup->lat_offsets);
        setup->lat_offsets = NULL;
    }
    if (setup->mu != NULL) {
        PyMem_Free(setup->mu);
        setup->mu = NULL;
    }
    if (setup->weights != NULL) {
        PyMem_Free(setup->weights);
        setup->weights = NULL;
    }
    setup->ndgl = 0;
    setup->ngptot = 0;
    setup->max_nlon = 0;
}

static void setup_capsule_destructor(PyObject *capsule) {
    EctransSetup *setup = (EctransSetup *) PyCapsule_GetPointer(
        capsule,
        SKYBORN_ECTRANS_SETUP_CAPSULE
    );
    if (setup == NULL) {
        PyErr_Clear();
        return;
    }
    free_setup_contents(setup);
    PyMem_Free(setup);
}

static int validate_nloen_array(PyArrayObject *nloen_arr, int *ndgl, int *ngptot) {
    int ierror = 0;
    if (PyArray_NDIM(nloen_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "nloen must be a 1D int32 array");
        return 0;
    }
    *ndgl = (int) PyArray_DIM(nloen_arr, 0);
    validate_nloen(
        *ndgl,
        (const int *) PyArray_DATA(nloen_arr),
        ngptot,
        &ierror
    );
    if (ierror != 0) {
        PyErr_Format(PyExc_ValueError, "invalid reduced Gaussian nloen metadata (ierror=%d)", ierror);
        return 0;
    }
    return 1;
}

static int get_setup_from_capsule(PyObject *handle_obj, EctransSetup **setup_out) {
    EctransSetup *setup = NULL;

    if (!PyCapsule_IsValid(handle_obj, SKYBORN_ECTRANS_SETUP_CAPSULE)) {
        PyErr_SetString(
            PyExc_TypeError,
            "setup handle must be a valid ectrans setup capsule"
        );
        return 0;
    }

    setup = (EctransSetup *) PyCapsule_GetPointer(
        handle_obj,
        SKYBORN_ECTRANS_SETUP_CAPSULE
    );
    if (setup == NULL) {
        return 0;
    }
    if (setup->nloen == NULL || setup->ndgl <= 0 || setup->ngptot <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "ectrans setup handle is closed");
        return 0;
    }

    *setup_out = setup;
    return 1;
}

static int flatten_grid_array(
    PyArrayObject *grid_arr,
    int ngptot,
    int *nt_out,
    int *rank_out,
    npy_intp **dims_out
) {
    npy_intp total_size = PyArray_SIZE(grid_arr);
    if (PyArray_NDIM(grid_arr) < 1) {
        PyErr_SetString(PyExc_ValueError, "grid array must be at least rank 1");
        return 0;
    }
    if (PyArray_DIM(grid_arr, 0) != (npy_intp) ngptot) {
        PyErr_Format(
            PyExc_ValueError,
            "grid leading dimension must be %d, got %lld",
            ngptot,
            (long long) PyArray_DIM(grid_arr, 0)
        );
        return 0;
    }
    if (total_size % ngptot != 0) {
        PyErr_SetString(PyExc_ValueError, "grid array size is not divisible by ngptot");
        return 0;
    }
    *nt_out = (int) (total_size / ngptot);
    *rank_out = PyArray_NDIM(grid_arr);
    *dims_out = PyArray_DIMS(grid_arr);
    return 1;
}

static int flatten_spectral_array(
    PyArrayObject *spec_arr,
    int max_ntrunc,
    int *ntrunc_out,
    int *nt_out,
    int *rank_out,
    npy_intp **dims_out
) {
    npy_intp total_size = PyArray_SIZE(spec_arr);
    int ntrunc;

    if (PyArray_NDIM(spec_arr) < 1) {
        PyErr_SetString(PyExc_ValueError, "spectral array must be at least rank 1");
        return 0;
    }

    ntrunc = infer_ntrunc_from_ncoeff(PyArray_DIM(spec_arr, 0));
    if (ntrunc < 0 || ntrunc > max_ntrunc) {
        PyErr_Format(
            PyExc_ValueError,
            "spectral leading dimension does not define a valid ntrunc <= %d",
            max_ntrunc
        );
        return 0;
    }

    *ntrunc_out = ntrunc;
    *nt_out = (int) (total_size / PyArray_DIM(spec_arr, 0));
    *rank_out = PyArray_NDIM(spec_arr);
    *dims_out = PyArray_DIMS(spec_arr);
    return 1;
}

static PyObject *backend_validate_nloen(PyObject *self, PyObject *args) {
    PyObject *nloen_obj = NULL;
    PyArrayObject *nloen_arr = NULL;
    int ndgl = 0;
    int ngptot = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "O", &nloen_obj)) {
        return NULL;
    }

    nloen_arr = require_array(nloen_obj, NPY_INT32);
    if (nloen_arr == NULL) {
        return NULL;
    }
    if (!validate_nloen_array(nloen_arr, &ndgl, &ngptot)) {
        Py_DECREF(nloen_arr);
        return NULL;
    }

    Py_DECREF(nloen_arr);
    return Py_BuildValue("(ii)", ndgl, ngptot);
}

static PyObject *backend_create_setup(PyObject *self, PyObject *args) {
    PyObject *nloen_obj = NULL;
    PyObject *mu_obj = NULL;
    PyObject *weights_obj = NULL;
    PyArrayObject *nloen_arr = NULL;
    PyArrayObject *mu_arr = NULL;
    PyArrayObject *weights_arr = NULL;
    EctransSetup *setup = NULL;
    PyObject *capsule = NULL;
    int ndgl = 0;
    int ngptot = 0;
    int ierror = 0;
    npy_intp ilat = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOO", &nloen_obj, &mu_obj, &weights_obj)) {
        return NULL;
    }

    nloen_arr = require_array(nloen_obj, NPY_INT32);
    mu_arr = require_array(mu_obj, NPY_FLOAT64);
    weights_arr = require_array(weights_obj, NPY_FLOAT64);
    if (nloen_arr == NULL || mu_arr == NULL || weights_arr == NULL) {
        return NULL;
    }
    if (!validate_nloen_array(nloen_arr, &ndgl, &ngptot)) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }
    if (PyArray_NDIM(mu_arr) != 1 || PyArray_DIM(mu_arr, 0) != (npy_intp) ndgl) {
        PyErr_SetString(PyExc_ValueError, "mu must be a 1D float64 array with length ndgl");
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }
    if (PyArray_NDIM(weights_arr) != 1 || PyArray_DIM(weights_arr, 0) != (npy_intp) ndgl) {
        PyErr_SetString(
            PyExc_ValueError,
            "weights must be a 1D float64 array with length ndgl"
        );
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        return NULL;
    }

    setup = (EctransSetup *) PyMem_Malloc(sizeof(EctransSetup));
    if (setup == NULL) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        return PyErr_NoMemory();
    }
    setup->ndgl = ndgl;
    setup->ngptot = ngptot;
    setup->max_nlon = 0;
    setup->lat_offsets = NULL;
    setup->mu = NULL;
    setup->weights = NULL;
    setup->nloen = (int *) PyMem_Malloc((size_t) ndgl * sizeof(int));
    if (setup->nloen == NULL) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        PyMem_Free(setup);
        return PyErr_NoMemory();
    }

    for (ilat = 0; ilat < (npy_intp) ndgl; ++ilat) {
        setup->nloen[ilat] = *(int *) PyArray_GETPTR1(nloen_arr, ilat);
    }
    setup->lat_offsets = (int *) PyMem_Malloc((size_t) ndgl * sizeof(int));
    if (setup->lat_offsets == NULL) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        free_setup_contents(setup);
        PyMem_Free(setup);
        return PyErr_NoMemory();
    }
    setup->mu = (double *) PyMem_Malloc((size_t) ndgl * sizeof(double));
    if (setup->mu == NULL) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        free_setup_contents(setup);
        PyMem_Free(setup);
        return PyErr_NoMemory();
    }
    setup->weights = (double *) PyMem_Malloc((size_t) ndgl * sizeof(double));
    if (setup->weights == NULL) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        free_setup_contents(setup);
        PyMem_Free(setup);
        return PyErr_NoMemory();
    }
    build_grid_layout(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        setup->lat_offsets,
        &setup->max_nlon,
        &ierror
    );
    if (ierror != 0) {
        Py_DECREF(nloen_arr);
        Py_DECREF(mu_arr);
        Py_DECREF(weights_arr);
        free_setup_contents(setup);
        PyMem_Free(setup);
        PyErr_Format(
            PyExc_ValueError,
            "failed to build reduced-grid layout (ierror=%d)",
            ierror
        );
        return NULL;
    }
    for (ilat = 0; ilat < (npy_intp) ndgl; ++ilat) {
        setup->mu[ilat] = *(double *) PyArray_GETPTR1(mu_arr, ilat);
        setup->weights[ilat] = *(double *) PyArray_GETPTR1(weights_arr, ilat);
    }

    capsule = PyCapsule_New(
        (void *) setup,
        SKYBORN_ECTRANS_SETUP_CAPSULE,
        setup_capsule_destructor
    );
    Py_DECREF(nloen_arr);
    Py_DECREF(mu_arr);
    Py_DECREF(weights_arr);
    if (capsule == NULL) {
        free_setup_contents(setup);
        PyMem_Free(setup);
        return NULL;
    }

    return capsule;
}

static PyObject *backend_describe_setup(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    EctransSetup *setup = NULL;
    npy_intp dims[1];
    PyArrayObject *nloen_arr = NULL;
    PyArrayObject *lat_offsets_arr = NULL;
    PyArrayObject *mu_arr = NULL;
    PyArrayObject *weights_arr = NULL;
    npy_intp ilat = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "O", &handle_obj)) {
        return NULL;
    }
    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }

    dims[0] = (npy_intp) setup->ndgl;
    nloen_arr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT32);
    lat_offsets_arr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT32);
    mu_arr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    weights_arr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_FLOAT64);
    if (nloen_arr == NULL || lat_offsets_arr == NULL || mu_arr == NULL || weights_arr == NULL) {
        Py_XDECREF(nloen_arr);
        Py_XDECREF(lat_offsets_arr);
        Py_XDECREF(mu_arr);
        Py_XDECREF(weights_arr);
        return NULL;
    }

    for (ilat = 0; ilat < dims[0]; ++ilat) {
        *(int *) PyArray_GETPTR1(nloen_arr, ilat) = setup->nloen[ilat];
        *(int *) PyArray_GETPTR1(lat_offsets_arr, ilat) = setup->lat_offsets[ilat];
        *(double *) PyArray_GETPTR1(mu_arr, ilat) = setup->mu[ilat];
        *(double *) PyArray_GETPTR1(weights_arr, ilat) = setup->weights[ilat];
    }

    return Py_BuildValue(
        "{s:i,s:i,s:i,s:N,s:N,s:N,s:N}",
        "ndgl",
        setup->ndgl,
        "ngptot",
        setup->ngptot,
        "max_nlon",
        setup->max_nlon,
        "nloen",
        nloen_arr,
        "lat_offsets",
        lat_offsets_arr,
        "mu",
        mu_arr,
        "weights",
        weights_arr
    );
}

static PyObject *backend_destroy_setup(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    EctransSetup *setup = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "O", &handle_obj)) {
        return NULL;
    }

    if (!PyCapsule_IsValid(handle_obj, SKYBORN_ECTRANS_SETUP_CAPSULE)) {
        PyErr_SetString(
            PyExc_TypeError,
            "setup handle must be a valid ectrans setup capsule"
        );
        return NULL;
    }

    setup = (EctransSetup *) PyCapsule_GetPointer(
        handle_obj,
        SKYBORN_ECTRANS_SETUP_CAPSULE
    );
    if (setup == NULL) {
        return NULL;
    }

    free_setup_contents(setup);
    Py_RETURN_NONE;
}

static PyObject *backend_scalar_analysis(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *grid_obj = NULL;
    int ntrunc = -1;
    EctransSetup *setup = NULL;
    PyArrayObject *grid_arr = NULL;
    PyArrayObject *out_arr = NULL;
    npy_intp *grid_dims = NULL;
    npy_intp *out_dims = NULL;
    int nt = 0, grid_rank = 0, ierror = 0;
    npy_intp nspec2 = 0;
    npy_intp idx;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOi", &handle_obj, &grid_obj, &ntrunc)) {
        return NULL;
    }

    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }
    grid_arr = require_array(grid_obj, NPY_FLOAT64);
    if (grid_arr == NULL) {
        goto fail;
    }
    if (!flatten_grid_array(grid_arr, setup->ngptot, &nt, &grid_rank, &grid_dims)) {
        goto fail;
    }
    if (ntrunc < 0 || ntrunc > setup->ndgl - 1) {
        PyErr_Format(
            PyExc_ValueError,
            "ntrunc must be between 0 and %d",
            setup->ndgl - 1
        );
        goto fail;
    }

    nspec2 = (npy_intp) (ntrunc + 1) * (npy_intp) (ntrunc + 2);
    out_dims = PyMem_Malloc((size_t) grid_rank * sizeof(npy_intp));
    if (out_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    out_dims[0] = nspec2;
    for (idx = 1; idx < grid_rank; ++idx) {
        out_dims[idx] = grid_dims[idx];
    }

    out_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_FLOAT64, 0);
    if (out_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    scalar_analysis_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        ntrunc,
        nt,
        (const double *) PyArray_DATA(grid_arr),
        (double *) PyArray_DATA(out_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    PyMem_Free(out_dims);
    Py_DECREF(grid_arr);
    return Py_BuildValue("(Ni)", out_arr, ierror);

fail:
    if (out_dims != NULL) {
        PyMem_Free(out_dims);
    }
    Py_XDECREF(grid_arr);
    Py_XDECREF(out_arr);
    return NULL;
}

static PyObject *backend_scalar_synthesis(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *spec_obj = NULL;
    EctransSetup *setup = NULL;
    PyArrayObject *spec_arr = NULL;
    PyArrayObject *grid_arr = NULL;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp *grid_dims = NULL;
    npy_intp idx;

    (void) self;

    if (!PyArg_ParseTuple(args, "OO", &handle_obj, &spec_obj)) {
        return NULL;
    }

    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }
    spec_arr = require_array(spec_obj, NPY_FLOAT64);
    if (spec_arr == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(spec_arr) < 1) {
        PyErr_SetString(PyExc_ValueError, "dataspec must be at least rank 1");
        goto fail;
    }
    if ((PyArray_DIM(spec_arr, 0) % 2) != 0) {
        PyErr_SetString(PyExc_ValueError, "dataspec leading dimension must be even");
        goto fail;
    }
    ntrunc = infer_ntrunc_from_ncoeff(PyArray_DIM(spec_arr, 0) / 2);
    if (ntrunc < 0 || ntrunc > setup->ndgl - 1) {
        PyErr_SetString(PyExc_ValueError, "dataspec leading dimension does not define a valid packed ntrunc");
        goto fail;
    }
    spec_rank = PyArray_NDIM(spec_arr);
    spec_dims = PyArray_DIMS(spec_arr);
    nt = (int) (PyArray_SIZE(spec_arr) / PyArray_DIM(spec_arr, 0));

    grid_dims = PyMem_Malloc((size_t) spec_rank * sizeof(npy_intp));
    if (grid_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    grid_dims[0] = setup->ngptot;
    for (idx = 1; idx < spec_rank; ++idx) {
        grid_dims[idx] = spec_dims[idx];
    }

    grid_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, grid_dims, NPY_FLOAT64, 0);
    if (grid_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    scalar_synthesis_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        ntrunc,
        nt,
        (const double *) PyArray_DATA(spec_arr),
        (double *) PyArray_DATA(grid_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    PyMem_Free(grid_dims);
    Py_DECREF(spec_arr);
    return Py_BuildValue("(Ni)", grid_arr, ierror);

fail:
    if (grid_dims != NULL) {
        PyMem_Free(grid_dims);
    }
    Py_XDECREF(spec_arr);
    Py_XDECREF(grid_arr);
    return NULL;
}

static PyObject *backend_vrtdiv_analysis(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *ugrid_obj = NULL;
    PyObject *vgrid_obj = NULL;
    int ntrunc = -1;
    EctransSetup *setup = NULL;
    PyArrayObject *ugrid_arr = NULL;
    PyArrayObject *vgrid_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *vrt_out_arr = NULL;
    PyArrayObject *div_out_arr = NULL;
    npy_intp *grid_dims = NULL;
    npy_intp *out_dims = NULL;
    int nt = 0, grid_rank = 0, ierror = 0;
    npy_intp ncoeff = 0;
    npy_intp idx;
    npy_complex128 *vrt_out_data = NULL;
    npy_complex128 *div_out_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOi", &handle_obj, &ugrid_obj, &vgrid_obj, &ntrunc)) {
        return NULL;
    }

    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }
    ugrid_arr = require_array(ugrid_obj, NPY_FLOAT64);
    vgrid_arr = require_array(vgrid_obj, NPY_FLOAT64);
    if (ugrid_arr == NULL || vgrid_arr == NULL) {
        goto fail;
    }
    if (!flatten_grid_array(ugrid_arr, setup->ngptot, &nt, &grid_rank, &grid_dims)) {
        goto fail;
    }
    if (!flatten_grid_array(vgrid_arr, setup->ngptot, &nt, &grid_rank, &grid_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(ugrid_arr) != PyArray_NDIM(vgrid_arr)) {
        PyErr_SetString(PyExc_ValueError, "ugrid and vgrid must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(ugrid_arr), PyArray_DIMS(vgrid_arr), PyArray_NDIM(ugrid_arr))) {
        PyErr_SetString(PyExc_ValueError, "ugrid and vgrid must have the same shape");
        goto fail;
    }
    if (ntrunc < 0 || ntrunc > setup->ndgl - 1) {
        PyErr_Format(
            PyExc_ValueError,
            "ntrunc must be between 0 and %d",
            setup->ndgl - 1
        );
        goto fail;
    }

    ncoeff = ((npy_intp) (ntrunc + 1) * (npy_intp) (ntrunc + 2)) / 2;
    out_dims = PyMem_Malloc((size_t) grid_rank * sizeof(npy_intp));
    if (out_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    out_dims[0] = ncoeff;
    for (idx = 1; idx < grid_rank; ++idx) {
        out_dims[idx] = grid_dims[idx];
    }

    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_FLOAT64, 0);
    vrt_out_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_COMPLEX128, 0);
    div_out_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank, out_dims, NPY_COMPLEX128, 0);
    if (
        vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        vrt_out_arr == NULL || div_out_arr == NULL
    ) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    vrtdiv_analysis_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        ntrunc,
        nt,
        (const double *) PyArray_DATA(ugrid_arr),
        (const double *) PyArray_DATA(vgrid_arr),
        (double *) PyArray_DATA(vrt_real_arr),
        (double *) PyArray_DATA(vrt_imag_arr),
        (double *) PyArray_DATA(div_real_arr),
        (double *) PyArray_DATA(div_imag_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    vrt_out_data = (npy_complex128 *) PyArray_DATA(vrt_out_arr);
    div_out_data = (npy_complex128 *) PyArray_DATA(div_out_arr);
    vrt_real_data = (double *) PyArray_DATA(vrt_real_arr);
    vrt_imag_data = (double *) PyArray_DATA(vrt_imag_arr);
    div_real_data = (double *) PyArray_DATA(div_real_arr);
    div_imag_data = (double *) PyArray_DATA(div_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(vrt_out_arr); ++idx) {
        vrt_out_data[idx] = vrt_real_data[idx] + vrt_imag_data[idx] * I;
        div_out_data[idx] = div_real_data[idx] + div_imag_data[idx] * I;
    }

    PyMem_Free(out_dims);
    Py_DECREF(ugrid_arr);
    Py_DECREF(vgrid_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    return Py_BuildValue("(NNi)", vrt_out_arr, div_out_arr, ierror);

fail:
    if (out_dims != NULL) {
        PyMem_Free(out_dims);
    }
    Py_XDECREF(ugrid_arr);
    Py_XDECREF(vgrid_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(vrt_out_arr);
    Py_XDECREF(div_out_arr);
    return NULL;
}

static PyObject *backend_uv_synthesis(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *vrt_obj = NULL;
    PyObject *div_obj = NULL;
    EctransSetup *setup = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *div_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *ugrid_arr = NULL;
    PyArrayObject *vgrid_arr = NULL;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, div_ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp *grid_dims = NULL;
    npy_intp idx;
    npy_complex128 *vrt_data = NULL;
    npy_complex128 *div_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOO", &handle_obj, &vrt_obj, &div_obj)) {
        return NULL;
    }

    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }
    vrt_arr = require_array(vrt_obj, NPY_COMPLEX128);
    div_arr = require_array(div_obj, NPY_COMPLEX128);
    if (vrt_arr == NULL || div_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(vrt_arr, setup->ndgl - 1, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(div_arr, setup->ndgl - 1, &div_ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(vrt_arr) != PyArray_NDIM(div_arr)) {
        PyErr_SetString(PyExc_ValueError, "vrtspec and divspec must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(vrt_arr), PyArray_DIMS(div_arr), PyArray_NDIM(vrt_arr))) {
        PyErr_SetString(PyExc_ValueError, "vrtspec and divspec must have the same shape");
        goto fail;
    }
    if (ntrunc != div_ntrunc) {
        PyErr_SetString(PyExc_ValueError, "vrtspec and divspec must use the same ntrunc");
        goto fail;
    }

    grid_dims = PyMem_Malloc((size_t) spec_rank * sizeof(npy_intp));
    if (grid_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    grid_dims[0] = setup->ngptot;
    for (idx = 1; idx < spec_rank; ++idx) {
        grid_dims[idx] = spec_dims[idx];
    }

    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    ugrid_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, grid_dims, NPY_FLOAT64, 0);
    vgrid_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, grid_dims, NPY_FLOAT64, 0);
    if (
        vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        ugrid_arr == NULL || vgrid_arr == NULL
    ) {
        goto fail;
    }

    vrt_data = (npy_complex128 *) PyArray_DATA(vrt_arr);
    div_data = (npy_complex128 *) PyArray_DATA(div_arr);
    vrt_real_data = (double *) PyArray_DATA(vrt_real_arr);
    vrt_imag_data = (double *) PyArray_DATA(vrt_imag_arr);
    div_real_data = (double *) PyArray_DATA(div_real_arr);
    div_imag_data = (double *) PyArray_DATA(div_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(vrt_arr); ++idx) {
        vrt_real_data[idx] = creal(vrt_data[idx]);
        vrt_imag_data[idx] = cimag(vrt_data[idx]);
        div_real_data[idx] = creal(div_data[idx]);
        div_imag_data[idx] = cimag(div_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    uv_synthesis_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        ntrunc,
        nt,
        (const double *) PyArray_DATA(vrt_real_arr),
        (const double *) PyArray_DATA(vrt_imag_arr),
        (const double *) PyArray_DATA(div_real_arr),
        (const double *) PyArray_DATA(div_imag_arr),
        (double *) PyArray_DATA(ugrid_arr),
        (double *) PyArray_DATA(vgrid_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    PyMem_Free(grid_dims);
    Py_DECREF(vrt_arr);
    Py_DECREF(div_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    return Py_BuildValue("(NNi)", ugrid_arr, vgrid_arr, ierror);

fail:
    if (grid_dims != NULL) {
        PyMem_Free(grid_dims);
    }
    Py_XDECREF(vrt_arr);
    Py_XDECREF(div_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(ugrid_arr);
    Py_XDECREF(vgrid_arr);
    return NULL;
}

static PyObject *backend_gradient_synthesis(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *chi_obj = NULL;
    EctransSetup *setup = NULL;
    PyArrayObject *chi_arr = NULL;
    PyArrayObject *chi_real_arr = NULL;
    PyArrayObject *chi_imag_arr = NULL;
    PyArrayObject *ugrad_arr = NULL;
    PyArrayObject *vgrad_arr = NULL;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp *grid_dims = NULL;
    npy_intp idx;
    npy_complex128 *chi_data = NULL;
    double *chi_real_data = NULL;
    double *chi_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OO", &handle_obj, &chi_obj)) {
        return NULL;
    }

    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }
    chi_arr = require_array(chi_obj, NPY_COMPLEX128);
    if (chi_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(chi_arr, setup->ndgl - 1, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }

    grid_dims = PyMem_Malloc((size_t) spec_rank * sizeof(npy_intp));
    if (grid_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    grid_dims[0] = setup->ngptot;
    for (idx = 1; idx < spec_rank; ++idx) {
        grid_dims[idx] = spec_dims[idx];
    }

    chi_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    chi_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    ugrad_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, grid_dims, NPY_FLOAT64, 0);
    vgrad_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, grid_dims, NPY_FLOAT64, 0);
    if (chi_real_arr == NULL || chi_imag_arr == NULL || ugrad_arr == NULL || vgrad_arr == NULL) {
        goto fail;
    }

    chi_data = (npy_complex128 *) PyArray_DATA(chi_arr);
    chi_real_data = (double *) PyArray_DATA(chi_real_arr);
    chi_imag_data = (double *) PyArray_DATA(chi_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(chi_arr); ++idx) {
        chi_real_data[idx] = creal(chi_data[idx]);
        chi_imag_data[idx] = cimag(chi_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    gradient_synthesis_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        ntrunc,
        nt,
        (const double *) PyArray_DATA(chi_real_arr),
        (const double *) PyArray_DATA(chi_imag_arr),
        (double *) PyArray_DATA(ugrad_arr),
        (double *) PyArray_DATA(vgrad_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    PyMem_Free(grid_dims);
    Py_DECREF(chi_arr);
    Py_DECREF(chi_real_arr);
    Py_DECREF(chi_imag_arr);
    return Py_BuildValue("(NNi)", ugrad_arr, vgrad_arr, ierror);

fail:
    if (grid_dims != NULL) {
        PyMem_Free(grid_dims);
    }
    Py_XDECREF(chi_arr);
    Py_XDECREF(chi_real_arr);
    Py_XDECREF(chi_imag_arr);
    Py_XDECREF(ugrad_arr);
    Py_XDECREF(vgrad_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"validate_nloen", backend_validate_nloen, METH_VARARGS, "Validate reduced-grid longitude counts."},
    {"create_setup", backend_create_setup, METH_VARARGS, "Create an experimental ectrans setup handle."},
    {"describe_setup", backend_describe_setup, METH_VARARGS, "Describe an experimental ectrans setup handle."},
    {"destroy_setup", backend_destroy_setup, METH_VARARGS, "Destroy an experimental ectrans setup handle."},
    {"scalar_analysis_stub", backend_scalar_analysis, METH_VARARGS, "Call the native scalar-analysis stub."},
    {"scalar_synthesis_stub", backend_scalar_synthesis, METH_VARARGS, "Call the native scalar-synthesis stub."},
    {"vrtdiv_analysis_stub", backend_vrtdiv_analysis, METH_VARARGS, "Call the native vector-analysis stub."},
    {"uv_synthesis_stub", backend_uv_synthesis, METH_VARARGS, "Call the native wind-synthesis stub."},
    {"gradient_synthesis_stub", backend_gradient_synthesis, METH_VARARGS, "Call the native gradient-synthesis stub."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_ectrans_backend",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit__ectrans_backend(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
