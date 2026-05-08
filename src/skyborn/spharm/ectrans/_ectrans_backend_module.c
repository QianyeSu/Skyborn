#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#define SKYBORN_ECTRANS_SETUP_CAPSULE "skyborn.spharm._ectrans_setup"

typedef struct {
    int ndgl;
    int ngptot;
    int max_nlon;
    int *nloen;
    int *lat_offsets;
    double *mu;
    double *weights;
    double rsphere;
} EctransSetup;

#define SKYBORN_UV_TO_VORDIV_BLOCK_CAPSULE "skyborn.spharm._ectrans_uv_to_vordiv_block_setup"
#define SKYBORN_UV_TO_VORDIV_BLOCK_CACHE_MAX_BYTES (96ULL * 1024ULL * 1024ULL)

typedef struct {
    int ntrunc;
    double rsphere;
    size_t basis_used_entries;
    size_t active_used_entries;
    size_t *basis_offsets;
    int *active_counts;
    size_t *active_offsets;
    int *active_columns;
    double *basis_real;
    double *basis_imag;
} UvToVordivBlockSetup;


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
    const double *mu,
    const double *weights,
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
    const double *mu,
    int ngptot,
    int ntrunc,
    int nt,
    const double *dataspec_packed,
    double *datagrid,
    int *ierror
);

void scalar_fourier_stub(
    int ndgl,
    const int *nloen,
    int ngptot,
    int mmax,
    int nt,
    const double *datagrid,
    double *fourier_real,
    double *fourier_imag,
    int *ierror
);

void scalar_block_solve_stub(
    int ndgl,
    int nblock,
    int nt,
    const double *weights,
    const double *basis_real,
    const double *basis_imag,
    const double *observed_real,
    const double *observed_imag,
    double *solution_real,
    double *solution_imag,
    int *ierror
);

void weighted_block_solve_stub(
    int nrow,
    int nblock,
    int nt,
    const double *weights,
    const double *basis_real,
    const double *basis_imag,
    const double *observed_real,
    const double *observed_imag,
    double *solution_real,
    double *solution_imag,
    int *ierror
);

void vrtdiv_analysis_stub(
    int ndgl,
    const int *nloen,
    const double *weights,
    int ngptot,
    double rsphere,
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
    double rsphere,
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
    double rsphere,
    int ntrunc,
    int nt,
    const double *chispec_r,
    const double *chispec_i,
    double *ugrad,
    double *vgrad,
    int *ierror
);

void vordiv_to_uv(
    int ntrunc,
    int nt,
    double rsphere,
    const double *vrtspec_r,
    const double *vrtspec_i,
    const double *divspec_r,
    const double *divspec_i,
    double *uspec_r,
    double *uspec_i,
    double *vspec_r,
    double *vspec_i,
    int *ierror
);

void uv_to_vordiv(
    int ntrunc,
    int nt,
    double rsphere,
    const double *uspec_r,
    const double *uspec_i,
    const double *vspec_r,
    const double *vspec_i,
    double *vrtspec_r,
    double *vrtspec_i,
    double *divspec_r,
    double *divspec_i,
    int *ierror
);

void ldfou2_uv_scaling(
    int ntrunc,
    int km,
    int kf_uv,
    double rsphere,
    const double *paia_in,
    const double *psia_in,
    double *paia_out,
    double *psia_out,
    int *ierror
);

void ledir_dgemm(
    int ntrunc,
    int km,
    int kfc,
    int kdglu,
    const double *paia,
    const double *psia,
    const double *rpnma,
    const double *rpnms,
    const double *pw,
    double *poa1,
    int *ierror
);

void prfi1b_uv_block(
    int ntrunc,
    int km,
    int nt,
    double rsphere,
    const double *uspec_r,
    const double *uspec_i,
    const double *vspec_r,
    const double *vspec_i,
    double *poa1_out,
    int *ierror
);

void vd2uv_uv_block(
    int ntrunc,
    int km,
    int nt,
    double rsphere,
    const double *vrtspec_r,
    const double *vrtspec_i,
    const double *divspec_r,
    const double *divspec_i,
    double *poa1_out,
    int *ierror
);

void poa1_to_vordiv(
    int ntrunc,
    int km,
    int kf_uv,
    double rsphere,
    const double *poa1_in,
    double *vrtspec_r,
    double *vrtspec_i,
    double *divspec_r,
    double *divspec_i,
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
    setup->rsphere = 0.0;
    setup->ndgl = 0;
    setup->ngptot = 0;
    setup->max_nlon = 0;
}

static void free_uv_to_vordiv_block_setup(UvToVordivBlockSetup *setup) {
    if (setup == NULL) {
        return;
    }
    if (setup->basis_offsets != NULL) {
        PyMem_Free(setup->basis_offsets);
        setup->basis_offsets = NULL;
    }
    if (setup->active_counts != NULL) {
        PyMem_Free(setup->active_counts);
        setup->active_counts = NULL;
    }
    if (setup->active_offsets != NULL) {
        PyMem_Free(setup->active_offsets);
        setup->active_offsets = NULL;
    }
    if (setup->active_columns != NULL) {
        PyMem_Free(setup->active_columns);
        setup->active_columns = NULL;
    }
    if (setup->basis_real != NULL) {
        PyMem_Free(setup->basis_real);
        setup->basis_real = NULL;
    }
    if (setup->basis_imag != NULL) {
        PyMem_Free(setup->basis_imag);
        setup->basis_imag = NULL;
    }
    setup->ntrunc = -1;
    setup->rsphere = 0.0;
    setup->basis_used_entries = 0;
    setup->active_used_entries = 0;
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

static void uv_to_vordiv_block_setup_capsule_destructor(PyObject *capsule) {
    UvToVordivBlockSetup *setup = (UvToVordivBlockSetup *) PyCapsule_GetPointer(
        capsule,
        SKYBORN_UV_TO_VORDIV_BLOCK_CAPSULE
    );
    if (setup == NULL) {
        PyErr_Clear();
        return;
    }
    free_uv_to_vordiv_block_setup(setup);
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

static int get_uv_to_vordiv_block_setup_from_capsule(
    PyObject *handle_obj,
    UvToVordivBlockSetup **setup_out
) {
    UvToVordivBlockSetup *setup = NULL;

    if (!PyCapsule_IsValid(handle_obj, SKYBORN_UV_TO_VORDIV_BLOCK_CAPSULE)) {
        PyErr_SetString(
            PyExc_TypeError,
            "uv_to_vordiv block handle must be a valid setup capsule"
        );
        return 0;
    }

    setup = (UvToVordivBlockSetup *) PyCapsule_GetPointer(
        handle_obj,
        SKYBORN_UV_TO_VORDIV_BLOCK_CAPSULE
    );
    if (setup == NULL) {
        return 0;
    }
    if (
        setup->basis_offsets == NULL || setup->active_counts == NULL ||
        setup->active_offsets == NULL || setup->active_columns == NULL ||
        setup->basis_real == NULL || setup->basis_imag == NULL ||
        setup->ntrunc < 0
    ) {
        PyErr_SetString(PyExc_RuntimeError, "uv_to_vordiv block setup handle is closed");
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
    double rsphere = 0.0;
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

    if (!PyArg_ParseTuple(args, "OOOd", &nloen_obj, &mu_obj, &weights_obj, &rsphere)) {
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
    setup->rsphere = rsphere;
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
        weights_arr,
        "rsphere",
        setup->rsphere
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

static PyObject *backend_create_uv_to_vordiv_block_setup(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    UvToVordivBlockSetup *setup = NULL;
    PyObject *capsule = NULL;
    double rsphere = 0.0;
    int ignored_nt = 0, ignored_rank = 0, ntrunc = -1, v_ntrunc = -1, ierror = 0;
    npy_intp *ignored_dims = NULL;
    size_t ncoeff = 0, total_dense_entries = 0, total_active_columns = 0;
    size_t cache_bytes = 0, running_basis_offset = 0, running_active_offset = 0;
    double *basis_vrt_real = NULL, *basis_vrt_imag = NULL;
    double *basis_div_real = NULL, *basis_div_imag = NULL;
    double *u_basis_real = NULL, *u_basis_imag = NULL;
    double *v_basis_real = NULL, *v_basis_imag = NULL;
    double *basis_real = NULL, *basis_imag = NULL;
    size_t start_idx = 0;
    int m = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOd", &u_obj, &v_obj, &rsphere)) {
        return NULL;
    }

    u_arr = require_array(u_obj, NPY_COMPLEX128);
    v_arr = require_array(v_obj, NPY_COMPLEX128);
    if (u_arr == NULL || v_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(u_arr, INT_MAX, &ntrunc, &ignored_nt, &ignored_rank, &ignored_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(v_arr, INT_MAX, &v_ntrunc, &ignored_nt, &ignored_rank, &ignored_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(u_arr) != PyArray_NDIM(v_arr)) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(u_arr), PyArray_DIMS(v_arr), PyArray_NDIM(u_arr))) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same shape");
        goto fail;
    }
    if (ntrunc != v_ntrunc) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must use the same ntrunc");
        goto fail;
    }
    ncoeff = (size_t) PyArray_DIM(u_arr, 0);
    for (m = 0; m <= ntrunc; ++m) {
        size_t block_size = (size_t) (ntrunc - m + 1);
        size_t nrow = 2U * block_size;
        total_dense_entries += nrow * nrow;
        total_active_columns += nrow;
    }
    cache_bytes =
        total_dense_entries * sizeof(double) * 2U +
        total_active_columns * sizeof(int) +
        (size_t) (ntrunc + 1) * (sizeof(size_t) * 2U + sizeof(int));
    if (cache_bytes > SKYBORN_UV_TO_VORDIV_BLOCK_CACHE_MAX_BYTES) {
        PyErr_Format(
            PyExc_MemoryError,
            "uv_to_vordiv block cache would require about %.2f MiB, above the %.2f MiB limit",
            (double) cache_bytes / (1024.0 * 1024.0),
            (double) SKYBORN_UV_TO_VORDIV_BLOCK_CACHE_MAX_BYTES / (1024.0 * 1024.0)
        );
        goto fail;
    }

    setup = (UvToVordivBlockSetup *) PyMem_Malloc(sizeof(UvToVordivBlockSetup));
    if (setup == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    memset(setup, 0, sizeof(UvToVordivBlockSetup));
    setup->ntrunc = ntrunc;
    setup->rsphere = rsphere;
    setup->basis_offsets = (size_t *) PyMem_Malloc((size_t) (ntrunc + 1) * sizeof(size_t));
    setup->active_counts = (int *) PyMem_Malloc((size_t) (ntrunc + 1) * sizeof(int));
    setup->active_offsets = (size_t *) PyMem_Malloc((size_t) (ntrunc + 1) * sizeof(size_t));
    setup->active_columns = (int *) PyMem_Malloc(total_active_columns * sizeof(int));
    setup->basis_real = (double *) PyMem_Calloc(total_dense_entries, sizeof(double));
    setup->basis_imag = (double *) PyMem_Calloc(total_dense_entries, sizeof(double));
    if (
        setup->basis_offsets == NULL || setup->active_counts == NULL ||
        setup->active_offsets == NULL || setup->active_columns == NULL ||
        setup->basis_real == NULL || setup->basis_imag == NULL
    ) {
        PyErr_NoMemory();
        goto fail;
    }

    basis_vrt_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_vrt_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_div_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_div_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    u_basis_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    u_basis_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    v_basis_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    v_basis_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_real = (double *) PyMem_Calloc((size_t) (2 * (ntrunc + 1)) * (size_t) (2 * (ntrunc + 1)), sizeof(double));
    basis_imag = (double *) PyMem_Calloc((size_t) (2 * (ntrunc + 1)) * (size_t) (2 * (ntrunc + 1)), sizeof(double));
    if (
        basis_vrt_real == NULL || basis_vrt_imag == NULL ||
        basis_div_real == NULL || basis_div_imag == NULL ||
        u_basis_real == NULL || u_basis_imag == NULL ||
        v_basis_real == NULL || v_basis_imag == NULL ||
        basis_real == NULL || basis_imag == NULL
    ) {
        PyErr_NoMemory();
        goto fail;
    }

    start_idx = 0;
    for (m = 0; m <= ntrunc; ++m) {
        int block_size = ntrunc - m + 1;
        int nrow = 2 * block_size;
        int ncol = 2 * block_size;
        int local_idx = 0;
        int row = 0;
        int active_count = 0;

        memset(basis_real, 0, (size_t) nrow * (size_t) ncol * sizeof(double));
        memset(basis_imag, 0, (size_t) nrow * (size_t) ncol * sizeof(double));

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;
            size_t idx_local = 0;

            basis_vrt_real[global_idx] = 1.0;
            vordiv_to_uv(
                ntrunc, 1, rsphere,
                basis_vrt_real, basis_vrt_imag,
                basis_div_real, basis_div_imag,
                u_basis_real, u_basis_imag,
                v_basis_real, v_basis_imag,
                &ierror
            );
            basis_vrt_real[global_idx] = 0.0;
            if (ierror != 0) {
                PyErr_Format(PyExc_RuntimeError, "vordiv_to_uv failed while building vrt basis (ierror=%d)", ierror);
                goto fail;
            }
            for (idx_local = 0; idx_local < (size_t) block_size; ++idx_local) {
                size_t coeff_idx = start_idx + idx_local;
                size_t u_row = idx_local * (size_t) ncol + (size_t) local_idx;
                size_t v_row = (size_t) (block_size + (int) idx_local) * (size_t) ncol + (size_t) local_idx;
                basis_real[u_row] = u_basis_real[coeff_idx];
                basis_imag[u_row] = u_basis_imag[coeff_idx];
                basis_real[v_row] = v_basis_real[coeff_idx];
                basis_imag[v_row] = v_basis_imag[coeff_idx];
            }

            basis_div_real[global_idx] = 1.0;
            vordiv_to_uv(
                ntrunc, 1, rsphere,
                basis_vrt_real, basis_vrt_imag,
                basis_div_real, basis_div_imag,
                u_basis_real, u_basis_imag,
                v_basis_real, v_basis_imag,
                &ierror
            );
            basis_div_real[global_idx] = 0.0;
            if (ierror != 0) {
                PyErr_Format(PyExc_RuntimeError, "vordiv_to_uv failed while building div basis (ierror=%d)", ierror);
                goto fail;
            }
            for (idx_local = 0; idx_local < (size_t) block_size; ++idx_local) {
                size_t coeff_idx = start_idx + idx_local;
                size_t u_row = idx_local * (size_t) ncol + (size_t) (block_size + local_idx);
                size_t v_row = (size_t) (block_size + (int) idx_local) * (size_t) ncol + (size_t) (block_size + local_idx);
                basis_real[u_row] = u_basis_real[coeff_idx];
                basis_imag[u_row] = u_basis_imag[coeff_idx];
                basis_real[v_row] = v_basis_real[coeff_idx];
                basis_imag[v_row] = v_basis_imag[coeff_idx];
            }
        }

        for (local_idx = 0; local_idx < ncol; ++local_idx) {
            double column_norm = 0.0;
            for (row = 0; row < nrow; ++row) {
                size_t basis_idx = (size_t) row * (size_t) ncol + (size_t) local_idx;
                double column_value = fabs(basis_real[basis_idx]) + fabs(basis_imag[basis_idx]);
                if (column_value > column_norm) {
                    column_norm = column_value;
                }
            }
            if (column_norm > 1.0e-18) {
                setup->active_columns[running_active_offset + (size_t) active_count] = local_idx;
                active_count += 1;
            }
        }

        setup->basis_offsets[m] = running_basis_offset;
        setup->active_counts[m] = active_count;
        setup->active_offsets[m] = running_active_offset;
        memcpy(
            setup->basis_real + running_basis_offset,
            basis_real,
            (size_t) nrow * (size_t) ncol * sizeof(double)
        );
        memcpy(
            setup->basis_imag + running_basis_offset,
            basis_imag,
            (size_t) nrow * (size_t) ncol * sizeof(double)
        );

        running_basis_offset += (size_t) nrow * (size_t) ncol;
        running_active_offset += (size_t) active_count;
        start_idx += (size_t) block_size;
    }

    setup->basis_used_entries = running_basis_offset;
    setup->active_used_entries = running_active_offset;

    capsule = PyCapsule_New(
        (void *) setup,
        SKYBORN_UV_TO_VORDIV_BLOCK_CAPSULE,
        uv_to_vordiv_block_setup_capsule_destructor
    );
    if (capsule == NULL) {
        goto fail;
    }

    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    if (basis_vrt_real != NULL) PyMem_Free(basis_vrt_real);
    if (basis_vrt_imag != NULL) PyMem_Free(basis_vrt_imag);
    if (basis_div_real != NULL) PyMem_Free(basis_div_real);
    if (basis_div_imag != NULL) PyMem_Free(basis_div_imag);
    if (u_basis_real != NULL) PyMem_Free(u_basis_real);
    if (u_basis_imag != NULL) PyMem_Free(u_basis_imag);
    if (v_basis_real != NULL) PyMem_Free(v_basis_real);
    if (v_basis_imag != NULL) PyMem_Free(v_basis_imag);
    if (basis_real != NULL) PyMem_Free(basis_real);
    if (basis_imag != NULL) PyMem_Free(basis_imag);
    return capsule;

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    if (basis_vrt_real != NULL) PyMem_Free(basis_vrt_real);
    if (basis_vrt_imag != NULL) PyMem_Free(basis_vrt_imag);
    if (basis_div_real != NULL) PyMem_Free(basis_div_real);
    if (basis_div_imag != NULL) PyMem_Free(basis_div_imag);
    if (u_basis_real != NULL) PyMem_Free(u_basis_real);
    if (u_basis_imag != NULL) PyMem_Free(u_basis_imag);
    if (v_basis_real != NULL) PyMem_Free(v_basis_real);
    if (v_basis_imag != NULL) PyMem_Free(v_basis_imag);
    if (basis_real != NULL) PyMem_Free(basis_real);
    if (basis_imag != NULL) PyMem_Free(basis_imag);
    if (setup != NULL) {
        free_uv_to_vordiv_block_setup(setup);
        PyMem_Free(setup);
    }
    Py_XDECREF(capsule);
    return NULL;
}

static PyObject *backend_destroy_uv_to_vordiv_block_setup(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    UvToVordivBlockSetup *setup = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "O", &handle_obj)) {
        return NULL;
    }
    if (!get_uv_to_vordiv_block_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }

    free_uv_to_vordiv_block_setup(setup);
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
        (const double *) setup->mu,
        (const double *) setup->weights,
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
        (const double *) setup->mu,
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

static PyObject *backend_scalar_fourier(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *grid_obj = NULL;
    int mmax = -1;
    EctransSetup *setup = NULL;
    PyArrayObject *grid_arr = NULL;
    PyArrayObject *real_arr = NULL;
    PyArrayObject *imag_arr = NULL;
    PyArrayObject *out_arr = NULL;
    int nt = 0, grid_rank = 0, ierror = 0;
    npy_intp *grid_dims = NULL;
    npy_intp *out_dims = NULL;
    npy_intp idx;
    npy_complex128 *out_data = NULL;
    double *real_data = NULL;
    double *imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOi", &handle_obj, &grid_obj, &mmax)) {
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
    if (mmax < 0 || mmax > setup->max_nlon - 1) {
        PyErr_Format(
            PyExc_ValueError,
            "mmax must be between 0 and %d",
            setup->max_nlon - 1
        );
        goto fail;
    }

    out_dims = PyMem_Malloc((size_t) (grid_rank + 1) * sizeof(npy_intp));
    if (out_dims == NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    out_dims[0] = setup->ndgl;
    out_dims[1] = (npy_intp) mmax + 1;
    for (idx = 1; idx < grid_rank; ++idx) {
        out_dims[idx + 1] = grid_dims[idx];
    }

    real_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank + 1, out_dims, NPY_FLOAT64, 0);
    imag_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank + 1, out_dims, NPY_FLOAT64, 0);
    out_arr = (PyArrayObject *) PyArray_ZEROS(grid_rank + 1, out_dims, NPY_COMPLEX128, 0);
    if (real_arr == NULL || imag_arr == NULL || out_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    scalar_fourier_stub(
        setup->ndgl,
        (const int *) setup->nloen,
        setup->ngptot,
        mmax,
        nt,
        (const double *) PyArray_DATA(grid_arr),
        (double *) PyArray_DATA(real_arr),
        (double *) PyArray_DATA(imag_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    out_data = (npy_complex128 *) PyArray_DATA(out_arr);
    real_data = (double *) PyArray_DATA(real_arr);
    imag_data = (double *) PyArray_DATA(imag_arr);
    for (idx = 0; idx < PyArray_SIZE(out_arr); ++idx) {
        out_data[idx] = real_data[idx] + imag_data[idx] * I;
    }

    PyMem_Free(out_dims);
    Py_DECREF(grid_arr);
    Py_DECREF(real_arr);
    Py_DECREF(imag_arr);
    return Py_BuildValue("(Ni)", out_arr, ierror);

fail:
    if (out_dims != NULL) {
        PyMem_Free(out_dims);
    }
    Py_XDECREF(grid_arr);
    Py_XDECREF(real_arr);
    Py_XDECREF(imag_arr);
    Py_XDECREF(out_arr);
    return NULL;
}

static PyObject *backend_scalar_block_solve(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *basis_obj = NULL;
    PyObject *observed_obj = NULL;
    EctransSetup *setup = NULL;
    PyArrayObject *basis_arr = NULL;
    PyArrayObject *observed_arr = NULL;
    PyArrayObject *basis_real_arr = NULL;
    PyArrayObject *basis_imag_arr = NULL;
    PyArrayObject *observed_real_arr = NULL;
    PyArrayObject *observed_imag_arr = NULL;
    PyArrayObject *solution_real_arr = NULL;
    PyArrayObject *solution_imag_arr = NULL;
    PyArrayObject *solution_out_arr = NULL;
    int nblock = 0, nt = 0, ierror = 0;
    npy_intp out_dims[2];
    npy_intp idx;
    npy_complex128 *basis_data = NULL;
    npy_complex128 *observed_data = NULL;
    npy_complex128 *solution_out_data = NULL;
    double *basis_real_data = NULL;
    double *basis_imag_data = NULL;
    double *observed_real_data = NULL;
    double *observed_imag_data = NULL;
    double *solution_real_data = NULL;
    double *solution_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOO", &handle_obj, &basis_obj, &observed_obj)) {
        return NULL;
    }
    if (!get_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }

    basis_arr = require_array(basis_obj, NPY_COMPLEX128);
    observed_arr = require_array(observed_obj, NPY_COMPLEX128);
    if (basis_arr == NULL || observed_arr == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(basis_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "basis must be a rank-2 complex128 array");
        goto fail;
    }
    if (PyArray_NDIM(observed_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "observed must be a rank-2 complex128 array");
        goto fail;
    }
    if (PyArray_DIM(basis_arr, 0) != (npy_intp) setup->ndgl) {
        PyErr_Format(PyExc_ValueError, "basis leading dimension must be %d", setup->ndgl);
        goto fail;
    }
    if (PyArray_DIM(observed_arr, 0) != (npy_intp) setup->ndgl) {
        PyErr_Format(PyExc_ValueError, "observed leading dimension must be %d", setup->ndgl);
        goto fail;
    }

    nblock = (int) PyArray_DIM(basis_arr, 1);
    nt = (int) PyArray_DIM(observed_arr, 1);
    out_dims[0] = (npy_intp) nblock;
    out_dims[1] = (npy_intp) nt;

    basis_real_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(basis_arr), NPY_FLOAT64, 0);
    basis_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(basis_arr), NPY_FLOAT64, 0);
    observed_real_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(observed_arr), NPY_FLOAT64, 0);
    observed_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(observed_arr), NPY_FLOAT64, 0);
    solution_real_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_FLOAT64, 0);
    solution_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_FLOAT64, 0);
    solution_out_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_COMPLEX128, 0);
    if (
        basis_real_arr == NULL || basis_imag_arr == NULL ||
        observed_real_arr == NULL || observed_imag_arr == NULL ||
        solution_real_arr == NULL || solution_imag_arr == NULL ||
        solution_out_arr == NULL
    ) {
        goto fail;
    }

    basis_data = (npy_complex128 *) PyArray_DATA(basis_arr);
    observed_data = (npy_complex128 *) PyArray_DATA(observed_arr);
    basis_real_data = (double *) PyArray_DATA(basis_real_arr);
    basis_imag_data = (double *) PyArray_DATA(basis_imag_arr);
    observed_real_data = (double *) PyArray_DATA(observed_real_arr);
    observed_imag_data = (double *) PyArray_DATA(observed_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(basis_arr); ++idx) {
        basis_real_data[idx] = creal(basis_data[idx]);
        basis_imag_data[idx] = cimag(basis_data[idx]);
    }
    for (idx = 0; idx < PyArray_SIZE(observed_arr); ++idx) {
        observed_real_data[idx] = creal(observed_data[idx]);
        observed_imag_data[idx] = cimag(observed_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    scalar_block_solve_stub(
        setup->ndgl,
        nblock,
        nt,
        (const double *) setup->weights,
        (const double *) PyArray_DATA(basis_real_arr),
        (const double *) PyArray_DATA(basis_imag_arr),
        (const double *) PyArray_DATA(observed_real_arr),
        (const double *) PyArray_DATA(observed_imag_arr),
        (double *) PyArray_DATA(solution_real_arr),
        (double *) PyArray_DATA(solution_imag_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    solution_out_data = (npy_complex128 *) PyArray_DATA(solution_out_arr);
    solution_real_data = (double *) PyArray_DATA(solution_real_arr);
    solution_imag_data = (double *) PyArray_DATA(solution_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(solution_out_arr); ++idx) {
        solution_out_data[idx] = solution_real_data[idx] + solution_imag_data[idx] * I;
    }

    Py_DECREF(basis_arr);
    Py_DECREF(observed_arr);
    Py_DECREF(basis_real_arr);
    Py_DECREF(basis_imag_arr);
    Py_DECREF(observed_real_arr);
    Py_DECREF(observed_imag_arr);
    Py_DECREF(solution_real_arr);
    Py_DECREF(solution_imag_arr);
    return Py_BuildValue("(Ni)", solution_out_arr, ierror);

fail:
    Py_XDECREF(basis_arr);
    Py_XDECREF(observed_arr);
    Py_XDECREF(basis_real_arr);
    Py_XDECREF(basis_imag_arr);
    Py_XDECREF(observed_real_arr);
    Py_XDECREF(observed_imag_arr);
    Py_XDECREF(solution_real_arr);
    Py_XDECREF(solution_imag_arr);
    Py_XDECREF(solution_out_arr);
    return NULL;
}

static PyObject *backend_weighted_block_solve(PyObject *self, PyObject *args) {
    PyObject *weights_obj = NULL;
    PyObject *basis_obj = NULL;
    PyObject *observed_obj = NULL;
    PyArrayObject *weights_arr = NULL;
    PyArrayObject *basis_arr = NULL;
    PyArrayObject *observed_arr = NULL;
    PyArrayObject *basis_real_arr = NULL;
    PyArrayObject *basis_imag_arr = NULL;
    PyArrayObject *observed_real_arr = NULL;
    PyArrayObject *observed_imag_arr = NULL;
    PyArrayObject *solution_real_arr = NULL;
    PyArrayObject *solution_imag_arr = NULL;
    PyArrayObject *solution_out_arr = NULL;
    int nrow = 0, nblock = 0, nt = 0, ierror = 0;
    npy_intp out_dims[2];
    npy_intp idx;
    npy_complex128 *basis_data = NULL;
    npy_complex128 *observed_data = NULL;
    npy_complex128 *solution_out_data = NULL;
    double *basis_real_data = NULL;
    double *basis_imag_data = NULL;
    double *observed_real_data = NULL;
    double *observed_imag_data = NULL;
    double *solution_real_data = NULL;
    double *solution_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOO", &weights_obj, &basis_obj, &observed_obj)) {
        return NULL;
    }

    weights_arr = require_array(weights_obj, NPY_FLOAT64);
    basis_arr = require_array(basis_obj, NPY_COMPLEX128);
    observed_arr = require_array(observed_obj, NPY_COMPLEX128);
    if (weights_arr == NULL || basis_arr == NULL || observed_arr == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(weights_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "weights must be a 1D float64 array");
        goto fail;
    }
    if (PyArray_NDIM(basis_arr) != 2 || PyArray_NDIM(observed_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "basis and observed must both be 2D complex128 arrays");
        goto fail;
    }

    nrow = (int) PyArray_DIM(weights_arr, 0);
    if (PyArray_DIM(basis_arr, 0) != (npy_intp) nrow || PyArray_DIM(observed_arr, 0) != (npy_intp) nrow) {
        PyErr_SetString(PyExc_ValueError, "weights, basis, and observed must use the same leading dimension");
        goto fail;
    }

    nblock = (int) PyArray_DIM(basis_arr, 1);
    nt = (int) PyArray_DIM(observed_arr, 1);
    out_dims[0] = (npy_intp) nblock;
    out_dims[1] = (npy_intp) nt;

    basis_real_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(basis_arr), NPY_FLOAT64, 0);
    basis_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(basis_arr), NPY_FLOAT64, 0);
    observed_real_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(observed_arr), NPY_FLOAT64, 0);
    observed_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, PyArray_DIMS(observed_arr), NPY_FLOAT64, 0);
    solution_real_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_FLOAT64, 0);
    solution_imag_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_FLOAT64, 0);
    solution_out_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_COMPLEX128, 0);
    if (
        basis_real_arr == NULL || basis_imag_arr == NULL ||
        observed_real_arr == NULL || observed_imag_arr == NULL ||
        solution_real_arr == NULL || solution_imag_arr == NULL ||
        solution_out_arr == NULL
    ) {
        goto fail;
    }

    basis_data = (npy_complex128 *) PyArray_DATA(basis_arr);
    observed_data = (npy_complex128 *) PyArray_DATA(observed_arr);
    basis_real_data = (double *) PyArray_DATA(basis_real_arr);
    basis_imag_data = (double *) PyArray_DATA(basis_imag_arr);
    observed_real_data = (double *) PyArray_DATA(observed_real_arr);
    observed_imag_data = (double *) PyArray_DATA(observed_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(basis_arr); ++idx) {
        basis_real_data[idx] = creal(basis_data[idx]);
        basis_imag_data[idx] = cimag(basis_data[idx]);
    }
    for (idx = 0; idx < PyArray_SIZE(observed_arr); ++idx) {
        observed_real_data[idx] = creal(observed_data[idx]);
        observed_imag_data[idx] = cimag(observed_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    weighted_block_solve_stub(
        nrow,
        nblock,
        nt,
        (const double *) PyArray_DATA(weights_arr),
        (const double *) PyArray_DATA(basis_real_arr),
        (const double *) PyArray_DATA(basis_imag_arr),
        (const double *) PyArray_DATA(observed_real_arr),
        (const double *) PyArray_DATA(observed_imag_arr),
        (double *) PyArray_DATA(solution_real_arr),
        (double *) PyArray_DATA(solution_imag_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    solution_out_data = (npy_complex128 *) PyArray_DATA(solution_out_arr);
    solution_real_data = (double *) PyArray_DATA(solution_real_arr);
    solution_imag_data = (double *) PyArray_DATA(solution_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(solution_out_arr); ++idx) {
        solution_out_data[idx] = solution_real_data[idx] + solution_imag_data[idx] * I;
    }

    Py_DECREF(weights_arr);
    Py_DECREF(basis_arr);
    Py_DECREF(observed_arr);
    Py_DECREF(basis_real_arr);
    Py_DECREF(basis_imag_arr);
    Py_DECREF(observed_real_arr);
    Py_DECREF(observed_imag_arr);
    Py_DECREF(solution_real_arr);
    Py_DECREF(solution_imag_arr);
    return Py_BuildValue("(Ni)", solution_out_arr, ierror);

fail:
    Py_XDECREF(weights_arr);
    Py_XDECREF(basis_arr);
    Py_XDECREF(observed_arr);
    Py_XDECREF(basis_real_arr);
    Py_XDECREF(basis_imag_arr);
    Py_XDECREF(observed_real_arr);
    Py_XDECREF(observed_imag_arr);
    Py_XDECREF(solution_real_arr);
    Py_XDECREF(solution_imag_arr);
    Py_XDECREF(solution_out_arr);
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
        (const double *) setup->weights,
        setup->ngptot,
        setup->rsphere,
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
        setup->rsphere,
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
        setup->rsphere,
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

static PyObject *backend_vordiv_to_uv(PyObject *self, PyObject *args) {
    PyObject *vrt_obj = NULL;
    PyObject *div_obj = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *div_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *u_real_arr = NULL;
    PyArrayObject *u_imag_arr = NULL;
    PyArrayObject *v_real_arr = NULL;
    PyArrayObject *v_imag_arr = NULL;
    PyArrayObject *u_out_arr = NULL;
    PyArrayObject *v_out_arr = NULL;
    double rsphere = 0.0;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, div_ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp idx;
    npy_complex128 *vrt_data = NULL;
    npy_complex128 *div_data = NULL;
    npy_complex128 *u_out_data = NULL;
    npy_complex128 *v_out_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;
    double *u_real_data = NULL;
    double *u_imag_data = NULL;
    double *v_real_data = NULL;
    double *v_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOd", &vrt_obj, &div_obj, &rsphere)) {
        return NULL;
    }

    vrt_arr = require_array(vrt_obj, NPY_COMPLEX128);
    div_arr = require_array(div_obj, NPY_COMPLEX128);
    if (vrt_arr == NULL || div_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(vrt_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(div_arr, INT_MAX, &div_ntrunc, &nt, &spec_rank, &spec_dims)) {
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

    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    u_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    u_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    u_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    v_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    if (
        vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        u_real_arr == NULL || u_imag_arr == NULL ||
        v_real_arr == NULL || v_imag_arr == NULL ||
        u_out_arr == NULL || v_out_arr == NULL
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
    vordiv_to_uv(
        ntrunc,
        nt,
        rsphere,
        (const double *) PyArray_DATA(vrt_real_arr),
        (const double *) PyArray_DATA(vrt_imag_arr),
        (const double *) PyArray_DATA(div_real_arr),
        (const double *) PyArray_DATA(div_imag_arr),
        (double *) PyArray_DATA(u_real_arr),
        (double *) PyArray_DATA(u_imag_arr),
        (double *) PyArray_DATA(v_real_arr),
        (double *) PyArray_DATA(v_imag_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    u_out_data = (npy_complex128 *) PyArray_DATA(u_out_arr);
    v_out_data = (npy_complex128 *) PyArray_DATA(v_out_arr);
    u_real_data = (double *) PyArray_DATA(u_real_arr);
    u_imag_data = (double *) PyArray_DATA(u_imag_arr);
    v_real_data = (double *) PyArray_DATA(v_real_arr);
    v_imag_data = (double *) PyArray_DATA(v_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(u_out_arr); ++idx) {
        u_out_data[idx] = u_real_data[idx] + u_imag_data[idx] * I;
        v_out_data[idx] = v_real_data[idx] + v_imag_data[idx] * I;
    }

    Py_DECREF(vrt_arr);
    Py_DECREF(div_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    Py_DECREF(u_real_arr);
    Py_DECREF(u_imag_arr);
    Py_DECREF(v_real_arr);
    Py_DECREF(v_imag_arr);
    return Py_BuildValue("(NNi)", u_out_arr, v_out_arr, ierror);

fail:
    Py_XDECREF(vrt_arr);
    Py_XDECREF(div_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(u_real_arr);
    Py_XDECREF(u_imag_arr);
    Py_XDECREF(v_real_arr);
    Py_XDECREF(v_imag_arr);
    Py_XDECREF(u_out_arr);
    Py_XDECREF(v_out_arr);
    return NULL;
}

static PyObject *backend_prfi1b_uv_block(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *u_real_arr = NULL;
    PyArrayObject *u_imag_arr = NULL;
    PyArrayObject *v_real_arr = NULL;
    PyArrayObject *v_imag_arr = NULL;
    PyArrayObject *poa1_arr = NULL;
    double rsphere = 0.0;
    int km = -1, nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, v_ntrunc = -1;
    int nlei1 = -1;
    npy_intp *spec_dims = NULL;
    npy_intp poa1_dims[2];
    npy_intp idx;
    npy_complex128 *u_data = NULL;
    npy_complex128 *v_data = NULL;
    double *u_real_data = NULL;
    double *u_imag_data = NULL;
    double *v_real_data = NULL;
    double *v_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "iOOd", &km, &u_obj, &v_obj, &rsphere)) {
        return NULL;
    }

    u_arr = require_array(u_obj, NPY_COMPLEX128);
    v_arr = require_array(v_obj, NPY_COMPLEX128);
    if (u_arr == NULL || v_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(u_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(v_arr, INT_MAX, &v_ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(u_arr) != PyArray_NDIM(v_arr)) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(u_arr), PyArray_DIMS(v_arr), PyArray_NDIM(u_arr))) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same shape");
        goto fail;
    }
    if (ntrunc != v_ntrunc) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must use the same ntrunc");
        goto fail;
    }
    if (km < 0 || km > ntrunc) {
        PyErr_SetString(PyExc_ValueError, "require 0 <= km <= ntrunc");
        goto fail;
    }

    u_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    u_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2);
    poa1_dims[0] = (npy_intp) nlei1;
    poa1_dims[1] = (npy_intp) (4 * nt);
    poa1_arr = (PyArrayObject *) PyArray_ZEROS(2, poa1_dims, NPY_FLOAT64, 0);
    if (
        u_real_arr == NULL || u_imag_arr == NULL ||
        v_real_arr == NULL || v_imag_arr == NULL ||
        poa1_arr == NULL
    ) {
        goto fail;
    }

    u_data = (npy_complex128 *) PyArray_DATA(u_arr);
    v_data = (npy_complex128 *) PyArray_DATA(v_arr);
    u_real_data = (double *) PyArray_DATA(u_real_arr);
    u_imag_data = (double *) PyArray_DATA(u_imag_arr);
    v_real_data = (double *) PyArray_DATA(v_real_arr);
    v_imag_data = (double *) PyArray_DATA(v_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(u_arr); ++idx) {
        u_real_data[idx] = creal(u_data[idx]);
        u_imag_data[idx] = cimag(u_data[idx]);
        v_real_data[idx] = creal(v_data[idx]);
        v_imag_data[idx] = cimag(v_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    prfi1b_uv_block(
        ntrunc,
        km,
        nt,
        rsphere,
        (const double *) PyArray_DATA(u_real_arr),
        (const double *) PyArray_DATA(u_imag_arr),
        (const double *) PyArray_DATA(v_real_arr),
        (const double *) PyArray_DATA(v_imag_arr),
        (double *) PyArray_DATA(poa1_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    Py_DECREF(u_real_arr);
    Py_DECREF(u_imag_arr);
    Py_DECREF(v_real_arr);
    Py_DECREF(v_imag_arr);
    return Py_BuildValue("(Ni)", poa1_arr, ierror);

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(u_real_arr);
    Py_XDECREF(u_imag_arr);
    Py_XDECREF(v_real_arr);
    Py_XDECREF(v_imag_arr);
    Py_XDECREF(poa1_arr);
    return NULL;
}

static PyObject *backend_vd2uv_uv_block(PyObject *self, PyObject *args) {
    PyObject *vrt_obj = NULL;
    PyObject *div_obj = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *div_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *poa1_arr = NULL;
    double rsphere = 0.0;
    int km = -1, nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, div_ntrunc = -1;
    int nlei1 = -1;
    npy_intp *spec_dims = NULL;
    npy_intp poa1_dims[2];
    npy_intp idx;
    npy_complex128 *vrt_data = NULL;
    npy_complex128 *div_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "iOOd", &km, &vrt_obj, &div_obj, &rsphere)) {
        return NULL;
    }

    vrt_arr = require_array(vrt_obj, NPY_COMPLEX128);
    div_arr = require_array(div_obj, NPY_COMPLEX128);
    if (vrt_arr == NULL || div_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(vrt_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(div_arr, INT_MAX, &div_ntrunc, &nt, &spec_rank, &spec_dims)) {
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
    if (km < 0 || km > ntrunc) {
        PyErr_SetString(PyExc_ValueError, "require 0 <= km <= ntrunc");
        goto fail;
    }

    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2);
    poa1_dims[0] = (npy_intp) nlei1;
    poa1_dims[1] = (npy_intp) (4 * nt);
    poa1_arr = (PyArrayObject *) PyArray_ZEROS(2, poa1_dims, NPY_FLOAT64, 0);
    if (
        vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        poa1_arr == NULL
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
    vd2uv_uv_block(
        ntrunc,
        km,
        nt,
        rsphere,
        (const double *) PyArray_DATA(vrt_real_arr),
        (const double *) PyArray_DATA(vrt_imag_arr),
        (const double *) PyArray_DATA(div_real_arr),
        (const double *) PyArray_DATA(div_imag_arr),
        (double *) PyArray_DATA(poa1_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(vrt_arr);
    Py_DECREF(div_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    return Py_BuildValue("(Ni)", poa1_arr, ierror);

fail:
    Py_XDECREF(vrt_arr);
    Py_XDECREF(div_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(poa1_arr);
    return NULL;
}

static PyObject *backend_uv_to_vordiv(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *u_real_arr = NULL;
    PyArrayObject *u_imag_arr = NULL;
    PyArrayObject *v_real_arr = NULL;
    PyArrayObject *v_imag_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *vrt_out_arr = NULL;
    PyArrayObject *div_out_arr = NULL;
    double rsphere = 0.0;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, v_ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp idx;
    npy_complex128 *u_data = NULL;
    npy_complex128 *v_data = NULL;
    npy_complex128 *vrt_out_data = NULL;
    npy_complex128 *div_out_data = NULL;
    double *u_real_data = NULL;
    double *u_imag_data = NULL;
    double *v_real_data = NULL;
    double *v_imag_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOd", &u_obj, &v_obj, &rsphere)) {
        return NULL;
    }

    u_arr = require_array(u_obj, NPY_COMPLEX128);
    v_arr = require_array(v_obj, NPY_COMPLEX128);
    if (u_arr == NULL || v_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(u_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(v_arr, INT_MAX, &v_ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(u_arr) != PyArray_NDIM(v_arr)) {
      PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same rank");
      goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(u_arr), PyArray_DIMS(v_arr), PyArray_NDIM(u_arr))) {
      PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same shape");
      goto fail;
    }
    if (ntrunc != v_ntrunc) {
      PyErr_SetString(PyExc_ValueError, "uspec and vspec must use the same ntrunc");
      goto fail;
    }

    u_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    u_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    v_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_FLOAT64, 0);
    vrt_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    div_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    if (
        u_real_arr == NULL || u_imag_arr == NULL ||
        v_real_arr == NULL || v_imag_arr == NULL ||
        vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        vrt_out_arr == NULL || div_out_arr == NULL
    ) {
        goto fail;
    }

    u_data = (npy_complex128 *) PyArray_DATA(u_arr);
    v_data = (npy_complex128 *) PyArray_DATA(v_arr);
    u_real_data = (double *) PyArray_DATA(u_real_arr);
    u_imag_data = (double *) PyArray_DATA(u_imag_arr);
    v_real_data = (double *) PyArray_DATA(v_real_arr);
    v_imag_data = (double *) PyArray_DATA(v_imag_arr);
    for (idx = 0; idx < PyArray_SIZE(u_arr); ++idx) {
        u_real_data[idx] = creal(u_data[idx]);
        u_imag_data[idx] = cimag(u_data[idx]);
        v_real_data[idx] = creal(v_data[idx]);
        v_imag_data[idx] = cimag(v_data[idx]);
    }

    Py_BEGIN_ALLOW_THREADS
    uv_to_vordiv(
        ntrunc,
        nt,
        rsphere,
        (const double *) PyArray_DATA(u_real_arr),
        (const double *) PyArray_DATA(u_imag_arr),
        (const double *) PyArray_DATA(v_real_arr),
        (const double *) PyArray_DATA(v_imag_arr),
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

    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    Py_DECREF(u_real_arr);
    Py_DECREF(u_imag_arr);
    Py_DECREF(v_real_arr);
    Py_DECREF(v_imag_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    return Py_BuildValue("(NNi)", vrt_out_arr, div_out_arr, ierror);

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(u_real_arr);
    Py_XDECREF(u_imag_arr);
    Py_XDECREF(v_real_arr);
    Py_XDECREF(v_imag_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(vrt_out_arr);
    Py_XDECREF(div_out_arr);
    return NULL;
}

static PyObject *backend_uv_to_vordiv_block_native(PyObject *self, PyObject *args) {
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *vrt_out_arr = NULL;
    PyArrayObject *div_out_arr = NULL;
    double rsphere = 0.0;
    int nt = 0, spec_rank = 0, ierror = 0, ntrunc = -1, v_ntrunc = -1;
    npy_intp *spec_dims = NULL;
    npy_intp idx;
    npy_complex128 *u_data = NULL;
    npy_complex128 *v_data = NULL;
    npy_complex128 *vrt_out_data = NULL;
    npy_complex128 *div_out_data = NULL;
    double *basis_vrt_real = NULL, *basis_vrt_imag = NULL;
    double *basis_div_real = NULL, *basis_div_imag = NULL;
    double *u_basis_real = NULL, *u_basis_imag = NULL;
    double *v_basis_real = NULL, *v_basis_imag = NULL;
    double *weights = NULL;
    double *basis_real = NULL, *basis_imag = NULL;
    double *reduced_basis_real = NULL, *reduced_basis_imag = NULL;
    double *observed_real = NULL, *observed_imag = NULL;
    double *solution_real = NULL, *solution_imag = NULL;
    double *reduced_solution_real = NULL, *reduced_solution_imag = NULL;
    int *active_columns = NULL;
    size_t ncoeff = 0, max_block = 0, max_nrow = 0;
    size_t start_idx = 0;
    int m = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOd", &u_obj, &v_obj, &rsphere)) {
        return NULL;
    }

    u_arr = require_array(u_obj, NPY_COMPLEX128);
    v_arr = require_array(v_obj, NPY_COMPLEX128);
    if (u_arr == NULL || v_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(u_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(v_arr, INT_MAX, &v_ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(u_arr) != PyArray_NDIM(v_arr)) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(u_arr), PyArray_DIMS(v_arr), PyArray_NDIM(u_arr))) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same shape");
        goto fail;
    }
    if (ntrunc != v_ntrunc) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must use the same ntrunc");
        goto fail;
    }

    vrt_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    div_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    if (vrt_out_arr == NULL || div_out_arr == NULL) {
        goto fail;
    }
    vrt_out_data = (npy_complex128 *) PyArray_DATA(vrt_out_arr);
    div_out_data = (npy_complex128 *) PyArray_DATA(div_out_arr);

    ncoeff = (size_t) PyArray_DIM(u_arr, 0);
    max_block = (size_t) (ntrunc + 1);
    max_nrow = 2U * max_block;

    weights = (double *) PyMem_Malloc(max_nrow * sizeof(double));
    reduced_basis_real = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    reduced_basis_imag = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    observed_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    observed_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    solution_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    solution_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    reduced_solution_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    reduced_solution_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    active_columns = (int *) PyMem_Malloc(max_nrow * sizeof(int));
    if (
        weights == NULL ||
        reduced_basis_real == NULL || reduced_basis_imag == NULL ||
        observed_real == NULL || observed_imag == NULL ||
        solution_real == NULL || solution_imag == NULL ||
        reduced_solution_real == NULL || reduced_solution_imag == NULL ||
        active_columns == NULL
    ) {
        PyErr_NoMemory();
        goto fail;
    }

    for (idx = 0; idx < (npy_intp) max_nrow; ++idx) {
        weights[idx] = 1.0;
    }

    u_data = (npy_complex128 *) PyArray_DATA(u_arr);
    v_data = (npy_complex128 *) PyArray_DATA(v_arr);
    vrt_out_data = (npy_complex128 *) PyArray_DATA(vrt_out_arr);
    div_out_data = (npy_complex128 *) PyArray_DATA(div_out_arr);

    basis_vrt_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_vrt_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_div_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_div_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    u_basis_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    u_basis_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    v_basis_real = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    v_basis_imag = (double *) PyMem_Calloc(ncoeff, sizeof(double));
    basis_real = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    basis_imag = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    if (
        basis_vrt_real == NULL || basis_vrt_imag == NULL ||
        basis_div_real == NULL || basis_div_imag == NULL ||
        u_basis_real == NULL || u_basis_imag == NULL ||
        v_basis_real == NULL || v_basis_imag == NULL ||
        basis_real == NULL || basis_imag == NULL
    ) {
        PyErr_NoMemory();
        goto fail;
    }

    start_idx = 0;
    for (m = 0; m <= ntrunc; ++m) {
        int block_size = ntrunc - m + 1;
        int nrow = 2 * block_size;
        int ncol = 2 * block_size;
        int local_idx = 0;
        int it = 0;
        int active_count = 0;
        int row = 0;
        size_t basis_bytes = (size_t) nrow * (size_t) ncol * sizeof(double);
        size_t rhs_bytes = (size_t) nrow * (size_t) nt * sizeof(double);

        memset(basis_real, 0, basis_bytes);
        memset(basis_imag, 0, basis_bytes);
        memset(reduced_basis_real, 0, basis_bytes);
        memset(reduced_basis_imag, 0, basis_bytes);
        memset(observed_real, 0, rhs_bytes);
        memset(observed_imag, 0, rhs_bytes);
        memset(solution_real, 0, rhs_bytes);
        memset(solution_imag, 0, rhs_bytes);
        memset(reduced_solution_real, 0, rhs_bytes);
        memset(reduced_solution_imag, 0, rhs_bytes);

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;
            for (it = 0; it < nt; ++it) {
                size_t input_idx = global_idx * (size_t) nt + (size_t) it;
                size_t u_obs_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                size_t v_obs_idx = (size_t) (block_size + local_idx) * (size_t) nt + (size_t) it;
                observed_real[u_obs_idx] = creal(u_data[input_idx]);
                observed_imag[u_obs_idx] = cimag(u_data[input_idx]);
                observed_real[v_obs_idx] = creal(v_data[input_idx]);
                observed_imag[v_obs_idx] = cimag(v_data[input_idx]);
            }
        }

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;

            basis_vrt_real[global_idx] = 1.0;
            vordiv_to_uv(
                ntrunc, 1, rsphere,
                basis_vrt_real, basis_vrt_imag,
                basis_div_real, basis_div_imag,
                u_basis_real, u_basis_imag,
                v_basis_real, v_basis_imag,
                &ierror
            );
            basis_vrt_real[global_idx] = 0.0;
            if (ierror != 0) {
                goto done;
            }
            for (idx = 0; idx < block_size; ++idx) {
                size_t coeff_idx = start_idx + idx;
                size_t u_row = idx * (size_t) ncol + (size_t) local_idx;
                size_t v_row = (size_t) (block_size + (int) idx) * (size_t) ncol + (size_t) local_idx;
                basis_real[u_row] = u_basis_real[coeff_idx];
                basis_imag[u_row] = u_basis_imag[coeff_idx];
                basis_real[v_row] = v_basis_real[coeff_idx];
                basis_imag[v_row] = v_basis_imag[coeff_idx];
            }

            basis_div_real[global_idx] = 1.0;
            vordiv_to_uv(
                ntrunc, 1, rsphere,
                basis_vrt_real, basis_vrt_imag,
                basis_div_real, basis_div_imag,
                u_basis_real, u_basis_imag,
                v_basis_real, v_basis_imag,
                &ierror
            );
            basis_div_real[global_idx] = 0.0;
            if (ierror != 0) {
                goto done;
            }
            for (idx = 0; idx < block_size; ++idx) {
                size_t coeff_idx = start_idx + idx;
                size_t u_row = idx * (size_t) ncol + (size_t) (block_size + local_idx);
                size_t v_row = (size_t) (block_size + (int) idx) * (size_t) ncol + (size_t) (block_size + local_idx);
                basis_real[u_row] = u_basis_real[coeff_idx];
                basis_imag[u_row] = u_basis_imag[coeff_idx];
                basis_real[v_row] = v_basis_real[coeff_idx];
                basis_imag[v_row] = v_basis_imag[coeff_idx];
            }
        }

        active_count = 0;
        for (local_idx = 0; local_idx < ncol; ++local_idx) {
            double column_norm = 0.0;
            for (row = 0; row < nrow; ++row) {
                size_t basis_idx = (size_t) row * (size_t) ncol + (size_t) local_idx;
                double column_value = fabs(basis_real[basis_idx]) + fabs(basis_imag[basis_idx]);
                if (column_value > column_norm) {
                    column_norm = column_value;
                }
            }
            if (column_norm > 1.0e-18) {
                active_columns[active_count++] = local_idx;
            }
        }

        if (active_count > 0) {
            for (row = 0; row < nrow; ++row) {
                int active_idx = 0;
                for (active_idx = 0; active_idx < active_count; ++active_idx) {
                    int col = active_columns[active_idx];
                    size_t source_idx = (size_t) row * (size_t) ncol + (size_t) col;
                    size_t target_idx = (size_t) row * (size_t) active_count + (size_t) active_idx;
                    reduced_basis_real[target_idx] = basis_real[source_idx];
                    reduced_basis_imag[target_idx] = basis_imag[source_idx];
                }
            }

            weighted_block_solve_stub(
                nrow, active_count, nt, weights,
                reduced_basis_real, reduced_basis_imag,
                observed_real, observed_imag,
                reduced_solution_real, reduced_solution_imag,
                &ierror
            );
            if (ierror != 0) {
                goto done;
            }

            for (local_idx = 0; local_idx < active_count; ++local_idx) {
                int col = active_columns[local_idx];
                for (it = 0; it < nt; ++it) {
                    size_t source_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                    size_t target_idx = (size_t) col * (size_t) nt + (size_t) it;
                    solution_real[target_idx] = reduced_solution_real[source_idx];
                    solution_imag[target_idx] = reduced_solution_imag[source_idx];
                }
            }
        }

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;
            for (it = 0; it < nt; ++it) {
                size_t output_idx = global_idx * (size_t) nt + (size_t) it;
                size_t vrt_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                size_t div_idx = (size_t) (block_size + local_idx) * (size_t) nt + (size_t) it;
                vrt_out_data[output_idx] = solution_real[vrt_idx] + solution_imag[vrt_idx] * I;
                div_out_data[output_idx] = solution_real[div_idx] + solution_imag[div_idx] * I;
            }
        }

        start_idx += (size_t) block_size;
    }

done:
    PyMem_Free(basis_vrt_real); PyMem_Free(basis_vrt_imag);
    PyMem_Free(basis_div_real); PyMem_Free(basis_div_imag);
    PyMem_Free(u_basis_real); PyMem_Free(u_basis_imag);
    PyMem_Free(v_basis_real); PyMem_Free(v_basis_imag);
    PyMem_Free(weights);
    PyMem_Free(basis_real); PyMem_Free(basis_imag);
    PyMem_Free(reduced_basis_real); PyMem_Free(reduced_basis_imag);
    PyMem_Free(observed_real); PyMem_Free(observed_imag);
    PyMem_Free(solution_real); PyMem_Free(solution_imag);
    PyMem_Free(reduced_solution_real); PyMem_Free(reduced_solution_imag);
    PyMem_Free(active_columns);

    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    return Py_BuildValue("(NNi)", vrt_out_arr, div_out_arr, ierror);

fail:
    if (basis_vrt_real != NULL) PyMem_Free(basis_vrt_real);
    if (basis_vrt_imag != NULL) PyMem_Free(basis_vrt_imag);
    if (basis_div_real != NULL) PyMem_Free(basis_div_real);
    if (basis_div_imag != NULL) PyMem_Free(basis_div_imag);
    if (u_basis_real != NULL) PyMem_Free(u_basis_real);
    if (u_basis_imag != NULL) PyMem_Free(u_basis_imag);
    if (v_basis_real != NULL) PyMem_Free(v_basis_real);
    if (v_basis_imag != NULL) PyMem_Free(v_basis_imag);
    if (weights != NULL) PyMem_Free(weights);
    if (basis_real != NULL) PyMem_Free(basis_real);
    if (basis_imag != NULL) PyMem_Free(basis_imag);
    if (reduced_basis_real != NULL) PyMem_Free(reduced_basis_real);
    if (reduced_basis_imag != NULL) PyMem_Free(reduced_basis_imag);
    if (observed_real != NULL) PyMem_Free(observed_real);
    if (observed_imag != NULL) PyMem_Free(observed_imag);
    if (solution_real != NULL) PyMem_Free(solution_real);
    if (solution_imag != NULL) PyMem_Free(solution_imag);
    if (reduced_solution_real != NULL) PyMem_Free(reduced_solution_real);
    if (reduced_solution_imag != NULL) PyMem_Free(reduced_solution_imag);
    if (active_columns != NULL) PyMem_Free(active_columns);
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(vrt_out_arr);
    Py_XDECREF(div_out_arr);
    return NULL;
}

static PyObject *backend_uv_to_vordiv_block_native_with_setup(PyObject *self, PyObject *args) {
    PyObject *handle_obj = NULL;
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    UvToVordivBlockSetup *setup = NULL;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *vrt_out_arr = NULL;
    PyArrayObject *div_out_arr = NULL;
    npy_intp *spec_dims = NULL;
    int nt = 0, spec_rank = 0, ntrunc = -1, v_ntrunc = -1, ierror = 0;
    npy_complex128 *vrt_out_data = NULL, *div_out_data = NULL;
    double *weights = NULL;
    double *reduced_basis_real = NULL, *reduced_basis_imag = NULL;
    double *observed_real = NULL, *observed_imag = NULL;
    double *solution_real = NULL, *solution_imag = NULL;
    double *reduced_solution_real = NULL, *reduced_solution_imag = NULL;
    int *active_columns = NULL;
    size_t max_block = 0, max_nrow = 0, start_idx = 0;
    int m = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOO", &handle_obj, &u_obj, &v_obj)) {
        return NULL;
    }
    if (!get_uv_to_vordiv_block_setup_from_capsule(handle_obj, &setup)) {
        return NULL;
    }

    u_arr = require_array(u_obj, NPY_COMPLEX128);
    v_arr = require_array(v_obj, NPY_COMPLEX128);
    if (u_arr == NULL || v_arr == NULL) {
        goto fail;
    }
    if (!flatten_spectral_array(u_arr, INT_MAX, &ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (!flatten_spectral_array(v_arr, INT_MAX, &v_ntrunc, &nt, &spec_rank, &spec_dims)) {
        goto fail;
    }
    if (PyArray_NDIM(u_arr) != PyArray_NDIM(v_arr)) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same rank");
        goto fail;
    }
    if (!PyArray_CompareLists(
            PyArray_DIMS(u_arr), PyArray_DIMS(v_arr), PyArray_NDIM(u_arr))) {
        PyErr_SetString(PyExc_ValueError, "uspec and vspec must have the same shape");
        goto fail;
    }
    if (ntrunc != v_ntrunc || ntrunc != setup->ntrunc) {
        PyErr_SetString(PyExc_ValueError, "uspec/vspec ntrunc must match the setup handle");
        goto fail;
    }
    vrt_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    div_out_arr = (PyArrayObject *) PyArray_ZEROS(spec_rank, spec_dims, NPY_COMPLEX128, 0);
    if (vrt_out_arr == NULL || div_out_arr == NULL) {
        goto fail;
    }
    vrt_out_data = (npy_complex128 *) PyArray_DATA(vrt_out_arr);
    div_out_data = (npy_complex128 *) PyArray_DATA(div_out_arr);

    max_block = (size_t) (ntrunc + 1);
    max_nrow = 2U * max_block;
    weights = (double *) PyMem_Malloc(max_nrow * sizeof(double));
    reduced_basis_real = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    reduced_basis_imag = (double *) PyMem_Calloc(max_nrow * max_nrow, sizeof(double));
    observed_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    observed_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    solution_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    solution_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    reduced_solution_real = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    reduced_solution_imag = (double *) PyMem_Calloc(max_nrow * (size_t) nt, sizeof(double));
    active_columns = (int *) PyMem_Malloc(max_nrow * sizeof(int));
    if (
        weights == NULL || reduced_basis_real == NULL || reduced_basis_imag == NULL ||
        observed_real == NULL || observed_imag == NULL ||
        solution_real == NULL || solution_imag == NULL ||
        reduced_solution_real == NULL || reduced_solution_imag == NULL ||
        active_columns == NULL
    ) {
        PyErr_NoMemory();
        goto fail;
    }

    for (m = 0; m <= ntrunc; ++m) {
        int block_size = ntrunc - m + 1;
        int nrow = 2 * block_size;
        int ncol = 2 * block_size;
        int local_idx = 0;
        int active_count = setup->active_counts[m];
        int it = 0;
        size_t basis_offset = setup->basis_offsets[m];
        size_t active_offset = setup->active_offsets[m];
        size_t rhs_bytes = (size_t) nrow * (size_t) nt * sizeof(double);

        memset(reduced_basis_real, 0, (size_t) nrow * (size_t) nrow * sizeof(double));
        memset(reduced_basis_imag, 0, (size_t) nrow * (size_t) nrow * sizeof(double));
        memset(observed_real, 0, rhs_bytes);
        memset(observed_imag, 0, rhs_bytes);
        memset(solution_real, 0, rhs_bytes);
        memset(solution_imag, 0, rhs_bytes);
        memset(reduced_solution_real, 0, rhs_bytes);
        memset(reduced_solution_imag, 0, rhs_bytes);

        for (local_idx = 0; local_idx < nrow; ++local_idx) {
            weights[local_idx] = 1.0;
        }
        for (local_idx = 0; local_idx < active_count; ++local_idx) {
            active_columns[local_idx] = setup->active_columns[active_offset + (size_t) local_idx];
        }

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;
            for (it = 0; it < nt; ++it) {
                size_t input_idx = global_idx * (size_t) nt + (size_t) it;
                size_t u_obs_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                size_t v_obs_idx = (size_t) (block_size + local_idx) * (size_t) nt + (size_t) it;
                observed_real[u_obs_idx] = creal(((npy_complex128 *) PyArray_DATA(u_arr))[input_idx]);
                observed_imag[u_obs_idx] = cimag(((npy_complex128 *) PyArray_DATA(u_arr))[input_idx]);
                observed_real[v_obs_idx] = creal(((npy_complex128 *) PyArray_DATA(v_arr))[input_idx]);
                observed_imag[v_obs_idx] = cimag(((npy_complex128 *) PyArray_DATA(v_arr))[input_idx]);
            }
        }

        for (local_idx = 0; local_idx < nrow; ++local_idx) {
            int active_idx = 0;
            for (active_idx = 0; active_idx < active_count; ++active_idx) {
                int col = active_columns[active_idx];
                size_t source_idx = basis_offset + (size_t) local_idx * (size_t) ncol + (size_t) col;
                size_t target_idx = (size_t) local_idx * (size_t) active_count + (size_t) active_idx;
                reduced_basis_real[target_idx] = setup->basis_real[source_idx];
                reduced_basis_imag[target_idx] = setup->basis_imag[source_idx];
            }
        }

        if (active_count > 0) {
            weighted_block_solve_stub(
                nrow, active_count, nt, weights,
                reduced_basis_real, reduced_basis_imag,
                observed_real, observed_imag,
                reduced_solution_real, reduced_solution_imag,
                &ierror
            );
            if (ierror != 0) {
                PyErr_Format(
                    PyExc_RuntimeError,
                    "uv_to_vordiv_block_native_with_setup solve failed (ierror=%d, m=%d)",
                    ierror,
                    m
                );
                goto fail;
            }

            for (local_idx = 0; local_idx < active_count; ++local_idx) {
                int col = active_columns[local_idx];
                for (it = 0; it < nt; ++it) {
                    size_t source_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                    size_t target_idx = (size_t) col * (size_t) nt + (size_t) it;
                    solution_real[target_idx] = reduced_solution_real[source_idx];
                    solution_imag[target_idx] = reduced_solution_imag[source_idx];
                }
            }
        }

        for (local_idx = 0; local_idx < block_size; ++local_idx) {
            size_t global_idx = start_idx + (size_t) local_idx;
            for (it = 0; it < nt; ++it) {
                size_t output_idx = global_idx * (size_t) nt + (size_t) it;
                size_t vrt_idx = (size_t) local_idx * (size_t) nt + (size_t) it;
                size_t div_idx = (size_t) (block_size + local_idx) * (size_t) nt + (size_t) it;
                vrt_out_data[output_idx] = solution_real[vrt_idx] + solution_imag[vrt_idx] * I;
                div_out_data[output_idx] = solution_real[div_idx] + solution_imag[div_idx] * I;
            }
        }

        start_idx += (size_t) block_size;
    }

    Py_DECREF(u_arr);
    Py_DECREF(v_arr);
    PyMem_Free(weights);
    PyMem_Free(reduced_basis_real);
    PyMem_Free(reduced_basis_imag);
    PyMem_Free(observed_real);
    PyMem_Free(observed_imag);
    PyMem_Free(solution_real);
    PyMem_Free(solution_imag);
    PyMem_Free(reduced_solution_real);
    PyMem_Free(reduced_solution_imag);
    PyMem_Free(active_columns);
    return Py_BuildValue("(NNi)", vrt_out_arr, div_out_arr, 0);

fail:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    if (weights != NULL) PyMem_Free(weights);
    if (reduced_basis_real != NULL) PyMem_Free(reduced_basis_real);
    if (reduced_basis_imag != NULL) PyMem_Free(reduced_basis_imag);
    if (observed_real != NULL) PyMem_Free(observed_real);
    if (observed_imag != NULL) PyMem_Free(observed_imag);
    if (solution_real != NULL) PyMem_Free(solution_real);
    if (solution_imag != NULL) PyMem_Free(solution_imag);
    if (reduced_solution_real != NULL) PyMem_Free(reduced_solution_real);
    if (reduced_solution_imag != NULL) PyMem_Free(reduced_solution_imag);
    if (active_columns != NULL) PyMem_Free(active_columns);
    Py_XDECREF(vrt_out_arr);
    Py_XDECREF(div_out_arr);
    return NULL;
}

static PyObject *backend_ldfou2_uv_scaling(PyObject *self, PyObject *args) {
    PyObject *paia_obj = NULL;
    PyObject *psia_obj = NULL;
    PyArrayObject *paia_arr = NULL;
    PyArrayObject *psia_arr = NULL;
    PyArrayObject *paia_out_arr = NULL;
    PyArrayObject *psia_out_arr = NULL;
    double rsphere = 0.0;
    int ntrunc = -1, km = -1, kf_uv = -1, ierror = 0;
    npy_intp dims[2];

    (void) self;

    if (!PyArg_ParseTuple(args, "iiOdO", &ntrunc, &km, &paia_obj, &rsphere, &psia_obj)) {
        return NULL;
    }

    paia_arr = require_array(paia_obj, NPY_FLOAT64);
    psia_arr = require_array(psia_obj, NPY_FLOAT64);
    if (paia_arr == NULL || psia_arr == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(paia_arr) != 2 || PyArray_NDIM(psia_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "paia and psia must be rank-2 float64 arrays");
        goto fail;
    }
    if (!PyArray_CompareLists(PyArray_DIMS(paia_arr), PyArray_DIMS(psia_arr), 2)) {
        PyErr_SetString(PyExc_ValueError, "paia and psia must have the same shape");
        goto fail;
    }
    if (ntrunc < 0 || km < 0 || km > ntrunc) {
        PyErr_SetString(PyExc_ValueError, "require 0 <= km <= ntrunc");
        goto fail;
    }

    dims[0] = PyArray_DIM(paia_arr, 0);
    dims[1] = PyArray_DIM(paia_arr, 1);
    if (dims[0] % 4 != 0) {
        PyErr_SetString(PyExc_ValueError, "paia/psia first dimension must be divisible by 4");
        goto fail;
    }
    kf_uv = (int) (dims[0] / 4);
    if (dims[1] != (npy_intp) ((ntrunc + 2) / 2)) {
        PyErr_SetString(PyExc_ValueError, "second dimension must equal (ntrunc + 2) // 2");
        goto fail;
    }

    paia_out_arr = (PyArrayObject *) PyArray_ZEROS(2, dims, NPY_FLOAT64, 0);
    psia_out_arr = (PyArrayObject *) PyArray_ZEROS(2, dims, NPY_FLOAT64, 0);
    if (paia_out_arr == NULL || psia_out_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    ldfou2_uv_scaling(
        ntrunc,
        km,
        kf_uv,
        rsphere,
        (const double *) PyArray_DATA(paia_arr),
        (const double *) PyArray_DATA(psia_arr),
        (double *) PyArray_DATA(paia_out_arr),
        (double *) PyArray_DATA(psia_out_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(paia_arr);
    Py_DECREF(psia_arr);
    return Py_BuildValue("(NNi)", paia_out_arr, psia_out_arr, ierror);

fail:
    Py_XDECREF(paia_arr);
    Py_XDECREF(psia_arr);
    Py_XDECREF(paia_out_arr);
    Py_XDECREF(psia_out_arr);
    return NULL;
}

static PyObject *backend_ledir_dgemm(PyObject *self, PyObject *args) {
    PyObject *paia_obj = NULL;
    PyObject *psia_obj = NULL;
    PyObject *rpnma_obj = NULL;
    PyObject *rpnms_obj = NULL;
    PyObject *pw_obj = NULL;
    PyArrayObject *paia_arr = NULL;
    PyArrayObject *psia_arr = NULL;
    PyArrayObject *rpnma_arr = NULL;
    PyArrayObject *rpnms_arr = NULL;
    PyArrayObject *pw_arr = NULL;
    PyArrayObject *poa1_arr = NULL;
    int ntrunc = -1, km = -1, ierror = 0;
    int kfc = -1, kdglu = -1, ila = -1, ils = -1, nlei1 = -1;
    npy_intp out_dims[2];

    (void) self;

    if (!PyArg_ParseTuple(
            args,
            "iiOOOOO",
            &ntrunc,
            &km,
            &paia_obj,
            &psia_obj,
            &rpnma_obj,
            &rpnms_obj,
            &pw_obj)) {
        return NULL;
    }

    paia_arr = require_array(paia_obj, NPY_FLOAT64);
    psia_arr = require_array(psia_obj, NPY_FLOAT64);
    rpnma_arr = require_array(rpnma_obj, NPY_FLOAT64);
    rpnms_arr = require_array(rpnms_obj, NPY_FLOAT64);
    pw_arr = require_array(pw_obj, NPY_FLOAT64);
    if (paia_arr == NULL || psia_arr == NULL || rpnma_arr == NULL ||
        rpnms_arr == NULL || pw_arr == NULL) {
        goto fail;
    }

    if (ntrunc < 0 || km < 0 || km > ntrunc) {
        PyErr_SetString(PyExc_ValueError, "require 0 <= km <= ntrunc");
        goto fail;
    }
    if (PyArray_NDIM(paia_arr) != 2 || PyArray_NDIM(psia_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "paia and psia must be rank-2 float64 arrays");
        goto fail;
    }
    if (!PyArray_CompareLists(PyArray_DIMS(paia_arr), PyArray_DIMS(psia_arr), 2)) {
        PyErr_SetString(PyExc_ValueError, "paia and psia must have the same shape");
        goto fail;
    }
    if (PyArray_NDIM(rpnma_arr) != 2 || PyArray_NDIM(rpnms_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "rpnma and rpnms must be rank-2 float64 arrays");
        goto fail;
    }
    if (PyArray_NDIM(pw_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "pw must be rank-1 float64 array");
        goto fail;
    }

    kfc = (int) PyArray_DIM(paia_arr, 0);
    kdglu = (int) PyArray_DIM(paia_arr, 1);
    ila = (ntrunc - km + 2) / 2;
    ils = (ntrunc - km + 3) / 2;
    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2);

    if ((int) PyArray_DIM(rpnma_arr, 0) != kdglu || (int) PyArray_DIM(rpnma_arr, 1) != ila) {
        PyErr_SetString(PyExc_ValueError, "rpnma shape must be (kdglu, (ntrunc-km+2)//2)");
        goto fail;
    }
    if ((int) PyArray_DIM(rpnms_arr, 0) != kdglu || (int) PyArray_DIM(rpnms_arr, 1) != ils) {
        PyErr_SetString(PyExc_ValueError, "rpnms shape must be (kdglu, (ntrunc-km+3)//2)");
        goto fail;
    }
    if ((int) PyArray_DIM(pw_arr, 0) != kdglu) {
        PyErr_SetString(PyExc_ValueError, "pw size must equal kdglu");
        goto fail;
    }

    out_dims[0] = nlei1;
    out_dims[1] = kfc;
    poa1_arr = (PyArrayObject *) PyArray_ZEROS(2, out_dims, NPY_FLOAT64, 0);
    if (poa1_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    ledir_dgemm(
        ntrunc,
        km,
        kfc,
        kdglu,
        (const double *) PyArray_DATA(paia_arr),
        (const double *) PyArray_DATA(psia_arr),
        (const double *) PyArray_DATA(rpnma_arr),
        (const double *) PyArray_DATA(rpnms_arr),
        (const double *) PyArray_DATA(pw_arr),
        (double *) PyArray_DATA(poa1_arr),
        &ierror
    );
    Py_END_ALLOW_THREADS

    Py_DECREF(paia_arr);
    Py_DECREF(psia_arr);
    Py_DECREF(rpnma_arr);
    Py_DECREF(rpnms_arr);
    Py_DECREF(pw_arr);
    return Py_BuildValue("(Ni)", poa1_arr, ierror);

fail:
    Py_XDECREF(paia_arr);
    Py_XDECREF(psia_arr);
    Py_XDECREF(rpnma_arr);
    Py_XDECREF(rpnms_arr);
    Py_XDECREF(pw_arr);
    Py_XDECREF(poa1_arr);
    return NULL;
}

static PyObject *backend_poa1_to_vordiv(PyObject *self, PyObject *args) {
    PyObject *poa1_obj = NULL;
    PyArrayObject *poa1_arr = NULL;
    PyArrayObject *vrt_real_arr = NULL;
    PyArrayObject *vrt_imag_arr = NULL;
    PyArrayObject *div_real_arr = NULL;
    PyArrayObject *div_imag_arr = NULL;
    PyArrayObject *vrt_out_arr = NULL;
    PyArrayObject *div_out_arr = NULL;
    double rsphere = 0.0;
    int ntrunc = -1, km = -1, ierror = 0, kf_uv = -1, ncoeff = -1, nlei1 = -1;
    npy_intp spec_dims[1];
    npy_intp idx;
    npy_complex128 *vrt_out_data = NULL;
    npy_complex128 *div_out_data = NULL;
    double *vrt_real_data = NULL;
    double *vrt_imag_data = NULL;
    double *div_real_data = NULL;
    double *div_imag_data = NULL;

    (void) self;

    if (!PyArg_ParseTuple(args, "iiOd", &ntrunc, &km, &poa1_obj, &rsphere)) {
        return NULL;
    }

    poa1_arr = require_array(poa1_obj, NPY_FLOAT64);
    if (poa1_arr == NULL) {
        goto fail;
    }
    if (ntrunc < 0 || km < 0 || km > ntrunc) {
        PyErr_SetString(PyExc_ValueError, "require 0 <= km <= ntrunc");
        goto fail;
    }
    if (PyArray_NDIM(poa1_arr) != 2) {
        PyErr_SetString(PyExc_ValueError, "poa1 must be rank-2 float64 array");
        goto fail;
    }

    nlei1 = ntrunc + 4 + ((ntrunc + 5) % 2);
    if ((int) PyArray_DIM(poa1_arr, 0) != nlei1) {
        PyErr_SetString(PyExc_ValueError, "poa1 first dimension must equal ntrunc + 4 + ((ntrunc + 5) % 2)");
        goto fail;
    }
    if (((int) PyArray_DIM(poa1_arr, 1)) % 4 != 0) {
        PyErr_SetString(PyExc_ValueError, "poa1 second dimension must be divisible by 4");
        goto fail;
    }
    kf_uv = (int) PyArray_DIM(poa1_arr, 1) / 4;

    ncoeff = ((ntrunc + 1) * (ntrunc + 2)) / 2;
    spec_dims[0] = ncoeff;
    vrt_real_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_FLOAT64, 0);
    vrt_imag_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_FLOAT64, 0);
    div_real_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_FLOAT64, 0);
    div_imag_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_FLOAT64, 0);
    vrt_out_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_COMPLEX128, 0);
    div_out_arr = (PyArrayObject *) PyArray_ZEROS(1, spec_dims, NPY_COMPLEX128, 0);
    if (vrt_real_arr == NULL || vrt_imag_arr == NULL ||
        div_real_arr == NULL || div_imag_arr == NULL ||
        vrt_out_arr == NULL || div_out_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    poa1_to_vordiv(
        ntrunc,
        km,
        kf_uv,
        rsphere,
        (const double *) PyArray_DATA(poa1_arr),
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

    Py_DECREF(poa1_arr);
    Py_DECREF(vrt_real_arr);
    Py_DECREF(vrt_imag_arr);
    Py_DECREF(div_real_arr);
    Py_DECREF(div_imag_arr);
    return Py_BuildValue("(NNi)", vrt_out_arr, div_out_arr, ierror);

fail:
    Py_XDECREF(poa1_arr);
    Py_XDECREF(vrt_real_arr);
    Py_XDECREF(vrt_imag_arr);
    Py_XDECREF(div_real_arr);
    Py_XDECREF(div_imag_arr);
    Py_XDECREF(vrt_out_arr);
    Py_XDECREF(div_out_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"validate_nloen", backend_validate_nloen, METH_VARARGS, "Validate reduced-grid longitude counts."},
    {"create_setup", backend_create_setup, METH_VARARGS, "Create an experimental ectrans setup handle."},
    {"describe_setup", backend_describe_setup, METH_VARARGS, "Describe an experimental ectrans setup handle."},
    {"destroy_setup", backend_destroy_setup, METH_VARARGS, "Destroy an experimental ectrans setup handle."},
    {"create_uv_to_vordiv_block_setup", backend_create_uv_to_vordiv_block_setup, METH_VARARGS, "Create a cached blockwise uv-to-vordiv setup handle."},
    {"destroy_uv_to_vordiv_block_setup", backend_destroy_uv_to_vordiv_block_setup, METH_VARARGS, "Destroy a cached blockwise uv-to-vordiv setup handle."},
    {"scalar_analysis_stub", backend_scalar_analysis, METH_VARARGS, "Call the native scalar-analysis stub."},
    {"scalar_synthesis_stub", backend_scalar_synthesis, METH_VARARGS, "Call the native scalar-synthesis stub."},
    {"scalar_fourier_stub", backend_scalar_fourier, METH_VARARGS, "Call the native reduced-grid Fourier helper."},
    {"scalar_block_solve_stub", backend_scalar_block_solve, METH_VARARGS, "Call the native reduced-grid weighted block solver."},
    {"weighted_block_solve_stub", backend_weighted_block_solve, METH_VARARGS, "Call the native general weighted complex block solver."},
    {"vrtdiv_analysis_stub", backend_vrtdiv_analysis, METH_VARARGS, "Call the native vector-analysis stub."},
    {"uv_synthesis_stub", backend_uv_synthesis, METH_VARARGS, "Call the native wind-synthesis stub."},
    {"gradient_synthesis_stub", backend_gradient_synthesis, METH_VARARGS, "Call the native gradient-synthesis stub."},
    {"vordiv_to_uv", backend_vordiv_to_uv, METH_VARARGS, "Call the experimental ectrans vordiv-to-uv spectral helper."},
    {"prfi1b_uv_block", backend_prfi1b_uv_block, METH_VARARGS, "Build a PRFI1B-style UV POA1 block from public spectral u/v."},
    {"vd2uv_uv_block", backend_vd2uv_uv_block, METH_VARARGS, "Build a VD2UV internal UV POA1 block from public spectral vrt/div."},
    {"uv_to_vordiv", backend_uv_to_vordiv, METH_VARARGS, "Call the experimental ectrans uv-to-vordiv spectral helper."},
    {"uv_to_vordiv_block_native", backend_uv_to_vordiv_block_native, METH_VARARGS, "Call the experimental ectrans blockwise native uv-to-vordiv spectral helper."},
    {"uv_to_vordiv_block_native_with_setup", backend_uv_to_vordiv_block_native_with_setup, METH_VARARGS, "Call the cached experimental ectrans blockwise native uv-to-vordiv helper."},
    {"ldfou2_uv_scaling", backend_ldfou2_uv_scaling, METH_VARARGS, "Apply the experimental ectrans LDFOU2 uv scaling helper."},
    {"ledir_dgemm", backend_ledir_dgemm, METH_VARARGS, "Apply the experimental DGEMM-only OpenIFS LEDIR kernel."},
    {"poa1_to_vordiv", backend_poa1_to_vordiv, METH_VARARGS, "Apply the experimental OpenIFS UVTVD+UPDSPB stage to POA1."},
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
