#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <numpy/arrayobject.h>

void gaqd_c(int nlat, void *theta, void *wts, int *ierror);
void getlegfunc_c(void *legfunc, float lat, int ntrunc);
void specintrp_c(float rlon, int ntrunc, void *datnm, void *pnm, float *ob);
void lap_c(void *dataspec, void *dataspec_lap, int nmdim, int nt, float rsphere);
void invlap_c(void *dataspec, void *dataspec_ilap, int nmdim, int nt, float rsphere);
void multsmoothfact_c(void *dataspec, void *dataspec_smooth, void *smooth, int nlat, int nmdim, int nt);
void onedtotwod_c(void *dataspec, void *a, void *b, int nlat, int nmdim, int nt);
void onedtotwod_vrtdiv_c(void *vrtspec, void *divspec, void *br, void *bi, void *cr, void *ci, int nlat, int nmdim, int nt, float rsphere);
void onedtotwod_vrt_c(void *vrtspec, void *cr, void *ci, int nlat, int nmdim, int nt, float rsphere);
void onedtotwod_div_c(void *divspec, void *br, void *bi, int nlat, int nmdim, int nt, float rsphere);
void twodtooned_c(void *dataspec, void *a, void *b, int nlat, int ntrunc, int nt);
void twodtooned_vrtdiv_c(void *vrtspec, void *divspec, void *br, void *bi, void *cr, void *ci, int nlat, int ntrunc, int nt, float rsphere);
void twodtooned_vrt_c(void *vrtspec, void *cr, void *ci, int nlat, int ntrunc, int nt, float rsphere);
void twodtooned_div_c(void *divspec, void *br, void *bi, int nlat, int ntrunc, int nt, float rsphere);
void ihgeod_c(int m, void *x, void *y, void *z);
void shaesi_c(int nlat, int nlon, int lshaes, int lwork, int ldwork, void *wshaes, int *ierror);
void shsesi_c(int nlat, int nlon, int lshses, int lwork, int ldwork, void *wshses, int *ierror);
void shagsi_c(int nlat, int nlon, int lshags, int lwork, int ldwork, void *wshags, int *ierror);
void shsgsi_c(int nlat, int nlon, int lshsgs, int lwork, int ldwork, void *wshsgs, int *ierror);
void shaeci_c(int nlat, int nlon, int lshaec, int ldwork, void *wshaec, int *ierror);
void shseci_c(int nlat, int nlon, int lshsec, int ldwork, void *wshsec, int *ierror);
void shagci_c(int nlat, int nlon, int lshagc, int ldwork, void *wshagc, int *ierror);
void shsgci_c(int nlat, int nlon, int lshsgc, int ldwork, void *wshsgc, int *ierror);
void shaes_c(void *g, int nlat, int nlon, int nt, void *wshaes, int lshaes, int lwork, void *a, void *b, int *ierror);
void shags_c(void *g, int nlat, int nlon, int nt, void *wshags, int lshags, int lwork, void *a, void *b, int *ierror);
void shaec_c(void *g, int nlat, int nlon, int nt, void *wshaec, int lshaec, int lwork, void *a, void *b, int *ierror);
void shagc_c(void *g, int nlat, int nlon, int nt, void *wshagc, int lshagc, int lwork, void *a, void *b, int *ierror);
void shses_c(int nlon, void *a, void *b, int nlat, int nt, void *wshses, int lshses, int lwork, void *g, int *ierror);
void shsgs_c(int nlon, void *a, void *b, int nlat, int nt, void *wshsgs, int lshsgs, int lwork, void *g, int *ierror);
void shsec_c(int nlon, void *a, void *b, int nlat, int nt, void *wshsec, int lshsec, int lwork, void *g, int *ierror);
void shsgc_c(int nlon, void *a, void *b, int nlat, int nt, void *wshsgc, int lshsgc, int lwork, void *g, int *ierror);
void vhaesi_c(int nlat, int nlon, int lvhaes, int lwork, int ldwork, void *wvhaes, int *ierror);
void vhagsi_c(int nlat, int nlon, int lvhags, int ldwork, void *wvhags, int *ierror);
void vhsesi_c(int nlat, int nlon, int lvhses, int lwork, int ldwork, void *wvhses, int *ierror);
void vhsgsi_c(int nlat, int nlon, int lvhsgs, int ldwork, void *wvhsgs, int *ierror);
void vhaeci_c(int nlat, int nlon, int lvhaec, int ldwork, void *wvhaec, int *ierror);
void vhagci_c(int nlat, int nlon, int lvhagc, int ldwork, void *wvhagc, int *ierror);
void vhseci_c(int nlat, int nlon, int lvhsec, int ldwork, void *wvhsec, int *ierror);
void vhsgci_c(int nlat, int nlon, int lvhsgc, int ldwork, void *wvhsgc, int *ierror);
void vhaes_c(void *v, void *w, int nlat, int nlon, int nt, int ityp, void *wvhaes, int lvhaes, int lwork, void *br, void *bi, void *cr, void *ci, int *ierror);
void vhags_c(void *v, void *w, int nlat, int nlon, int nt, int ityp, void *wvhags, int lvhags, int lwork, void *br, void *bi, void *cr, void *ci, int *ierror);
void vhaec_c(void *v, void *w, int nlat, int nlon, int nt, int ityp, void *wvhaec, int lvhaec, int lwork, void *br, void *bi, void *cr, void *ci, int *ierror);
void vhagc_c(void *v, void *w, int nlat, int nlon, int nt, int ityp, void *wvhagc, int lvhagc, int lwork, void *br, void *bi, void *cr, void *ci, int *ierror);
void vhses_c(int nlon, void *br, void *bi, void *cr, void *ci, int nlat, int nt, int ityp, void *wvhses, int lvhses, int lwork, void *v, void *w, int *ierror);
void vhsgs_c(int nlon, void *br, void *bi, void *cr, void *ci, int nlat, int nt, int ityp, void *wvhsgs, int lvhsgs, int lwork, void *v, void *w, int *ierror);
void vhsec_c(int nlon, void *br, void *bi, void *cr, void *ci, int nlat, int nt, int ityp, void *wvhsec, int lvhsec, int lwork, void *v, void *w, int *ierror);
void vhsgc_c(int nlon, void *br, void *bi, void *cr, void *ci, int nlat, int nt, int ityp, void *wvhsgc, int lvhsgc, int lwork, void *v, void *w, int *ierror);

static PyArrayObject *as_real32_1d_aligned(PyObject *obj) {
    return (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED);
}

static int as_real32_3d_fortran(PyObject *obj, PyArrayObject **arr, int *d0, int *d1, int *d2) {
    *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (*arr == NULL) {
        return -1;
    }
    if (PyArray_NDIM(*arr) != 3) {
        PyErr_SetString(PyExc_ValueError, "expected a 3D float32 Fortran-contiguous array");
        Py_DECREF(*arr);
        *arr = NULL;
        return -1;
    }
    *d0 = (int)PyArray_DIM(*arr, 0);
    *d1 = (int)PyArray_DIM(*arr, 1);
    *d2 = (int)PyArray_DIM(*arr, 2);
    return 0;
}

static PyObject *py_gaqd(PyObject *self, PyObject *args) {
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

    (void)dwork_arr;
    (void)ldwork;
    gaqd_c(nlat, PyArray_DATA(theta_arr), PyArray_DATA(wts_arr), &ierror);

    result = Py_BuildValue("NNi", theta_arr, wts_arr, ierror);
    Py_DECREF(dwork_arr);
    return result;

fail:
    Py_XDECREF(theta_arr);
    Py_XDECREF(wts_arr);
    Py_XDECREF(dwork_arr);
    return NULL;
}

static PyObject *py_getlegfunc(PyObject *self, PyObject *args) {
    float lat;
    int ntrunc;
    int nmdim;
    npy_intp dims[1];
    PyArrayObject *out_arr = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "fi", &lat, &ntrunc)) {
        return NULL;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nmdim;
    out_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (out_arr == NULL) {
        return NULL;
    }

    getlegfunc_c(PyArray_DATA(out_arr), lat, ntrunc);
    return (PyObject *)out_arr;
}

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

static PyObject *py_lap(PyObject *self, PyObject *args) {
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
    lap_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), nmdim, nt, rsphere);
    Py_DECREF(spec_arr);
    return (PyObject *)out_arr;
}

static PyObject *py_invlap(PyObject *self, PyObject *args) {
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
    invlap_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), nmdim, nt, rsphere);
    Py_DECREF(spec_arr);
    return (PyObject *)out_arr;
}

static PyObject *py_multsmoothfact(PyObject *self, PyObject *args) {
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
    smooth_arr = (PyArrayObject *)PyArray_FROM_OTF(smooth_obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (smooth_arr == NULL) {
        Py_DECREF(spec_arr);
        return NULL;
    }
    if (PyArray_NDIM(smooth_arr) != 1) {
        PyErr_SetString(PyExc_ValueError, "smooth must be 1D float32");
        goto fail;
    }
    nlat = (int)PyArray_DIM(smooth_arr, 0);
    out_arr = (PyArrayObject *)PyArray_EMPTY(PyArray_NDIM(spec_arr), PyArray_DIMS(spec_arr), NPY_COMPLEX64, 1);
    if (out_arr == NULL) {
        goto fail;
    }
    multsmoothfact_c(PyArray_DATA(spec_arr), PyArray_DATA(out_arr), PyArray_DATA(smooth_arr), nlat, nmdim, nt);
    Py_DECREF(spec_arr);
    Py_DECREF(smooth_arr);
    return (PyObject *)out_arr;

fail:
    Py_XDECREF(spec_arr);
    Py_XDECREF(smooth_arr);
    Py_XDECREF(out_arr);
    return NULL;
}

static PyObject *py_specintrp(PyObject *self, PyObject *args) {
    float rlon;
    int ntrunc;
    float ob = 0.0f;
    PyObject *datnm_obj = NULL;
    PyObject *pnm_obj = NULL;
    PyArrayObject *datnm_arr = NULL;
    PyArrayObject *pnm_arr = NULL;
    PyArrayObject *scrm_arr = NULL;
    npy_intp dims[1];
    (void)self;

    if (!PyArg_ParseTuple(args, "fiOO", &rlon, &ntrunc, &datnm_obj, &pnm_obj)) {
        return NULL;
    }
    datnm_arr = (PyArrayObject *)PyArray_FROM_OTF(datnm_obj, NPY_COMPLEX64, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    pnm_arr = (PyArrayObject *)PyArray_FROM_OTF(pnm_obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS);
    if (datnm_arr == NULL || pnm_arr == NULL) {
        goto fail;
    }
    dims[0] = ntrunc + 1;
    scrm_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_COMPLEX64, 1);
    if (scrm_arr == NULL) {
        goto fail;
    }
    (void)scrm_arr;
    specintrp_c(rlon, ntrunc, PyArray_DATA(datnm_arr), PyArray_DATA(pnm_arr), &ob);
    Py_DECREF(datnm_arr);
    Py_DECREF(pnm_arr);
    Py_DECREF(scrm_arr);
    return PyFloat_FromDouble((double)ob);

fail:
    Py_XDECREF(datnm_arr);
    Py_XDECREF(pnm_arr);
    Py_XDECREF(scrm_arr);
    return NULL;
}

static PyObject *py_onedtotwod(PyObject *self, PyObject *args) {
    PyObject *spec_obj = NULL;
    PyArrayObject *spec_arr = NULL;
    PyArrayObject *a_arr = NULL;
    PyArrayObject *b_arr = NULL;
    int nlat;
    int nmdim;
    int nt;
    npy_intp dims[3];
    (void)self;

    if (!PyArg_ParseTuple(args, "Oi", &spec_obj, &nlat)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(spec_obj, &spec_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    dims[0] = nlat;
    dims[1] = nlat;
    dims[2] = nt;
    a_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    b_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (a_arr == NULL || b_arr == NULL) {
        goto fail;
    }
    onedtotwod_c(PyArray_DATA(spec_arr), PyArray_DATA(a_arr), PyArray_DATA(b_arr), nlat, nmdim, nt);
    Py_DECREF(spec_arr);
    return Py_BuildValue("NN", a_arr, b_arr);

fail:
    Py_XDECREF(spec_arr);
    Py_XDECREF(a_arr);
    Py_XDECREF(b_arr);
    return NULL;
}

static PyObject *py_onedtotwod_vrtdiv(PyObject *self, PyObject *args) {
    PyObject *vrt_obj = NULL;
    PyObject *div_obj = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *div_arr = NULL;
    PyArrayObject *br_arr = NULL;
    PyArrayObject *bi_arr = NULL;
    PyArrayObject *cr_arr = NULL;
    PyArrayObject *ci_arr = NULL;
    int nlat;
    int nmdim_vrt;
    int nt_vrt;
    int nmdim_div;
    int nt_div;
    float rsphere;
    npy_intp dims[3];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOif", &vrt_obj, &div_obj, &nlat, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(vrt_obj, &vrt_arr, &nmdim_vrt, &nt_vrt) != 0) {
        return NULL;
    }
    if (as_complex64_2d_fortran(div_obj, &div_arr, &nmdim_div, &nt_div) != 0) {
        Py_DECREF(vrt_arr);
        return NULL;
    }
    if (nmdim_vrt != nmdim_div || nt_vrt != nt_div) {
        PyErr_SetString(PyExc_ValueError, "vrtspec and divspec must have the same 2D shape");
        goto fail_vrtdiv;
    }
    dims[0] = nlat;
    dims[1] = nlat;
    dims[2] = nt_vrt;
    br_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    bi_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    cr_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    ci_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (br_arr == NULL || bi_arr == NULL || cr_arr == NULL || ci_arr == NULL) {
        goto fail_vrtdiv;
    }
    onedtotwod_vrtdiv_c(PyArray_DATA(vrt_arr), PyArray_DATA(div_arr), PyArray_DATA(br_arr), PyArray_DATA(bi_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), nlat, nmdim_vrt, nt_vrt, rsphere);
    Py_DECREF(vrt_arr);
    Py_DECREF(div_arr);
    return Py_BuildValue("NNNN", br_arr, bi_arr, cr_arr, ci_arr);

fail_vrtdiv:
    Py_XDECREF(vrt_arr);
    Py_XDECREF(div_arr);
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    return NULL;
}

static PyObject *py_onedtotwod_vrt(PyObject *self, PyObject *args) {
    PyObject *vrt_obj = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *cr_arr = NULL;
    PyArrayObject *ci_arr = NULL;
    int nlat;
    int nmdim;
    int nt;
    float rsphere;
    npy_intp dims[3];
    (void)self;

    if (!PyArg_ParseTuple(args, "Oif", &vrt_obj, &nlat, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(vrt_obj, &vrt_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    dims[0] = nlat;
    dims[1] = nlat;
    dims[2] = nt;
    cr_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    ci_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (cr_arr == NULL || ci_arr == NULL) {
        goto fail_vrt;
    }
    onedtotwod_vrt_c(PyArray_DATA(vrt_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), nlat, nmdim, nt, rsphere);
    Py_DECREF(vrt_arr);
    return Py_BuildValue("NN", cr_arr, ci_arr);

fail_vrt:
    Py_XDECREF(vrt_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    return NULL;
}

static PyObject *py_onedtotwod_div(PyObject *self, PyObject *args) {
    PyObject *div_obj = NULL;
    PyArrayObject *div_arr = NULL;
    PyArrayObject *br_arr = NULL;
    PyArrayObject *bi_arr = NULL;
    int nlat;
    int nmdim;
    int nt;
    float rsphere;
    npy_intp dims[3];
    (void)self;

    if (!PyArg_ParseTuple(args, "Oif", &div_obj, &nlat, &rsphere)) {
        return NULL;
    }
    if (as_complex64_2d_fortran(div_obj, &div_arr, &nmdim, &nt) != 0) {
        return NULL;
    }
    dims[0] = nlat;
    dims[1] = nlat;
    dims[2] = nt;
    br_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    bi_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (br_arr == NULL || bi_arr == NULL) {
        goto fail_div;
    }
    onedtotwod_div_c(PyArray_DATA(div_arr), PyArray_DATA(br_arr), PyArray_DATA(bi_arr), nlat, nmdim, nt, rsphere);
    Py_DECREF(div_arr);
    return Py_BuildValue("NN", br_arr, bi_arr);

fail_div:
    Py_XDECREF(div_arr);
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    return NULL;
}

static PyObject *py_twodtooned(PyObject *self, PyObject *args) {
    PyObject *a_obj = NULL;
    PyObject *b_obj = NULL;
    PyArrayObject *a_arr = NULL;
    PyArrayObject *b_arr = NULL;
    PyArrayObject *out_arr = NULL;
    int nlat_a;
    int nlat2_a;
    int nt_a;
    int nlat_b;
    int nlat2_b;
    int nt_b;
    int ntrunc;
    npy_intp dims[2];
    int nmdim;
    (void)self;

    if (!PyArg_ParseTuple(args, "OOi", &a_obj, &b_obj, &ntrunc)) {
        return NULL;
    }
    if (as_real32_3d_fortran(a_obj, &a_arr, &nlat_a, &nlat2_a, &nt_a) != 0) {
        return NULL;
    }
    if (as_real32_3d_fortran(b_obj, &b_arr, &nlat_b, &nlat2_b, &nt_b) != 0) {
        Py_DECREF(a_arr);
        return NULL;
    }
    if (nlat_a != nlat_b || nlat2_a != nlat2_b || nt_a != nt_b) {
        PyErr_SetString(PyExc_ValueError, "a and b must have the same 3D shape");
        goto fail_twod;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nmdim;
    dims[1] = nt_a;
    out_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (out_arr == NULL) {
        goto fail_twod;
    }
    twodtooned_c(PyArray_DATA(out_arr), PyArray_DATA(a_arr), PyArray_DATA(b_arr), nlat_a, ntrunc, nt_a);
    Py_DECREF(a_arr);
    Py_DECREF(b_arr);
    return (PyObject *)out_arr;

fail_twod:
    Py_XDECREF(a_arr);
    Py_XDECREF(b_arr);
    Py_XDECREF(out_arr);
    return NULL;
}

static PyObject *py_twodtooned_vrtdiv(PyObject *self, PyObject *args) {
    PyObject *br_obj = NULL;
    PyObject *bi_obj = NULL;
    PyObject *cr_obj = NULL;
    PyObject *ci_obj = NULL;
    PyArrayObject *br_arr = NULL;
    PyArrayObject *bi_arr = NULL;
    PyArrayObject *cr_arr = NULL;
    PyArrayObject *ci_arr = NULL;
    PyArrayObject *vrt_arr = NULL;
    PyArrayObject *div_arr = NULL;
    int nlat_br, nlat2_br, nt_br;
    int nlat_bi, nlat2_bi, nt_bi;
    int nlat_cr, nlat2_cr, nt_cr;
    int nlat_ci, nlat2_ci, nt_ci;
    int ntrunc;
    int nmdim;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOif", &br_obj, &bi_obj, &cr_obj, &ci_obj, &ntrunc, &rsphere)) {
        return NULL;
    }
    if (as_real32_3d_fortran(br_obj, &br_arr, &nlat_br, &nlat2_br, &nt_br) != 0) return NULL;
    if (as_real32_3d_fortran(bi_obj, &bi_arr, &nlat_bi, &nlat2_bi, &nt_bi) != 0) goto fail_twod_vrtdiv;
    if (as_real32_3d_fortran(cr_obj, &cr_arr, &nlat_cr, &nlat2_cr, &nt_cr) != 0) goto fail_twod_vrtdiv;
    if (as_real32_3d_fortran(ci_obj, &ci_arr, &nlat_ci, &nlat2_ci, &nt_ci) != 0) goto fail_twod_vrtdiv;
    if (nlat_br != nlat_bi || nlat_br != nlat_cr || nlat_br != nlat_ci || nlat2_br != nlat2_bi || nlat2_br != nlat2_cr || nlat2_br != nlat2_ci || nt_br != nt_bi || nt_br != nt_cr || nt_br != nt_ci) {
        PyErr_SetString(PyExc_ValueError, "br, bi, cr, and ci must have the same 3D shape");
        goto fail_twod_vrtdiv;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nmdim;
    dims[1] = nt_br;
    vrt_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    div_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (vrt_arr == NULL || div_arr == NULL) goto fail_twod_vrtdiv;
    twodtooned_vrtdiv_c(PyArray_DATA(vrt_arr), PyArray_DATA(div_arr), PyArray_DATA(br_arr), PyArray_DATA(bi_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), nlat_br, ntrunc, nt_br, rsphere);
    Py_DECREF(br_arr); Py_DECREF(bi_arr); Py_DECREF(cr_arr); Py_DECREF(ci_arr);
    return Py_BuildValue("NN", vrt_arr, div_arr);

fail_twod_vrtdiv:
    Py_XDECREF(br_arr); Py_XDECREF(bi_arr); Py_XDECREF(cr_arr); Py_XDECREF(ci_arr);
    Py_XDECREF(vrt_arr); Py_XDECREF(div_arr);
    return NULL;
}

static PyObject *py_twodtooned_vrt(PyObject *self, PyObject *args) {
    PyObject *cr_obj = NULL;
    PyObject *ci_obj = NULL;
    PyArrayObject *cr_arr = NULL;
    PyArrayObject *ci_arr = NULL;
    PyArrayObject *vrt_arr = NULL;
    int nlat_cr, nlat2_cr, nt_cr;
    int nlat_ci, nlat2_ci, nt_ci;
    int ntrunc;
    int nmdim;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOif", &cr_obj, &ci_obj, &ntrunc, &rsphere)) return NULL;
    if (as_real32_3d_fortran(cr_obj, &cr_arr, &nlat_cr, &nlat2_cr, &nt_cr) != 0) return NULL;
    if (as_real32_3d_fortran(ci_obj, &ci_arr, &nlat_ci, &nlat2_ci, &nt_ci) != 0) goto fail_twod_vrt;
    if (nlat_cr != nlat_ci || nlat2_cr != nlat2_ci || nt_cr != nt_ci) {
        PyErr_SetString(PyExc_ValueError, "cr and ci must have the same 3D shape");
        goto fail_twod_vrt;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nmdim;
    dims[1] = nt_cr;
    vrt_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (vrt_arr == NULL) goto fail_twod_vrt;
    twodtooned_vrt_c(PyArray_DATA(vrt_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), nlat_cr, ntrunc, nt_cr, rsphere);
    Py_DECREF(cr_arr); Py_DECREF(ci_arr);
    return (PyObject *)vrt_arr;

fail_twod_vrt:
    Py_XDECREF(cr_arr); Py_XDECREF(ci_arr); Py_XDECREF(vrt_arr);
    return NULL;
}

static PyObject *py_twodtooned_div(PyObject *self, PyObject *args) {
    PyObject *br_obj = NULL;
    PyObject *bi_obj = NULL;
    PyArrayObject *br_arr = NULL;
    PyArrayObject *bi_arr = NULL;
    PyArrayObject *div_arr = NULL;
    int nlat_br, nlat2_br, nt_br;
    int nlat_bi, nlat2_bi, nt_bi;
    int ntrunc;
    int nmdim;
    float rsphere;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOif", &br_obj, &bi_obj, &ntrunc, &rsphere)) return NULL;
    if (as_real32_3d_fortran(br_obj, &br_arr, &nlat_br, &nlat2_br, &nt_br) != 0) return NULL;
    if (as_real32_3d_fortran(bi_obj, &bi_arr, &nlat_bi, &nlat2_bi, &nt_bi) != 0) goto fail_twod_div;
    if (nlat_br != nlat_bi || nlat2_br != nlat2_bi || nt_br != nt_bi) {
        PyErr_SetString(PyExc_ValueError, "br and bi must have the same 3D shape");
        goto fail_twod_div;
    }
    nmdim = (ntrunc + 1) * (ntrunc + 2) / 2;
    dims[0] = nmdim;
    dims[1] = nt_br;
    div_arr = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_COMPLEX64, 1);
    if (div_arr == NULL) goto fail_twod_div;
    twodtooned_div_c(PyArray_DATA(div_arr), PyArray_DATA(br_arr), PyArray_DATA(bi_arr), nlat_br, ntrunc, nt_br, rsphere);
    Py_DECREF(br_arr); Py_DECREF(bi_arr);
    return (PyObject *)div_arr;

fail_twod_div:
    Py_XDECREF(br_arr); Py_XDECREF(bi_arr); Py_XDECREF(div_arr);
    return NULL;
}

static PyObject *py_ihgeod(PyObject *self, PyObject *args) {
    int m;
    int idp;
    int jdp;
    npy_intp dims[3];
    PyArrayObject *x_arr = NULL;
    PyArrayObject *y_arr = NULL;
    PyArrayObject *z_arr = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "i", &m)) {
        return NULL;
    }
    idp = 2 * m - 1;
    jdp = m;
    dims[0] = idp;
    dims[1] = jdp;
    dims[2] = 5;
    x_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    y_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    z_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (x_arr == NULL || y_arr == NULL || z_arr == NULL) {
        goto fail_ihgeod;
    }
    (void)idp;
    (void)jdp;
    ihgeod_c(m, PyArray_DATA(x_arr), PyArray_DATA(y_arr), PyArray_DATA(z_arr));
    return Py_BuildValue("NNN", x_arr, y_arr, z_arr);

fail_ihgeod:
    Py_XDECREF(x_arr);
    Py_XDECREF(y_arr);
    Py_XDECREF(z_arr);
    return NULL;
}

static PyObject *build_init_array_result(int ierror, PyArrayObject *arr) {
    return Py_BuildValue("Ni", arr, ierror);
}

static PyObject *py_shaesi(PyObject *self, PyObject *args) {
    int nlat, nlon, lshaes, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;

    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lshaes, &lwork, &ldwork)) return NULL;
    dims[0] = lshaes;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shaesi_c(nlat, nlon, lshaes, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shsesi(PyObject *self, PyObject *args) {
    int nlat, nlon, lshses, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lshses, &lwork, &ldwork)) return NULL;
    dims[0] = lshses;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shsesi_c(nlat, nlon, lshses, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shagsi(PyObject *self, PyObject *args) {
    int nlat, nlon, lshags, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lshags, &lwork, &ldwork)) return NULL;
    dims[0] = lshags;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shagsi_c(nlat, nlon, lshags, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shsgsi(PyObject *self, PyObject *args) {
    int nlat, nlon, lshsgs, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lshsgs, &lwork, &ldwork)) return NULL;
    dims[0] = lshsgs;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shsgsi_c(nlat, nlon, lshsgs, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shaeci(PyObject *self, PyObject *args) {
    int nlat, nlon, lshaec, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lshaec, &ldwork)) return NULL;
    dims[0] = lshaec;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shaeci_c(nlat, nlon, lshaec, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shseci(PyObject *self, PyObject *args) {
    int nlat, nlon, lshsec, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lshsec, &ldwork)) return NULL;
    dims[0] = lshsec;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shseci_c(nlat, nlon, lshsec, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shagci(PyObject *self, PyObject *args) {
    int nlat, nlon, lshagc, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lshagc, &ldwork)) return NULL;
    dims[0] = lshagc;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shagci_c(nlat, nlon, lshagc, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shsgci(PyObject *self, PyObject *args) {
    int nlat, nlon, lshsgc, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lshsgc, &ldwork)) return NULL;
    dims[0] = lshsgc;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    shsgci_c(nlat, nlon, lshsgc, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_shaes_like(PyObject *g_obj, PyObject *w_obj, int lwork, void (*func)(void *, int, int, int, void *, int, int, void *, void *, int *)) {
    PyArrayObject *g_arr = NULL;
    PyArrayObject *w_arr = NULL;
    PyArrayObject *a_arr = NULL;
    PyArrayObject *b_arr = NULL;
    int nlat, nlon, nt;
    npy_intp dims[3];
    int ierror = 0;

    g_arr = (PyArrayObject *)PyArray_FROM_OTF(g_obj, NPY_FLOAT32, NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS);
    if (g_arr == NULL) return NULL;
    if (PyArray_NDIM(g_arr) != 3) {
        PyErr_SetString(PyExc_ValueError, "expected a 3D float32 Fortran-contiguous grid array");
        goto fail;
    }
    nlat = (int)PyArray_DIM(g_arr, 0);
    nlon = (int)PyArray_DIM(g_arr, 1);
    nt = (int)PyArray_DIM(g_arr, 2);
    w_arr = as_real32_1d_aligned(w_obj);
    if (w_arr == NULL) goto fail;
    dims[0] = nlat;
    dims[1] = nlat;
    dims[2] = nt;
    a_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    b_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (a_arr == NULL || b_arr == NULL) goto fail;
    func(PyArray_DATA(g_arr), nlat, nlon, nt, PyArray_DATA(w_arr), (int)PyArray_DIM(w_arr, 0), lwork, PyArray_DATA(a_arr), PyArray_DATA(b_arr), &ierror);
    Py_DECREF(g_arr);
    Py_DECREF(w_arr);
    return Py_BuildValue("NNi", a_arr, b_arr, ierror);

fail:
    Py_XDECREF(g_arr);
    Py_XDECREF(w_arr);
    Py_XDECREF(a_arr);
    Py_XDECREF(b_arr);
    return NULL;
}

static PyObject *py_shaes(PyObject *self, PyObject *args) {
    PyObject *g_obj = NULL, *w_obj = NULL;
    int lwork;
    (void)self;
    if (!PyArg_ParseTuple(args, "OOi", &g_obj, &w_obj, &lwork)) return NULL;
    return py_shaes_like(g_obj, w_obj, lwork, shaes_c);
}

static PyObject *py_shags(PyObject *self, PyObject *args) {
    PyObject *g_obj = NULL, *w_obj = NULL;
    int lwork;
    (void)self;
    if (!PyArg_ParseTuple(args, "OOi", &g_obj, &w_obj, &lwork)) return NULL;
    return py_shaes_like(g_obj, w_obj, lwork, shags_c);
}

static PyObject *py_shaec(PyObject *self, PyObject *args) {
    PyObject *g_obj = NULL, *w_obj = NULL;
    int lwork;
    (void)self;
    if (!PyArg_ParseTuple(args, "OOi", &g_obj, &w_obj, &lwork)) return NULL;
    return py_shaes_like(g_obj, w_obj, lwork, shaec_c);
}

static PyObject *py_shagc(PyObject *self, PyObject *args) {
    PyObject *g_obj = NULL, *w_obj = NULL;
    int lwork;
    (void)self;
    if (!PyArg_ParseTuple(args, "OOi", &g_obj, &w_obj, &lwork)) return NULL;
    return py_shaes_like(g_obj, w_obj, lwork, shagc_c);
}

static PyObject *py_shses_like(int nlon, PyObject *a_obj, PyObject *b_obj, PyObject *w_obj, int lwork, void (*func)(int, void *, void *, int, int, void *, int, int, void *, int *)) {
    PyArrayObject *a_arr = NULL;
    PyArrayObject *b_arr = NULL;
    PyArrayObject *w_arr = NULL;
    PyArrayObject *g_arr = NULL;
    int nlat_a, nlat2_a, nt_a;
    int nlat_b, nlat2_b, nt_b;
    int ierror = 0;
    npy_intp dims[3];

    if (as_real32_3d_fortran(a_obj, &a_arr, &nlat_a, &nlat2_a, &nt_a) != 0) return NULL;
    if (as_real32_3d_fortran(b_obj, &b_arr, &nlat_b, &nlat2_b, &nt_b) != 0) goto fail;
    if (nlat_a != nlat_b || nlat2_a != nlat2_b || nt_a != nt_b) {
        PyErr_SetString(PyExc_ValueError, "a and b must have the same 3D shape");
        goto fail;
    }
    w_arr = as_real32_1d_aligned(w_obj);
    if (w_arr == NULL) goto fail;
    dims[0] = nlat_a;
    dims[1] = nlon;
    dims[2] = nt_a;
    g_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (g_arr == NULL) goto fail;
    func(nlon, PyArray_DATA(a_arr), PyArray_DATA(b_arr), nlat_a, nt_a, PyArray_DATA(w_arr), (int)PyArray_DIM(w_arr, 0), lwork, PyArray_DATA(g_arr), &ierror);
    Py_DECREF(a_arr);
    Py_DECREF(b_arr);
    Py_DECREF(w_arr);
    return Py_BuildValue("Ni", g_arr, ierror);

fail:
    Py_XDECREF(a_arr);
    Py_XDECREF(b_arr);
    Py_XDECREF(w_arr);
    Py_XDECREF(g_arr);
    return NULL;
}

static PyObject *py_shses(PyObject *self, PyObject *args) {
    int nlon, lwork;
    PyObject *a_obj = NULL, *b_obj = NULL, *w_obj = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iOOOi", &nlon, &a_obj, &b_obj, &w_obj, &lwork)) return NULL;
    return py_shses_like(nlon, a_obj, b_obj, w_obj, lwork, shses_c);
}

static PyObject *py_shsgs(PyObject *self, PyObject *args) {
    int nlon, lwork;
    PyObject *a_obj = NULL, *b_obj = NULL, *w_obj = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iOOOi", &nlon, &a_obj, &b_obj, &w_obj, &lwork)) return NULL;
    return py_shses_like(nlon, a_obj, b_obj, w_obj, lwork, shsgs_c);
}

static PyObject *py_shsec(PyObject *self, PyObject *args) {
    int nlon, lwork;
    PyObject *a_obj = NULL, *b_obj = NULL, *w_obj = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iOOOi", &nlon, &a_obj, &b_obj, &w_obj, &lwork)) return NULL;
    return py_shses_like(nlon, a_obj, b_obj, w_obj, lwork, shsec_c);
}

static PyObject *py_shsgc(PyObject *self, PyObject *args) {
    int nlon, lwork;
    PyObject *a_obj = NULL, *b_obj = NULL, *w_obj = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iOOOi", &nlon, &a_obj, &b_obj, &w_obj, &lwork)) return NULL;
    return py_shses_like(nlon, a_obj, b_obj, w_obj, lwork, shsgc_c);
}

static PyObject *py_vhaesi(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhaes, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lvhaes, &lwork, &ldwork)) return NULL;
    dims[0] = lvhaes;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhaesi_c(nlat, nlon, lvhaes, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhagsi(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhags, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhags, &ldwork)) return NULL;
    dims[0] = lvhags;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhagsi_c(nlat, nlon, lvhags, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhsesi(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhses, lwork, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiiii", &nlat, &nlon, &lvhses, &lwork, &ldwork)) return NULL;
    dims[0] = lvhses;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhsesi_c(nlat, nlon, lvhses, lwork, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhsgsi(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhsgs, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhsgs, &ldwork)) return NULL;
    dims[0] = lvhsgs;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhsgsi_c(nlat, nlon, lvhsgs, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhaeci(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhaec, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhaec, &ldwork)) return NULL;
    dims[0] = lvhaec;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhaeci_c(nlat, nlon, lvhaec, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhagci(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhagc, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhagc, &ldwork)) return NULL;
    dims[0] = lvhagc;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhagci_c(nlat, nlon, lvhagc, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhseci(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhsec, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhsec, &ldwork)) return NULL;
    dims[0] = lvhsec;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhseci_c(nlat, nlon, lvhsec, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vhsgci(PyObject *self, PyObject *args) {
    int nlat, nlon, lvhsgc, ldwork, ierror = 0;
    npy_intp dims[1];
    PyArrayObject *w_arr = NULL;
    (void)self;
    if (!PyArg_ParseTuple(args, "iiii", &nlat, &nlon, &lvhsgc, &ldwork)) return NULL;
    dims[0] = lvhsgc;
    w_arr = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_FLOAT32, 1);
    if (w_arr == NULL) return NULL;
    vhsgci_c(nlat, nlon, lvhsgc, ldwork, PyArray_DATA(w_arr), &ierror);
    return build_init_array_result(ierror, w_arr);
}

static PyObject *py_vha_like(PyObject *args, int default_ityp,
                             void (*func)(void *, void *, int, int, int, int, void *, int, int, void *, void *, void *, void *, int *)) {
    PyObject *v_obj = NULL, *w_obj = NULL, *wvha_obj = NULL;
    PyArrayObject *v_arr = NULL, *w_arr = NULL, *wvha_arr = NULL;
    PyArrayObject *br_arr = NULL, *bi_arr = NULL, *cr_arr = NULL, *ci_arr = NULL;
    int nlat_v, nlon_v, nt_v;
    int nlat_w, nlon_w, nt_w;
    int lwork, ityp = default_ityp, ierror = 0;
    npy_intp dims[3];

    if (!PyArg_ParseTuple(args, "OOOii", &v_obj, &w_obj, &wvha_obj, &lwork, &ityp)) return NULL;
    if (as_real32_3d_fortran(v_obj, &v_arr, &nlat_v, &nlon_v, &nt_v) != 0) return NULL;
    if (as_real32_3d_fortran(w_obj, &w_arr, &nlat_w, &nlon_w, &nt_w) != 0) goto fail;
    if (nlat_v != nlat_w || nlon_v != nlon_w || nt_v != nt_w) {
        PyErr_SetString(PyExc_ValueError, "v and w must have the same 3D shape");
        goto fail;
    }
    wvha_arr = as_real32_1d_aligned(wvha_obj);
    if (wvha_arr == NULL) goto fail;
    dims[0] = nlat_v;
    dims[1] = nlat_v;
    dims[2] = nt_v;
    br_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    bi_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    cr_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    ci_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (br_arr == NULL || bi_arr == NULL || cr_arr == NULL || ci_arr == NULL) goto fail;
    func(PyArray_DATA(v_arr), PyArray_DATA(w_arr), nlat_v, nlon_v, nt_v, ityp, PyArray_DATA(wvha_arr),
         (int)PyArray_DIM(wvha_arr, 0), lwork, PyArray_DATA(br_arr), PyArray_DATA(bi_arr),
         PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), &ierror);
    Py_DECREF(v_arr);
    Py_DECREF(w_arr);
    Py_DECREF(wvha_arr);
    return Py_BuildValue("NNNNi", br_arr, bi_arr, cr_arr, ci_arr, ierror);

fail:
    Py_XDECREF(v_arr);
    Py_XDECREF(w_arr);
    Py_XDECREF(wvha_arr);
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    return NULL;
}

static PyObject *py_vhaes(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_like(args, 0, vhaes_c);
}

static PyObject *py_vhags(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_like(args, 0, vhags_c);
}

static PyObject *py_vhaec(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_like(args, 0, vhaec_c);
}

static PyObject *py_vhagc(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_like(args, 0, vhagc_c);
}

static PyObject *py_vha_component_like(PyObject *args, int ityp, int needs_special_lwork,
                                       void (*func)(void *, void *, int, int, int, int, void *, int, int, void *, void *, void *, void *, int *),
                                       int return_div) {
    PyObject *v_obj = NULL, *w_obj = NULL, *wvha_obj = NULL;
    PyArrayObject *v_arr = NULL, *w_arr = NULL, *wvha_arr = NULL;
    PyArrayObject *br_arr = NULL, *bi_arr = NULL, *cr_arr = NULL, *ci_arr = NULL;
    int nlat_v, nlon_v, nt_v;
    int nlat_w, nlon_w, nt_w;
    int lwork = 0, ierror = 0;
    npy_intp dims_full[3] = {1, 1, 1};
    npy_intp dims_dummy[3] = {1, 1, 1};

    if (needs_special_lwork) {
        if (!PyArg_ParseTuple(args, "OOO", &v_obj, &w_obj, &wvha_obj)) return NULL;
    } else {
        if (!PyArg_ParseTuple(args, "OOOi", &v_obj, &w_obj, &wvha_obj, &lwork)) return NULL;
    }
    if (as_real32_3d_fortran(v_obj, &v_arr, &nlat_v, &nlon_v, &nt_v) != 0) return NULL;
    if (as_real32_3d_fortran(w_obj, &w_arr, &nlat_w, &nlon_w, &nt_w) != 0) goto fail;
    if (nlat_v != nlat_w || nlon_v != nlon_w || nt_v != nt_w) {
        PyErr_SetString(PyExc_ValueError, "v and w must have the same 3D shape");
        goto fail;
    }
    wvha_arr = as_real32_1d_aligned(wvha_obj);
    if (wvha_arr == NULL) goto fail;
    if (needs_special_lwork) {
        const int nbatch = 24;
        lwork = (2 * nbatch + 1) * nlat_v * nlon_v;
    }
    dims_full[0] = nlat_v;
    dims_full[1] = nlat_v;
    dims_full[2] = nt_v;
    if (return_div) {
        br_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_full, NPY_FLOAT32, 1);
        bi_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_full, NPY_FLOAT32, 1);
        cr_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_dummy, NPY_FLOAT32, 1);
        ci_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_dummy, NPY_FLOAT32, 1);
    } else {
        br_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_dummy, NPY_FLOAT32, 1);
        bi_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_dummy, NPY_FLOAT32, 1);
        cr_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_full, NPY_FLOAT32, 1);
        ci_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_full, NPY_FLOAT32, 1);
    }
    if (br_arr == NULL || bi_arr == NULL || cr_arr == NULL || ci_arr == NULL) goto fail;
    func(PyArray_DATA(v_arr), PyArray_DATA(w_arr), nlat_v, nlon_v, nt_v, ityp, PyArray_DATA(wvha_arr),
         (int)PyArray_DIM(wvha_arr, 0), lwork, PyArray_DATA(br_arr), PyArray_DATA(bi_arr),
         PyArray_DATA(cr_arr), PyArray_DATA(ci_arr), &ierror);
    Py_DECREF(v_arr);
    Py_DECREF(w_arr);
    Py_DECREF(wvha_arr);
    if (return_div) {
        Py_DECREF(cr_arr);
        Py_DECREF(ci_arr);
        return Py_BuildValue("NNi", br_arr, bi_arr, ierror);
    }
    Py_DECREF(br_arr);
    Py_DECREF(bi_arr);
    return Py_BuildValue("NNi", cr_arr, ci_arr, ierror);

fail:
    Py_XDECREF(v_arr);
    Py_XDECREF(w_arr);
    Py_XDECREF(wvha_arr);
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    return NULL;
}

static PyObject *py_vhaesdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 1, 1, vhaes_c, 1);
}

static PyObject *py_vhaesvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 2, 0, vhaes_c, 0);
}

static PyObject *py_vhagsdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 1, 0, vhags_c, 1);
}

static PyObject *py_vhagsvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 2, 0, vhags_c, 0);
}

static PyObject *py_vhaecdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 1, 0, vhaec_c, 1);
}

static PyObject *py_vhaecvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 2, 0, vhaec_c, 0);
}

static PyObject *py_vhagcdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 1, 0, vhagc_c, 1);
}

static PyObject *py_vhagcvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vha_component_like(args, 2, 0, vhagc_c, 0);
}

static PyObject *py_vhs_like(PyObject *args, int default_ityp,
                             void (*func)(int, void *, void *, void *, void *, int, int, int, void *, int, int, void *, void *, int *)) {
    int nlon, lwork, ityp = default_ityp;
    PyObject *br_obj = NULL, *bi_obj = NULL, *cr_obj = NULL, *ci_obj = NULL, *wvh_obj = NULL;
    PyArrayObject *br_arr = NULL, *bi_arr = NULL, *cr_arr = NULL, *ci_arr = NULL, *wvh_arr = NULL;
    PyArrayObject *v_arr = NULL, *w_arr = NULL;
    int nlat_br, nlat2_br, nt_br;
    int nlat_bi, nlat2_bi, nt_bi;
    int nlat_cr, nlat2_cr, nt_cr;
    int nlat_ci, nlat2_ci, nt_ci;
    int ierror = 0;
    npy_intp dims[3];

    if (!PyArg_ParseTuple(args, "iOOOOOii", &nlon, &br_obj, &bi_obj, &cr_obj, &ci_obj, &wvh_obj, &lwork, &ityp)) return NULL;
    if (as_real32_3d_fortran(br_obj, &br_arr, &nlat_br, &nlat2_br, &nt_br) != 0) return NULL;
    if (as_real32_3d_fortran(bi_obj, &bi_arr, &nlat_bi, &nlat2_bi, &nt_bi) != 0) goto fail;
    if (as_real32_3d_fortran(cr_obj, &cr_arr, &nlat_cr, &nlat2_cr, &nt_cr) != 0) goto fail;
    if (as_real32_3d_fortran(ci_obj, &ci_arr, &nlat_ci, &nlat2_ci, &nt_ci) != 0) goto fail;
    if (nlat_br != nlat_bi || nlat_br != nlat_cr || nlat_br != nlat_ci ||
        nlat2_br != nlat2_bi || nlat2_br != nlat2_cr || nlat2_br != nlat2_ci ||
        nt_br != nt_bi || nt_br != nt_cr || nt_br != nt_ci) {
        PyErr_SetString(PyExc_ValueError, "br, bi, cr, and ci must have the same 3D shape");
        goto fail;
    }
    wvh_arr = as_real32_1d_aligned(wvh_obj);
    if (wvh_arr == NULL) goto fail;
    dims[0] = nlat_br;
    dims[1] = nlon;
    dims[2] = nt_br;
    v_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    w_arr = (PyArrayObject *)PyArray_EMPTY(3, dims, NPY_FLOAT32, 1);
    if (v_arr == NULL || w_arr == NULL) goto fail;
    func(nlon, PyArray_DATA(br_arr), PyArray_DATA(bi_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr),
         nlat_br, nt_br, ityp, PyArray_DATA(wvh_arr), (int)PyArray_DIM(wvh_arr, 0), lwork,
         PyArray_DATA(v_arr), PyArray_DATA(w_arr), &ierror);
    Py_DECREF(br_arr);
    Py_DECREF(bi_arr);
    Py_DECREF(cr_arr);
    Py_DECREF(ci_arr);
    Py_DECREF(wvh_arr);
    return Py_BuildValue("NNi", v_arr, w_arr, ierror);

fail:
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    Py_XDECREF(wvh_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(w_arr);
    return NULL;
}

static PyObject *py_vhses(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_like(args, 0, vhses_c);
}

static PyObject *py_vhsgs(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_like(args, 0, vhsgs_c);
}

static PyObject *py_vhsec(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_like(args, 0, vhsec_c);
}

static PyObject *py_vhsgc(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_like(args, 0, vhsgc_c);
}

static PyObject *py_vhs_component_like(PyObject *args, int ityp,
                                       void (*func)(int, void *, void *, void *, void *, int, int, int, void *, int, int, void *, void *, int *),
                                       int is_div) {
    int nlon, lwork;
    PyObject *a_obj = NULL, *b_obj = NULL, *wvh_obj = NULL;
    PyArrayObject *a_arr = NULL, *b_arr = NULL, *wvh_arr = NULL;
    PyArrayObject *br_arr = NULL, *bi_arr = NULL, *cr_arr = NULL, *ci_arr = NULL;
    PyArrayObject *v_arr = NULL, *w_arr = NULL;
    int nlat_a, nlat2_a, nt_a;
    int nlat_b, nlat2_b, nt_b;
    int ierror = 0;
    npy_intp dims_dummy[3] = {1, 1, 1};
    npy_intp dims_out[3];

    if (!PyArg_ParseTuple(args, "iOOOi", &nlon, &a_obj, &b_obj, &wvh_obj, &lwork)) return NULL;
    if (as_real32_3d_fortran(a_obj, &a_arr, &nlat_a, &nlat2_a, &nt_a) != 0) return NULL;
    if (as_real32_3d_fortran(b_obj, &b_arr, &nlat_b, &nlat2_b, &nt_b) != 0) goto fail;
    if (nlat_a != nlat_b || nlat2_a != nlat2_b || nt_a != nt_b) {
        PyErr_SetString(PyExc_ValueError, "component arrays must have the same 3D shape");
        goto fail;
    }
    wvh_arr = as_real32_1d_aligned(wvh_obj);
    if (wvh_arr == NULL) goto fail;
    dims_out[0] = nlat_a;
    dims_out[1] = nlon;
    dims_out[2] = nt_a;
    if (is_div) {
        br_arr = a_arr;
        bi_arr = b_arr;
        Py_INCREF(br_arr);
        Py_INCREF(bi_arr);
        cr_arr = (PyArrayObject *)PyArray_ZEROS(3, dims_dummy, NPY_FLOAT32, 1);
        ci_arr = (PyArrayObject *)PyArray_ZEROS(3, dims_dummy, NPY_FLOAT32, 1);
    } else {
        cr_arr = a_arr;
        ci_arr = b_arr;
        Py_INCREF(cr_arr);
        Py_INCREF(ci_arr);
        br_arr = (PyArrayObject *)PyArray_ZEROS(3, dims_dummy, NPY_FLOAT32, 1);
        bi_arr = (PyArrayObject *)PyArray_ZEROS(3, dims_dummy, NPY_FLOAT32, 1);
    }
    if (br_arr == NULL || bi_arr == NULL || cr_arr == NULL || ci_arr == NULL) goto fail;
    v_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_out, NPY_FLOAT32, 1);
    w_arr = (PyArrayObject *)PyArray_EMPTY(3, dims_out, NPY_FLOAT32, 1);
    if (v_arr == NULL || w_arr == NULL) goto fail;
    func(nlon, PyArray_DATA(br_arr), PyArray_DATA(bi_arr), PyArray_DATA(cr_arr), PyArray_DATA(ci_arr),
         nlat_a, nt_a, ityp, PyArray_DATA(wvh_arr), (int)PyArray_DIM(wvh_arr, 0), lwork,
         PyArray_DATA(v_arr), PyArray_DATA(w_arr), &ierror);
    Py_DECREF(a_arr);
    Py_DECREF(b_arr);
    Py_DECREF(wvh_arr);
    Py_DECREF(br_arr);
    Py_DECREF(bi_arr);
    Py_DECREF(cr_arr);
    Py_DECREF(ci_arr);
    return Py_BuildValue("NNi", v_arr, w_arr, ierror);

fail:
    Py_XDECREF(a_arr);
    Py_XDECREF(b_arr);
    Py_XDECREF(wvh_arr);
    Py_XDECREF(br_arr);
    Py_XDECREF(bi_arr);
    Py_XDECREF(cr_arr);
    Py_XDECREF(ci_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(w_arr);
    return NULL;
}

static PyObject *py_vhsesdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 1, vhses_c, 1);
}

static PyObject *py_vhsesvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 2, vhses_c, 0);
}

static PyObject *py_vhsgsdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 1, vhsgs_c, 1);
}

static PyObject *py_vhsgsvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 2, vhsgs_c, 0);
}

static PyObject *py_vhsecdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 1, vhsec_c, 1);
}

static PyObject *py_vhsecvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 2, vhsec_c, 0);
}

static PyObject *py_vhsgcdiv(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 1, vhsgc_c, 1);
}

static PyObject *py_vhsgcvrt(PyObject *self, PyObject *args) {
    (void)self;
    return py_vhs_component_like(args, 2, vhsgc_c, 0);
}

static PyMethodDef module_methods[] = {
    {"gaqd", py_gaqd, METH_VARARGS, "Compute Gaussian latitudes and weights."},
    {"getlegfunc", py_getlegfunc, METH_VARARGS, "Compute associated Legendre functions."},
    {"lap", py_lap, METH_VARARGS, "Apply Laplacian in spectral space."},
    {"invlap", py_invlap, METH_VARARGS, "Apply inverse Laplacian in spectral space."},
    {"multsmoothfact", py_multsmoothfact, METH_VARARGS, "Apply spectral smoothing factors."},
    {"specintrp", py_specintrp, METH_VARARGS, "Interpolate one spectral field to one longitude."},
    {"shaesi", py_shaesi, METH_VARARGS, "Initialize stored regular scalar analysis workspace."},
    {"shsesi", py_shsesi, METH_VARARGS, "Initialize stored regular scalar synthesis workspace."},
    {"shagsi", py_shagsi, METH_VARARGS, "Initialize stored gaussian scalar analysis workspace."},
    {"shsgsi", py_shsgsi, METH_VARARGS, "Initialize stored gaussian scalar synthesis workspace."},
    {"vhaesi", py_vhaesi, METH_VARARGS, "Initialize stored regular vector analysis workspace."},
    {"vhagsi", py_vhagsi, METH_VARARGS, "Initialize stored gaussian vector analysis workspace."},
    {"vhsesi", py_vhsesi, METH_VARARGS, "Initialize stored regular vector synthesis workspace."},
    {"vhsgsi", py_vhsgsi, METH_VARARGS, "Initialize stored gaussian vector synthesis workspace."},
    {"shaeci", py_shaeci, METH_VARARGS, "Initialize computed regular scalar analysis workspace."},
    {"shseci", py_shseci, METH_VARARGS, "Initialize computed regular scalar synthesis workspace."},
    {"shagci", py_shagci, METH_VARARGS, "Initialize computed gaussian scalar analysis workspace."},
    {"shsgci", py_shsgci, METH_VARARGS, "Initialize computed gaussian scalar synthesis workspace."},
    {"vhaeci", py_vhaeci, METH_VARARGS, "Initialize computed regular vector analysis workspace."},
    {"vhagci", py_vhagci, METH_VARARGS, "Initialize computed gaussian vector analysis workspace."},
    {"vhseci", py_vhseci, METH_VARARGS, "Initialize computed regular vector synthesis workspace."},
    {"vhsgci", py_vhsgci, METH_VARARGS, "Initialize computed gaussian vector synthesis workspace."},
    {"shaes", py_shaes, METH_VARARGS, "Regular stored scalar harmonic analysis."},
    {"shags", py_shags, METH_VARARGS, "Gaussian stored scalar harmonic analysis."},
    {"shaec", py_shaec, METH_VARARGS, "Regular computed scalar harmonic analysis."},
    {"shagc", py_shagc, METH_VARARGS, "Gaussian computed scalar harmonic analysis."},
    {"vhaes", py_vhaes, METH_VARARGS, "Regular stored vector harmonic analysis."},
    {"vhags", py_vhags, METH_VARARGS, "Gaussian stored vector harmonic analysis."},
    {"vhaec", py_vhaec, METH_VARARGS, "Regular computed vector harmonic analysis."},
    {"vhagc", py_vhagc, METH_VARARGS, "Gaussian computed vector harmonic analysis."},
    {"vhaesdiv", py_vhaesdiv, METH_VARARGS, "Regular stored vector divergence-only analysis."},
    {"vhaesvrt", py_vhaesvrt, METH_VARARGS, "Regular stored vector vorticity-only analysis."},
    {"vhagsdiv", py_vhagsdiv, METH_VARARGS, "Gaussian stored vector divergence-only analysis."},
    {"vhagsvrt", py_vhagsvrt, METH_VARARGS, "Gaussian stored vector vorticity-only analysis."},
    {"vhaecdiv", py_vhaecdiv, METH_VARARGS, "Regular computed vector divergence-only analysis."},
    {"vhaecvrt", py_vhaecvrt, METH_VARARGS, "Regular computed vector vorticity-only analysis."},
    {"vhagcdiv", py_vhagcdiv, METH_VARARGS, "Gaussian computed vector divergence-only analysis."},
    {"vhagcvrt", py_vhagcvrt, METH_VARARGS, "Gaussian computed vector vorticity-only analysis."},
    {"shses", py_shses, METH_VARARGS, "Regular stored scalar harmonic synthesis."},
    {"shsgs", py_shsgs, METH_VARARGS, "Gaussian stored scalar harmonic synthesis."},
    {"shsec", py_shsec, METH_VARARGS, "Regular computed scalar harmonic synthesis."},
    {"shsgc", py_shsgc, METH_VARARGS, "Gaussian computed scalar harmonic synthesis."},
    {"vhses", py_vhses, METH_VARARGS, "Regular stored vector harmonic synthesis."},
    {"vhsgs", py_vhsgs, METH_VARARGS, "Gaussian stored vector harmonic synthesis."},
    {"vhsec", py_vhsec, METH_VARARGS, "Regular computed vector harmonic synthesis."},
    {"vhsgc", py_vhsgc, METH_VARARGS, "Gaussian computed vector harmonic synthesis."},
    {"vhsesdiv", py_vhsesdiv, METH_VARARGS, "Regular stored vector divergence-only synthesis."},
    {"vhsesvrt", py_vhsesvrt, METH_VARARGS, "Regular stored vector vorticity-only synthesis."},
    {"vhsgsdiv", py_vhsgsdiv, METH_VARARGS, "Gaussian stored vector divergence-only synthesis."},
    {"vhsgsvrt", py_vhsgsvrt, METH_VARARGS, "Gaussian stored vector vorticity-only synthesis."},
    {"vhsecdiv", py_vhsecdiv, METH_VARARGS, "Regular computed vector divergence-only synthesis."},
    {"vhsecvrt", py_vhsecvrt, METH_VARARGS, "Regular computed vector vorticity-only synthesis."},
    {"vhsgcdiv", py_vhsgcdiv, METH_VARARGS, "Gaussian computed vector divergence-only synthesis."},
    {"vhsgcvrt", py_vhsgcvrt, METH_VARARGS, "Gaussian computed vector vorticity-only synthesis."},
    {"onedtotwod", py_onedtotwod, METH_VARARGS, "Convert 1D spectral coefficients to 2D harmonic arrays."},
    {"onedtotwod_vrtdiv", py_onedtotwod_vrtdiv, METH_VARARGS, "Convert vorticity/divergence spectra to vector harmonic arrays."},
    {"onedtotwod_vrt", py_onedtotwod_vrt, METH_VARARGS, "Convert vorticity spectra to vector harmonic arrays."},
    {"onedtotwod_div", py_onedtotwod_div, METH_VARARGS, "Convert divergence spectra to vector harmonic arrays."},
    {"twodtooned", py_twodtooned, METH_VARARGS, "Convert 2D harmonic arrays to 1D spectral coefficients."},
    {"twodtooned_vrtdiv", py_twodtooned_vrtdiv, METH_VARARGS, "Convert vector harmonic arrays to vorticity/divergence spectra."},
    {"twodtooned_vrt", py_twodtooned_vrt, METH_VARARGS, "Convert vector harmonic arrays to vorticity spectra."},
    {"twodtooned_div", py_twodtooned_div, METH_VARARGS, "Convert vector harmonic arrays to divergence spectra."},
    {"ihgeod", py_ihgeod, METH_VARARGS, "Compute icosahedral geodesic points."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_spherepack_direct",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit__spherepack_direct(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
