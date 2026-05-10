#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>

void dvinth2p_nodes_pa_c(
    void *dati,
    void *dato,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nlevi,
    int ncol,
    int nlevo
);
void dvinth2p_nodes_pa_into_c(
    void *dati,
    void *dato,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nlevi,
    int ncol,
    int nlevo
);
void dvinth2p_ecmwf_nodes_pa_c(
    void *dati,
    void *dato,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nlevi,
    int ncol,
    int nlevo,
    int varflg,
    void *tbot,
    void *phis
);
void dvinth2p_ecmwf_nodes_pa_into_c(
    void *dati,
    void *dato,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nlevi,
    int ncol,
    int nlevo,
    int varflg,
    void *tbot,
    void *phis
);
void ddelta_pressure_hybrid_pa_c(
    void *psfc,
    void *dph,
    void *hbcofa,
    void *hbcofb,
    double p0,
    int ncol,
    int nlev,
    int nlevo
);
void ddelta_pressure_hybrid_pa_into_c(
    void *psfc,
    void *dph,
    void *hbcofa,
    void *hbcofb,
    double p0,
    int ncol,
    int nlev,
    int nlevo
);
void dpressure_at_hybrid_levels_pa_c(
    void *psfc,
    void *pressure,
    void *hbcofa,
    void *hbcofb,
    double p0,
    int ncol,
    int nlev
);
void dpressure_at_hybrid_levels_pa_into_c(
    void *psfc,
    void *pressure,
    void *hbcofa,
    void *hbcofb,
    double p0,
    int ncol,
    int nlev
);
void dgeopotential_height_hybrid_corder_pa_into_c(
    void *temp_flat,
    void *q_flat,
    void *z3_flat,
    void *psfc,
    void *phis,
    void *hyai,
    void *hybi,
    double p0,
    int nouter,
    int nlev,
    int ninner
);
void dsigma2hybrid_nodes_c(
    void *dati,
    void *dato,
    void *sigmai,
    void *sigmao,
    int intyp,
    double spvl,
    int nlevi,
    int ncol,
    int nlevo
);
void dsigma2hybrid_nodes_into_c(
    void *dati,
    void *dato,
    void *sigmai,
    void *sigmao,
    int intyp,
    double spvl,
    int nlevi,
    int ncol,
    int nlevo
);
void dsigma2hybrid_nodes_corder_into_c(
    void *dati_flat,
    void *dato_flat,
    void *sigmai,
    void *hyam,
    void *hybm,
    double p0,
    void *psfc,
    int intyp,
    double spvl,
    int nouter,
    int nlevi,
    int ninner,
    int nlevo
);
void dvinth2p_nodes_corder_pa_into_c(
    void *dati_flat,
    void *dato_flat,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nouter,
    int nlevi,
    int ninner,
    int nlevo
);
void dvinth2p_ecmwf_nodes_corder_pa_into_c(
    void *dati_flat,
    void *dato_flat,
    void *hbcofa,
    void *hbcofb,
    double p0,
    void *plevo,
    int intyp,
    void *psfc,
    double spvl,
    int kxtrp,
    int nouter,
    int nlevi,
    int ninner,
    int nlevo,
    int varflg,
    void *tbot,
    void *phis
);

static PyArrayObject *to_double_1d(PyObject *obj, int writable) {
    int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS;
    if (writable) {
        flags |= NPY_ARRAY_WRITEABLE;
    }
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, flags);
}

static PyArrayObject *to_double_2d_fortran(PyObject *obj, int writable) {
    int flags = NPY_ARRAY_ALIGNED | NPY_ARRAY_F_CONTIGUOUS;
    if (writable) {
        flags |= NPY_ARRAY_WRITEABLE;
    }
    return (PyArrayObject *) PyArray_FROM_OTF(obj, NPY_FLOAT64, flags);
}

static int require_2d_fortran(PyArrayObject *arr, const char *name) {
    if (arr == NULL || PyArray_NDIM(arr) != 2) {
        PyErr_Format(PyExc_ValueError, "%s must be a 2D array", name);
        return 0;
    }
    return 1;
}

static int require_1d(PyArrayObject *arr, const char *name) {
    if (arr == NULL || PyArray_NDIM(arr) != 1) {
        PyErr_Format(PyExc_ValueError, "%s must be a 1D array", name);
        return 0;
    }
    return 1;
}

static PyObject *py_dvinth2p_nodes_pa(PyObject *self, PyObject *args) {
    PyObject *dati_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    npy_intp dato_dims[2];
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int nlevi = 0;
    int ncol = 0;
    int nlevo = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOdOiOdi", &dati_obj, &hbcofa_obj, &hbcofb_obj, &p0, &plevo_obj, &intyp, &psfc_obj, &spvl, &kxtrp)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_1d(hbcofa_arr, "hbcofa") ||
        !require_1d(hbcofb_arr, "hbcofb") || !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc")) {
        goto fail;
    }

    nlevi = (int) PyArray_DIM(dati_arr, 0);
    ncol = (int) PyArray_DIM(dati_arr, 1);
    nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(hbcofa_arr, 0) != nlevi || (int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(psfc_arr, 0) != ncol) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_nodes_pa input shapes are inconsistent");
        goto fail;
    }

    dato_dims[0] = nlevo;
    dato_dims[1] = ncol;
    dato_arr = (PyArrayObject *) PyArray_EMPTY(2, dato_dims, NPY_FLOAT64, 1);
    if (dato_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_nodes_pa_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                      p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nlevi, ncol, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    return PyArray_Return(dato_arr);

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    Py_XDECREF(dato_arr);
    return NULL;
}

static PyObject *py_dvinth2p_nodes_pa_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"dati", "dato", "hbcofa", "hbcofb", "p0", "plevo", "intyp", "psfc", "spvl", "kxtrp", "nlevi", "ncol", "nlevo", NULL};
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int nlevi = 0;
    int ncol = 0;
    int nlevo = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOdOiOdi|iii", kwlist, &dati_obj, &dato_obj, &hbcofa_obj, &hbcofb_obj, &p0,
                                     &plevo_obj, &intyp, &psfc_obj, &spvl, &kxtrp, &nlevi, &ncol, &nlevo)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    dato_arr = to_double_2d_fortran(dato_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_2d_fortran(dato_arr, "dato") || !require_1d(hbcofa_arr, "hbcofa") ||
        !require_1d(hbcofb_arr, "hbcofb") || !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc")) {
        goto fail;
    }

    if (nlevi == 0) nlevi = (int) PyArray_DIM(dati_arr, 0);
    if (ncol == 0) ncol = (int) PyArray_DIM(dati_arr, 1);
    if (nlevo == 0) nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(dati_arr, 0) != nlevi || (int) PyArray_DIM(dati_arr, 1) != ncol ||
        (int) PyArray_DIM(dato_arr, 0) != nlevo || (int) PyArray_DIM(dato_arr, 1) != ncol ||
        (int) PyArray_DIM(hbcofa_arr, 0) != nlevi || (int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(plevo_arr, 0) != nlevo || (int) PyArray_DIM(psfc_arr, 0) != ncol) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_nodes_pa_into input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_nodes_pa_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                           p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nlevi, ncol, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    return NULL;
}

static PyObject *py_dvinth2p_ecmwf_nodes_pa(PyObject *self, PyObject *args) {
    PyObject *dati_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyObject *tbot_obj = NULL;
    PyObject *phis_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *tbot_arr = NULL;
    PyArrayObject *phis_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    npy_intp dato_dims[2];
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int varflg = 0;
    int nlevi = 0;
    int ncol = 0;
    int nlevo = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOdOiOdiiOO", &dati_obj, &hbcofa_obj, &hbcofb_obj, &p0, &plevo_obj, &intyp,
                          &psfc_obj, &spvl, &kxtrp, &varflg, &tbot_obj, &phis_obj)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    tbot_arr = to_double_1d(tbot_obj, 0);
    phis_arr = to_double_1d(phis_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_1d(hbcofa_arr, "hbcofa") || !require_1d(hbcofb_arr, "hbcofb") ||
        !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc") || !require_1d(tbot_arr, "tbot") ||
        !require_1d(phis_arr, "phis")) {
        goto fail;
    }

    nlevi = (int) PyArray_DIM(dati_arr, 0);
    ncol = (int) PyArray_DIM(dati_arr, 1);
    nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(hbcofa_arr, 0) != nlevi || (int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(psfc_arr, 0) != ncol || (int) PyArray_DIM(tbot_arr, 0) != ncol ||
        (int) PyArray_DIM(phis_arr, 0) != ncol) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_ecmwf_nodes_pa input shapes are inconsistent");
        goto fail;
    }

    dato_dims[0] = nlevo;
    dato_dims[1] = ncol;
    dato_arr = (PyArrayObject *) PyArray_EMPTY(2, dato_dims, NPY_FLOAT64, 1);
    if (dato_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_ecmwf_nodes_pa_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                            p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nlevi, ncol, nlevo,
                            varflg, PyArray_DATA(tbot_arr), PyArray_DATA(phis_arr));
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    Py_DECREF(tbot_arr);
    Py_DECREF(phis_arr);
    return PyArray_Return(dato_arr);

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    Py_XDECREF(tbot_arr);
    Py_XDECREF(phis_arr);
    Py_XDECREF(dato_arr);
    return NULL;
}

static PyObject *py_dvinth2p_ecmwf_nodes_pa_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"dati", "dato", "hbcofa", "hbcofb", "p0", "plevo", "intyp", "psfc", "spvl", "kxtrp", "varflg", "tbot", "phis", "nlevi", "ncol", "nlevo", NULL};
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyObject *tbot_obj = NULL;
    PyObject *phis_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *tbot_arr = NULL;
    PyArrayObject *phis_arr = NULL;
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int varflg = 0;
    int nlevi = 0;
    int ncol = 0;
    int nlevo = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOdOiOdiiOO|iii", kwlist, &dati_obj, &dato_obj, &hbcofa_obj, &hbcofb_obj, &p0,
                                     &plevo_obj, &intyp, &psfc_obj, &spvl, &kxtrp, &varflg, &tbot_obj, &phis_obj,
                                     &nlevi, &ncol, &nlevo)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    dato_arr = to_double_2d_fortran(dato_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    tbot_arr = to_double_1d(tbot_obj, 0);
    phis_arr = to_double_1d(phis_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_2d_fortran(dato_arr, "dato") || !require_1d(hbcofa_arr, "hbcofa") ||
        !require_1d(hbcofb_arr, "hbcofb") || !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc") ||
        !require_1d(tbot_arr, "tbot") || !require_1d(phis_arr, "phis")) {
        goto fail;
    }

    if (nlevi == 0) nlevi = (int) PyArray_DIM(dati_arr, 0);
    if (ncol == 0) ncol = (int) PyArray_DIM(dati_arr, 1);
    if (nlevo == 0) nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(dati_arr, 0) != nlevi || (int) PyArray_DIM(dati_arr, 1) != ncol ||
        (int) PyArray_DIM(dato_arr, 0) != nlevo || (int) PyArray_DIM(dato_arr, 1) != ncol ||
        (int) PyArray_DIM(hbcofa_arr, 0) != nlevi || (int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(plevo_arr, 0) != nlevo || (int) PyArray_DIM(psfc_arr, 0) != ncol ||
        (int) PyArray_DIM(tbot_arr, 0) != ncol || (int) PyArray_DIM(phis_arr, 0) != ncol) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_ecmwf_nodes_pa_into input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_ecmwf_nodes_pa_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                                 p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nlevi, ncol, nlevo,
                                 varflg, PyArray_DATA(tbot_arr), PyArray_DATA(phis_arr));
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    Py_DECREF(tbot_arr);
    Py_DECREF(phis_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    Py_XDECREF(tbot_arr);
    Py_XDECREF(phis_arr);
    return NULL;
}

static PyObject *py_ddelta_pressure_hybrid_pa(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"psfc", "hbcofa", "hbcofb", "p0", "nlevo", "ncol", "nlev", NULL};
    PyObject *psfc_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *dph_arr = NULL;
    npy_intp dph_dims[2];
    double p0 = 0.0;
    int nlevo = 0;
    int ncol = 0;
    int nlev = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOd|iii", kwlist, &psfc_obj, &hbcofa_obj, &hbcofb_obj, &p0, &nlevo, &ncol, &nlev)) {
        return NULL;
    }

    psfc_arr = to_double_1d(psfc_obj, 0);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    if (!require_1d(psfc_arr, "psfc") || !require_1d(hbcofa_arr, "hbcofa") || !require_1d(hbcofb_arr, "hbcofb")) {
        goto fail;
    }

    if (ncol == 0) ncol = (int) PyArray_DIM(psfc_arr, 0);
    if (nlev == 0) nlev = (int) PyArray_DIM(hbcofa_arr, 0);
    if (nlevo == 0) nlevo = nlev - 1;
    if ((int) PyArray_DIM(psfc_arr, 0) != ncol || (int) PyArray_DIM(hbcofa_arr, 0) != nlev ||
        (int) PyArray_DIM(hbcofb_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "delta_pressure_hybrid input shapes are inconsistent");
        goto fail;
    }

    dph_dims[0] = nlevo;
    dph_dims[1] = ncol;
    dph_arr = (PyArrayObject *) PyArray_EMPTY(2, dph_dims, NPY_FLOAT64, 1);
    if (dph_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    ddelta_pressure_hybrid_pa_c(PyArray_DATA(psfc_arr), PyArray_DATA(dph_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr), p0, ncol, nlev, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(psfc_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    return PyArray_Return(dph_arr);

fail:
    Py_XDECREF(psfc_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(dph_arr);
    return NULL;
}

static PyObject *py_ddelta_pressure_hybrid_pa_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"psfc", "dph", "hbcofa", "hbcofb", "p0", "nlevo", "ncol", "nlev", NULL};
    PyObject *psfc_obj = NULL;
    PyObject *dph_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *dph_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    double p0 = 0.0;
    int nlevo = 0;
    int ncol = 0;
    int nlev = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOd|iii", kwlist, &psfc_obj, &dph_obj, &hbcofa_obj, &hbcofb_obj, &p0, &nlevo, &ncol, &nlev)) {
        return NULL;
    }

    psfc_arr = to_double_1d(psfc_obj, 0);
    dph_arr = to_double_2d_fortran(dph_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    if (!require_1d(psfc_arr, "psfc") || !require_2d_fortran(dph_arr, "dph") || !require_1d(hbcofa_arr, "hbcofa") || !require_1d(hbcofb_arr, "hbcofb")) {
        goto fail;
    }

    if (ncol == 0) ncol = (int) PyArray_DIM(psfc_arr, 0);
    if (nlev == 0) nlev = (int) PyArray_DIM(hbcofa_arr, 0);
    if (nlevo == 0) nlevo = (int) PyArray_DIM(dph_arr, 0);
    if ((int) PyArray_DIM(psfc_arr, 0) != ncol || (int) PyArray_DIM(dph_arr, 0) != nlevo ||
        (int) PyArray_DIM(dph_arr, 1) != ncol || (int) PyArray_DIM(hbcofa_arr, 0) != nlev ||
        (int) PyArray_DIM(hbcofb_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "delta_pressure_hybrid_into input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    ddelta_pressure_hybrid_pa_into_c(PyArray_DATA(psfc_arr), PyArray_DATA(dph_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr), p0, ncol, nlev, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(psfc_arr);
    Py_DECREF(dph_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(psfc_arr);
    Py_XDECREF(dph_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    return NULL;
}

static PyObject *py_dpressure_at_hybrid_levels_pa(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"psfc", "hbcofa", "hbcofb", "p0", "ncol", "nlev", NULL};
    PyObject *psfc_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    npy_intp pressure_dims[2];
    double p0 = 0.0;
    int ncol = 0;
    int nlev = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOd|ii", kwlist, &psfc_obj, &hbcofa_obj, &hbcofb_obj, &p0, &ncol, &nlev)) {
        return NULL;
    }

    psfc_arr = to_double_1d(psfc_obj, 0);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    if (!require_1d(psfc_arr, "psfc") || !require_1d(hbcofa_arr, "hbcofa") || !require_1d(hbcofb_arr, "hbcofb")) {
        goto fail;
    }

    if (ncol == 0) ncol = (int) PyArray_DIM(psfc_arr, 0);
    if (nlev == 0) nlev = (int) PyArray_DIM(hbcofa_arr, 0);
    if ((int) PyArray_DIM(psfc_arr, 0) != ncol || (int) PyArray_DIM(hbcofa_arr, 0) != nlev || (int) PyArray_DIM(hbcofb_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pressure_at_hybrid_levels input shapes are inconsistent");
        goto fail;
    }

    pressure_dims[0] = nlev;
    pressure_dims[1] = ncol;
    pressure_arr = (PyArrayObject *) PyArray_EMPTY(2, pressure_dims, NPY_FLOAT64, 1);
    if (pressure_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dpressure_at_hybrid_levels_pa_c(PyArray_DATA(psfc_arr), PyArray_DATA(pressure_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr), p0, ncol, nlev);
    Py_END_ALLOW_THREADS

    Py_DECREF(psfc_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    return PyArray_Return(pressure_arr);

fail:
    Py_XDECREF(psfc_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(pressure_arr);
    return NULL;
}

static PyObject *py_dpressure_at_hybrid_levels_pa_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"psfc", "pressure", "hbcofa", "hbcofb", "p0", "ncol", "nlev", NULL};
    PyObject *psfc_obj = NULL;
    PyObject *pressure_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *pressure_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    double p0 = 0.0;
    int ncol = 0;
    int nlev = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOd|ii", kwlist, &psfc_obj, &pressure_obj, &hbcofa_obj, &hbcofb_obj, &p0, &ncol, &nlev)) {
        return NULL;
    }

    psfc_arr = to_double_1d(psfc_obj, 0);
    pressure_arr = to_double_2d_fortran(pressure_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    if (!require_1d(psfc_arr, "psfc") || !require_2d_fortran(pressure_arr, "pressure") || !require_1d(hbcofa_arr, "hbcofa") || !require_1d(hbcofb_arr, "hbcofb")) {
        goto fail;
    }

    if (ncol == 0) ncol = (int) PyArray_DIM(psfc_arr, 0);
    if (nlev == 0) nlev = (int) PyArray_DIM(hbcofa_arr, 0);
    if ((int) PyArray_DIM(psfc_arr, 0) != ncol || (int) PyArray_DIM(pressure_arr, 0) != nlev ||
        (int) PyArray_DIM(pressure_arr, 1) != ncol || (int) PyArray_DIM(hbcofa_arr, 0) != nlev ||
        (int) PyArray_DIM(hbcofb_arr, 0) != nlev) {
        PyErr_SetString(PyExc_ValueError, "pressure_at_hybrid_levels_into input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dpressure_at_hybrid_levels_pa_into_c(PyArray_DATA(psfc_arr), PyArray_DATA(pressure_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr), p0, ncol, nlev);
    Py_END_ALLOW_THREADS

    Py_DECREF(psfc_arr);
    Py_DECREF(pressure_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(psfc_arr);
    Py_XDECREF(pressure_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    return NULL;
}

static PyObject *py_dgeopotential_height_hybrid_corder_pa_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"temp_flat", "q_flat", "z3_flat", "psfc", "phis", "hyai", "hybi", "p0", "nouter", "nlev", "ninner", NULL};
    PyObject *temp_obj = NULL;
    PyObject *q_obj = NULL;
    PyObject *z3_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyObject *phis_obj = NULL;
    PyObject *hyai_obj = NULL;
    PyObject *hybi_obj = NULL;
    PyArrayObject *temp_arr = NULL;
    PyArrayObject *q_arr = NULL;
    PyArrayObject *z3_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *phis_arr = NULL;
    PyArrayObject *hyai_arr = NULL;
    PyArrayObject *hybi_arr = NULL;
    double p0 = 0.0;
    int nouter = 0;
    int nlev = 0;
    int ninner = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOOOOd|iii",
            kwlist,
            &temp_obj,
            &q_obj,
            &z3_obj,
            &psfc_obj,
            &phis_obj,
            &hyai_obj,
            &hybi_obj,
            &p0,
            &nouter,
            &nlev,
            &ninner
        )) {
        return NULL;
    }

    temp_arr = to_double_1d(temp_obj, 0);
    q_arr = to_double_1d(q_obj, 0);
    z3_arr = to_double_1d(z3_obj, 1);
    psfc_arr = to_double_1d(psfc_obj, 0);
    phis_arr = to_double_1d(phis_obj, 0);
    hyai_arr = to_double_1d(hyai_obj, 0);
    hybi_arr = to_double_1d(hybi_obj, 0);
    if (!require_1d(temp_arr, "temp_flat") || !require_1d(q_arr, "q_flat") || !require_1d(z3_arr, "z3_flat") ||
        !require_1d(psfc_arr, "psfc") || !require_1d(phis_arr, "phis") || !require_1d(hyai_arr, "hyai") || !require_1d(hybi_arr, "hybi")) {
        goto fail;
    }

    if ((int) PyArray_DIM(temp_arr, 0) != nouter * nlev * ninner || (int) PyArray_DIM(q_arr, 0) != nouter * nlev * ninner ||
        (int) PyArray_DIM(z3_arr, 0) != nouter * nlev * ninner || (int) PyArray_DIM(psfc_arr, 0) != nouter * ninner ||
        (int) PyArray_DIM(phis_arr, 0) != nouter * ninner || (int) PyArray_DIM(hyai_arr, 0) != nlev + 1 ||
        (int) PyArray_DIM(hybi_arr, 0) != nlev + 1) {
        PyErr_SetString(PyExc_ValueError, "geopotential_height_hybrid_corder input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dgeopotential_height_hybrid_corder_pa_into_c(PyArray_DATA(temp_arr), PyArray_DATA(q_arr), PyArray_DATA(z3_arr), PyArray_DATA(psfc_arr),
                                               PyArray_DATA(phis_arr), PyArray_DATA(hyai_arr), PyArray_DATA(hybi_arr), p0, nouter, nlev, ninner);
    Py_END_ALLOW_THREADS

    Py_DECREF(temp_arr);
    Py_DECREF(q_arr);
    Py_DECREF(z3_arr);
    Py_DECREF(psfc_arr);
    Py_DECREF(phis_arr);
    Py_DECREF(hyai_arr);
    Py_DECREF(hybi_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(temp_arr);
    Py_XDECREF(q_arr);
    Py_XDECREF(z3_arr);
    Py_XDECREF(psfc_arr);
    Py_XDECREF(phis_arr);
    Py_XDECREF(hyai_arr);
    Py_XDECREF(hybi_arr);
    return NULL;
}

static PyObject *py_dsigma2hybrid_nodes(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"dati", "sigmai", "sigmao", "intyp", "spvl", "nlevo", "nlevi", "ncol", NULL};
    PyObject *dati_obj = NULL;
    PyObject *sigmai_obj = NULL;
    PyObject *sigmao_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *sigmai_arr = NULL;
    PyArrayObject *sigmao_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    double spvl = 0.0;
    int intyp = 0;
    int nlevo = 0;
    int nlevi = 0;
    int ncol = 0;
    npy_intp dato_dims[2];

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOidi|ii", kwlist, &dati_obj, &sigmai_obj, &sigmao_obj, &intyp, &spvl, &nlevo, &nlevi, &ncol)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    sigmai_arr = to_double_1d(sigmai_obj, 0);
    sigmao_arr = to_double_2d_fortran(sigmao_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_1d(sigmai_arr, "sigmai") || !require_2d_fortran(sigmao_arr, "sigmao")) {
        goto fail;
    }

    if (nlevi == 0) nlevi = (int) PyArray_DIM(dati_arr, 0);
    if (ncol == 0) ncol = (int) PyArray_DIM(dati_arr, 1);
    if (nlevo == 0) nlevo = (int) PyArray_DIM(sigmao_arr, 0);
    if ((int) PyArray_DIM(dati_arr, 0) != nlevi || (int) PyArray_DIM(dati_arr, 1) != ncol ||
        (int) PyArray_DIM(sigmai_arr, 0) != nlevi || (int) PyArray_DIM(sigmao_arr, 0) != nlevo ||
        (int) PyArray_DIM(sigmao_arr, 1) != ncol) {
        PyErr_SetString(PyExc_ValueError, "sigma2hybrid input shapes are inconsistent");
        goto fail;
    }

    dato_dims[0] = nlevo;
    dato_dims[1] = ncol;
    dato_arr = (PyArrayObject *) PyArray_EMPTY(2, dato_dims, NPY_FLOAT64, 1);
    if (dato_arr == NULL) {
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dsigma2hybrid_nodes_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(sigmai_arr), PyArray_DATA(sigmao_arr),
                        intyp, spvl, nlevi, ncol, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(sigmai_arr);
    Py_DECREF(sigmao_arr);
    return PyArray_Return(dato_arr);

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(sigmai_arr);
    Py_XDECREF(sigmao_arr);
    Py_XDECREF(dato_arr);
    return NULL;
}

static PyObject *py_dsigma2hybrid_nodes_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"dati", "dato", "sigmai", "sigmao", "intyp", "spvl", "nlevo", "nlevi", "ncol", NULL};
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *sigmai_obj = NULL;
    PyObject *sigmao_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *sigmai_arr = NULL;
    PyArrayObject *sigmao_arr = NULL;
    double spvl = 0.0;
    int intyp = 0;
    int nlevo = 0;
    int nlevi = 0;
    int ncol = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOidi|ii", kwlist, &dati_obj, &dato_obj, &sigmai_obj, &sigmao_obj, &intyp, &spvl, &nlevo, &nlevi, &ncol)) {
        return NULL;
    }

    dati_arr = to_double_2d_fortran(dati_obj, 0);
    dato_arr = to_double_2d_fortran(dato_obj, 1);
    sigmai_arr = to_double_1d(sigmai_obj, 0);
    sigmao_arr = to_double_2d_fortran(sigmao_obj, 0);
    if (!require_2d_fortran(dati_arr, "dati") || !require_2d_fortran(dato_arr, "dato") || !require_1d(sigmai_arr, "sigmai") ||
        !require_2d_fortran(sigmao_arr, "sigmao")) {
        goto fail;
    }

    if (nlevi == 0) nlevi = (int) PyArray_DIM(dati_arr, 0);
    if (ncol == 0) ncol = (int) PyArray_DIM(dati_arr, 1);
    if (nlevo == 0) nlevo = (int) PyArray_DIM(dato_arr, 0);
    if ((int) PyArray_DIM(dati_arr, 0) != nlevi || (int) PyArray_DIM(dati_arr, 1) != ncol ||
        (int) PyArray_DIM(dato_arr, 0) != nlevo || (int) PyArray_DIM(dato_arr, 1) != ncol ||
        (int) PyArray_DIM(sigmai_arr, 0) != nlevi || (int) PyArray_DIM(sigmao_arr, 0) != nlevo ||
        (int) PyArray_DIM(sigmao_arr, 1) != ncol) {
        PyErr_SetString(PyExc_ValueError, "sigma2hybrid_into input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dsigma2hybrid_nodes_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(sigmai_arr), PyArray_DATA(sigmao_arr),
                             intyp, spvl, nlevi, ncol, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(sigmai_arr);
    Py_DECREF(sigmao_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(sigmai_arr);
    Py_XDECREF(sigmao_arr);
    return NULL;
}

static PyObject *py_dsigma2hybrid_nodes_corder_into(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"dati_flat", "dato_flat", "sigmai", "hyam", "hybm", "p0", "psfc", "intyp", "spvl", "nouter", "ninner", "nlevo", "nlevi", NULL};
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *sigmai_obj = NULL;
    PyObject *hyam_obj = NULL;
    PyObject *hybm_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *sigmai_arr = NULL;
    PyArrayObject *hyam_arr = NULL;
    PyArrayObject *hybm_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int nouter = 0;
    int ninner = 0;
    int nlevo = 0;
    int nlevi = 0;

    (void) self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOdOidiii|i", kwlist, &dati_obj, &dato_obj, &sigmai_obj, &hyam_obj, &hybm_obj,
                                     &p0, &psfc_obj, &intyp, &spvl, &nouter, &ninner, &nlevo, &nlevi)) {
        return NULL;
    }

    dati_arr = to_double_1d(dati_obj, 0);
    dato_arr = to_double_1d(dato_obj, 1);
    sigmai_arr = to_double_1d(sigmai_obj, 0);
    hyam_arr = to_double_1d(hyam_obj, 0);
    hybm_arr = to_double_1d(hybm_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    if (!require_1d(dati_arr, "dati_flat") || !require_1d(dato_arr, "dato_flat") || !require_1d(sigmai_arr, "sigmai") ||
        !require_1d(hyam_arr, "hyam") || !require_1d(hybm_arr, "hybm") || !require_1d(psfc_arr, "psfc")) {
        goto fail;
    }

    if (nlevi == 0) nlevi = (int) PyArray_DIM(sigmai_arr, 0);
    if ((int) PyArray_DIM(hyam_arr, 0) != nlevo || (int) PyArray_DIM(hybm_arr, 0) != nlevo ||
        (int) PyArray_DIM(psfc_arr, 0) != nouter * ninner || (int) PyArray_DIM(dati_arr, 0) != nouter * nlevi * ninner ||
        (int) PyArray_DIM(dato_arr, 0) != nouter * nlevo * ninner) {
        PyErr_SetString(PyExc_ValueError, "sigma2hybrid_corder input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dsigma2hybrid_nodes_corder_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(sigmai_arr), PyArray_DATA(hyam_arr),
                                    PyArray_DATA(hybm_arr), p0, PyArray_DATA(psfc_arr), intyp, spvl, nouter, nlevi, ninner, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(sigmai_arr);
    Py_DECREF(hyam_arr);
    Py_DECREF(hybm_arr);
    Py_DECREF(psfc_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(sigmai_arr);
    Py_XDECREF(hyam_arr);
    Py_XDECREF(hybm_arr);
    Py_XDECREF(psfc_arr);
    return NULL;
}

static PyObject *py_dvinth2p_nodes_corder_pa_into(PyObject *self, PyObject *args) {
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int nouter = 0;
    int ninner = 0;
    int nlevi = 0;
    int nlevo = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOdOiOdiii", &dati_obj, &dato_obj, &hbcofa_obj, &hbcofb_obj, &p0, &plevo_obj,
                          &intyp, &psfc_obj, &spvl, &kxtrp, &nouter, &ninner)) {
        return NULL;
    }

    dati_arr = to_double_1d(dati_obj, 0);
    dato_arr = to_double_1d(dato_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    if (!require_1d(dati_arr, "dati_flat") || !require_1d(dato_arr, "dato_flat") || !require_1d(hbcofa_arr, "hbcofa") ||
        !require_1d(hbcofb_arr, "hbcofb") || !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc")) {
        goto fail;
    }

    nlevi = (int) PyArray_DIM(hbcofa_arr, 0);
    nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(psfc_arr, 0) != nouter * ninner || (int) PyArray_DIM(dati_arr, 0) != nouter * nlevi * ninner ||
        (int) PyArray_DIM(dato_arr, 0) != nouter * nlevo * ninner) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_nodes_corder input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_nodes_corder_pa_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                                  p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nouter, nlevi, ninner, nlevo);
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    return NULL;
}

static PyObject *py_dvinth2p_ecmwf_nodes_corder_pa_into(PyObject *self, PyObject *args) {
    PyObject *dati_obj = NULL;
    PyObject *dato_obj = NULL;
    PyObject *hbcofa_obj = NULL;
    PyObject *hbcofb_obj = NULL;
    PyObject *plevo_obj = NULL;
    PyObject *psfc_obj = NULL;
    PyObject *tbot_obj = NULL;
    PyObject *phis_obj = NULL;
    PyArrayObject *dati_arr = NULL;
    PyArrayObject *dato_arr = NULL;
    PyArrayObject *hbcofa_arr = NULL;
    PyArrayObject *hbcofb_arr = NULL;
    PyArrayObject *plevo_arr = NULL;
    PyArrayObject *psfc_arr = NULL;
    PyArrayObject *tbot_arr = NULL;
    PyArrayObject *phis_arr = NULL;
    double p0 = 0.0;
    double spvl = 0.0;
    int intyp = 0;
    int kxtrp = 0;
    int nouter = 0;
    int ninner = 0;
    int nlevo = 0;
    int varflg = 0;
    int nlevi = 0;

    (void) self;

    if (!PyArg_ParseTuple(args, "OOOOdOiOdiiiiOO", &dati_obj, &dato_obj, &hbcofa_obj, &hbcofb_obj, &p0, &plevo_obj,
                          &intyp, &psfc_obj, &spvl, &kxtrp, &nouter, &ninner, &varflg, &tbot_obj, &phis_obj)) {
        return NULL;
    }

    dati_arr = to_double_1d(dati_obj, 0);
    dato_arr = to_double_1d(dato_obj, 1);
    hbcofa_arr = to_double_1d(hbcofa_obj, 0);
    hbcofb_arr = to_double_1d(hbcofb_obj, 0);
    plevo_arr = to_double_1d(plevo_obj, 0);
    psfc_arr = to_double_1d(psfc_obj, 0);
    tbot_arr = to_double_1d(tbot_obj, 0);
    phis_arr = to_double_1d(phis_obj, 0);
    if (!require_1d(dati_arr, "dati_flat") || !require_1d(dato_arr, "dato_flat") || !require_1d(hbcofa_arr, "hbcofa") ||
        !require_1d(hbcofb_arr, "hbcofb") || !require_1d(plevo_arr, "plevo") || !require_1d(psfc_arr, "psfc") ||
        !require_1d(tbot_arr, "tbot") || !require_1d(phis_arr, "phis")) {
        goto fail;
    }

    nlevi = (int) PyArray_DIM(hbcofa_arr, 0);
    nlevo = (int) PyArray_DIM(plevo_arr, 0);
    if ((int) PyArray_DIM(hbcofb_arr, 0) != nlevi ||
        (int) PyArray_DIM(psfc_arr, 0) != nouter * ninner || (int) PyArray_DIM(tbot_arr, 0) != nouter * ninner ||
        (int) PyArray_DIM(phis_arr, 0) != nouter * ninner || (int) PyArray_DIM(dati_arr, 0) != nouter * nlevi * ninner ||
        (int) PyArray_DIM(dato_arr, 0) != nouter * nlevo * ninner) {
        PyErr_SetString(PyExc_ValueError, "vinth2p_ecmwf_nodes_corder input shapes are inconsistent");
        goto fail;
    }

    Py_BEGIN_ALLOW_THREADS
    dvinth2p_ecmwf_nodes_corder_pa_into_c(PyArray_DATA(dati_arr), PyArray_DATA(dato_arr), PyArray_DATA(hbcofa_arr), PyArray_DATA(hbcofb_arr),
                                        p0, PyArray_DATA(plevo_arr), intyp, PyArray_DATA(psfc_arr), spvl, kxtrp, nouter, nlevi,
                                        ninner, nlevo, varflg, PyArray_DATA(tbot_arr), PyArray_DATA(phis_arr));
    Py_END_ALLOW_THREADS

    Py_DECREF(dati_arr);
    Py_DECREF(dato_arr);
    Py_DECREF(hbcofa_arr);
    Py_DECREF(hbcofb_arr);
    Py_DECREF(plevo_arr);
    Py_DECREF(psfc_arr);
    Py_DECREF(tbot_arr);
    Py_DECREF(phis_arr);
    Py_RETURN_NONE;

fail:
    Py_XDECREF(dati_arr);
    Py_XDECREF(dato_arr);
    Py_XDECREF(hbcofa_arr);
    Py_XDECREF(hbcofb_arr);
    Py_XDECREF(plevo_arr);
    Py_XDECREF(psfc_arr);
    Py_XDECREF(tbot_arr);
    Py_XDECREF(phis_arr);
    return NULL;
}

static PyMethodDef module_methods[] = {
    {"dvinth2p_nodes_pa", py_dvinth2p_nodes_pa, METH_VARARGS, "Return-allocating hybrid-to-pressure interpolation."},
    {"dvinth2p_nodes_pa_into", (PyCFunction) (void (*)(void)) py_dvinth2p_nodes_pa_into, METH_VARARGS | METH_KEYWORDS, "In-place hybrid-to-pressure interpolation."},
    {"dvinth2p_ecmwf_nodes_pa", py_dvinth2p_ecmwf_nodes_pa, METH_VARARGS, "Return-allocating ECMWF hybrid-to-pressure interpolation."},
    {"dvinth2p_ecmwf_nodes_pa_into", (PyCFunction) (void (*)(void)) py_dvinth2p_ecmwf_nodes_pa_into, METH_VARARGS | METH_KEYWORDS, "In-place ECMWF hybrid-to-pressure interpolation."},
    {"ddelta_pressure_hybrid_pa", (PyCFunction) (void (*)(void)) py_ddelta_pressure_hybrid_pa, METH_VARARGS | METH_KEYWORDS, "Return hybrid layer thicknesses."},
    {"ddelta_pressure_hybrid_pa_into", (PyCFunction) (void (*)(void)) py_ddelta_pressure_hybrid_pa_into, METH_VARARGS | METH_KEYWORDS, "Write hybrid layer thicknesses into a provided buffer."},
    {"dpressure_at_hybrid_levels_pa", (PyCFunction) (void (*)(void)) py_dpressure_at_hybrid_levels_pa, METH_VARARGS | METH_KEYWORDS, "Return hybrid-level pressures."},
    {"dpressure_at_hybrid_levels_pa_into", (PyCFunction) (void (*)(void)) py_dpressure_at_hybrid_levels_pa_into, METH_VARARGS | METH_KEYWORDS, "Write hybrid-level pressures into a provided buffer."},
    {"dgeopotential_height_hybrid_corder_pa_into", (PyCFunction) (void (*)(void)) py_dgeopotential_height_hybrid_corder_pa_into, METH_VARARGS | METH_KEYWORDS, "In-place C-order geopotential-height kernel."},
    {"dsigma2hybrid_nodes", (PyCFunction) (void (*)(void)) py_dsigma2hybrid_nodes, METH_VARARGS | METH_KEYWORDS, "Return-allocating sigma-to-hybrid interpolation."},
    {"dsigma2hybrid_nodes_into", (PyCFunction) (void (*)(void)) py_dsigma2hybrid_nodes_into, METH_VARARGS | METH_KEYWORDS, "In-place sigma-to-hybrid interpolation."},
    {"dsigma2hybrid_nodes_corder_into", (PyCFunction) (void (*)(void)) py_dsigma2hybrid_nodes_corder_into, METH_VARARGS | METH_KEYWORDS, "In-place C-order sigma-to-hybrid interpolation."},
    {"dvinth2p_nodes_corder_pa_into", py_dvinth2p_nodes_corder_pa_into, METH_VARARGS, "In-place C-order hybrid-to-pressure interpolation."},
    {"dvinth2p_ecmwf_nodes_corder_pa_into", py_dvinth2p_ecmwf_nodes_corder_pa_into, METH_VARARGS, "In-place C-order ECMWF hybrid-to-pressure interpolation."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "vinth2p_kernels",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC PyInit_vinth2p_kernels(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
