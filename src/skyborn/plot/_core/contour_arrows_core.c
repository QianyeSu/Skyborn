/*
 * contour_arrows_core.c
 *
 * C-accelerated core functions for arrow_contour.
 * Provides point_at_distance, local_straightness_score, and select_arrow_end_distances.
 */

#define PY_SSIZE_T_CLEAN
#define _USE_MATH_DEFINES  /* Enable M_PI on MSVC */
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Helper: compute 2D Euclidean distance */
static inline double hypot2(double x, double y) {
    return sqrt(x * x + y * y);
}

/* Helper: clamp value to [lower, upper] */
static inline double clip_double(double value, double lower, double upper) {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

/*
 * point_at_distance(vertices, distance) -> ndarray or None
 *
 * Find the point at given arc-length distance along the path.
 */
static PyObject* point_at_distance(PyObject* self, PyObject* args) {
    PyArrayObject *vertices_obj;
    PyArrayObject *vertices_contiguous = NULL;
    double distance;

    if (!PyArg_ParseTuple(args, "O!d", &PyArray_Type, &vertices_obj, &distance)) {
        return NULL;
    }

    if (PyArray_NDIM(vertices_obj) != 2 || PyArray_DIM(vertices_obj, 1) != 2) {
        PyErr_SetString(PyExc_ValueError, "vertices must be (N, 2) array");
        return NULL;
    }

    /* Check dtype - must be float64 */
    if (PyArray_TYPE(vertices_obj) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "vertices must be float64 (np.float64) array");
        return NULL;
    }

    /* Ensure array is contiguous for safe stride calculation */
    if (!PyArray_ISCARRAY_RO(vertices_obj) && !PyArray_ISFARRAY_RO(vertices_obj)) {
        vertices_contiguous = (PyArrayObject*)PyArray_GETCONTIGUOUS(vertices_obj);
        if (!vertices_contiguous) {
            return NULL;
        }
        vertices_obj = vertices_contiguous;
    }

    npy_intp n = PyArray_DIM(vertices_obj, 0);
    if (n < 2) {
        Py_XDECREF(vertices_contiguous);
        Py_RETURN_NONE;
    }

    double *vertices = (double*)PyArray_DATA(vertices_obj);
    npy_intp stride0 = PyArray_STRIDE(vertices_obj, 0) / sizeof(double);
    npy_intp stride1 = PyArray_STRIDE(vertices_obj, 1) / sizeof(double);

    /* Compute total length */
    double total = 0.0;
    for (npy_intp i = 0; i < n - 1; i++) {
        double x0 = vertices[i * stride0];
        double y0 = vertices[i * stride0 + stride1];
        double x1 = vertices[(i + 1) * stride0];
        double y1 = vertices[(i + 1) * stride0 + stride1];
        double dx = x1 - x0;
        double dy = y1 - y0;
        total += hypot2(dx, dy);
    }

    if (total <= 0.0) {
        Py_XDECREF(vertices_contiguous);
        Py_RETURN_NONE;
    }

    distance = clip_double(distance, 0.0, total);

    /* Find segment containing target distance */
    double cumulative = 0.0;
    for (npy_intp i = 0; i < n - 1; i++) {
        double x0 = vertices[i * stride0];
        double y0 = vertices[i * stride0 + stride1];
        double x1 = vertices[(i + 1) * stride0];
        double y1 = vertices[(i + 1) * stride0 + stride1];
        double dx = x1 - x0;
        double dy = y1 - y0;
        double segment_length = hypot2(dx, dy);

        if (cumulative + segment_length >= distance) {
            double overshoot = distance - cumulative;
            double ratio = (segment_length > 0.0) ? (overshoot / segment_length) : 0.0;

            npy_intp dims[1] = {2};
            PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
            if (!result) return NULL;

            double *result_data = (double*)PyArray_DATA(result);
            result_data[0] = x0 + ratio * dx;
            result_data[1] = y0 + ratio * dy;

            Py_XDECREF(vertices_contiguous);
            return (PyObject*)result;
        }

        cumulative += segment_length;
    }

    /* Return last vertex */
    npy_intp dims[1] = {2};
    PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!result) return NULL;

    double *result_data = (double*)PyArray_DATA(result);
    result_data[0] = vertices[(n - 1) * stride0];
    result_data[1] = vertices[(n - 1) * stride0 + stride1];

    Py_XDECREF(vertices_contiguous);
    return (PyObject*)result;
}

/*
 * local_straightness_score(vertices, distance, arrow_length) -> float
 *
 * Score the straightness at given distance for arrow placement.
 * Uses three-point method: start, middle, end along arrow span.
 */
static PyObject* local_straightness_score(PyObject* self, PyObject* args) {
    PyArrayObject *vertices_obj;
    double distance, arrow_length;

    if (!PyArg_ParseTuple(args, "O!dd", &PyArray_Type, &vertices_obj, &distance, &arrow_length)) {
        return NULL;
    }

    /* Call point_at_distance for three points */
    PyObject *start_args = Py_BuildValue("(Od)", vertices_obj, fmax(0.0, distance - arrow_length));
    PyObject *middle_args = Py_BuildValue("(Od)", vertices_obj, fmax(0.0, distance - arrow_length * 0.5));
    PyObject *end_args = Py_BuildValue("(Od)", vertices_obj, distance);

    if (!start_args || !middle_args || !end_args) {
        Py_XDECREF(start_args);
        Py_XDECREF(middle_args);
        Py_XDECREF(end_args);
        return NULL;
    }

    PyObject *start_obj = point_at_distance(self, start_args);
    PyObject *middle_obj = point_at_distance(self, middle_args);
    PyObject *end_obj = point_at_distance(self, end_args);

    Py_DECREF(start_args);
    Py_DECREF(middle_args);
    Py_DECREF(end_args);

    if (!start_obj || !middle_obj || !end_obj) {
        Py_XDECREF(start_obj);
        Py_XDECREF(middle_obj);
        Py_XDECREF(end_obj);
        return NULL;
    }

    /* Check for None returns */
    if (start_obj == Py_None || middle_obj == Py_None || end_obj == Py_None) {
        Py_DECREF(start_obj);
        Py_DECREF(middle_obj);
        Py_DECREF(end_obj);
        return PyFloat_FromDouble(-INFINITY);
    }

    PyArrayObject *start = (PyArrayObject*)start_obj;
    PyArrayObject *middle = (PyArrayObject*)middle_obj;
    PyArrayObject *end = (PyArrayObject*)end_obj;

    double *start_data = (double*)PyArray_DATA(start);
    double *middle_data = (double*)PyArray_DATA(middle);
    double *end_data = (double*)PyArray_DATA(end);

    /* Compute vectors: first = middle - start, second = end - middle */
    double first_x = middle_data[0] - start_data[0];
    double first_y = middle_data[1] - start_data[1];
    double second_x = end_data[0] - middle_data[0];
    double second_y = end_data[1] - middle_data[1];

    double first_length = hypot2(first_x, first_y);
    double second_length = hypot2(second_x, second_y);

    /* Chord from start to end */
    double chord_x = end_data[0] - start_data[0];
    double chord_y = end_data[1] - start_data[1];
    double chord_length = hypot2(chord_x, chord_y);

    Py_DECREF(start_obj);
    Py_DECREF(middle_obj);
    Py_DECREF(end_obj);

    if (first_length <= 0.0 || second_length <= 0.0 || chord_length <= 0.0) {
        return PyFloat_FromDouble(-INFINITY);
    }

    /* Compute angle between vectors */
    double dot_product = first_x * second_x + first_y * second_y;
    double cosine = clip_double(dot_product / (first_length * second_length), -1.0, 1.0);
    double turn_angle = acos(cosine);

    if (!isfinite(turn_angle)) {
        return PyFloat_FromDouble(-INFINITY);
    }

    /* Score formula matching Python version */
    double straight_ratio = fmin(chord_length / fmax(arrow_length, 1e-12), 1.0);
    double score = straight_ratio - 0.65 * (turn_angle / M_PI);

    return PyFloat_FromDouble(score);
}

/*
 * _compute_local_straightness_score_inline
 *
 * Internal C-only version that avoids Python function call overhead.
 */
static double _compute_local_straightness_score_inline(
    double *vertices_data,
    npy_intp n_vertices,
    npy_intp stride0,
    npy_intp stride1,
    double distance,
    double arrow_length
) {
    if (n_vertices < 2) return -INFINITY;

    /* Compute total length */
    double total = 0.0;
    for (npy_intp i = 0; i < n_vertices - 1; i++) {
        double x0 = vertices_data[i * stride0];
        double y0 = vertices_data[i * stride0 + stride1];
        double x1 = vertices_data[(i + 1) * stride0];
        double y1 = vertices_data[(i + 1) * stride0 + stride1];
        double dx = x1 - x0;
        double dy = y1 - y0;
        total += hypot2(dx, dy);
    }

    if (total <= 0.0) return -INFINITY;

    /* Find three points: start, middle, end */
    double start_dist = fmax(0.0, distance - arrow_length);
    double middle_dist = fmax(0.0, distance - arrow_length * 0.5);
    double end_dist = clip_double(distance, 0.0, total);

    double start_x, start_y, middle_x, middle_y, end_x, end_y;
    bool found_start = false, found_middle = false, found_end = false;

    double cumulative = 0.0;
    for (npy_intp i = 0; i < n_vertices - 1; i++) {
        double x0 = vertices_data[i * stride0];
        double y0 = vertices_data[i * stride0 + stride1];
        double x1 = vertices_data[(i + 1) * stride0];
        double y1 = vertices_data[(i + 1) * stride0 + stride1];
        double dx = x1 - x0;
        double dy = y1 - y0;
        double segment_length = hypot2(dx, dy);

        if (!found_start && cumulative + segment_length >= start_dist) {
            double overshoot = start_dist - cumulative;
            double ratio = (segment_length > 0.0) ? (overshoot / segment_length) : 0.0;
            start_x = x0 + ratio * dx;
            start_y = y0 + ratio * dy;
            found_start = true;
        }

        if (!found_middle && cumulative + segment_length >= middle_dist) {
            double overshoot = middle_dist - cumulative;
            double ratio = (segment_length > 0.0) ? (overshoot / segment_length) : 0.0;
            middle_x = x0 + ratio * dx;
            middle_y = y0 + ratio * dy;
            found_middle = true;
        }

        if (!found_end && cumulative + segment_length >= end_dist) {
            double overshoot = end_dist - cumulative;
            double ratio = (segment_length > 0.0) ? (overshoot / segment_length) : 0.0;
            end_x = x0 + ratio * dx;
            end_y = y0 + ratio * dy;
            found_end = true;
            break;
        }

        cumulative += segment_length;
    }

    if (!found_start || !found_middle || !found_end) return -INFINITY;

    /* Compute vectors: first = middle - start, second = end - middle */
    double first_x = middle_x - start_x;
    double first_y = middle_y - start_y;
    double second_x = end_x - middle_x;
    double second_y = end_y - middle_y;

    double first_length = hypot2(first_x, first_y);
    double second_length = hypot2(second_x, second_y);

    /* Chord from start to end */
    double chord_x = end_x - start_x;
    double chord_y = end_y - start_y;
    double chord_length = hypot2(chord_x, chord_y);

    if (first_length <= 0.0 || second_length <= 0.0 || chord_length <= 0.0) {
        return -INFINITY;
    }

    /* Compute angle between vectors */
    double dot_product = first_x * second_x + first_y * second_y;
    double cosine = clip_double(dot_product / (first_length * second_length), -1.0, 1.0);
    double turn_angle = acos(cosine);

    if (!isfinite(turn_angle)) return -INFINITY;

    /* Score formula matching Python version */
    double straight_ratio = fmin(chord_length / fmax(arrow_length, 1e-12), 1.0);
    double score = straight_ratio - 0.65 * (turn_angle / M_PI);

    return score;
}
typedef struct {
    double score;
    double distance;
} ScoredCandidate;

static int compare_scored_desc(const void *a, const void *b) {
    double diff = ((ScoredCandidate*)b)->score - ((ScoredCandidate*)a)->score;
    if (diff > 0.0) return 1;
    if (diff < 0.0) return -1;
    return 0;
}

/*
 * select_arrow_end_distances(vertices, total_length, arrow_count, arrow_length, min_spacing_override) -> ndarray
 *
 * Select optimal arrow placement distances using straightness scoring.
 */
static PyObject* select_arrow_end_distances(PyObject* self, PyObject* args) {
    PyArrayObject *vertices_obj;
    PyArrayObject *vertices_contiguous = NULL;
    double total_length, arrow_length;
    int arrow_count;
    PyObject *min_spacing_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!didO", &PyArray_Type, &vertices_obj, &total_length,
                          &arrow_count, &arrow_length, &min_spacing_obj)) {
        return NULL;
    }

    /* Check dtype - must be float64 */
    if (PyArray_TYPE(vertices_obj) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "vertices must be float64 (np.float64) array");
        return NULL;
    }

    /* Ensure array is contiguous for safe stride calculation */
    if (!PyArray_ISCARRAY_RO(vertices_obj) && !PyArray_ISFARRAY_RO(vertices_obj)) {
        vertices_contiguous = (PyArrayObject*)PyArray_GETCONTIGUOUS(vertices_obj);
        if (!vertices_contiguous) {
            return NULL;
        }
        vertices_obj = vertices_contiguous;
    }

    if (arrow_count <= 1) {
        /* Simple uniform spacing */
        npy_intp dims[1] = {arrow_count};
        PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (!result) {
            Py_XDECREF(vertices_contiguous);
            return NULL;
        }

        double *result_data = (double*)PyArray_DATA(result);
        if (arrow_count == 1) {
            result_data[0] = total_length / 2.0;
        }
        Py_XDECREF(vertices_contiguous);
        return (PyObject*)result;
    }

    double lower = fmin(total_length, arrow_length * 1.1);
    double upper = fmax(lower, total_length - arrow_length * 0.5);

    if (upper <= lower) {
        /* Fallback to uniform */
        npy_intp dims[1] = {arrow_count};
        PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (!result) {
            Py_XDECREF(vertices_contiguous);
            return NULL;
        }

        double *result_data = (double*)PyArray_DATA(result);
        double spacing = total_length / (arrow_count + 1.0);
        for (int i = 0; i < arrow_count; i++) {
            result_data[i] = spacing * (i + 1);
        }
        Py_XDECREF(vertices_contiguous);
        return (PyObject*)result;
    }

    /* Sample and score candidates */
    int sample_count = (arrow_count * 28 > 64) ? (arrow_count * 28) : 64;
    ScoredCandidate *scored = malloc(sizeof(ScoredCandidate) * sample_count);
    if (!scored) {
        Py_XDECREF(vertices_contiguous);
        PyErr_NoMemory();
        return NULL;
    }

    /* Get vertices data pointers for inline scoring */
    npy_intp n_vertices = PyArray_DIM(vertices_obj, 0);
    double *vertices_data = (double*)PyArray_DATA(vertices_obj);
    npy_intp stride0 = PyArray_STRIDE(vertices_obj, 0) / sizeof(double);
    npy_intp stride1 = PyArray_STRIDE(vertices_obj, 1) / sizeof(double);

    int valid_count = 0;
    for (int i = 0; i < sample_count; i++) {
        double distance = lower + (upper - lower) * i / (sample_count - 1.0);

        /* Use inline C version to avoid Python call overhead */
        double score = _compute_local_straightness_score_inline(
            vertices_data, n_vertices, stride0, stride1, distance, arrow_length
        );

        if (isfinite(score)) {
            scored[valid_count].score = score;
            scored[valid_count].distance = distance;
            valid_count++;
        }
    }

    if (valid_count == 0) {
        free(scored);
        Py_XDECREF(vertices_contiguous);
        /* Fallback to uniform */
        npy_intp dims[1] = {arrow_count};
        PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        if (!result) return NULL;

        double *result_data = (double*)PyArray_DATA(result);
        double spacing = total_length / (arrow_count + 1.0);
        for (int i = 0; i < arrow_count; i++) {
            result_data[i] = spacing * (i + 1);
        }
        return (PyObject*)result;
    }

    /* Sort by score descending */
    qsort(scored, valid_count, sizeof(ScoredCandidate), compare_scored_desc);

    /* Determine minimum spacing */
    double min_spacing;
    if (min_spacing_obj != NULL && min_spacing_obj != Py_None) {
        min_spacing = PyFloat_AsDouble(min_spacing_obj);
        if (PyErr_Occurred()) {
            free(scored);
            Py_XDECREF(vertices_contiguous);
            return NULL;
        }
    } else {
        double adaptive = total_length / fmax(arrow_count * 2.2, 1.0);
        min_spacing = fmax(arrow_length * 2.5, adaptive);
    }

    /* Greedy selection with spacing constraint */
    double *selected = malloc(sizeof(double) * arrow_count);
    if (!selected) {
        free(scored);
        Py_XDECREF(vertices_contiguous);
        PyErr_NoMemory();
        return NULL;
    }

    int selected_count = 0;
    for (int i = 0; i < valid_count && selected_count < arrow_count; i++) {
        double candidate = scored[i].distance;
        bool can_add = true;

        for (int j = 0; j < selected_count; j++) {
            if (fabs(candidate - selected[j]) < min_spacing) {
                can_add = false;
                break;
            }
        }

        if (can_add) {
            selected[selected_count++] = candidate;
        }
    }

    free(scored);

    /* Fill remaining with uniform spacing if needed */
    if (selected_count < arrow_count) {
        double spacing = total_length / (arrow_count + 1.0);
        for (int i = 1; i <= arrow_count && selected_count < arrow_count; i++) {
            double candidate = spacing * i;
            bool can_add = true;

            for (int j = 0; j < selected_count; j++) {
                if (fabs(candidate - selected[j]) < arrow_length) {
                    can_add = false;
                    break;
                }
            }

            if (can_add) {
                selected[selected_count++] = candidate;
            }
        }
    }

    /* Sort selected distances */
    for (int i = 0; i < selected_count - 1; i++) {
        for (int j = i + 1; j < selected_count; j++) {
            if (selected[i] > selected[j]) {
                double temp = selected[i];
                selected[i] = selected[j];
                selected[j] = temp;
            }
        }
    }

    /* Create result array */
    npy_intp dims[1] = {selected_count};
    PyArrayObject *result = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!result) {
        free(selected);
        Py_XDECREF(vertices_contiguous);
        return NULL;
    }

    double *result_data = (double*)PyArray_DATA(result);
    for (int i = 0; i < selected_count; i++) {
        result_data[i] = selected[i];
    }

    free(selected);
    Py_XDECREF(vertices_contiguous);
    return (PyObject*)result;
}

/* Module method table */
static PyMethodDef ContourArrowsCoreMethods[] = {
    {"point_at_distance", point_at_distance, METH_VARARGS,
     "Find point at given arc-length distance along path"},
    {"local_straightness_score", local_straightness_score, METH_VARARGS,
     "Score straightness at given distance for arrow placement"},
    {"select_arrow_end_distances", select_arrow_end_distances, METH_VARARGS,
     "Select optimal arrow placement distances"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef contour_arrows_core_module = {
    PyModuleDef_HEAD_INIT,
    "contour_arrows_core",
    "C-accelerated core functions for arrow_contour",
    -1,
    ContourArrowsCoreMethods
};

/* Module initialization */
PyMODINIT_FUNC PyInit_contour_arrows_core(void) {
    import_array();
    return PyModule_Create(&contour_arrows_core_module);
}
