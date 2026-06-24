#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double score;
    double distance;
} Candidate;

static double hypot2(double x, double y) {
    return hypot(x, y);
}

static int point_at_distance(
    const double *vertices,
    npy_intp nvertices,
    double distance,
    double *out_x,
    double *out_y
) {
    double total = 0.0;
    npy_intp i;

    if (nvertices < 2) {
        return 0;
    }

    for (i = 0; i < nvertices - 1; ++i) {
        const double dx = vertices[(i + 1) * 2] - vertices[i * 2];
        const double dy = vertices[(i + 1) * 2 + 1] - vertices[i * 2 + 1];
        total += hypot2(dx, dy);
    }
    if (total <= 0.0) {
        return 0;
    }
    if (distance < 0.0) {
        distance = 0.0;
    } else if (distance > total) {
        distance = total;
    }

    double cumulative = 0.0;
    for (i = 0; i < nvertices - 1; ++i) {
        const double x0 = vertices[i * 2];
        const double y0 = vertices[i * 2 + 1];
        const double dx = vertices[(i + 1) * 2] - x0;
        const double dy = vertices[(i + 1) * 2 + 1] - y0;
        const double length = hypot2(dx, dy);
        if (length <= 0.0) {
            continue;
        }
        if (cumulative + length >= distance) {
            const double fraction = (distance - cumulative) / length;
            *out_x = x0 + fraction * dx;
            *out_y = y0 + fraction * dy;
            return 1;
        }
        cumulative += length;
    }

    *out_x = vertices[(nvertices - 1) * 2];
    *out_y = vertices[(nvertices - 1) * 2 + 1];
    return 1;
}

static double total_length(const double *vertices, npy_intp nvertices) {
    double total = 0.0;
    npy_intp i;
    for (i = 0; i < nvertices - 1; ++i) {
        const double dx = vertices[(i + 1) * 2] - vertices[i * 2];
        const double dy = vertices[(i + 1) * 2 + 1] - vertices[i * 2 + 1];
        total += hypot2(dx, dy);
    }
    return total;
}

static double local_straightness_score(
    const double *vertices,
    npy_intp nvertices,
    double distance,
    double arrow_length
) {
    double sx, sy, mx, my, ex, ey;
    if (!point_at_distance(vertices, nvertices, fmax(0.0, distance - arrow_length), &sx, &sy) ||
        !point_at_distance(vertices, nvertices, fmax(0.0, distance - arrow_length * 0.5), &mx, &my) ||
        !point_at_distance(vertices, nvertices, distance, &ex, &ey)) {
        return -INFINITY;
    }

    const double fx = mx - sx;
    const double fy = my - sy;
    const double sx2 = ex - mx;
    const double sy2 = ey - my;
    const double first_length = hypot2(fx, fy);
    const double second_length = hypot2(sx2, sy2);
    const double chord_length = hypot2(ex - sx, ey - sy);
    if (first_length <= 0.0 || second_length <= 0.0 || chord_length <= 0.0) {
        return -INFINITY;
    }

    double cosine = (fx * sx2 + fy * sy2) / (first_length * second_length);
    if (cosine < -1.0) {
        cosine = -1.0;
    } else if (cosine > 1.0) {
        cosine = 1.0;
    }
    const double turn_angle = acos(cosine);
    const double straight_ratio = fmin(chord_length / fmax(arrow_length, 1e-12), 1.0);
    return straight_ratio - 0.65 * (turn_angle / M_PI);
}

static int compare_candidates_desc(const void *a, const void *b) {
    const Candidate *ca = (const Candidate *)a;
    const Candidate *cb = (const Candidate *)b;
    if (ca->score < cb->score) {
        return 1;
    }
    if (ca->score > cb->score) {
        return -1;
    }
    return 0;
}

static int select_arrow_distances(
    const double *vertices,
    npy_intp nvertices,
    double path_length,
    int arrow_count,
    double arrow_length,
    double *selected
) {
    int selected_count = 0;
    int i, j;

    if (arrow_count <= 0 || path_length <= 0.0) {
        return 0;
    }

    if (arrow_count <= 1) {
        selected[0] = path_length * 0.5;
        return 1;
    }

    const double lower = fmin(path_length, arrow_length * 1.1);
    const double upper = fmax(lower, path_length - arrow_length * 0.5);
    if (upper <= lower) {
        const double spacing = path_length / (double)arrow_count;
        for (i = 0; i < arrow_count; ++i) {
            selected[i] = spacing * ((double)i + 0.5);
        }
        return arrow_count;
    }

    const int sample_count = arrow_count * 28 > 64 ? arrow_count * 28 : 64;
    Candidate *candidates = (Candidate *)malloc((size_t)sample_count * sizeof(Candidate));
    if (candidates == NULL) {
        PyErr_NoMemory();
        return -1;
    }

    int candidate_count = 0;
    for (i = 0; i < sample_count; ++i) {
        const double fraction = sample_count == 1 ? 0.0 : (double)i / (double)(sample_count - 1);
        const double distance = lower + fraction * (upper - lower);
        const double score = local_straightness_score(vertices, nvertices, distance, arrow_length);
        if (isfinite(score)) {
            candidates[candidate_count].score = score;
            candidates[candidate_count].distance = distance;
            candidate_count++;
        }
    }

    if (candidate_count > 0) {
        qsort(candidates, (size_t)candidate_count, sizeof(Candidate), compare_candidates_desc);
        const double min_spacing = fmax(arrow_length * 2.5, path_length / fmax((double)arrow_count * 2.2, 1.0));
        for (i = 0; i < candidate_count && selected_count < arrow_count; ++i) {
            int ok = 1;
            for (j = 0; j < selected_count; ++j) {
                if (fabs(candidates[i].distance - selected[j]) < min_spacing) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                selected[selected_count++] = candidates[i].distance;
            }
        }
    }
    free(candidates);

    if (selected_count < arrow_count) {
        const double spacing = path_length / (double)arrow_count;
        for (i = 0; i < arrow_count && selected_count < arrow_count; ++i) {
            const double distance = spacing * ((double)i + 0.5);
            int ok = 1;
            for (j = 0; j < selected_count; ++j) {
                if (fabs(distance - selected[j]) < arrow_length) {
                    ok = 0;
                    break;
                }
            }
            if (ok) {
                selected[selected_count++] = distance;
            }
        }
    }

    for (i = 0; i < selected_count - 1; ++i) {
        for (j = i + 1; j < selected_count; ++j) {
            if (selected[j] < selected[i]) {
                const double temp = selected[i];
                selected[i] = selected[j];
                selected[j] = temp;
            }
        }
    }
    return selected_count;
}

static PyObject *build_arrow_segments(PyObject *self, PyObject *args) {
    PyObject *vertices_obj = NULL;
    PyArrayObject *vertices_array = NULL;
    int arrow_count = 1;
    double arrow_length = 0.0;
    double arrow_size = 0.45;

    (void)self;
    if (!PyArg_ParseTuple(
            args,
            "Oid|d:build_arrow_segments",
            &vertices_obj,
            &arrow_count,
            &arrow_length,
            &arrow_size
        )) {
        return NULL;
    }

    vertices_array = (PyArrayObject *)PyArray_FROM_OTF(vertices_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (vertices_array == NULL) {
        return NULL;
    }
    if (PyArray_NDIM(vertices_array) != 2 || PyArray_DIM(vertices_array, 1) != 2) {
        Py_DECREF(vertices_array);
        PyErr_SetString(PyExc_ValueError, "vertices must be a 2D array with shape (n, 2)");
        return NULL;
    }
    const npy_intp nvertices = PyArray_DIM(vertices_array, 0);
    if (nvertices < 2 || arrow_count < 1 || arrow_length <= 0.0 || arrow_size <= 0.0) {
        Py_DECREF(vertices_array);
        npy_intp head_dims[3] = {0, 2, 2};
        npy_intp meta_dims[3] = {0, 2, 2};
        PyObject *empty_heads = PyArray_SimpleNew(3, head_dims, NPY_DOUBLE);
        PyObject *empty_meta = PyArray_SimpleNew(3, meta_dims, NPY_DOUBLE);
        if (empty_heads == NULL || empty_meta == NULL) {
            Py_XDECREF(empty_heads);
            Py_XDECREF(empty_meta);
            return NULL;
        }
        return Py_BuildValue("NN", empty_heads, empty_meta);
    }

    const double *vertices = (const double *)PyArray_DATA(vertices_array);
    const double path_length = total_length(vertices, nvertices);
    double *distances = (double *)malloc((size_t)arrow_count * sizeof(double));
    if (distances == NULL) {
        Py_DECREF(vertices_array);
        return PyErr_NoMemory();
    }

    const int selected_count = select_arrow_distances(
        vertices,
        nvertices,
        path_length,
        arrow_count,
        arrow_length,
        distances
    );
    if (selected_count < 0) {
        free(distances);
        Py_DECREF(vertices_array);
        return NULL;
    }

    npy_intp head_dims[3] = {selected_count * 2, 2, 2};
    npy_intp meta_dims[3] = {selected_count, 2, 2};
    PyArrayObject *head_array = (PyArrayObject *)PyArray_SimpleNew(3, head_dims, NPY_DOUBLE);
    PyArrayObject *meta_array = (PyArrayObject *)PyArray_SimpleNew(3, meta_dims, NPY_DOUBLE);
    if (head_array == NULL || meta_array == NULL) {
        Py_XDECREF(head_array);
        Py_XDECREF(meta_array);
        free(distances);
        Py_DECREF(vertices_array);
        return NULL;
    }

    double *heads = (double *)PyArray_DATA(head_array);
    double *meta = (double *)PyArray_DATA(meta_array);
    int written = 0;
    int i;
    for (i = 0; i < selected_count; ++i) {
        double sx, sy, ex, ey;
        if (!point_at_distance(vertices, nvertices, fmax(0.0, distances[i] - arrow_length), &sx, &sy) ||
            !point_at_distance(vertices, nvertices, distances[i], &ex, &ey)) {
            continue;
        }
        const double vx = ex - sx;
        const double vy = ey - sy;
        const double vector_length = hypot2(vx, vy);
        if (vector_length <= 0.0) {
            continue;
        }
        const double tx = vx / vector_length;
        const double ty = vy / vector_length;
        const double nx = -ty;
        const double ny = tx;
        const double width = vector_length * arrow_size;
        const double bx = sx;
        const double by = sy;

        double *h0 = heads + (written * 2) * 4;
        h0[0] = ex;
        h0[1] = ey;
        h0[2] = bx + nx * width * 0.5;
        h0[3] = by + ny * width * 0.5;
        h0[4] = ex;
        h0[5] = ey;
        h0[6] = bx - nx * width * 0.5;
        h0[7] = by - ny * width * 0.5;

        double *m0 = meta + written * 4;
        m0[0] = sx;
        m0[1] = sy;
        m0[2] = ex;
        m0[3] = ey;
        written++;
    }

    free(distances);
    Py_DECREF(vertices_array);

    return Py_BuildValue("NN", (PyObject *)head_array, (PyObject *)meta_array);
}

static PyMethodDef ContourCoreMethods[] = {
    {
        "build_arrow_segments",
        (PyCFunction)build_arrow_segments,
        METH_VARARGS,
        "Build display-space line-arrow head segments for arrow_contour."
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_contour_core",
    "Native contour plotting helpers.",
    -1,
    ContourCoreMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__contour_core(void) {
    PyObject *module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }
    import_array();
    return module;
}
