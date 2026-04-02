/*
 * Native C core for Skyborn NCL-like curly-vector tracing.
 *
 * Author: Qianye Su <suqianye2000@gmail.com>
 * Copyright (c) 2025-2026 Qianye Su
 * Created: 2026-03-01 14:58:56
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_19_API_VERSION
#include <numpy/arrayobject.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    npy_intp nx;
    npy_intp ny;
    double x_origin;
    double y_origin;
    double dx_safe;
    double dy_safe;
    double inv_dx;
    double inv_dy;
    double max_x_index;
    double max_y_index;
    double max_cell_x;
    double max_cell_y;
} GridInfo;

/* Interpolated flow state and local data->display Jacobian at one sample point. */
typedef struct
{
    double display_x;
    double display_y;
    double j00;
    double j01;
    double j10;
    double j11;
    double dir_x;
    double dir_y;
    double speed;
} SampleState;

/* Bucketed display-space candidate index used by the thinning pass. */
typedef struct
{
    long bx;
    long by;
    npy_intp idx;
} BucketEntry;

/*
 * Sample one vector state on the regular data grid.
 *
 * This performs bilinear interpolation of u/v in data space, then uses the
 * local display-grid Jacobian to estimate the display-space tangent that drives
 * the native NCL-like stepping loop.
 */
static int sample_state(
    const GridInfo *grid,
    const double *u_data,
    const double *v_data,
    const double *display_data,
    const npy_bool *cell_valid,
    double x,
    double y,
    SampleState *out)
{
    double xi;
    double yi;
    npy_intp ix;
    npy_intp iy;
    npy_intp nx;
    npy_intp idx00;
    npy_intp idx01;
    npy_intp idx10;
    npy_intp idx11;
    double sx;
    double sy;
    double one_minus_sx;
    double one_minus_sy;
    double u00;
    double u01;
    double u10;
    double u11;
    double v00;
    double v01;
    double v10;
    double v11;
    double u_value;
    double v_value;
    double p00x;
    double p00y;
    double p01x;
    double p01y;
    double p10x;
    double p10y;
    double p11x;
    double p11y;
    double det;
    double vec_x;
    double vec_y;
    double vec_norm;

    if (grid->nx < 2 || grid->ny < 2)
    {
        return 0;
    }

    xi = (x - grid->x_origin) * grid->inv_dx;
    yi = (y - grid->y_origin) * grid->inv_dy;
    if (!(0.0 <= xi && xi <= grid->max_x_index && 0.0 <= yi && yi <= grid->max_y_index))
    {
        return 0;
    }

    ix = (npy_intp)floor(fmin(fmax(xi, 0.0), grid->max_cell_x));
    iy = (npy_intp)floor(fmin(fmax(yi, 0.0), grid->max_cell_y));
    if (!cell_valid[iy * (grid->nx - 1) + ix])
    {
        return 0;
    }

    sx = xi - (double)ix;
    sy = yi - (double)iy;
    one_minus_sx = 1.0 - sx;
    one_minus_sy = 1.0 - sy;
    nx = grid->nx;

    idx00 = iy * nx + ix;
    idx01 = iy * nx + (ix + 1);
    idx10 = (iy + 1) * nx + ix;
    idx11 = (iy + 1) * nx + (ix + 1);

    u00 = u_data[idx00];
    u01 = u_data[idx01];
    u10 = u_data[idx10];
    u11 = u_data[idx11];
    v00 = v_data[idx00];
    v01 = v_data[idx01];
    v10 = v_data[idx10];
    v11 = v_data[idx11];

    u_value = (u00 * one_minus_sx + u01 * sx) * one_minus_sy + (u10 * one_minus_sx + u11 * sx) * sy;
    v_value = (v00 * one_minus_sx + v01 * sx) * one_minus_sy + (v10 * one_minus_sx + v11 * sx) * sy;
    if (!isfinite(u_value) || !isfinite(v_value))
    {
        return 0;
    }

    idx00 *= 2;
    idx01 *= 2;
    idx10 *= 2;
    idx11 *= 2;

    p00x = display_data[idx00];
    p00y = display_data[idx00 + 1];
    p01x = display_data[idx01];
    p01y = display_data[idx01 + 1];
    p10x = display_data[idx10];
    p10y = display_data[idx10 + 1];
    p11x = display_data[idx11];
    p11y = display_data[idx11 + 1];

    out->display_x = p00x * one_minus_sx * one_minus_sy + p01x * sx * one_minus_sy + p10x * one_minus_sx * sy + p11x * sx * sy;
    out->display_y = p00y * one_minus_sx * one_minus_sy + p01y * sx * one_minus_sy + p10y * one_minus_sx * sy + p11y * sx * sy;
    out->j00 = ((p01x - p00x) * one_minus_sy + (p11x - p10x) * sy) / grid->dx_safe;
    out->j10 = ((p01y - p00y) * one_minus_sy + (p11y - p10y) * sy) / grid->dx_safe;
    out->j01 = ((p10x - p00x) * one_minus_sx + (p11x - p01x) * sx) / grid->dy_safe;
    out->j11 = ((p10y - p00y) * one_minus_sx + (p11y - p01y) * sx) / grid->dy_safe;

    det = out->j00 * out->j11 - out->j01 * out->j10;
    if (!isfinite(out->display_x) || !isfinite(out->display_y) || !isfinite(det) || fabs(det) <= 1e-12)
    {
        return 0;
    }

    vec_x = out->j00 * u_value + out->j01 * v_value;
    vec_y = out->j10 * u_value + out->j11 * v_value;
    vec_norm = hypot(vec_x, vec_y);
    if (!isfinite(vec_norm) || vec_norm <= 1e-12)
    {
        return 0;
    }

    out->dir_x = vec_x / vec_norm;
    out->dir_y = vec_y / vec_norm;
    out->speed = hypot(u_value, v_value);
    return isfinite(out->speed) ? 1 : 0;
}

/* Convert one display-space step back into a data-space increment. */
static int display_step_to_data(
    double j00,
    double j01,
    double j10,
    double j11,
    double step_x,
    double step_y,
    double *data_x,
    double *data_y)
{
    double det = j00 * j11 - j01 * j10;
    if (!isfinite(det) || fabs(det) <= 1e-12)
    {
        return 0;
    }

    *data_x = (j11 * step_x - j01 * step_y) / det;
    *data_y = (-j10 * step_x + j00 * step_y) / det;
    return isfinite(*data_x) && isfinite(*data_y);
}

/* Match the Python renderer's NCL-like speed-dependent pixel step law. */
static double ncl_step_length_px(double base_step_px, double local_speed, double speed_scale)
{
    double speed_fraction;
    if (speed_scale <= 1e-12)
    {
        speed_scale = 1e-12;
    }
    speed_fraction = local_speed / speed_scale;
    if (speed_fraction < 0.0)
    {
        speed_fraction = 0.0;
    }
    else if (speed_fraction > 1.0)
    {
        speed_fraction = 1.0;
    }
    {
        double step = base_step_px * speed_fraction * speed_fraction;
        return (step > 0.35) ? step : 0.35;
    }
}

static int point_within_grid_data(const GridInfo *grid, double x, double y)
{
    double xi = (x - grid->x_origin) * grid->inv_dx;
    double yi = (y - grid->y_origin) * grid->inv_dy;
    return 0.0 <= xi && xi <= grid->max_x_index && 0.0 <= yi && yi <= grid->max_y_index;
}

/* Bilinear sampling helper for scalar fields on the same regular grid. */
static int sample_scalar_field(
    const GridInfo *grid,
    const double *field_data,
    double x,
    double y,
    double *out)
{
    double xi;
    double yi;
    npy_intp ix;
    npy_intp iy;
    npy_intp ix_next;
    npy_intp iy_next;
    double sx;
    double sy;
    double a00;
    double a01;
    double a10;
    double a11;
    double value;

    if (grid->nx < 1 || grid->ny < 1)
    {
        return 0;
    }

    xi = (x - grid->x_origin) * grid->inv_dx;
    yi = (y - grid->y_origin) * grid->inv_dy;
    if (!(0.0 <= xi && xi <= grid->max_x_index && 0.0 <= yi && yi <= grid->max_y_index))
    {
        return 0;
    }

    ix = (npy_intp)xi;
    iy = (npy_intp)yi;
    ix_next = (ix == grid->nx - 1) ? ix : (ix + 1);
    iy_next = (iy == grid->ny - 1) ? iy : (iy + 1);

    sx = xi - (double)ix;
    sy = yi - (double)iy;

    a00 = field_data[iy * grid->nx + ix];
    a01 = field_data[iy * grid->nx + ix_next];
    a10 = field_data[iy_next * grid->nx + ix];
    a11 = field_data[iy_next * grid->nx + ix_next];

    value = ((a00 * (1.0 - sx) + a01 * sx) * (1.0 - sy) + (a10 * (1.0 - sx) + a11 * sx) * sy);
    if (!isfinite(value))
    {
        return 0;
    }

    *out = value;
    return 1;
}

static int compare_bucket_entries(const void *left_ptr, const void *right_ptr)
{
    const BucketEntry *left = (const BucketEntry *)left_ptr;
    const BucketEntry *right = (const BucketEntry *)right_ptr;

    if (left->bx < right->bx)
    {
        return -1;
    }
    if (left->bx > right->bx)
    {
        return 1;
    }
    if (left->by < right->by)
    {
        return -1;
    }
    if (left->by > right->by)
    {
        return 1;
    }
    if (left->idx < right->idx)
    {
        return -1;
    }
    if (left->idx > right->idx)
    {
        return 1;
    }
    return 0;
}

static npy_intp lower_bound_bucket(
    const BucketEntry *entries,
    npy_intp count,
    long bx,
    long by)
{
    npy_intp left = 0;
    npy_intp right = count;

    while (left < right)
    {
        npy_intp middle = left + (right - left) / 2;
        const BucketEntry *entry = &entries[middle];
        if (entry->bx < bx || (entry->bx == bx && entry->by < by))
        {
            left = middle + 1;
        }
        else
        {
            right = middle;
        }
    }

    return left;
}

static int viewport_contains(double x0, double y0, double x1, double y1, double x, double y)
{
    return isfinite(x) && isfinite(y) && x0 <= x && x <= x1 && y0 <= y && y <= y1;
}

/*
 * Clip a candidate display-space segment against the active viewport.
 *
 * The tracer still records the shortened step, then terminates, which keeps the
 * last arrow tail inside the map boundary instead of overshooting the panel.
 */
static int clip_display_step_to_viewport(
    double start_x,
    double start_y,
    double *end_x,
    double *end_y,
    double viewport_x0,
    double viewport_y0,
    double viewport_x1,
    double viewport_y1)
{
    double delta_x;
    double delta_y;
    double factor;

    if (viewport_contains(viewport_x0, viewport_y0, viewport_x1, viewport_y1, *end_x, *end_y))
    {
        return 0;
    }

    delta_x = *end_x - start_x;
    delta_y = *end_y - start_y;
    if (!isfinite(delta_x) || !isfinite(delta_y))
    {
        *end_x = start_x;
        *end_y = start_y;
        return 1;
    }

    factor = 1.0;
    if (delta_x < -1e-12)
    {
        double candidate = (viewport_x0 - start_x) / delta_x;
        if (candidate < factor)
        {
            factor = candidate;
        }
    }
    else if (delta_x > 1e-12)
    {
        double candidate = (viewport_x1 - start_x) / delta_x;
        if (candidate < factor)
        {
            factor = candidate;
        }
    }

    if (delta_y < -1e-12)
    {
        double candidate = (viewport_y0 - start_y) / delta_y;
        if (candidate < factor)
        {
            factor = candidate;
        }
    }
    else if (delta_y > 1e-12)
    {
        double candidate = (viewport_y1 - start_y) / delta_y;
        if (candidate < factor)
        {
            factor = candidate;
        }
    }

    if (factor < 0.0)
    {
        factor = 0.0;
    }
    else if (factor > 1.0)
    {
        factor = 1.0;
    }

    *end_x = start_x + delta_x * factor;
    *end_y = start_y + delta_y * factor;
    return 1;
}

/*
 * Trace a single forward or backward curly-vector branch.
 *
 * The caller supplies data-space u/v, a precomputed display grid, and the valid
 * cell mask. The loop works in display space, applies the one-third backward
 * correction used by the Python implementation, converts the chosen step back
 * into data coordinates, and returns a compact (N, 2) curve.
 */
static PyObject *trace_ncl_direction(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *u_obj;
    PyObject *v_obj;
    PyObject *display_grid_obj;
    PyObject *cell_valid_obj;
    PyArrayObject *u_arr = NULL;
    PyArrayObject *v_arr = NULL;
    PyArrayObject *display_arr = NULL;
    PyArrayObject *cell_valid_arr = NULL;
    PyArrayObject *full_output = NULL;
    PyObject *result = NULL;
    double x_origin;
    double y_origin;
    double dx;
    double dy;
    double start_x;
    double start_y;
    double max_length_px;
    double direction_sign;
    double step_px;
    double speed_scale;
    double viewport_x0;
    double viewport_y0;
    double viewport_x1;
    double viewport_y1;
    int max_steps = 512;
    static char *kwlist[] = {
        "u",
        "v",
        "display_grid",
        "cell_valid",
        "x_origin",
        "y_origin",
        "dx",
        "dy",
        "start_x",
        "start_y",
        "max_length_px",
        "direction_sign",
        "step_px",
        "speed_scale",
        "viewport_x0",
        "viewport_y0",
        "viewport_x1",
        "viewport_y1",
        "max_steps",
        NULL,
    };
    GridInfo grid;
    const double *u_data;
    const double *v_data;
    const double *display_data;
    const npy_bool *cell_valid_data;
    npy_intp nx;
    npy_intp ny;
    npy_intp point_count = 1;
    npy_intp output_dims[2];
    double *output_data;
    double current_data_x;
    double current_data_y;
    double current_display_x;
    double current_display_y;
    double previous_display_x = 0.0;
    double previous_display_y = 0.0;
    int has_previous_display = 0;
    double travelled = 0.0;
    SampleState initial_state;
    int step_index;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOdddddddddddddd|i:trace_ncl_direction",
            kwlist,
            &u_obj,
            &v_obj,
            &display_grid_obj,
            &cell_valid_obj,
            &x_origin,
            &y_origin,
            &dx,
            &dy,
            &start_x,
            &start_y,
            &max_length_px,
            &direction_sign,
            &step_px,
            &speed_scale,
            &viewport_x0,
            &viewport_y0,
            &viewport_x1,
            &viewport_y1,
            &max_steps))
    {
        return NULL;
    }

    if (max_steps < 1)
    {
        max_steps = 1;
    }
    else if (max_steps > 4096)
    {
        max_steps = 4096;
    }

    u_arr = (PyArrayObject *)PyArray_FROM_OTF(u_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    v_arr = (PyArrayObject *)PyArray_FROM_OTF(v_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    display_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_grid_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    cell_valid_arr = (PyArrayObject *)PyArray_FROM_OTF(
        cell_valid_obj, NPY_BOOL, NPY_ARRAY_IN_ARRAY);
    if (u_arr == NULL || v_arr == NULL || display_arr == NULL || cell_valid_arr == NULL)
    {
        goto cleanup;
    }

    if (PyArray_NDIM(u_arr) != 2 || PyArray_NDIM(v_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "u and v must be 2D float64 arrays");
        goto cleanup;
    }
    if (PyArray_NDIM(display_arr) != 3 || PyArray_DIM(display_arr, 2) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_grid must have shape (ny, nx, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(cell_valid_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "cell_valid must be a 2D boolean array");
        goto cleanup;
    }

    ny = PyArray_DIM(u_arr, 0);
    nx = PyArray_DIM(u_arr, 1);
    if (PyArray_DIM(v_arr, 0) != ny || PyArray_DIM(v_arr, 1) != nx)
    {
        PyErr_SetString(PyExc_ValueError, "u and v must share the same shape");
        goto cleanup;
    }
    if (PyArray_DIM(display_arr, 0) != ny || PyArray_DIM(display_arr, 1) != nx)
    {
        PyErr_SetString(PyExc_ValueError, "display_grid must match the u/v grid shape");
        goto cleanup;
    }
    if (ny < 2 || nx < 2)
    {
        Py_INCREF(Py_None);
        result = Py_None;
        goto cleanup;
    }
    if (PyArray_DIM(cell_valid_arr, 0) != ny - 1 || PyArray_DIM(cell_valid_arr, 1) != nx - 1)
    {
        PyErr_SetString(PyExc_ValueError, "cell_valid must have shape (ny-1, nx-1)");
        goto cleanup;
    }

    grid.nx = nx;
    grid.ny = ny;
    grid.x_origin = x_origin;
    grid.y_origin = y_origin;
    grid.dx_safe = fabs(dx) > 1e-12 ? fabs(dx) : 1e-12;
    grid.dy_safe = fabs(dy) > 1e-12 ? fabs(dy) : 1e-12;
    grid.inv_dx = 1.0 / grid.dx_safe;
    grid.inv_dy = 1.0 / grid.dy_safe;
    grid.max_x_index = (double)(nx - 1);
    grid.max_y_index = (double)(ny - 1);
    grid.max_cell_x = (double)(nx - 2);
    grid.max_cell_y = (double)(ny - 2);

    u_data = (const double *)PyArray_DATA(u_arr);
    v_data = (const double *)PyArray_DATA(v_arr);
    display_data = (const double *)PyArray_DATA(display_arr);
    cell_valid_data = (const npy_bool *)PyArray_DATA(cell_valid_arr);

    if (!sample_state(
            &grid,
            u_data,
            v_data,
            display_data,
            cell_valid_data,
            start_x,
            start_y,
            &initial_state))
    {
        Py_INCREF(Py_None);
        result = Py_None;
        goto cleanup;
    }

    output_dims[0] = max_steps + 1;
    output_dims[1] = 2;
    full_output = (PyArrayObject *)PyArray_SimpleNew(2, output_dims, NPY_DOUBLE);
    if (full_output == NULL)
    {
        goto cleanup;
    }

    output_data = (double *)PyArray_DATA(full_output);
    output_data[0] = start_x;
    output_data[1] = start_y;
    current_data_x = start_x;
    current_data_y = start_y;
    current_display_x = initial_state.display_x;
    current_display_y = initial_state.display_y;

    Py_BEGIN_ALLOW_THREADS

    for (step_index = 0; step_index < max_steps; ++step_index)
    {
        SampleState state;
        double remaining = max_length_px - travelled;
        double step_length;
        double corrected_display_x;
        double corrected_display_y;
        double candidate_display_x;
        double candidate_display_y;
        double display_step_x;
        double display_step_y;
        double data_step_x;
        double data_step_y;
        double candidate_x;
        double candidate_y;
        int clipped;
        double actual_step;

        if (remaining <= 1e-6)
        {
            break;
        }

        /* Re-sample the local state after the first accepted step. */
        if (step_index == 0)
        {
            state = initial_state;
        }
        else if (!sample_state(
                     &grid,
                     u_data,
                     v_data,
                     display_data,
                     cell_valid_data,
                     current_data_x,
                     current_data_y,
                     &state))
        {
            break;
        }

        /* Build the nominal NCL-like step length in display pixels. */
        current_display_x = state.display_x;
        current_display_y = state.display_y;
        step_length = ncl_step_length_px(step_px, state.speed, speed_scale);
        if (step_length > remaining)
        {
            step_length = remaining;
        }
        if (step_length <= 1e-6)
        {
            break;
        }

        /* Apply the one-third backward correction before stepping away again. */
        if (has_previous_display)
        {
            corrected_display_x = current_display_x - (current_display_x - previous_display_x) / 3.0;
            corrected_display_y = current_display_y - (current_display_y - previous_display_y) / 3.0;
        }
        else
        {
            corrected_display_x = current_display_x;
            corrected_display_y = current_display_y;
        }

        /* First candidate uses the local tangent at the current sample point. */
        candidate_display_x = corrected_display_x + direction_sign * state.dir_x * step_length;
        candidate_display_y = corrected_display_y + direction_sign * state.dir_y * step_length;
        clipped = clip_display_step_to_viewport(
            corrected_display_x,
            corrected_display_y,
            &candidate_display_x,
            &candidate_display_y,
            viewport_x0,
            viewport_y0,
            viewport_x1,
            viewport_y1);

        display_step_x = candidate_display_x - current_display_x;
        display_step_y = candidate_display_y - current_display_y;
        if (!display_step_to_data(
                state.j00,
                state.j01,
                state.j10,
                state.j11,
                display_step_x,
                display_step_y,
                &data_step_x,
                &data_step_y))
        {
            break;
        }
        candidate_x = current_data_x + data_step_x;
        candidate_y = current_data_y + data_step_y;
        if (!isfinite(candidate_x) || !isfinite(candidate_y) ||
            !point_within_grid_data(&grid, candidate_x, candidate_y))
        {
            break;
        }

        /*
         * If the next sample is valid, blend the two local tangents and local
         * speeds to reduce sharp kinks. This mirrors the smoother second-pass
         * behavior used by the Python renderer.
         */
        {
            SampleState next_state;
            if (sample_state(
                    &grid,
                    u_data,
                    v_data,
                    display_data,
                    cell_valid_data,
                    candidate_x,
                    candidate_y,
                    &next_state))
            {
                double average_direction_x = direction_sign * (state.dir_x + next_state.dir_x);
                double average_direction_y = direction_sign * (state.dir_y + next_state.dir_y);
                double average_norm = hypot(average_direction_x, average_direction_y);
                if (average_norm > 1e-12)
                {
                    double average_speed = 0.5 * (state.speed + next_state.speed);
                    double average_step = ncl_step_length_px(step_px, average_speed, speed_scale);
                    double avg_j00;
                    double avg_j01;
                    double avg_j10;
                    double avg_j11;
                    if (average_step > remaining)
                    {
                        average_step = remaining;
                    }
                    candidate_display_x = corrected_display_x + average_direction_x / average_norm * average_step;
                    candidate_display_y = corrected_display_y + average_direction_y / average_norm * average_step;
                    clipped = clip_display_step_to_viewport(
                        corrected_display_x,
                        corrected_display_y,
                        &candidate_display_x,
                        &candidate_display_y,
                        viewport_x0,
                        viewport_y0,
                        viewport_x1,
                        viewport_y1);

                    avg_j00 = 0.5 * (state.j00 + next_state.j00);
                    avg_j01 = 0.5 * (state.j01 + next_state.j01);
                    avg_j10 = 0.5 * (state.j10 + next_state.j10);
                    avg_j11 = 0.5 * (state.j11 + next_state.j11);
                    display_step_x = candidate_display_x - current_display_x;
                    display_step_y = candidate_display_y - current_display_y;
                    if (!display_step_to_data(
                            avg_j00,
                            avg_j01,
                            avg_j10,
                            avg_j11,
                            display_step_x,
                            display_step_y,
                            &data_step_x,
                            &data_step_y))
                    {
                        break;
                    }
                    candidate_x = current_data_x + data_step_x;
                    candidate_y = current_data_y + data_step_y;
                    if (!isfinite(candidate_x) || !isfinite(candidate_y) ||
                        !point_within_grid_data(&grid, candidate_x, candidate_y))
                    {
                        break;
                    }
                }
            }
        }

        /* Ignore effectively zero-length steps after clipping or inversion. */
        actual_step = hypot(candidate_display_x - current_display_x, candidate_display_y - current_display_y);
        if (actual_step <= 0.2)
        {
            break;
        }

        output_data[point_count * 2] = candidate_x;
        output_data[point_count * 2 + 1] = candidate_y;
        point_count += 1;

        previous_display_x = current_display_x;
        previous_display_y = current_display_y;
        has_previous_display = 1;
        current_data_x = candidate_x;
        current_data_y = candidate_y;
        travelled += actual_step;

        if (clipped)
        {
            break;
        }
    }

    Py_END_ALLOW_THREADS

    if (point_count < 2)
    {
        Py_INCREF(Py_None);
        result = Py_None;
        goto cleanup;
    }

    if (point_count == output_dims[0])
    {
        result = (PyObject *)full_output;
        full_output = NULL;
        goto cleanup;
    }

    {
        npy_intp trimmed_dims[2];
        PyArrayObject *trimmed_output;

        trimmed_dims[0] = point_count;
        trimmed_dims[1] = 2;
        trimmed_output = (PyArrayObject *)PyArray_SimpleNew(2, trimmed_dims, NPY_DOUBLE);
        if (trimmed_output == NULL)
        {
            goto cleanup;
        }
        memcpy(
            PyArray_DATA(trimmed_output),
            PyArray_DATA(full_output),
            (size_t)(point_count * 2) * sizeof(double));
        result = (PyObject *)trimmed_output;
    }

cleanup:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(display_arr);
    Py_XDECREF(cell_valid_arr);
    Py_XDECREF(full_output);
    return result;
}

/* Python wrapper for scalar sampling at one point. */
static PyObject *sample_grid_field(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *field_obj;
    PyArrayObject *field_arr = NULL;
    PyObject *result = NULL;
    double x_origin;
    double y_origin;
    double dx;
    double dy;
    double x;
    double y;
    GridInfo grid;
    double value;
    static char *kwlist[] = {
        "field",
        "x_origin",
        "y_origin",
        "dx",
        "dy",
        "x",
        "y",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "Odddddd:sample_grid_field",
            kwlist,
            &field_obj,
            &x_origin,
            &y_origin,
            &dx,
            &dy,
            &x,
            &y))
    {
        return NULL;
    }

    field_arr = (PyArrayObject *)PyArray_FROM_OTF(field_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (field_arr == NULL)
    {
        goto cleanup;
    }
    if (PyArray_NDIM(field_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "field must be a 2D float64 array");
        goto cleanup;
    }

    grid.nx = PyArray_DIM(field_arr, 1);
    grid.ny = PyArray_DIM(field_arr, 0);
    grid.x_origin = x_origin;
    grid.y_origin = y_origin;
    grid.dx_safe = fabs(dx) > 1e-12 ? fabs(dx) : 1e-12;
    grid.dy_safe = fabs(dy) > 1e-12 ? fabs(dy) : 1e-12;
    grid.inv_dx = 1.0 / grid.dx_safe;
    grid.inv_dy = 1.0 / grid.dy_safe;
    grid.max_x_index = (double)(grid.nx - 1);
    grid.max_y_index = (double)(grid.ny - 1);
    grid.max_cell_x = (double)(grid.nx - 2);
    grid.max_cell_y = (double)(grid.ny - 2);

    if (!sample_scalar_field(
            &grid,
            (const double *)PyArray_DATA(field_arr),
            x,
            y,
            &value))
    {
        Py_INCREF(Py_None);
        result = Py_None;
        goto cleanup;
    }

    result = PyFloat_FromDouble(value);

cleanup:
    Py_XDECREF(field_arr);
    return result;
}

/* Python wrapper for vectorized scalar sampling at many points. */
static PyObject *sample_grid_field_array(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *field_obj;
    PyObject *points_obj;
    PyArrayObject *field_arr = NULL;
    PyArrayObject *points_arr = NULL;
    PyArrayObject *output_arr = NULL;
    PyObject *result = NULL;
    double x_origin;
    double y_origin;
    double dx;
    double dy;
    GridInfo grid;
    npy_intp point_count;
    double *points_data;
    double *output_data;
    npy_intp output_dims[1];
    npy_intp idx;
    static char *kwlist[] = {
        "field",
        "x_origin",
        "y_origin",
        "dx",
        "dy",
        "points",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OddddO:sample_grid_field_array",
            kwlist,
            &field_obj,
            &x_origin,
            &y_origin,
            &dx,
            &dy,
            &points_obj))
    {
        return NULL;
    }

    field_arr = (PyArrayObject *)PyArray_FROM_OTF(field_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    points_arr = (PyArrayObject *)PyArray_FROM_OTF(points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (field_arr == NULL || points_arr == NULL)
    {
        goto cleanup;
    }
    if (PyArray_NDIM(field_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "field must be a 2D float64 array");
        goto cleanup;
    }
    if (PyArray_NDIM(points_arr) == 1)
    {
        if (PyArray_DIM(points_arr, 0) != 2)
        {
            PyErr_SetString(PyExc_ValueError, "1D points must have length 2");
            goto cleanup;
        }
        point_count = 1;
    }
    else if (PyArray_NDIM(points_arr) == 2 && PyArray_DIM(points_arr, 1) == 2)
    {
        point_count = PyArray_DIM(points_arr, 0);
    }
    else
    {
        PyErr_SetString(PyExc_ValueError, "points must have shape (N, 2) or (2,)");
        goto cleanup;
    }

    grid.nx = PyArray_DIM(field_arr, 1);
    grid.ny = PyArray_DIM(field_arr, 0);
    grid.x_origin = x_origin;
    grid.y_origin = y_origin;
    grid.dx_safe = fabs(dx) > 1e-12 ? fabs(dx) : 1e-12;
    grid.dy_safe = fabs(dy) > 1e-12 ? fabs(dy) : 1e-12;
    grid.inv_dx = 1.0 / grid.dx_safe;
    grid.inv_dy = 1.0 / grid.dy_safe;
    grid.max_x_index = (double)(grid.nx - 1);
    grid.max_y_index = (double)(grid.ny - 1);
    grid.max_cell_x = (double)(grid.nx - 2);
    grid.max_cell_y = (double)(grid.ny - 2);

    output_dims[0] = point_count;
    output_arr = (PyArrayObject *)PyArray_SimpleNew(1, output_dims, NPY_DOUBLE);
    if (output_arr == NULL)
    {
        goto cleanup;
    }

    output_data = (double *)PyArray_DATA(output_arr);
    for (idx = 0; idx < point_count; ++idx)
    {
        output_data[idx] = Py_NAN;
    }

    points_data = (double *)PyArray_DATA(points_arr);
    for (idx = 0; idx < point_count; ++idx)
    {
        double value;
        double point_x = points_data[idx * 2];
        double point_y = points_data[idx * 2 + 1];
        if (sample_scalar_field(
                &grid,
                (const double *)PyArray_DATA(field_arr),
                point_x,
                point_y,
                &value))
        {
            output_data[idx] = value;
        }
    }

    result = (PyObject *)output_arr;
    output_arr = NULL;

cleanup:
    Py_XDECREF(field_arr);
    Py_XDECREF(points_arr);
    Py_XDECREF(output_arr);
    return result;
}

/*
 * Thin candidate centers in mapped/display space.
 *
 * Candidates are bucketed on a spacing-sized lattice so each accepted point
 * only checks its own and neighboring buckets instead of scanning all N points.
 */
static PyObject *thin_mapped_candidates(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *mapped_points_obj;
    PyArrayObject *mapped_points_arr = NULL;
    PyArrayObject *selected_arr = NULL;
    PyObject *result = NULL;
    double spacing_frac;
    double spacing_sq;
    double bucket_scale;
    npy_intp point_count;
    BucketEntry *entries = NULL;
    long *bucket_x = NULL;
    long *bucket_y = NULL;
    unsigned char *culled = NULL;
    npy_intp *selected_indices = NULL;
    npy_intp selected_count = 0;
    double *mapped_points_data;
    npy_intp idx;
    static char *kwlist[] = {
        "mapped_points",
        "spacing_frac",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "Od:thin_mapped_candidates",
            kwlist,
            &mapped_points_obj,
            &spacing_frac))
    {
        return NULL;
    }

    mapped_points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        mapped_points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (mapped_points_arr == NULL)
    {
        goto cleanup;
    }
    if (PyArray_NDIM(mapped_points_arr) != 2 || PyArray_DIM(mapped_points_arr, 1) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "mapped_points must have shape (N, 2)");
        goto cleanup;
    }

    point_count = PyArray_DIM(mapped_points_arr, 0);
    if (point_count == 0)
    {
        npy_intp dims[1] = {0};
        selected_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INTP);
        if (selected_arr == NULL)
        {
            goto cleanup;
        }
        result = (PyObject *)selected_arr;
        selected_arr = NULL;
        goto cleanup;
    }

    spacing_frac = fmax(spacing_frac, 1e-6);
    spacing_sq = spacing_frac * spacing_frac;
    bucket_scale = 1.0 / spacing_frac;
    mapped_points_data = (double *)PyArray_DATA(mapped_points_arr);

    entries = (BucketEntry *)PyMem_Malloc((size_t)point_count * sizeof(BucketEntry));
    bucket_x = (long *)PyMem_Malloc((size_t)point_count * sizeof(long));
    bucket_y = (long *)PyMem_Malloc((size_t)point_count * sizeof(long));
    culled = (unsigned char *)PyMem_Calloc((size_t)point_count, sizeof(unsigned char));
    selected_indices = (npy_intp *)PyMem_Malloc((size_t)point_count * sizeof(npy_intp));
    if (
        entries == NULL || bucket_x == NULL || bucket_y == NULL || culled == NULL || selected_indices == NULL)
    {
        PyErr_NoMemory();
        goto cleanup;
    }

    for (idx = 0; idx < point_count; ++idx)
    {
        long bx = (long)floor(mapped_points_data[idx * 2] * bucket_scale);
        long by = (long)floor(mapped_points_data[idx * 2 + 1] * bucket_scale);
        bucket_x[idx] = bx;
        bucket_y[idx] = by;
        entries[idx].bx = bx;
        entries[idx].by = by;
        entries[idx].idx = idx;
    }

    qsort(entries, (size_t)point_count, sizeof(BucketEntry), compare_bucket_entries);

    for (idx = 0; idx < point_count; ++idx)
    {
        npy_intp pos;
        double base_x;
        double base_y;
        long bx;
        long by;
        long search_x;

        if (culled[idx])
        {
            continue;
        }

        selected_indices[selected_count++] = idx;
        base_x = mapped_points_data[idx * 2];
        base_y = mapped_points_data[idx * 2 + 1];
        bx = bucket_x[idx];
        by = bucket_y[idx];

        for (search_x = bx - 1; search_x <= bx + 1; ++search_x)
        {
            long search_y;
            for (search_y = by - 1; search_y <= by + 1; ++search_y)
            {
                npy_intp lower = lower_bound_bucket(entries, point_count, search_x, search_y);
                for (pos = lower; pos < point_count; ++pos)
                {
                    npy_intp other_idx;
                    double dx_value;
                    double dy_value;
                    const BucketEntry *entry = &entries[pos];

                    if (entry->bx != search_x || entry->by != search_y)
                    {
                        break;
                    }

                    other_idx = entry->idx;
                    if (other_idx <= idx || culled[other_idx])
                    {
                        continue;
                    }

                    dx_value = base_x - mapped_points_data[other_idx * 2];
                    dy_value = base_y - mapped_points_data[other_idx * 2 + 1];
                    if (dx_value * dx_value + dy_value * dy_value < spacing_sq)
                    {
                        culled[other_idx] = 1;
                    }
                }
            }
        }
    }

    {
        npy_intp dims[1] = {selected_count};
        selected_arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INTP);
        if (selected_arr == NULL)
        {
            goto cleanup;
        }
        memcpy(
            PyArray_DATA(selected_arr),
            selected_indices,
            (size_t)selected_count * sizeof(npy_intp));
    }

    result = (PyObject *)selected_arr;
    selected_arr = NULL;

cleanup:
    Py_XDECREF(mapped_points_arr);
    Py_XDECREF(selected_arr);
    PyMem_Free(entries);
    PyMem_Free(bucket_x);
    PyMem_Free(bucket_y);
    PyMem_Free(culled);
    PyMem_Free(selected_indices);
    return result;
}

/* Native methods exported to Python. */
static PyMethodDef module_methods[] = {
    {
        "trace_ncl_direction",
        (PyCFunction)trace_ncl_direction,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Trace one NCL-like curved-vector direction on plain NumPy arrays."),
    },
    {
        "sample_grid_field",
        (PyCFunction)sample_grid_field,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Sample one scalar field value on a regular 2D grid."),
    },
    {
        "sample_grid_field_array",
        (PyCFunction)sample_grid_field_array,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Vectorized bilinear sampling of a scalar field on a regular 2D grid."),
    },
    {
        "thin_mapped_candidates",
        (PyCFunction)thin_mapped_candidates,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Cull later nearby mapped candidate points in scan order."),
    },
    {NULL, NULL, 0, NULL},
};

/* Standard CPython module definition. */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_ncl_curly_core",
    "Native C helpers for NCL-like curved-vector tracing.",
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

/* Module init: import NumPy C-API and create the extension module. */
PyMODINIT_FUNC PyInit__ncl_curly_core(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    if (module == NULL)
    {
        return NULL;
    }
    import_array();
    return module;
}
