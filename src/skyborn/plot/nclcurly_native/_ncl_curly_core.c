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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Regular rectilinear data-grid metadata shared by the sampling helpers.
 *
 * The tracer works in "data coordinates" (the original x/y indexing space of
 * the input arrays), but many hot loops need the same origin, spacing and
 * bounds checks over and over. The *_safe values guard against a zero spacing
 * input, while the reciprocals and precomputed maxima avoid repeating divisions
 * and shape arithmetic inside the inner loops.
 */
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

/*
 * Interpolated flow state at one sample position.
 *
 * display_x/display_y are the mapped panel coordinates of the sample.
 * j00..j11 are the local Jacobian of the data->display mapping:
 *
 *     [dx_display/dx_data  dx_display/dy_data]
 *     [dy_display/dx_data  dy_display/dy_data]
 *
 * dir_x/dir_y is the unit tangent used by the display-space stepping loop, and
 * speed is the original vector magnitude in data space. Keeping both lets the
 * tracer use NCL-like spacing rules while still following the mapped geometry.
 */
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

/*
 * One candidate point projected onto a spacing-sized bucket lattice.
 *
 * The thinning pass sorts these entries by bucket index so each accepted point
 * only has to inspect nearby buckets instead of scanning the whole candidate
 * cloud.
 */
typedef struct
{
    long bx;
    long by;
    npy_intp idx;
} BucketEntry;

static int compare_double_values(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return (a > b) - (a < b);
}

/*
 * Walk backward along one display-space curve and recover the tip direction.
 *
 * Both the open-line head and the filled triangular head reuse the same tip
 * tangent estimation logic. backoff_px controls how far from the tip the
 * tangent probe starts, which mirrors the Python helpers:
 *
 * - open head  : head_length_px * 1.35
 * - filled head: head_length_px * 1.25
 */
static int resolve_curve_tip_direction(
    const double *curve,
    npy_intp n_points,
    double backoff_px,
    double *tip_x,
    double *tip_y,
    double *unit_x,
    double *unit_y)
{
    npy_intp idx;
    double remaining;
    double tail_x;
    double tail_y;
    double dir_x;
    double dir_y;
    double dir_norm;

    if (n_points < 2 || !isfinite(backoff_px))
    {
        return 0;
    }

    for (idx = 0; idx < n_points * 2; ++idx)
    {
        if (!isfinite(curve[idx]))
        {
            return 0;
        }
    }

    *tip_x = curve[(n_points - 1) * 2];
    *tip_y = curve[(n_points - 1) * 2 + 1];
    tail_x = curve[0];
    tail_y = curve[1];
    remaining = fmax(backoff_px, 0.0);

    for (idx = n_points - 1; idx > 0; --idx)
    {
        double right_x = curve[idx * 2];
        double right_y = curve[idx * 2 + 1];
        double left_x = curve[(idx - 1) * 2];
        double left_y = curve[(idx - 1) * 2 + 1];
        double seg_x = right_x - left_x;
        double seg_y = right_y - left_y;
        double seg_length = hypot(seg_x, seg_y);

        if (seg_length <= 1e-12)
        {
            continue;
        }
        if (remaining <= seg_length)
        {
            double fraction = remaining / seg_length;
            tail_x = right_x - seg_x * fraction;
            tail_y = right_y - seg_y * fraction;
            break;
        }
        remaining -= seg_length;
    }

    dir_x = *tip_x - tail_x;
    dir_y = *tip_y - tail_y;
    dir_norm = hypot(dir_x, dir_y);
    if (!(isfinite(dir_norm)) || dir_norm <= 1e-12)
    {
        return 0;
    }

    *unit_x = dir_x / dir_norm;
    *unit_y = dir_y / dir_norm;
    return 1;
}

/*
 * Build one open-arrow head directly in display space.
 *
 * curve points are already projected into the final display coordinate system.
 * The helper walks backward along the polyline to estimate the tip tangent,
 * then emits left/tip/right display vertices for the open arrow head.
 */
static int build_open_arrow_vertices_for_curve(
    const double *curve,
    npy_intp n_points,
    double head_length_px,
    double head_width_px,
    double *out_vertices)
{
    double tip_x;
    double tip_y;
    double unit_x;
    double unit_y;
    double base_x;
    double base_y;
    double normal_x;
    double normal_y;

    if (!isfinite(head_length_px) || !isfinite(head_width_px))
    {
        return 0;
    }

    if (!resolve_curve_tip_direction(
            curve,
            n_points,
            head_length_px * 1.35,
            &tip_x,
            &tip_y,
            &unit_x,
            &unit_y))
    {
        return 0;
    }

    base_x = tip_x - unit_x * head_length_px;
    base_y = tip_y - unit_y * head_length_px;
    normal_x = -unit_y;
    normal_y = unit_x;

    out_vertices[0] = base_x + normal_x * head_width_px * 0.5;
    out_vertices[1] = base_y + normal_y * head_width_px * 0.5;
    out_vertices[2] = tip_x;
    out_vertices[3] = tip_y;
    out_vertices[4] = base_x - normal_x * head_width_px * 0.5;
    out_vertices[5] = base_y - normal_y * head_width_px * 0.5;
    return 1;
}

/*
 * Build one filled triangular arrow head directly in display space.
 *
 * The returned vertex order matches the Python polygon helper:
 * tip -> left base corner -> right base corner.
 */
static int build_filled_arrow_vertices_for_curve(
    const double *curve,
    npy_intp n_points,
    double head_length_px,
    double head_width_px,
    double *out_vertices)
{
    double tip_x;
    double tip_y;
    double unit_x;
    double unit_y;
    double base_x;
    double base_y;
    double normal_x;
    double normal_y;

    if (!isfinite(head_length_px) || !isfinite(head_width_px))
    {
        return 0;
    }

    if (!resolve_curve_tip_direction(
            curve,
            n_points,
            head_length_px * 1.25,
            &tip_x,
            &tip_y,
            &unit_x,
            &unit_y))
    {
        return 0;
    }

    normal_x = -unit_y;
    normal_y = unit_x;
    base_x = tip_x - unit_x * head_length_px;
    base_y = tip_y - unit_y * head_length_px;

    out_vertices[0] = tip_x;
    out_vertices[1] = tip_y;
    out_vertices[2] = base_x + normal_x * head_width_px * 0.5;
    out_vertices[3] = base_y + normal_y * head_width_px * 0.5;
    out_vertices[4] = base_x - normal_x * head_width_px * 0.5;
    out_vertices[5] = base_y - normal_y * head_width_px * 0.5;
    return 1;
}

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

    /* Convert the requested sample position into fractional grid coordinates. */
    xi = (x - grid->x_origin) * grid->inv_dx;
    yi = (y - grid->y_origin) * grid->inv_dy;
    if (!(0.0 <= xi && xi <= grid->max_x_index && 0.0 <= yi && yi <= grid->max_y_index))
    {
        return 0;
    }

    /*
     * Identify the enclosing cell. cell_valid is defined on cells rather than
     * nodes, so rejecting the cell here guarantees that the full bilinear
     * interpolation stencil is trustworthy.
     */
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

    /* Bilinearly interpolate the vector components in data space first. */
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

    /*
     * display_data stores the mapped x/y location of every grid node. We
     * interpolate those mapped coordinates and differentiate the bilinear patch
     * to obtain a local Jacobian for converting vectors and step increments
     * between the data grid and the rendered panel.
     */
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

    /* Singular local mappings cannot define a reliable display-space tangent. */
    det = out->j00 * out->j11 - out->j01 * out->j10;
    if (!isfinite(out->display_x) || !isfinite(out->display_y) || !isfinite(det) || fabs(det) <= 1e-12)
    {
        return 0;
    }

    /* Push the data-space vector through the local Jacobian into display space. */
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

/*
 * Convert one display-space increment back into the data grid's coordinate
 * system by inverting the local Jacobian.
 *
 * The tracer chooses its direction and step length in mapped/display space so
 * the geometry matches the projection. To advance the integration state we then
 * need the corresponding delta in the original rectilinear data grid.
 */
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

/*
 * Match the Python renderer's NCL-like speed-dependent pixel step law.
 *
 * Faster local wind gets longer steps, but the quadratic scaling keeps weak
 * flow from producing overly long, noisy wiggles. The hard floor prevents the
 * branch from stalling completely in very small but still valid flow regions.
 */
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

/*
 * Bilinear sampling helper for scalar fields on the same regular grid.
 *
 * Unlike the vector tracer this can safely clamp to the last node on the outer
 * edge, so ix_next/iy_next collapse to the edge index when the sample lands on
 * the boundary.
 */
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

/*
 * Shared implementation for display-space thinning.
 *
 * The input coordinates can already be normalized viewport coordinates
 * (origin 0, scale 1), or raw display/pixel coordinates paired with a viewport
 * origin and inverse width/height. Keeping the scaling here lets callers avoid
 * building an intermediate mapped_points array when they already have display
 * coordinates.
 */
static PyObject *thin_scaled_points(
    PyArrayObject *points_arr,
    double origin_x,
    double origin_y,
    double inv_width,
    double inv_height,
    double spacing_frac)
{
    PyArrayObject *selected_arr = NULL;
    PyObject *result = NULL;
    double spacing_sq;
    double bucket_scale;
    double bucket_scale_x;
    double bucket_scale_y;
    npy_intp point_count;
    BucketEntry *entries = NULL;
    unsigned char *culled = NULL;
    npy_intp *selected_indices = NULL;
    npy_intp selected_count = 0;
    const double *points_data;
    npy_intp idx;

    point_count = PyArray_DIM(points_arr, 0);
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
    bucket_scale_x = inv_width * bucket_scale;
    bucket_scale_y = inv_height * bucket_scale;
    points_data = (const double *)PyArray_DATA(points_arr);

    entries = (BucketEntry *)PyMem_Malloc((size_t)point_count * sizeof(BucketEntry));
    culled = (unsigned char *)PyMem_Calloc((size_t)point_count, sizeof(unsigned char));
    selected_indices = (npy_intp *)PyMem_Malloc((size_t)point_count * sizeof(npy_intp));
    if (entries == NULL || culled == NULL || selected_indices == NULL)
    {
        PyErr_NoMemory();
        goto cleanup;
    }

    Py_BEGIN_ALLOW_THREADS

    for (idx = 0; idx < point_count; ++idx)
    {
        const double px = points_data[idx * 2];
        const double py = points_data[idx * 2 + 1];
        entries[idx].bx = (long)floor((px - origin_x) * bucket_scale_x);
        entries[idx].by = (long)floor((py - origin_y) * bucket_scale_y);
        entries[idx].idx = idx;
    }

    /*
     * Sorting by bucket lets lower_bound_bucket jump straight to the first
     * entry in a bucket. The final keep/cull decision still follows the
     * original candidate order because idx drives the outer loop below.
     */
    qsort(entries, (size_t)point_count, sizeof(BucketEntry), compare_bucket_entries);

    for (idx = 0; idx < point_count; ++idx)
    {
        npy_intp pos;
        const double base_x = points_data[idx * 2];
        const double base_y = points_data[idx * 2 + 1];
        const long bx = (long)floor((base_x - origin_x) * bucket_scale_x);
        const long by = (long)floor((base_y - origin_y) * bucket_scale_y);
        long search_x;

        if (culled[idx])
        {
            continue;
        }

        selected_indices[selected_count++] = idx;

        /*
         * A point closer than spacing_frac must fall in the same bucket or one
         * of the eight neighbors because bucket width equals the target spacing.
         * That is why a 3x3 neighborhood is sufficient here.
         */
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

                    dx_value = (base_x - points_data[other_idx * 2]) * inv_width;
                    dy_value = (base_y - points_data[other_idx * 2 + 1]) * inv_height;
                    /* Cull any later candidate that violates the spacing rule. */
                    if (dx_value * dx_value + dy_value * dy_value < spacing_sq)
                    {
                        culled[other_idx] = 1;
                    }
                }
            }
        }
    }

    Py_END_ALLOW_THREADS

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
    Py_XDECREF(selected_arr);
    PyMem_Free(entries);
    PyMem_Free(culled);
    PyMem_Free(selected_indices);
    return result;
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

    /*
     * Find the earliest intersection between the proposed segment and the four
     * viewport sides. Because the segment starts inside the panel, the smallest
     * valid factor in [0, 1] is the amount we can keep.
     */
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
    PyArrayObject *full_display_output = NULL;
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
    int return_display = 0;
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
        "return_display",
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
    double *display_output_data = NULL;
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
            "OOOOdddddddddddddd|ii:trace_ncl_direction",
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
            &max_steps,
            &return_display))
    {
        return NULL;
    }

    /* Clamp the Python-provided step budget so preallocation stays bounded. */
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

    /* Validate the plain NumPy contract expected by the native kernel. */
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

    /*
     * Materialize a compact grid descriptor once so the stepping loop only
     * passes around a pointer instead of recomputing shape- and spacing-derived
     * constants every iteration.
     */
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

    /*
     * Preallocate the worst-case polyline length. If the branch terminates
     * early, the result is trimmed after tracing.
     */
    output_dims[0] = max_steps + 1;
    output_dims[1] = 2;
    full_output = (PyArrayObject *)PyArray_SimpleNew(2, output_dims, NPY_DOUBLE);
    if (full_output == NULL)
    {
        goto cleanup;
    }
    if (return_display)
    {
        full_display_output = (PyArrayObject *)PyArray_SimpleNew(2, output_dims, NPY_DOUBLE);
        if (full_display_output == NULL)
        {
            goto cleanup;
        }
        display_output_data = (double *)PyArray_DATA(full_display_output);
    }

    output_data = (double *)PyArray_DATA(full_output);
    output_data[0] = start_x;
    output_data[1] = start_y;
    current_data_x = start_x;
    current_data_y = start_y;
    current_display_x = initial_state.display_x;
    current_display_y = initial_state.display_y;
    if (return_display)
    {
        display_output_data[0] = current_display_x;
        display_output_data[1] = current_display_y;
    }

    /* The numeric stepping loop is CPU-bound and does not need the GIL. */
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
        if (return_display)
        {
            display_output_data[point_count * 2] = candidate_display_x;
            display_output_data[point_count * 2 + 1] = candidate_display_y;
        }
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
        if (return_display)
        {
            result = Py_BuildValue("NN", (PyObject *)full_output, (PyObject *)full_display_output);
            full_output = NULL;
            full_display_output = NULL;
        }
        else
        {
            result = (PyObject *)full_output;
            full_output = NULL;
        }
        goto cleanup;
    }

    {
        npy_intp trimmed_dims[2];
        PyArrayObject *trimmed_output = NULL;
        PyArrayObject *trimmed_display_output = NULL;

        /* Shrink the overallocated buffer to the exact accepted vertex count. */
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
        if (return_display)
        {
            trimmed_display_output = (PyArrayObject *)PyArray_SimpleNew(2, trimmed_dims, NPY_DOUBLE);
            if (trimmed_display_output == NULL)
            {
                Py_DECREF(trimmed_output);
                goto cleanup;
            }
            memcpy(
                PyArray_DATA(trimmed_display_output),
                PyArray_DATA(full_display_output),
                (size_t)(point_count * 2) * sizeof(double));
            result = Py_BuildValue("NN", (PyObject *)trimmed_output, (PyObject *)trimmed_display_output);
        }
        else
        {
            result = (PyObject *)trimmed_output;
        }
    }

cleanup:
    Py_XDECREF(u_arr);
    Py_XDECREF(v_arr);
    Py_XDECREF(display_arr);
    Py_XDECREF(cell_valid_arr);
    Py_XDECREF(full_output);
    Py_XDECREF(full_display_output);
    return result;
}

/*
 * Compute display-grid cell validity in C.
 *
 * This mirrors _NCLDisplaySampler._compute_cell_valid: a cell is valid only when
 * all four display nodes are finite and none of the four cell edges crosses a
 * projection seam according to the same robust jump threshold.
 */
static PyObject *compute_display_cell_valid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_grid_obj;
    PyArrayObject *display_arr = NULL;
    PyArrayObject *cell_valid_arr = NULL;
    PyObject *result = NULL;
    static char *kwlist[] = {
        "display_grid",
        NULL,
    };
    npy_intp ny;
    npy_intp nx;
    npy_intp n_edges;
    npy_intp edge_count = 0;
    double *edges = NULL;
    double jump_limit = 0.0;
    int has_jump_limit = 0;
    const double *display;
    npy_bool *cell_valid;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "O:compute_display_cell_valid",
            kwlist,
            &display_grid_obj))
    {
        return NULL;
    }

    display_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_grid_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (display_arr == NULL)
    {
        goto cleanup;
    }
    if (PyArray_NDIM(display_arr) != 3 ||
        PyArray_DIM(display_arr, 2) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_grid must have shape (ny, nx, 2)");
        goto cleanup;
    }

    ny = PyArray_DIM(display_arr, 0);
    nx = PyArray_DIM(display_arr, 1);
    if (ny < 2 || nx < 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_grid must contain at least one cell");
        goto cleanup;
    }

    {
        npy_intp output_dims[2] = {ny - 1, nx - 1};
        cell_valid_arr = (PyArrayObject *)PyArray_SimpleNew(2, output_dims, NPY_BOOL);
        if (cell_valid_arr == NULL)
        {
            goto cleanup;
        }
    }

    display = (const double *)PyArray_DATA(display_arr);
    cell_valid = (npy_bool *)PyArray_DATA(cell_valid_arr);
    n_edges = ny * (nx - 1) + (ny - 1) * nx;
    edges = (double *)malloc((size_t)n_edges * sizeof(double));
    if (edges == NULL)
    {
        PyErr_NoMemory();
        goto cleanup;
    }

    for (npy_intp iy = 0; iy < ny; ++iy)
    {
        for (npy_intp ix = 0; ix < nx - 1; ++ix)
        {
            npy_intp left_idx = (iy * nx + ix) * 2;
            npy_intp right_idx = (iy * nx + ix + 1) * 2;
            double length = hypot(
                display[right_idx] - display[left_idx],
                display[right_idx + 1] - display[left_idx + 1]);
            if (isfinite(length) && length > 1e-9)
            {
                edges[edge_count++] = length;
            }
        }
    }
    for (npy_intp iy = 0; iy < ny - 1; ++iy)
    {
        for (npy_intp ix = 0; ix < nx; ++ix)
        {
            npy_intp top_idx = (iy * nx + ix) * 2;
            npy_intp bottom_idx = ((iy + 1) * nx + ix) * 2;
            double length = hypot(
                display[bottom_idx] - display[top_idx],
                display[bottom_idx + 1] - display[top_idx + 1]);
            if (isfinite(length) && length > 1e-9)
            {
                edges[edge_count++] = length;
            }
        }
    }

    if (edge_count > 0)
    {
        npy_intp lower_count = edge_count / 2;
        if (lower_count < 1)
        {
            lower_count = 1;
        }
        qsort(edges, (size_t)edge_count, sizeof(double), compare_double_values);
        if (lower_count % 2 == 1)
        {
            jump_limit = edges[lower_count / 2];
        }
        else
        {
            jump_limit = 0.5 * (edges[lower_count / 2 - 1] + edges[lower_count / 2]);
        }
        jump_limit = fmax(jump_limit * 12.0, 1e-6);
        has_jump_limit = 1;
    }

    Py_BEGIN_ALLOW_THREADS
    for (npy_intp iy = 0; iy < ny - 1; ++iy)
    {
        for (npy_intp ix = 0; ix < nx - 1; ++ix)
        {
            npy_intp cell_idx = iy * (nx - 1) + ix;
            npy_intp idx00 = (iy * nx + ix) * 2;
            npy_intp idx01 = (iy * nx + ix + 1) * 2;
            npy_intp idx10 = ((iy + 1) * nx + ix) * 2;
            npy_intp idx11 = ((iy + 1) * nx + ix + 1) * 2;
            double top_edge;
            double bottom_edge;
            double left_edge;
            double right_edge;
            double max_edge;
            int finite_cell =
                isfinite(display[idx00]) &&
                isfinite(display[idx00 + 1]) &&
                isfinite(display[idx01]) &&
                isfinite(display[idx01 + 1]) &&
                isfinite(display[idx10]) &&
                isfinite(display[idx10 + 1]) &&
                isfinite(display[idx11]) &&
                isfinite(display[idx11 + 1]);

            if (!finite_cell)
            {
                cell_valid[cell_idx] = NPY_FALSE;
                continue;
            }
            if (!has_jump_limit)
            {
                cell_valid[cell_idx] = NPY_TRUE;
                continue;
            }

            top_edge = hypot(display[idx01] - display[idx00], display[idx01 + 1] - display[idx00 + 1]);
            bottom_edge = hypot(display[idx11] - display[idx10], display[idx11 + 1] - display[idx10 + 1]);
            left_edge = hypot(display[idx10] - display[idx00], display[idx10 + 1] - display[idx00 + 1]);
            right_edge = hypot(display[idx11] - display[idx01], display[idx11 + 1] - display[idx01 + 1]);
            max_edge = fmax(fmax(top_edge, bottom_edge), fmax(left_edge, right_edge));
            cell_valid[cell_idx] = (isfinite(max_edge) && max_edge <= jump_limit) ? NPY_TRUE : NPY_FALSE;
        }
    }
    Py_END_ALLOW_THREADS

    result = (PyObject *)cell_valid_arr;
    cell_valid_arr = NULL;

cleanup:
    free(edges);
    Py_XDECREF(display_arr);
    Py_XDECREF(cell_valid_arr);
    return result;
}

/*
 * Bilinearly sample display-grid positions and optional local Jacobians for a
 * batch of data-space points. Invalid rows remain NaN and false in the validity
 * mask, matching the previous Python vectorized implementation.
 */
static PyObject *sample_display_grid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_grid_obj;
    PyObject *cell_valid_obj;
    PyObject *points_obj;
    PyArrayObject *display_arr = NULL;
    PyArrayObject *cell_valid_arr = NULL;
    PyArrayObject *points_arr = NULL;
    PyArrayObject *display_points_arr = NULL;
    PyArrayObject *jacobians_arr = NULL;
    PyArrayObject *valid_arr = NULL;
    PyObject *result = NULL;
    double x_origin;
    double y_origin;
    double dx;
    double dy;
    int include_jacobian = 0;
    GridInfo grid;
    npy_intp point_count;
    static char *kwlist[] = {
        "display_grid",
        "cell_valid",
        "x_origin",
        "y_origin",
        "dx",
        "dy",
        "points",
        "include_jacobian",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOddddO|p:sample_display_grid",
            kwlist,
            &display_grid_obj,
            &cell_valid_obj,
            &x_origin,
            &y_origin,
            &dx,
            &dy,
            &points_obj,
            &include_jacobian))
    {
        return NULL;
    }

    display_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_grid_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    cell_valid_arr = (PyArrayObject *)PyArray_FROM_OTF(
        cell_valid_obj, NPY_BOOL, NPY_ARRAY_IN_ARRAY);
    points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (display_arr == NULL || cell_valid_arr == NULL || points_arr == NULL)
    {
        goto cleanup;
    }

    if (PyArray_NDIM(display_arr) != 3 ||
        PyArray_DIM(display_arr, 2) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_grid must have shape (ny, nx, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(cell_valid_arr) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "cell_valid must be a 2D boolean array");
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

    grid.ny = PyArray_DIM(display_arr, 0);
    grid.nx = PyArray_DIM(display_arr, 1);
    if (PyArray_DIM(cell_valid_arr, 0) != grid.ny - 1 ||
        PyArray_DIM(cell_valid_arr, 1) != grid.nx - 1)
    {
        PyErr_SetString(PyExc_ValueError, "cell_valid must have shape (ny-1, nx-1)");
        goto cleanup;
    }

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

    {
        npy_intp display_dims[2] = {point_count, 2};
        npy_intp valid_dims[1] = {point_count};
        display_points_arr = (PyArrayObject *)PyArray_SimpleNew(2, display_dims, NPY_DOUBLE);
        valid_arr = (PyArrayObject *)PyArray_SimpleNew(1, valid_dims, NPY_BOOL);
        if (display_points_arr == NULL || valid_arr == NULL)
        {
            goto cleanup;
        }
        if (include_jacobian)
        {
            npy_intp jacobian_dims[3] = {point_count, 2, 2};
            jacobians_arr = (PyArrayObject *)PyArray_SimpleNew(3, jacobian_dims, NPY_DOUBLE);
            if (jacobians_arr == NULL)
            {
                goto cleanup;
            }
        }
    }

    {
        const double *display = (const double *)PyArray_DATA(display_arr);
        const npy_bool *cell_valid = (const npy_bool *)PyArray_DATA(cell_valid_arr);
        const double *points = (const double *)PyArray_DATA(points_arr);
        double *display_points = (double *)PyArray_DATA(display_points_arr);
        double *jacobians = include_jacobian ? (double *)PyArray_DATA(jacobians_arr) : NULL;
        npy_bool *valid = (npy_bool *)PyArray_DATA(valid_arr);

        Py_BEGIN_ALLOW_THREADS
        for (npy_intp idx = 0; idx < point_count; ++idx)
        {
            display_points[idx * 2] = Py_NAN;
            display_points[idx * 2 + 1] = Py_NAN;
            valid[idx] = NPY_FALSE;
            if (include_jacobian)
            {
                for (int jidx = 0; jidx < 4; ++jidx)
                {
                    jacobians[idx * 4 + jidx] = Py_NAN;
                }
            }
        }

        if (grid.nx >= 2 && grid.ny >= 2)
        {
            for (npy_intp idx = 0; idx < point_count; ++idx)
            {
                double point_x = points[idx * 2];
                double point_y = points[idx * 2 + 1];
                double xi = (point_x - grid.x_origin) * grid.inv_dx;
                double yi = (point_y - grid.y_origin) * grid.inv_dy;
                npy_intp ix;
                npy_intp iy;
                npy_intp idx00;
                npy_intp idx01;
                npy_intp idx10;
                npy_intp idx11;
                double sx;
                double sy;
                double one_minus_sx;
                double one_minus_sy;
                double p00x;
                double p00y;
                double p01x;
                double p01y;
                double p10x;
                double p10y;
                double p11x;
                double p11y;
                double display_x;
                double display_y;
                double dfdx_x;
                double dfdx_y;
                double dfdy_x;
                double dfdy_y;
                double det;

                if (!(0.0 <= xi && xi <= grid.max_x_index && 0.0 <= yi && yi <= grid.max_y_index))
                {
                    continue;
                }
                ix = (npy_intp)floor(fmin(fmax(xi, 0.0), grid.max_cell_x));
                iy = (npy_intp)floor(fmin(fmax(yi, 0.0), grid.max_cell_y));
                if (!cell_valid[iy * (grid.nx - 1) + ix])
                {
                    continue;
                }

                sx = xi - (double)ix;
                sy = yi - (double)iy;
                one_minus_sx = 1.0 - sx;
                one_minus_sy = 1.0 - sy;
                idx00 = (iy * grid.nx + ix) * 2;
                idx01 = (iy * grid.nx + ix + 1) * 2;
                idx10 = ((iy + 1) * grid.nx + ix) * 2;
                idx11 = ((iy + 1) * grid.nx + ix + 1) * 2;

                p00x = display[idx00];
                p00y = display[idx00 + 1];
                p01x = display[idx01];
                p01y = display[idx01 + 1];
                p10x = display[idx10];
                p10y = display[idx10 + 1];
                p11x = display[idx11];
                p11y = display[idx11 + 1];

                display_x = p00x * one_minus_sx * one_minus_sy + p01x * sx * one_minus_sy + p10x * one_minus_sx * sy + p11x * sx * sy;
                display_y = p00y * one_minus_sx * one_minus_sy + p01y * sx * one_minus_sy + p10y * one_minus_sx * sy + p11y * sx * sy;
                if (!isfinite(display_x) || !isfinite(display_y))
                {
                    continue;
                }
                display_points[idx * 2] = display_x;
                display_points[idx * 2 + 1] = display_y;

                if (!include_jacobian)
                {
                    valid[idx] = NPY_TRUE;
                    continue;
                }

                dfdx_x = ((p01x - p00x) * one_minus_sy + (p11x - p10x) * sy) / grid.dx_safe;
                dfdx_y = ((p01y - p00y) * one_minus_sy + (p11y - p10y) * sy) / grid.dx_safe;
                dfdy_x = ((p10x - p00x) * one_minus_sx + (p11x - p01x) * sx) / grid.dy_safe;
                dfdy_y = ((p10y - p00y) * one_minus_sx + (p11y - p01y) * sx) / grid.dy_safe;
                det = dfdx_x * dfdy_y - dfdy_x * dfdx_y;
                if (!isfinite(dfdx_x) || !isfinite(dfdx_y) ||
                    !isfinite(dfdy_x) || !isfinite(dfdy_y) ||
                    !isfinite(det) || fabs(det) <= 1e-12)
                {
                    display_points[idx * 2] = Py_NAN;
                    display_points[idx * 2 + 1] = Py_NAN;
                    continue;
                }

                jacobians[idx * 4] = dfdx_x;
                jacobians[idx * 4 + 1] = dfdy_x;
                jacobians[idx * 4 + 2] = dfdx_y;
                jacobians[idx * 4 + 3] = dfdy_y;
                valid[idx] = NPY_TRUE;
            }
        }
        Py_END_ALLOW_THREADS
    }

    if (include_jacobian)
    {
        result = Py_BuildValue("NNN", (PyObject *)display_points_arr, (PyObject *)jacobians_arr, (PyObject *)valid_arr);
        display_points_arr = NULL;
        jacobians_arr = NULL;
        valid_arr = NULL;
    }
    else
    {
        result = Py_BuildValue("NON", (PyObject *)display_points_arr, Py_None, (PyObject *)valid_arr);
        display_points_arr = NULL;
        valid_arr = NULL;
    }

cleanup:
    Py_XDECREF(display_arr);
    Py_XDECREF(cell_valid_arr);
    Py_XDECREF(points_arr);
    Py_XDECREF(display_points_arr);
    Py_XDECREF(jacobians_arr);
    Py_XDECREF(valid_arr);
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

/*
 * Python wrapper for vectorized scalar sampling at many points.
 *
 * Invalid or out-of-bounds points are left as NaN so the Python caller can
 * preserve positional alignment with its original request array.
 */
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
    PyObject *result = NULL;
    double spacing_frac;
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

    result = thin_scaled_points(
        mapped_points_arr,
        0.0,
        0.0,
        1.0,
        1.0,
        spacing_frac);

cleanup:
    Py_XDECREF(mapped_points_arr);
    return result;
}

/*
 * Thin display/pixel-space candidate centers using the same viewport mapping
 * as _map_ncl_display_points_to_viewport, but without allocating mapped_points
 * on the Python side.
 */
static PyObject *thin_display_candidates(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_points_obj;
    PyArrayObject *display_points_arr = NULL;
    PyObject *result = NULL;
    double viewport_x0;
    double viewport_y0;
    double viewport_width;
    double viewport_height;
    double spacing_frac;
    static char *kwlist[] = {
        "display_points",
        "viewport_x0",
        "viewport_y0",
        "viewport_width",
        "viewport_height",
        "spacing_frac",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "Oddddd:thin_display_candidates",
            kwlist,
            &display_points_obj,
            &viewport_x0,
            &viewport_y0,
            &viewport_width,
            &viewport_height,
            &spacing_frac))
    {
        return NULL;
    }

    display_points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (display_points_arr == NULL)
    {
        goto cleanup;
    }
    if (PyArray_NDIM(display_points_arr) != 2 || PyArray_DIM(display_points_arr, 1) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_points must have shape (N, 2)");
        goto cleanup;
    }

    viewport_width = viewport_width > 1.0 ? viewport_width : 1.0;
    viewport_height = viewport_height > 1.0 ? viewport_height : 1.0;
    result = thin_scaled_points(
        display_points_arr,
        viewport_x0,
        viewport_y0,
        1.0 / viewport_width,
        1.0 / viewport_height,
        spacing_frac);

cleanup:
    Py_XDECREF(display_points_arr);
    return result;
}

/*
 * Expand selected scatter cells into interior stipple candidates.
 *
 * Python still owns coordinate extraction and Matplotlib/Cartopy transforms.
 * This native loop only consumes already-built data-space corners and
 * viewport-normalized display corners, then emits the same candidate lattice as
 * the previous Python implementation without per-cell np.linspace/meshgrid
 * allocations.
 */
static PyObject *generate_cell_candidates(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *corners_obj;
    PyObject *mapped_corners_obj;
    PyObject *source_points_obj;
    PyArrayObject *corners_arr = NULL;
    PyArrayObject *mapped_corners_arr = NULL;
    PyArrayObject *source_points_arr = NULL;
    PyArrayObject *candidate_points_arr = NULL;
    PyArrayObject *source_positions_arr = NULL;
    PyObject *result = NULL;
    double spacing_frac;
    double spacing;
    npy_intp n_cells;
    npy_intp total_count = 0;
    npy_intp output_dims_points[2];
    npy_intp output_dims_positions[1];
    const double *corners;
    const double *mapped_corners;
    const double *source_points;
    double *candidate_points;
    npy_intp *source_positions;
    static char *kwlist[] = {
        "corners",
        "mapped_corners",
        "source_points",
        "spacing_frac",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOd:generate_cell_candidates",
            kwlist,
            &corners_obj,
            &mapped_corners_obj,
            &source_points_obj,
            &spacing_frac))
    {
        return NULL;
    }

    corners_arr = (PyArrayObject *)PyArray_FROM_OTF(
        corners_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    mapped_corners_arr = (PyArrayObject *)PyArray_FROM_OTF(
        mapped_corners_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    source_points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        source_points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (corners_arr == NULL || mapped_corners_arr == NULL || source_points_arr == NULL)
    {
        goto cleanup;
    }

    if (PyArray_NDIM(corners_arr) != 3 ||
        PyArray_DIM(corners_arr, 1) != 4 ||
        PyArray_DIM(corners_arr, 2) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "corners must have shape (N, 4, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(mapped_corners_arr) != 3 ||
        PyArray_DIM(mapped_corners_arr, 1) != 4 ||
        PyArray_DIM(mapped_corners_arr, 2) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "mapped_corners must have shape (N, 4, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(source_points_arr) != 2 || PyArray_DIM(source_points_arr, 1) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "source_points must have shape (N, 2)");
        goto cleanup;
    }

    n_cells = PyArray_DIM(corners_arr, 0);
    if (PyArray_DIM(mapped_corners_arr, 0) != n_cells ||
        PyArray_DIM(source_points_arr, 0) != n_cells)
    {
        PyErr_SetString(PyExc_ValueError, "corners, mapped_corners, and source_points must share N");
        goto cleanup;
    }

    spacing = isfinite(spacing_frac) ? spacing_frac : 0.03;
    if (spacing < 1e-6)
    {
        spacing = 1e-6;
    }

    corners = (const double *)PyArray_DATA(corners_arr);
    mapped_corners = (const double *)PyArray_DATA(mapped_corners_arr);
    source_points = (const double *)PyArray_DATA(source_points_arr);

    /* First pass: count exactly how many candidate rows are needed. */
    for (npy_intp cell = 0; cell < n_cells; ++cell)
    {
        const double *mapped = mapped_corners + cell * 8;
        int valid = 1;
        double min_x = mapped[0];
        double max_x = mapped[0];
        double min_y = mapped[1];
        double max_y = mapped[1];
        int sub_x;
        int sub_y;

        for (int corner = 0; corner < 4; ++corner)
        {
            double x = mapped[corner * 2];
            double y = mapped[corner * 2 + 1];
            if (!isfinite(x) || !isfinite(y))
            {
                valid = 0;
                break;
            }
            if (x < min_x)
            {
                min_x = x;
            }
            if (x > max_x)
            {
                max_x = x;
            }
            if (y < min_y)
            {
                min_y = y;
            }
            if (y > max_y)
            {
                max_y = y;
            }
        }

        if (!valid)
        {
            total_count += 1;
            continue;
        }

        sub_x = (int)ceil((max_x - min_x) / spacing);
        sub_y = (int)ceil((max_y - min_y) / spacing);
        if (sub_x < 1)
        {
            sub_x = 1;
        }
        else if (sub_x > 12)
        {
            sub_x = 12;
        }
        if (sub_y < 1)
        {
            sub_y = 1;
        }
        else if (sub_y > 12)
        {
            sub_y = 12;
        }
        total_count += (npy_intp)sub_x * (npy_intp)sub_y;
    }

    output_dims_points[0] = total_count;
    output_dims_points[1] = 2;
    output_dims_positions[0] = total_count;
    candidate_points_arr = (PyArrayObject *)PyArray_SimpleNew(2, output_dims_points, NPY_DOUBLE);
    source_positions_arr = (PyArrayObject *)PyArray_SimpleNew(1, output_dims_positions, NPY_INTP);
    if (candidate_points_arr == NULL || source_positions_arr == NULL)
    {
        goto cleanup;
    }

    candidate_points = (double *)PyArray_DATA(candidate_points_arr);
    source_positions = (npy_intp *)PyArray_DATA(source_positions_arr);

    /* Second pass: fill the output arrays in the same scan order as Python. */
    {
        npy_intp out_pos = 0;
        for (npy_intp cell = 0; cell < n_cells; ++cell)
        {
            const double *corner_data = corners + cell * 8;
            const double *mapped = mapped_corners + cell * 8;
            int valid = 1;
            double min_x = mapped[0];
            double max_x = mapped[0];
            double min_y = mapped[1];
            double max_y = mapped[1];
            int sub_x;
            int sub_y;

            for (int corner = 0; corner < 4; ++corner)
            {
                double x = mapped[corner * 2];
                double y = mapped[corner * 2 + 1];
                if (!isfinite(x) || !isfinite(y))
                {
                    valid = 0;
                    break;
                }
                if (x < min_x)
                {
                    min_x = x;
                }
                if (x > max_x)
                {
                    max_x = x;
                }
                if (y < min_y)
                {
                    min_y = y;
                }
                if (y > max_y)
                {
                    max_y = y;
                }
            }

            if (!valid)
            {
                candidate_points[out_pos * 2] = source_points[cell * 2];
                candidate_points[out_pos * 2 + 1] = source_points[cell * 2 + 1];
                source_positions[out_pos] = cell;
                out_pos += 1;
                continue;
            }

            sub_x = (int)ceil((max_x - min_x) / spacing);
            sub_y = (int)ceil((max_y - min_y) / spacing);
            if (sub_x < 1)
            {
                sub_x = 1;
            }
            else if (sub_x > 12)
            {
                sub_x = 12;
            }
            if (sub_y < 1)
            {
                sub_y = 1;
            }
            else if (sub_y > 12)
            {
                sub_y = 12;
            }

            for (int row = 0; row < sub_y; ++row)
            {
                double vv = ((double)row + 0.5) / (double)sub_y;
                double one_minus_v = 1.0 - vv;
                for (int col = 0; col < sub_x; ++col)
                {
                    double uu = ((double)col + 0.5) / (double)sub_x;
                    double one_minus_u = 1.0 - uu;
                    double weight00 = one_minus_u * one_minus_v;
                    double weight10 = uu * one_minus_v;
                    double weight11 = uu * vv;
                    double weight01 = one_minus_u * vv;

                    candidate_points[out_pos * 2] =
                        weight00 * corner_data[0] +
                        weight10 * corner_data[2] +
                        weight11 * corner_data[4] +
                        weight01 * corner_data[6];
                    candidate_points[out_pos * 2 + 1] =
                        weight00 * corner_data[1] +
                        weight10 * corner_data[3] +
                        weight11 * corner_data[5] +
                        weight01 * corner_data[7];
                    source_positions[out_pos] = cell;
                    out_pos += 1;
                }
            }
        }
    }

    result = Py_BuildValue("NN", (PyObject *)candidate_points_arr, (PyObject *)source_positions_arr);
    candidate_points_arr = NULL;
    source_positions_arr = NULL;

cleanup:
    Py_XDECREF(corners_arr);
    Py_XDECREF(mapped_corners_arr);
    Py_XDECREF(source_points_arr);
    Py_XDECREF(candidate_points_arr);
    Py_XDECREF(source_positions_arr);
    return result;
}

/*
 * Validate display-space curly-vector geometry.
 *
 * This mirrors skyborn.plot._core.geometry._evaluate_ncl_display_curve after the
 * Python transform has already been applied. Keeping the thresholds identical is
 * important because this is a visual-quality guard, not a new algorithm.
 */
static PyObject *validate_display_curve(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_curve_obj;
    PyArrayObject *curve_arr = NULL;
    double viewport_width = 0.0;
    double viewport_height = 0.0;
    static char *kwlist[] = {
        "display_curve",
        "viewport_width",
        "viewport_height",
        NULL,
    };
    npy_intp n_points;
    const double *curve;
    double *dir_x = NULL;
    double *dir_y = NULL;
    npy_intp valid_count = 0;
    double arc_length = 0.0;
    double chord_length;
    double tortuosity;
    double max_turn = 0.0;
    double total_turn = 0.0;
    int previous_turn_sign = 0;
    int sign_changes = 0;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "O|dd:validate_display_curve",
            kwlist,
            &display_curve_obj,
            &viewport_width,
            &viewport_height))
    {
        return NULL;
    }

    curve_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_curve_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (curve_arr == NULL)
    {
        return NULL;
    }
    if (PyArray_NDIM(curve_arr) != 2 || PyArray_DIM(curve_arr, 1) != 2)
    {
        Py_DECREF(curve_arr);
        PyErr_SetString(PyExc_ValueError, "display_curve must have shape (N, 2)");
        return NULL;
    }

    n_points = PyArray_DIM(curve_arr, 0);
    if (n_points < 2)
    {
        Py_DECREF(curve_arr);
        Py_RETURN_FALSE;
    }

    curve = (const double *)PyArray_DATA(curve_arr);
    for (npy_intp idx = 0; idx < n_points * 2; ++idx)
    {
        if (!isfinite(curve[idx]))
        {
            Py_DECREF(curve_arr);
            Py_RETURN_FALSE;
        }
    }

    dir_x = (double *)calloc((size_t)(n_points - 1), sizeof(double));
    dir_y = (double *)calloc((size_t)(n_points - 1), sizeof(double));
    if (dir_x == NULL || dir_y == NULL)
    {
        free(dir_x);
        free(dir_y);
        Py_DECREF(curve_arr);
        return PyErr_NoMemory();
    }

    {
        double viewport_diag = hypot(viewport_width, viewport_height);
        double jump_limit = 0.0;
        int use_jump_limit = isfinite(viewport_diag) && viewport_diag > 0.0;
        if (use_jump_limit)
        {
            jump_limit = fmax(viewport_diag * 0.35, 1e-6);
        }

        for (npy_intp idx = 0; idx < n_points - 1; ++idx)
        {
            double dx = curve[(idx + 1) * 2] - curve[idx * 2];
            double dy = curve[(idx + 1) * 2 + 1] - curve[idx * 2 + 1];
            double length = hypot(dx, dy);
            if (use_jump_limit && length > jump_limit)
            {
                free(dir_x);
                free(dir_y);
                Py_DECREF(curve_arr);
                Py_RETURN_FALSE;
            }
            if (length > 1e-6)
            {
                dir_x[valid_count] = dx / length;
                dir_y[valid_count] = dy / length;
                arc_length += length;
                valid_count += 1;
            }
        }
    }

    if (valid_count < 1)
    {
        free(dir_x);
        free(dir_y);
        Py_DECREF(curve_arr);
        Py_RETURN_FALSE;
    }

    chord_length = hypot(
        curve[(n_points - 1) * 2] - curve[0],
        curve[(n_points - 1) * 2 + 1] - curve[1]);
    if (arc_length <= 1e-6 || chord_length <= 1e-6)
    {
        free(dir_x);
        free(dir_y);
        Py_DECREF(curve_arr);
        Py_RETURN_FALSE;
    }

    if (valid_count >= 2)
    {
        for (npy_intp idx = 0; idx < valid_count - 1; ++idx)
        {
            double dot = dir_x[idx] * dir_x[idx + 1] + dir_y[idx] * dir_y[idx + 1];
            double angle;
            double cross;
            int current_turn_sign = 0;

            if (dot < -1.0)
            {
                dot = -1.0;
            }
            else if (dot > 1.0)
            {
                dot = 1.0;
            }
            angle = acos(dot) * 180.0 / M_PI;
            if (angle > max_turn)
            {
                max_turn = angle;
            }
            total_turn += angle;

            cross = dir_x[idx] * dir_y[idx + 1] - dir_y[idx] * dir_x[idx + 1];
            if (fabs(cross) > 1e-6)
            {
                current_turn_sign = cross > 0.0 ? 1 : -1;
                if (previous_turn_sign != 0 && current_turn_sign * previous_turn_sign < 0)
                {
                    sign_changes += 1;
                }
                previous_turn_sign = current_turn_sign;
            }
        }
    }

    tortuosity = arc_length / chord_length;
    free(dir_x);
    free(dir_y);
    Py_DECREF(curve_arr);

    if (tortuosity > 1.28)
    {
        Py_RETURN_FALSE;
    }
    if (max_turn > 52.0)
    {
        Py_RETURN_FALSE;
    }
    if (total_turn > 120.0)
    {
        Py_RETURN_FALSE;
    }
    if (sign_changes > 0 && total_turn > 70.0)
    {
        Py_RETURN_FALSE;
    }
    Py_RETURN_TRUE;
}

/*
 * Build compact open-arrow display segments directly in native code.
 *
 * Output layout:
 * - segments: (K, 2, 2), where each valid curve contributes two segments
 * - source_positions: (K,), repeating the source curve index for both segments
 */
static PyObject *build_open_arrow_segments(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_points_obj;
    PyObject *curve_offsets_obj;
    PyObject *head_lengths_obj;
    PyObject *head_widths_obj;
    PyArrayObject *display_points_arr = NULL;
    PyArrayObject *curve_offsets_arr = NULL;
    PyArrayObject *head_lengths_arr = NULL;
    PyArrayObject *head_widths_arr = NULL;
    PyArrayObject *segments_arr = NULL;
    PyArrayObject *source_positions_arr = NULL;
    PyObject *result = NULL;
    npy_intp n_curves;
    npy_intp n_points;
    npy_intp valid_curves = 0;
    npy_intp output_dims_segments[3];
    npy_intp output_dims_positions[1];
    const double *display_points;
    const npy_intp *curve_offsets;
    const double *head_lengths;
    const double *head_widths;
    double *segments;
    npy_intp *source_positions;
    static char *kwlist[] = {
        "display_points",
        "curve_offsets",
        "head_lengths_px",
        "head_widths_px",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOO:build_open_arrow_segments",
            kwlist,
            &display_points_obj,
            &curve_offsets_obj,
            &head_lengths_obj,
            &head_widths_obj))
    {
        return NULL;
    }

    display_points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    curve_offsets_arr = (PyArrayObject *)PyArray_FROM_OTF(
        curve_offsets_obj, NPY_INTP, NPY_ARRAY_IN_ARRAY);
    head_lengths_arr = (PyArrayObject *)PyArray_FROM_OTF(
        head_lengths_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    head_widths_arr = (PyArrayObject *)PyArray_FROM_OTF(
        head_widths_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (display_points_arr == NULL || curve_offsets_arr == NULL ||
        head_lengths_arr == NULL || head_widths_arr == NULL)
    {
        goto cleanup;
    }

    if (PyArray_NDIM(display_points_arr) != 2 || PyArray_DIM(display_points_arr, 1) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_points must have shape (M, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(curve_offsets_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "curve_offsets must be a 1D array");
        goto cleanup;
    }
    if (PyArray_NDIM(head_lengths_arr) != 1 || PyArray_NDIM(head_widths_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "head_lengths_px and head_widths_px must be 1D arrays");
        goto cleanup;
    }

    n_curves = PyArray_DIM(head_lengths_arr, 0);
    if (PyArray_DIM(head_widths_arr, 0) != n_curves)
    {
        PyErr_SetString(PyExc_ValueError, "head_lengths_px and head_widths_px must have the same length");
        goto cleanup;
    }
    if (PyArray_DIM(curve_offsets_arr, 0) != n_curves + 1)
    {
        PyErr_SetString(PyExc_ValueError, "curve_offsets must have length len(head_lengths_px) + 1");
        goto cleanup;
    }

    display_points = (const double *)PyArray_DATA(display_points_arr);
    curve_offsets = (const npy_intp *)PyArray_DATA(curve_offsets_arr);
    head_lengths = (const double *)PyArray_DATA(head_lengths_arr);
    head_widths = (const double *)PyArray_DATA(head_widths_arr);
    n_points = PyArray_DIM(display_points_arr, 0);

    for (npy_intp idx = 0; idx < n_curves; ++idx)
    {
        if (curve_offsets[idx] < 0 || curve_offsets[idx] > curve_offsets[idx + 1] ||
            curve_offsets[idx + 1] > n_points)
        {
            PyErr_SetString(PyExc_ValueError, "curve_offsets must be monotonic and stay within display_points");
            goto cleanup;
        }
    }

    for (npy_intp curve_idx = 0; curve_idx < n_curves; ++curve_idx)
    {
        npy_intp start = curve_offsets[curve_idx];
        npy_intp end = curve_offsets[curve_idx + 1];
        double vertices[6];
        if (build_open_arrow_vertices_for_curve(
                display_points + start * 2,
                end - start,
                head_lengths[curve_idx],
                head_widths[curve_idx],
                vertices))
        {
            valid_curves += 1;
        }
    }

    output_dims_segments[0] = valid_curves * 2;
    output_dims_segments[1] = 2;
    output_dims_segments[2] = 2;
    output_dims_positions[0] = valid_curves * 2;
    segments_arr = (PyArrayObject *)PyArray_SimpleNew(3, output_dims_segments, NPY_DOUBLE);
    source_positions_arr = (PyArrayObject *)PyArray_SimpleNew(1, output_dims_positions, NPY_INTP);
    if (segments_arr == NULL || source_positions_arr == NULL)
    {
        goto cleanup;
    }

    segments = (double *)PyArray_DATA(segments_arr);
    source_positions = (npy_intp *)PyArray_DATA(source_positions_arr);

    {
        npy_intp out_idx = 0;
        for (npy_intp curve_idx = 0; curve_idx < n_curves; ++curve_idx)
        {
            npy_intp start = curve_offsets[curve_idx];
            npy_intp end = curve_offsets[curve_idx + 1];
            double vertices[6];
            if (!build_open_arrow_vertices_for_curve(
                    display_points + start * 2,
                    end - start,
                    head_lengths[curve_idx],
                    head_widths[curve_idx],
                    vertices))
            {
                continue;
            }

            segments[out_idx * 4 + 0] = vertices[0];
            segments[out_idx * 4 + 1] = vertices[1];
            segments[out_idx * 4 + 2] = vertices[2];
            segments[out_idx * 4 + 3] = vertices[3];
            source_positions[out_idx] = curve_idx;
            out_idx += 1;

            segments[out_idx * 4 + 0] = vertices[4];
            segments[out_idx * 4 + 1] = vertices[5];
            segments[out_idx * 4 + 2] = vertices[2];
            segments[out_idx * 4 + 3] = vertices[3];
            source_positions[out_idx] = curve_idx;
            out_idx += 1;
        }
    }

    result = Py_BuildValue("NN", (PyObject *)segments_arr, (PyObject *)source_positions_arr);
    segments_arr = NULL;
    source_positions_arr = NULL;

cleanup:
    Py_XDECREF(display_points_arr);
    Py_XDECREF(curve_offsets_arr);
    Py_XDECREF(head_lengths_arr);
    Py_XDECREF(head_widths_arr);
    Py_XDECREF(segments_arr);
    Py_XDECREF(source_positions_arr);
    return result;
}

/*
 * Build compact filled-arrow polygon vertices directly in native code.
 *
 * Output layout:
 * - polygons: (K, 3, 2), one triangle per valid curve
 * - source_positions: (K,), one source curve index per polygon
 */
static PyObject *build_filled_arrow_polygons(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *display_points_obj;
    PyObject *curve_offsets_obj;
    PyObject *head_lengths_obj;
    PyObject *head_widths_obj;
    PyArrayObject *display_points_arr = NULL;
    PyArrayObject *curve_offsets_arr = NULL;
    PyArrayObject *head_lengths_arr = NULL;
    PyArrayObject *head_widths_arr = NULL;
    PyArrayObject *polygons_arr = NULL;
    PyArrayObject *source_positions_arr = NULL;
    PyObject *result = NULL;
    npy_intp n_curves;
    npy_intp n_points;
    npy_intp valid_curves = 0;
    npy_intp output_dims_polygons[3];
    npy_intp output_dims_positions[1];
    const double *display_points;
    const npy_intp *curve_offsets;
    const double *head_lengths;
    const double *head_widths;
    double *polygons;
    npy_intp *source_positions;
    static char *kwlist[] = {
        "display_points",
        "curve_offsets",
        "head_lengths_px",
        "head_widths_px",
        NULL,
    };
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOO:build_filled_arrow_polygons",
            kwlist,
            &display_points_obj,
            &curve_offsets_obj,
            &head_lengths_obj,
            &head_widths_obj))
    {
        return NULL;
    }

    display_points_arr = (PyArrayObject *)PyArray_FROM_OTF(
        display_points_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    curve_offsets_arr = (PyArrayObject *)PyArray_FROM_OTF(
        curve_offsets_obj, NPY_INTP, NPY_ARRAY_IN_ARRAY);
    head_lengths_arr = (PyArrayObject *)PyArray_FROM_OTF(
        head_lengths_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    head_widths_arr = (PyArrayObject *)PyArray_FROM_OTF(
        head_widths_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (display_points_arr == NULL || curve_offsets_arr == NULL ||
        head_lengths_arr == NULL || head_widths_arr == NULL)
    {
        goto cleanup;
    }

    if (PyArray_NDIM(display_points_arr) != 2 || PyArray_DIM(display_points_arr, 1) != 2)
    {
        PyErr_SetString(PyExc_ValueError, "display_points must have shape (M, 2)");
        goto cleanup;
    }
    if (PyArray_NDIM(curve_offsets_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "curve_offsets must be a 1D array");
        goto cleanup;
    }
    if (PyArray_NDIM(head_lengths_arr) != 1 || PyArray_NDIM(head_widths_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "head_lengths_px and head_widths_px must be 1D arrays");
        goto cleanup;
    }

    n_curves = PyArray_DIM(head_lengths_arr, 0);
    if (PyArray_DIM(head_widths_arr, 0) != n_curves)
    {
        PyErr_SetString(PyExc_ValueError, "head_lengths_px and head_widths_px must have the same length");
        goto cleanup;
    }
    if (PyArray_DIM(curve_offsets_arr, 0) != n_curves + 1)
    {
        PyErr_SetString(PyExc_ValueError, "curve_offsets must have length len(head_lengths_px) + 1");
        goto cleanup;
    }

    display_points = (const double *)PyArray_DATA(display_points_arr);
    curve_offsets = (const npy_intp *)PyArray_DATA(curve_offsets_arr);
    head_lengths = (const double *)PyArray_DATA(head_lengths_arr);
    head_widths = (const double *)PyArray_DATA(head_widths_arr);
    n_points = PyArray_DIM(display_points_arr, 0);

    for (npy_intp idx = 0; idx < n_curves; ++idx)
    {
        if (curve_offsets[idx] < 0 || curve_offsets[idx] > curve_offsets[idx + 1] ||
            curve_offsets[idx + 1] > n_points)
        {
            PyErr_SetString(PyExc_ValueError, "curve_offsets must be monotonic and stay within display_points");
            goto cleanup;
        }
    }

    for (npy_intp curve_idx = 0; curve_idx < n_curves; ++curve_idx)
    {
        npy_intp start = curve_offsets[curve_idx];
        npy_intp end = curve_offsets[curve_idx + 1];
        double vertices[6];
        if (build_filled_arrow_vertices_for_curve(
                display_points + start * 2,
                end - start,
                head_lengths[curve_idx],
                head_widths[curve_idx],
                vertices))
        {
            valid_curves += 1;
        }
    }

    output_dims_polygons[0] = valid_curves;
    output_dims_polygons[1] = 3;
    output_dims_polygons[2] = 2;
    output_dims_positions[0] = valid_curves;
    polygons_arr = (PyArrayObject *)PyArray_SimpleNew(3, output_dims_polygons, NPY_DOUBLE);
    source_positions_arr = (PyArrayObject *)PyArray_SimpleNew(1, output_dims_positions, NPY_INTP);
    if (polygons_arr == NULL || source_positions_arr == NULL)
    {
        goto cleanup;
    }

    polygons = (double *)PyArray_DATA(polygons_arr);
    source_positions = (npy_intp *)PyArray_DATA(source_positions_arr);

    {
        npy_intp out_idx = 0;
        for (npy_intp curve_idx = 0; curve_idx < n_curves; ++curve_idx)
        {
            npy_intp start = curve_offsets[curve_idx];
            npy_intp end = curve_offsets[curve_idx + 1];
            double vertices[6];
            if (!build_filled_arrow_vertices_for_curve(
                    display_points + start * 2,
                    end - start,
                    head_lengths[curve_idx],
                    head_widths[curve_idx],
                    vertices))
            {
                continue;
            }

            memcpy(polygons + out_idx * 6, vertices, 6 * sizeof(double));
            source_positions[out_idx] = curve_idx;
            out_idx += 1;
        }
    }

    result = Py_BuildValue("NN", (PyObject *)polygons_arr, (PyObject *)source_positions_arr);
    polygons_arr = NULL;
    source_positions_arr = NULL;

cleanup:
    Py_XDECREF(display_points_arr);
    Py_XDECREF(curve_offsets_arr);
    Py_XDECREF(head_lengths_arr);
    Py_XDECREF(head_widths_arr);
    Py_XDECREF(polygons_arr);
    Py_XDECREF(source_positions_arr);
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
        "compute_display_cell_valid",
        (PyCFunction)compute_display_cell_valid,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Compute valid display-grid cells for NCL-like tracing."),
    },
    {
        "sample_display_grid",
        (PyCFunction)sample_display_grid,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Bilinearly sample display-grid positions and Jacobians."),
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
    {
        "thin_display_candidates",
        (PyCFunction)thin_display_candidates,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Cull later nearby display candidate points after viewport normalization."),
    },
    {
        "generate_cell_candidates",
        (PyCFunction)generate_cell_candidates,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Expand selected scatter cells into interior candidate points."),
    },
    {
        "validate_display_curve",
        (PyCFunction)validate_display_curve,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Validate display-space curly-vector geometry."),
    },
    {
        "build_open_arrow_segments",
        (PyCFunction)build_open_arrow_segments,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Build batched display-space open-arrow head line segments."),
    },
    {
        "build_filled_arrow_polygons",
        (PyCFunction)build_filled_arrow_polygons,
        METH_VARARGS | METH_KEYWORDS,
        PyDoc_STR("Build batched display-space filled-arrow triangle polygons."),
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
