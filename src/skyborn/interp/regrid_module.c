#include <Python.h>

#include <math.h>
#include <string.h>

#include <numpy/arrayobject.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int
is_strictly_increasing(const double *values, npy_intp n)
{
    if (n <= 1) return 1;

    for (npy_intp i = 1; i < n; ++i)
        if (!(values[i] > values[i - 1])) return 0;

    return 1;
}

static int
normalize_weight_rows(PyArrayObject *weights_arr)
{
    npy_intp nrows = PyArray_DIM(weights_arr, 0);
    npy_intp ncols = PyArray_DIM(weights_arr, 1);
    double *weights = (double *) PyArray_DATA(weights_arr);

    if (ncols <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "weight matrix must have at least one source column");
        return -1;
    }

    for (npy_intp i = 0; i < nrows; ++i)
    {
        double row_sum = 0.0;
        double *row = weights + i * ncols;

        for (npy_intp j = 0; j < ncols; ++j)
            row_sum += row[j];

        if (row_sum > 1e-15)
        {
            for (npy_intp j = 0; j < ncols; ++j)
                row[j] /= row_sum;
        }
        else
        {
            double uniform = 1.0 / (double) ncols;
            for (npy_intp j = 0; j < ncols; ++j)
                row[j] = uniform;
        }
    }

    return 0;
}

static PyObject *
conservative_latitude_weights(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *src_lat_obj = NULL;
    PyObject *tgt_lat_obj = NULL;

    if (!PyArg_ParseTuple(args, "OO", &src_lat_obj, &tgt_lat_obj)) return NULL;

    PyArrayObject *src_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if (src_lat_arr == NULL || tgt_lat_arr == NULL) goto fail;

    if (PyArray_NDIM(src_lat_arr) != 1 || PyArray_NDIM(tgt_lat_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "conservative_latitude_weights expects 1D latitude arrays");
        goto fail;
    }

    npy_intp nsrc = PyArray_DIM(src_lat_arr, 0);
    npy_intp ntgt = PyArray_DIM(tgt_lat_arr, 0);
    if (nsrc <= 0 || ntgt <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "latitude coordinate arrays must not be empty");
        goto fail;
    }

    const double *src = (const double *) PyArray_DATA(src_lat_arr);
    const double *tgt = (const double *) PyArray_DATA(tgt_lat_arr);
    if (!is_strictly_increasing(src, nsrc) || !is_strictly_increasing(tgt, ntgt))
    {
        PyErr_SetString(PyExc_ValueError, "conservative_latitude_weights requires strictly increasing coordinates");
        goto fail;
    }

    npy_intp src_bounds_len = nsrc + 1;
    npy_intp tgt_bounds_len = ntgt + 1;
    double *src_bounds = PyMem_Malloc((size_t) src_bounds_len * sizeof(double));
    double *tgt_bounds = PyMem_Malloc((size_t) tgt_bounds_len * sizeof(double));
    if (src_bounds == NULL || tgt_bounds == NULL)
    {
        PyErr_NoMemory();
        PyMem_Free(src_bounds);
        PyMem_Free(tgt_bounds);
        goto fail;
    }

    src_bounds[0] = -M_PI / 2.0;
    src_bounds[nsrc] = M_PI / 2.0;
    for (npy_intp i = 1; i < nsrc; ++i)
        src_bounds[i] = 0.5 * (src[i - 1] + src[i]);

    tgt_bounds[0] = -M_PI / 2.0;
    tgt_bounds[ntgt] = M_PI / 2.0;
    for (npy_intp i = 1; i < ntgt; ++i)
        tgt_bounds[i] = 0.5 * (tgt[i - 1] + tgt[i]);

    npy_intp dims[2] = { ntgt, nsrc };
    PyArrayObject *weights_arr = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (weights_arr == NULL)
    {
        PyMem_Free(src_bounds);
        PyMem_Free(tgt_bounds);
        goto fail;
    }

    double *weights = (double *) PyArray_DATA(weights_arr);
    for (npy_intp itgt = 0; itgt < ntgt; ++itgt)
    {
        for (npy_intp isrc = 0; isrc < nsrc; ++isrc)
        {
            double upper = fmin(tgt_bounds[itgt + 1], src_bounds[isrc + 1]);
            double lower = fmax(tgt_bounds[itgt], src_bounds[isrc]);
            double overlap = (upper > lower) ? (sin(upper) - sin(lower)) : 0.0;
            weights[itgt * nsrc + isrc] = overlap;
        }
    }

    PyMem_Free(src_bounds);
    PyMem_Free(tgt_bounds);

    if (normalize_weight_rows(weights_arr) != 0)
    {
        Py_DECREF(weights_arr);
        goto fail;
    }

    Py_DECREF(src_lat_arr);
    Py_DECREF(tgt_lat_arr);
    return PyArray_Return(weights_arr);

fail:
    Py_XDECREF(src_lat_arr);
    Py_XDECREF(tgt_lat_arr);
    return NULL;
}

static double
align_phase_with(double x, double target, double period)
{
    if (x > target + period / 2.0) x -= period;
    if (x < target - period / 2.0) x += period;
    return x;
}

static PyObject *
conservative_longitude_weights(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *src_lon_obj = NULL;
    PyObject *tgt_lon_obj = NULL;

    if (!PyArg_ParseTuple(args, "OO", &src_lon_obj, &tgt_lon_obj)) return NULL;

    PyArrayObject *src_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if (src_lon_arr == NULL || tgt_lon_arr == NULL) goto fail;

    if (PyArray_NDIM(src_lon_arr) != 1 || PyArray_NDIM(tgt_lon_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "conservative_longitude_weights expects 1D longitude arrays");
        goto fail;
    }

    npy_intp nsrc = PyArray_DIM(src_lon_arr, 0);
    npy_intp ntgt = PyArray_DIM(tgt_lon_arr, 0);
    if (nsrc <= 0 || ntgt <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "longitude coordinate arrays must not be empty");
        goto fail;
    }

    const double *src = (const double *) PyArray_DATA(src_lon_arr);
    const double *tgt = (const double *) PyArray_DATA(tgt_lon_arr);
    if (!is_strictly_increasing(src, nsrc) || !is_strictly_increasing(tgt, ntgt))
    {
        PyErr_SetString(PyExc_ValueError, "conservative_longitude_weights requires strictly increasing coordinates");
        goto fail;
    }

    double period = 2.0 * M_PI;
    double *src_mod = PyMem_Malloc((size_t) nsrc * sizeof(double));
    double *tgt_mod = PyMem_Malloc((size_t) ntgt * sizeof(double));
    double *src_lower = PyMem_Malloc((size_t) nsrc * sizeof(double));
    double *src_upper = PyMem_Malloc((size_t) nsrc * sizeof(double));
    double *tgt_lower = PyMem_Malloc((size_t) ntgt * sizeof(double));
    double *tgt_upper = PyMem_Malloc((size_t) ntgt * sizeof(double));
    if (src_mod == NULL || tgt_mod == NULL || src_lower == NULL || src_upper == NULL || tgt_lower == NULL
        || tgt_upper == NULL)
    {
        PyErr_NoMemory();
        PyMem_Free(src_mod);
        PyMem_Free(tgt_mod);
        PyMem_Free(src_lower);
        PyMem_Free(src_upper);
        PyMem_Free(tgt_lower);
        PyMem_Free(tgt_upper);
        goto fail;
    }

    for (npy_intp i = 0; i < nsrc; ++i)
    {
        src_mod[i] = fmod(src[i], period);
        if (src_mod[i] < 0.0) src_mod[i] += period;
    }
    for (npy_intp i = 0; i < ntgt; ++i)
    {
        tgt_mod[i] = fmod(tgt[i], period);
        if (tgt_mod[i] < 0.0) tgt_mod[i] += period;
    }

    for (npy_intp i = 0; i < nsrc; ++i)
    {
        double x_plus = align_phase_with(src_mod[(i + 1) % nsrc], src_mod[i], period);
        double x_minus = align_phase_with(src_mod[(i + nsrc - 1) % nsrc], src_mod[i], period);
        src_upper[i] = 0.5 * (src_mod[i] + x_plus);
        src_lower[i] = 0.5 * (x_minus + src_mod[i]);
    }

    for (npy_intp i = 0; i < ntgt; ++i)
    {
        double x_plus = align_phase_with(tgt_mod[(i + 1) % ntgt], tgt_mod[i], period);
        double x_minus = align_phase_with(tgt_mod[(i + ntgt - 1) % ntgt], tgt_mod[i], period);
        tgt_upper[i] = 0.5 * (tgt_mod[i] + x_plus);
        tgt_lower[i] = 0.5 * (x_minus + tgt_mod[i]);
    }

    npy_intp dims[2] = { ntgt, nsrc };
    PyArrayObject *weights_arr = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_FLOAT64);
    if (weights_arr == NULL)
    {
        PyMem_Free(src_mod);
        PyMem_Free(tgt_mod);
        PyMem_Free(src_lower);
        PyMem_Free(src_upper);
        PyMem_Free(tgt_lower);
        PyMem_Free(tgt_upper);
        goto fail;
    }

    double *weights = (double *) PyArray_DATA(weights_arr);
    for (npy_intp itgt = 0; itgt < ntgt; ++itgt)
    {
        for (npy_intp isrc = 0; isrc < nsrc; ++isrc)
        {
            double y0 = align_phase_with(src_lower[isrc], tgt_lower[itgt], period);
            double y1 = align_phase_with(src_upper[isrc], tgt_lower[itgt], period);
            double upper = fmin(tgt_upper[itgt], y1);
            double lower = fmax(tgt_lower[itgt], y0);
            double overlap = upper - lower;
            weights[itgt * nsrc + isrc] = (overlap > 0.0) ? overlap : 0.0;
        }
    }

    PyMem_Free(src_mod);
    PyMem_Free(tgt_mod);
    PyMem_Free(src_lower);
    PyMem_Free(src_upper);
    PyMem_Free(tgt_lower);
    PyMem_Free(tgt_upper);

    if (normalize_weight_rows(weights_arr) != 0)
    {
        Py_DECREF(weights_arr);
        goto fail;
    }

    Py_DECREF(src_lon_arr);
    Py_DECREF(tgt_lon_arr);
    return PyArray_Return(weights_arr);

fail:
    Py_XDECREF(src_lon_arr);
    Py_XDECREF(tgt_lon_arr);
    return NULL;
}

static PyObject *
conservative_regrid(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *field_obj = NULL;
    PyObject *lon_weights_obj = NULL;
    PyObject *lat_weights_obj = NULL;

    if (!PyArg_ParseTuple(args, "OOO", &field_obj, &lon_weights_obj, &lat_weights_obj)) return NULL;

    PyArrayObject *field_arr = (PyArrayObject *) PyArray_FROM_OTF(field_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *lon_weights_arr =
        (PyArrayObject *) PyArray_FROM_OTF(lon_weights_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *lat_weights_arr =
        (PyArrayObject *) PyArray_FROM_OTF(lat_weights_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *field_contig = NULL;
    PyArrayObject *lon_weights_contig = NULL;
    PyArrayObject *lat_weights_contig = NULL;

    if (field_arr == NULL || lon_weights_arr == NULL || lat_weights_arr == NULL) goto fail;

    if (PyArray_NDIM(field_arr) < 2 || PyArray_NDIM(lon_weights_arr) != 2 || PyArray_NDIM(lat_weights_arr) != 2)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "conservative_regrid expects field ndim >= 2 and 2D weight matrices");
        goto fail;
    }

    field_contig = (PyArrayObject *) PyArray_GETCONTIGUOUS(field_arr);
    lon_weights_contig = (PyArrayObject *) PyArray_GETCONTIGUOUS(lon_weights_arr);
    lat_weights_contig = (PyArrayObject *) PyArray_GETCONTIGUOUS(lat_weights_arr);
    if (field_contig == NULL || lon_weights_contig == NULL || lat_weights_contig == NULL) goto fail;

    int field_ndim = PyArray_NDIM(field_contig);
    npy_intp nsrc_lon = PyArray_DIM(field_contig, field_ndim - 2);
    npy_intp nsrc_lat = PyArray_DIM(field_contig, field_ndim - 1);
    npy_intp ntgt_lon = PyArray_DIM(lon_weights_contig, 0);
    npy_intp lon_weight_src = PyArray_DIM(lon_weights_contig, 1);
    npy_intp ntgt_lat = PyArray_DIM(lat_weights_contig, 0);
    npy_intp lat_weight_src = PyArray_DIM(lat_weights_contig, 1);

    if (nsrc_lon != lon_weight_src || nsrc_lat != lat_weight_src)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "conservative_regrid field shape must match source dimensions implied by weights");
        goto fail;
    }

    npy_intp out_dims[NPY_MAXDIMS];
    if (field_ndim > NPY_MAXDIMS)
    {
        PyErr_SetString(PyExc_ValueError, "field has too many dimensions");
        goto fail;
    }
    for (int i = 0; i < field_ndim - 2; ++i)
        out_dims[i] = PyArray_DIM(field_contig, i);
    out_dims[field_ndim - 2] = ntgt_lon;
    out_dims[field_ndim - 1] = ntgt_lat;

    PyArrayObject *out_arr = (PyArrayObject *) PyArray_SimpleNew(field_ndim, out_dims, NPY_FLOAT64);
    if (out_arr == NULL) goto fail;

    npy_intp outer_size = 1;
    for (int i = 0; i < field_ndim - 2; ++i)
        outer_size *= PyArray_DIM(field_contig, i);

    npy_intp src_plane = nsrc_lon * nsrc_lat;
    npy_intp tmp_plane = nsrc_lon * ntgt_lat;
    npy_intp tgt_plane = ntgt_lon * ntgt_lat;

    double *tmp_total = PyMem_Malloc((size_t) tmp_plane * sizeof(double));
    double *tmp_count = PyMem_Malloc((size_t) tmp_plane * sizeof(double));
    npy_intp *lon_start = PyMem_Malloc((size_t) ntgt_lon * sizeof(npy_intp));
    npy_intp *lon_stop = PyMem_Malloc((size_t) ntgt_lon * sizeof(npy_intp));
    npy_intp *lat_start = PyMem_Malloc((size_t) ntgt_lat * sizeof(npy_intp));
    npy_intp *lat_stop = PyMem_Malloc((size_t) ntgt_lat * sizeof(npy_intp));
    if (tmp_total == NULL || tmp_count == NULL || lon_start == NULL || lon_stop == NULL || lat_start == NULL
        || lat_stop == NULL)
    {
        PyErr_NoMemory();
        PyMem_Free(tmp_total);
        PyMem_Free(tmp_count);
        PyMem_Free(lon_start);
        PyMem_Free(lon_stop);
        PyMem_Free(lat_start);
        PyMem_Free(lat_stop);
        Py_DECREF(out_arr);
        goto fail;
    }

    const double *field = (const double *) PyArray_DATA(field_contig);
    const double *lon_weights = (const double *) PyArray_DATA(lon_weights_contig);
    const double *lat_weights = (const double *) PyArray_DATA(lat_weights_contig);
    double *out = (double *) PyArray_DATA(out_arr);

    for (npy_intp itgt_lon = 0; itgt_lon < ntgt_lon; ++itgt_lon)
    {
        const double *row = lon_weights + itgt_lon * nsrc_lon;
        npy_intp start = nsrc_lon;
        npy_intp stop = 0;

        for (npy_intp isrc_lon = 0; isrc_lon < nsrc_lon; ++isrc_lon)
        {
            if (row[isrc_lon] > 0.0)
            {
                if (start == nsrc_lon) start = isrc_lon;
                stop = isrc_lon + 1;
            }
        }

        if (start == nsrc_lon)
        {
            start = 0;
            stop = nsrc_lon;
        }
        lon_start[itgt_lon] = start;
        lon_stop[itgt_lon] = stop;
    }

    for (npy_intp itgt_lat = 0; itgt_lat < ntgt_lat; ++itgt_lat)
    {
        const double *row = lat_weights + itgt_lat * nsrc_lat;
        npy_intp start = nsrc_lat;
        npy_intp stop = 0;

        for (npy_intp isrc_lat = 0; isrc_lat < nsrc_lat; ++isrc_lat)
        {
            if (row[isrc_lat] > 0.0)
            {
                if (start == nsrc_lat) start = isrc_lat;
                stop = isrc_lat + 1;
            }
        }

        if (start == nsrc_lat)
        {
            start = 0;
            stop = nsrc_lat;
        }
        lat_start[itgt_lat] = start;
        lat_stop[itgt_lat] = stop;
    }

    for (npy_intp outer = 0; outer < outer_size; ++outer)
    {
        const double *field_slice = field + outer * src_plane;
        double *out_slice = out + outer * tgt_plane;

        for (npy_intp isrc_lon = 0; isrc_lon < nsrc_lon; ++isrc_lon)
        {
            const double *field_row = field_slice + isrc_lon * nsrc_lat;

            for (npy_intp itgt_lat = 0; itgt_lat < ntgt_lat; ++itgt_lat)
            {
                const double *lat_row = lat_weights + itgt_lat * nsrc_lat;
                npy_intp lat_begin = lat_start[itgt_lat];
                npy_intp lat_end = lat_stop[itgt_lat];
                double total = 0.0;
                double count = 0.0;

                for (npy_intp isrc_lat = lat_begin; isrc_lat < lat_end; ++isrc_lat)
                {
                    double value = field_row[isrc_lat];
                    double weight = lat_row[isrc_lat];
                    if (!isnan(value))
                    {
                        total += weight * value;
                        count += weight;
                    }
                }

                tmp_total[isrc_lon * ntgt_lat + itgt_lat] = total;
                tmp_count[isrc_lon * ntgt_lat + itgt_lat] = count;
            }
        }

        for (npy_intp itgt_lon = 0; itgt_lon < ntgt_lon; ++itgt_lon)
        {
            const double *lon_row = lon_weights + itgt_lon * nsrc_lon;
            npy_intp lon_begin = lon_start[itgt_lon];
            npy_intp lon_end = lon_stop[itgt_lon];

            for (npy_intp itgt_lat = 0; itgt_lat < ntgt_lat; ++itgt_lat)
            {
                double total = 0.0;
                double count = 0.0;

                for (npy_intp isrc_lon = lon_begin; isrc_lon < lon_end; ++isrc_lon)
                {
                    double weight = lon_row[isrc_lon];
                    npy_intp tmp_index = isrc_lon * ntgt_lat + itgt_lat;
                    total += weight * tmp_total[tmp_index];
                    count += weight * tmp_count[tmp_index];
                }

                out_slice[itgt_lon * ntgt_lat + itgt_lat] = (count > 0.0) ? (total / count) : Py_NAN;
            }
        }
    }

    PyMem_Free(tmp_total);
    PyMem_Free(tmp_count);
    PyMem_Free(lon_start);
    PyMem_Free(lon_stop);
    PyMem_Free(lat_start);
    PyMem_Free(lat_stop);
    Py_DECREF(field_arr);
    Py_DECREF(lon_weights_arr);
    Py_DECREF(lat_weights_arr);
    Py_DECREF(field_contig);
    Py_DECREF(lon_weights_contig);
    Py_DECREF(lat_weights_contig);
    return PyArray_Return(out_arr);

fail:
    Py_XDECREF(field_arr);
    Py_XDECREF(lon_weights_arr);
    Py_XDECREF(lat_weights_arr);
    Py_XDECREF(field_contig);
    Py_XDECREF(lon_weights_contig);
    Py_XDECREF(lat_weights_contig);
    return NULL;
}

static int
bracket_clamped(double x, const double *values, npy_intp n, npy_intp *left, double *frac)
{
    if (n <= 1)
    {
        *left = 0;
        *frac = 0.0;
        return 0;
    }

    if (x <= values[0])
    {
        *left = 0;
        *frac = 0.0;
        return 0;
    }

    if (x >= values[n - 1])
    {
        *left = n - 2;
        *frac = 1.0;
        return 0;
    }

    npy_intp lo = 0;
    npy_intp hi = n - 1;

    while (hi - lo > 1)
    {
        npy_intp mid = lo + (hi - lo) / 2;
        if (values[mid] <= x)
            lo = mid;
        else
            hi = mid;
    }

    double denom = values[lo + 1] - values[lo];
    if (denom <= 0.0)
    {
        *left = lo;
        *frac = 0.0;
        return -1;
    }

    *left = lo;
    *frac = (x - values[lo]) / denom;
    return 0;
}

static int
bilinear_regrid_2d_core(
    const double *field,
    npy_intp nsrc_lon,
    npy_intp nsrc_lat,
    const double *src_lon,
    const double *src_lat,
    npy_intp ntgt_lon,
    npy_intp ntgt_lat,
    const double *tgt_lon,
    const double *tgt_lat,
    double *out)
{
    for (npy_intp ilon = 0; ilon < ntgt_lon; ++ilon)
    {
        npy_intp lon_left = 0;
        double lon_frac = 0.0;
        if (bracket_clamped(tgt_lon[ilon], src_lon, nsrc_lon, &lon_left, &lon_frac) != 0)
        {
            return -1;
        }

        npy_intp lon_right = (nsrc_lon > 1) ? (lon_left + 1) : lon_left;
        for (npy_intp ilat = 0; ilat < ntgt_lat; ++ilat)
        {
            npy_intp lat_left = 0;
            double lat_frac = 0.0;
            if (bracket_clamped(tgt_lat[ilat], src_lat, nsrc_lat, &lat_left, &lat_frac) != 0)
            {
                return -1;
            }

            npy_intp lat_right = (nsrc_lat > 1) ? (lat_left + 1) : lat_left;

            double v00 = field[lon_left * nsrc_lat + lat_left];
            double v10 = field[lon_right * nsrc_lat + lat_left];
            double v01 = field[lon_left * nsrc_lat + lat_right];
            double v11 = field[lon_right * nsrc_lat + lat_right];

            double v0 = v00 + lon_frac * (v10 - v00);
            double v1 = v01 + lon_frac * (v11 - v01);
            out[ilon * ntgt_lat + ilat] = v0 + lat_frac * (v1 - v0);
        }
    }

    return 0;
}

static PyObject *
bilinear_regrid(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *field_obj = NULL;
    PyObject *src_lon_obj = NULL;
    PyObject *src_lat_obj = NULL;
    PyObject *tgt_lon_obj = NULL;
    PyObject *tgt_lat_obj = NULL;

    if (!PyArg_ParseTuple(
            args,
            "OOOOO",
            &field_obj,
            &src_lon_obj,
            &src_lat_obj,
            &tgt_lon_obj,
            &tgt_lat_obj))
        return NULL;

    PyArrayObject *field_arr = (PyArrayObject *) PyArray_FROM_OTF(field_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *src_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *src_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if (field_arr == NULL || src_lon_arr == NULL || src_lat_arr == NULL || tgt_lon_arr == NULL || tgt_lat_arr == NULL)
        goto fail;

    if (PyArray_NDIM(field_arr) != 2 || PyArray_NDIM(src_lon_arr) != 1 || PyArray_NDIM(src_lat_arr) != 1
        || PyArray_NDIM(tgt_lon_arr) != 1 || PyArray_NDIM(tgt_lat_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "bilinear_regrid expects one 2D field and four 1D coordinate arrays");
        goto fail;
    }

    npy_intp nsrc_lon = PyArray_DIM(src_lon_arr, 0);
    npy_intp nsrc_lat = PyArray_DIM(src_lat_arr, 0);
    npy_intp ntgt_lon = PyArray_DIM(tgt_lon_arr, 0);
    npy_intp ntgt_lat = PyArray_DIM(tgt_lat_arr, 0);

    if (PyArray_DIM(field_arr, 0) != nsrc_lon || PyArray_DIM(field_arr, 1) != nsrc_lat)
    {
        PyErr_SetString(PyExc_ValueError, "field shape must match source longitude/latitude lengths");
        goto fail;
    }

    if (nsrc_lon <= 0 || nsrc_lat <= 0 || ntgt_lon <= 0 || ntgt_lat <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "grid coordinate arrays must not be empty");
        goto fail;
    }

    const double *src_lon = (const double *) PyArray_DATA(src_lon_arr);
    const double *src_lat = (const double *) PyArray_DATA(src_lat_arr);
    const double *tgt_lon = (const double *) PyArray_DATA(tgt_lon_arr);
    const double *tgt_lat = (const double *) PyArray_DATA(tgt_lat_arr);
    const double *field = (const double *) PyArray_DATA(field_arr);

    if (!is_strictly_increasing(src_lon, nsrc_lon) || !is_strictly_increasing(src_lat, nsrc_lat)
        || !is_strictly_increasing(tgt_lon, ntgt_lon) || !is_strictly_increasing(tgt_lat, ntgt_lat))
    {
        PyErr_SetString(PyExc_ValueError, "bilinear_regrid requires strictly increasing source and target coordinates");
        goto fail;
    }

    npy_intp out_dims[2] = { ntgt_lon, ntgt_lat };
    PyArrayObject *out_arr = (PyArrayObject *) PyArray_SimpleNew(2, out_dims, NPY_FLOAT64);
    if (out_arr == NULL) goto fail;

    double *out = (double *) PyArray_DATA(out_arr);
    if (bilinear_regrid_2d_core(
            field, nsrc_lon, nsrc_lat, src_lon, src_lat, ntgt_lon, ntgt_lat, tgt_lon, tgt_lat, out)
        != 0)
    {
        PyErr_SetString(PyExc_ValueError, "bilinear_regrid received invalid source or target coordinates");
        goto fail_output;
    }

    Py_DECREF(field_arr);
    Py_DECREF(src_lon_arr);
    Py_DECREF(src_lat_arr);
    Py_DECREF(tgt_lon_arr);
    Py_DECREF(tgt_lat_arr);

    return PyArray_Return(out_arr);

fail_output:
    Py_DECREF(out_arr);
    goto fail;

fail:
    Py_XDECREF(field_arr);
    Py_XDECREF(src_lon_arr);
    Py_XDECREF(src_lat_arr);
    Py_XDECREF(tgt_lon_arr);
    Py_XDECREF(tgt_lat_arr);
    return NULL;
}

static PyObject *
bilinear_regrid_nd(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *field_obj = NULL;
    PyObject *src_lon_obj = NULL;
    PyObject *src_lat_obj = NULL;
    PyObject *tgt_lon_obj = NULL;
    PyObject *tgt_lat_obj = NULL;

    if (!PyArg_ParseTuple(
            args,
            "OOOOO",
            &field_obj,
            &src_lon_obj,
            &src_lat_obj,
            &tgt_lon_obj,
            &tgt_lat_obj))
        return NULL;

    PyArrayObject *field_arr = (PyArrayObject *) PyArray_FROM_OTF(field_obj, NPY_FLOAT64, NPY_ARRAY_CARRAY_RO);
    PyArrayObject *src_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *src_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if (field_arr == NULL || src_lon_arr == NULL || src_lat_arr == NULL || tgt_lon_arr == NULL || tgt_lat_arr == NULL)
        goto fail;

    if (PyArray_NDIM(field_arr) < 2 || PyArray_NDIM(src_lon_arr) != 1 || PyArray_NDIM(src_lat_arr) != 1
        || PyArray_NDIM(tgt_lon_arr) != 1 || PyArray_NDIM(tgt_lat_arr) != 1)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "bilinear_regrid_nd expects field ndim >= 2 and four 1D coordinate arrays");
        goto fail;
    }

    npy_intp nsrc_lon = PyArray_DIM(src_lon_arr, 0);
    npy_intp nsrc_lat = PyArray_DIM(src_lat_arr, 0);
    npy_intp ntgt_lon = PyArray_DIM(tgt_lon_arr, 0);
    npy_intp ntgt_lat = PyArray_DIM(tgt_lat_arr, 0);
    int field_ndim = PyArray_NDIM(field_arr);

    if (PyArray_DIM(field_arr, field_ndim - 2) != nsrc_lon || PyArray_DIM(field_arr, field_ndim - 1) != nsrc_lat)
    {
        PyErr_SetString(PyExc_ValueError, "field shape must match source longitude/latitude lengths");
        goto fail;
    }

    if (nsrc_lon <= 0 || nsrc_lat <= 0 || ntgt_lon <= 0 || ntgt_lat <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "grid coordinate arrays must not be empty");
        goto fail;
    }

    const double *src_lon = (const double *) PyArray_DATA(src_lon_arr);
    const double *src_lat = (const double *) PyArray_DATA(src_lat_arr);
    const double *tgt_lon = (const double *) PyArray_DATA(tgt_lon_arr);
    const double *tgt_lat = (const double *) PyArray_DATA(tgt_lat_arr);

    if (!is_strictly_increasing(src_lon, nsrc_lon) || !is_strictly_increasing(src_lat, nsrc_lat)
        || !is_strictly_increasing(tgt_lon, ntgt_lon) || !is_strictly_increasing(tgt_lat, ntgt_lat))
    {
        PyErr_SetString(
            PyExc_ValueError,
            "bilinear_regrid_nd requires strictly increasing source and target coordinates");
        goto fail;
    }

    npy_intp out_dims[NPY_MAXDIMS];
    if (field_ndim > NPY_MAXDIMS)
    {
        PyErr_SetString(PyExc_ValueError, "field has too many dimensions");
        goto fail;
    }
    for (int i = 0; i < field_ndim - 2; ++i)
        out_dims[i] = PyArray_DIM(field_arr, i);
    out_dims[field_ndim - 2] = ntgt_lon;
    out_dims[field_ndim - 1] = ntgt_lat;

    PyArrayObject *out_arr = (PyArrayObject *) PyArray_SimpleNew(field_ndim, out_dims, NPY_FLOAT64);
    if (out_arr == NULL) goto fail;

    npy_intp outer_size = 1;
    for (int i = 0; i < field_ndim - 2; ++i)
        outer_size *= PyArray_DIM(field_arr, i);

    npy_intp src_plane = nsrc_lon * nsrc_lat;
    npy_intp tgt_plane = ntgt_lon * ntgt_lat;
    const double *field = (const double *) PyArray_DATA(field_arr);
    double *out = (double *) PyArray_DATA(out_arr);
    int status = 0;

    NPY_BEGIN_ALLOW_THREADS;
    for (npy_intp outer = 0; outer < outer_size; ++outer)
    {
        const double *field_slice = field + outer * src_plane;
        double *out_slice = out + outer * tgt_plane;
        status = bilinear_regrid_2d_core(
            field_slice,
            nsrc_lon,
            nsrc_lat,
            src_lon,
            src_lat,
            ntgt_lon,
            ntgt_lat,
            tgt_lon,
            tgt_lat,
            out_slice);
        if (status != 0)
            break;
    }
    NPY_END_ALLOW_THREADS;

    if (status != 0)
    {
        PyErr_SetString(PyExc_ValueError, "bilinear_regrid received invalid source or target coordinates");
        Py_DECREF(out_arr);
        goto fail;
    }

    Py_DECREF(field_arr);
    Py_DECREF(src_lon_arr);
    Py_DECREF(src_lat_arr);
    Py_DECREF(tgt_lon_arr);
    Py_DECREF(tgt_lat_arr);
    return PyArray_Return(out_arr);

fail:
    Py_XDECREF(field_arr);
    Py_XDECREF(src_lon_arr);
    Py_XDECREF(src_lat_arr);
    Py_XDECREF(tgt_lon_arr);
    Py_XDECREF(tgt_lat_arr);
    return NULL;
}

static PyObject *
nearest_neighbor_indices(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *src_lon_obj = NULL;
    PyObject *src_lat_obj = NULL;
    PyObject *tgt_lon_obj = NULL;
    PyObject *tgt_lat_obj = NULL;

    if (!PyArg_ParseTuple(args, "OOOO", &src_lon_obj, &src_lat_obj, &tgt_lon_obj, &tgt_lat_obj)) return NULL;

    PyArrayObject *src_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *src_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(src_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lon_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lon_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *tgt_lat_arr = (PyArrayObject *) PyArray_FROM_OTF(tgt_lat_obj, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);

    if (src_lon_arr == NULL || src_lat_arr == NULL || tgt_lon_arr == NULL || tgt_lat_arr == NULL) goto fail;

    if (PyArray_NDIM(src_lon_arr) != 1 || PyArray_NDIM(src_lat_arr) != 1 || PyArray_NDIM(tgt_lon_arr) != 1
        || PyArray_NDIM(tgt_lat_arr) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "nearest_neighbor_indices expects 1D longitude and latitude arrays");
        goto fail;
    }

    npy_intp nsrc_lon = PyArray_DIM(src_lon_arr, 0);
    npy_intp nsrc_lat = PyArray_DIM(src_lat_arr, 0);
    npy_intp ntgt_lon = PyArray_DIM(tgt_lon_arr, 0);
    npy_intp ntgt_lat = PyArray_DIM(tgt_lat_arr, 0);

    if (nsrc_lon <= 0 || nsrc_lat <= 0 || ntgt_lon <= 0 || ntgt_lat <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "grid coordinate arrays must not be empty");
        goto fail;
    }

    npy_intp out_size = ntgt_lon * ntgt_lat;
    npy_intp dims[1] = { out_size };
    PyArrayObject *indices_arr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INTP);
    if (indices_arr == NULL) goto fail;

    const double *src_lon = (const double *) PyArray_DATA(src_lon_arr);
    const double *src_lat = (const double *) PyArray_DATA(src_lat_arr);
    const double *tgt_lon = (const double *) PyArray_DATA(tgt_lon_arr);
    const double *tgt_lat = (const double *) PyArray_DATA(tgt_lat_arr);
    npy_intp *indices = (npy_intp *) PyArray_DATA(indices_arr);
    double *src_sin_lat = NULL;
    double *src_cos_lat = NULL;
    double *src_sin_lon = NULL;
    double *src_cos_lon = NULL;
    double *tgt_sin_lat = NULL;
    double *tgt_cos_lat = NULL;
    double *tgt_sin_lon = NULL;
    double *tgt_cos_lon = NULL;

    src_sin_lat = PyMem_Malloc((size_t) nsrc_lat * sizeof(double));
    src_cos_lat = PyMem_Malloc((size_t) nsrc_lat * sizeof(double));
    src_sin_lon = PyMem_Malloc((size_t) nsrc_lon * sizeof(double));
    src_cos_lon = PyMem_Malloc((size_t) nsrc_lon * sizeof(double));
    tgt_sin_lat = PyMem_Malloc((size_t) ntgt_lat * sizeof(double));
    tgt_cos_lat = PyMem_Malloc((size_t) ntgt_lat * sizeof(double));
    tgt_sin_lon = PyMem_Malloc((size_t) ntgt_lon * sizeof(double));
    tgt_cos_lon = PyMem_Malloc((size_t) ntgt_lon * sizeof(double));
    if (src_sin_lat == NULL || src_cos_lat == NULL || src_sin_lon == NULL || src_cos_lon == NULL
        || tgt_sin_lat == NULL || tgt_cos_lat == NULL || tgt_sin_lon == NULL || tgt_cos_lon == NULL)
    {
        PyErr_NoMemory();
        Py_XDECREF(indices_arr);
        PyMem_Free(src_sin_lat);
        PyMem_Free(src_cos_lat);
        PyMem_Free(src_sin_lon);
        PyMem_Free(src_cos_lon);
        PyMem_Free(tgt_sin_lat);
        PyMem_Free(tgt_cos_lat);
        PyMem_Free(tgt_sin_lon);
        PyMem_Free(tgt_cos_lon);
        goto fail;
    }

    for (npy_intp ilat = 0; ilat < nsrc_lat; ++ilat)
    {
        src_sin_lat[ilat] = sin(src_lat[ilat]);
        src_cos_lat[ilat] = cos(src_lat[ilat]);
    }
    for (npy_intp ilon = 0; ilon < nsrc_lon; ++ilon)
    {
        src_sin_lon[ilon] = sin(src_lon[ilon]);
        src_cos_lon[ilon] = cos(src_lon[ilon]);
    }
    for (npy_intp itgt_lat = 0; itgt_lat < ntgt_lat; ++itgt_lat)
    {
        tgt_sin_lat[itgt_lat] = sin(tgt_lat[itgt_lat]);
        tgt_cos_lat[itgt_lat] = cos(tgt_lat[itgt_lat]);
    }
    for (npy_intp itgt_lon = 0; itgt_lon < ntgt_lon; ++itgt_lon)
    {
        tgt_sin_lon[itgt_lon] = sin(tgt_lon[itgt_lon]);
        tgt_cos_lon[itgt_lon] = cos(tgt_lon[itgt_lon]);
    }

    for (npy_intp ilon = 0; ilon < ntgt_lon; ++ilon)
    {
        double sin_tgt_lon = tgt_sin_lon[ilon];
        double cos_tgt_lon = tgt_cos_lon[ilon];

        for (npy_intp itgt_lat = 0; itgt_lat < ntgt_lat; ++itgt_lat)
        {
            double sin_tgt_phi = tgt_sin_lat[itgt_lat];
            double cos_tgt_phi = tgt_cos_lat[itgt_lat];
            double best_score = -2.0;
            npy_intp best_index = 0;

            for (npy_intp isrc_lon = 0; isrc_lon < nsrc_lon; ++isrc_lon)
            {
                double cos_dlon = cos_tgt_lon * src_cos_lon[isrc_lon] + sin_tgt_lon * src_sin_lon[isrc_lon];

                for (npy_intp isrc_lat = 0; isrc_lat < nsrc_lat; ++isrc_lat)
                {
                    double score = sin_tgt_phi * src_sin_lat[isrc_lat]
                                   + cos_tgt_phi * src_cos_lat[isrc_lat] * cos_dlon;
                    npy_intp current_index = isrc_lon * nsrc_lat + isrc_lat;

                    if (score > best_score
                        || (fabs(score - best_score) <= 1e-15 && current_index < best_index))
                    {
                        best_score = score;
                        best_index = current_index;
                    }
                }
            }

            indices[ilon * ntgt_lat + itgt_lat] = best_index;
        }
    }

    PyMem_Free(src_sin_lat);
    PyMem_Free(src_cos_lat);
    PyMem_Free(src_sin_lon);
    PyMem_Free(src_cos_lon);
    PyMem_Free(tgt_sin_lat);
    PyMem_Free(tgt_cos_lat);
    PyMem_Free(tgt_sin_lon);
    PyMem_Free(tgt_cos_lon);

    Py_DECREF(src_lon_arr);
    Py_DECREF(src_lat_arr);
    Py_DECREF(tgt_lon_arr);
    Py_DECREF(tgt_lat_arr);

    return PyArray_Return(indices_arr);

fail:
    Py_XDECREF(src_lon_arr);
    Py_XDECREF(src_lat_arr);
    Py_XDECREF(tgt_lon_arr);
    Py_XDECREF(tgt_lat_arr);
    return NULL;
}

static PyObject *
nearest_regrid_apply(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *field_obj = NULL;
    PyObject *indices_obj = NULL;
    Py_ssize_t ntgt_lon_ssize = 0;
    Py_ssize_t ntgt_lat_ssize = 0;

    if (!PyArg_ParseTuple(args, "OOnn", &field_obj, &indices_obj, &ntgt_lon_ssize, &ntgt_lat_ssize)) return NULL;

    if (ntgt_lon_ssize <= 0 || ntgt_lat_ssize <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "target longitude/latitude sizes must be positive");
        return NULL;
    }

    PyArrayObject *field_arr =
        (PyArrayObject *) PyArray_FromAny(field_obj, NULL, 0, 0, NPY_ARRAY_CARRAY_RO, NULL);
    PyArrayObject *indices_arr =
        (PyArrayObject *) PyArray_FROM_OTF(indices_obj, NPY_INTP, NPY_ARRAY_CARRAY_RO);

    if (field_arr == NULL || indices_arr == NULL) goto fail;

    if (PyArray_NDIM(field_arr) < 2 || PyArray_NDIM(indices_arr) != 1)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "nearest_regrid_apply expects field ndim >= 2 and a 1D indices array");
        goto fail;
    }

    int field_ndim = PyArray_NDIM(field_arr);
    npy_intp nsrc_lon = PyArray_DIM(field_arr, field_ndim - 2);
    npy_intp nsrc_lat = PyArray_DIM(field_arr, field_ndim - 1);
    npy_intp ntgt_lon = (npy_intp) ntgt_lon_ssize;
    npy_intp ntgt_lat = (npy_intp) ntgt_lat_ssize;
    npy_intp target_size = ntgt_lon * ntgt_lat;
    npy_intp source_size = nsrc_lon * nsrc_lat;
    npy_intp indices_size = PyArray_DIM(indices_arr, 0);

    if (indices_size != target_size)
    {
        PyErr_SetString(
            PyExc_ValueError,
            "nearest_regrid_apply indices length must match target longitude/latitude sizes");
        goto fail;
    }

    npy_intp out_dims[NPY_MAXDIMS];
    if (field_ndim > NPY_MAXDIMS)
    {
        PyErr_SetString(PyExc_ValueError, "field has too many dimensions");
        goto fail;
    }
    for (int i = 0; i < field_ndim - 2; ++i)
        out_dims[i] = PyArray_DIM(field_arr, i);
    out_dims[field_ndim - 2] = ntgt_lon;
    out_dims[field_ndim - 1] = ntgt_lat;

    PyArray_Descr *field_descr = PyArray_DESCR(field_arr);
    Py_INCREF(field_descr);
    PyArrayObject *out_arr = (PyArrayObject *) PyArray_SimpleNewFromDescr(field_ndim, out_dims, field_descr);
    if (out_arr == NULL) goto fail;

    npy_intp outer_size = 1;
    for (int i = 0; i < field_ndim - 2; ++i)
        outer_size *= PyArray_DIM(field_arr, i);

    char *field_data = PyArray_BYTES(field_arr);
    char *out_data = PyArray_BYTES(out_arr);
    const npy_intp *indices = (const npy_intp *) PyArray_DATA(indices_arr);
    npy_intp itemsize = PyArray_ITEMSIZE(field_arr);

    for (npy_intp i = 0; i < indices_size; ++i)
    {
        if (indices[i] < 0 || indices[i] >= source_size)
        {
            PyErr_SetString(PyExc_ValueError, "nearest_regrid_apply indices are out of bounds");
            Py_DECREF(out_arr);
            goto fail;
        }
    }

    NPY_BEGIN_ALLOW_THREADS;
    for (npy_intp outer = 0; outer < outer_size; ++outer)
    {
        char *field_slice = field_data + outer * source_size * itemsize;
        char *out_slice = out_data + outer * target_size * itemsize;

        for (npy_intp itgt = 0; itgt < target_size; ++itgt)
        {
            memcpy(
                out_slice + itgt * itemsize,
                field_slice + indices[itgt] * itemsize,
                (size_t) itemsize);
        }
    }
    NPY_END_ALLOW_THREADS;

    Py_DECREF(field_arr);
    Py_DECREF(indices_arr);
    return PyArray_Return(out_arr);

fail:
    Py_XDECREF(field_arr);
    Py_XDECREF(indices_arr);
    return NULL;
}

static PyMethodDef nearest_methods[] = {
    { "bilinear_regrid", bilinear_regrid, METH_VARARGS, "Bilinear interpolation on a regular lon/lat grid." },
    { "bilinear_regrid_nd", bilinear_regrid_nd, METH_VARARGS, "Bilinear interpolation on regular lon/lat grids with arbitrary leading dimensions." },
    { "conservative_latitude_weights", conservative_latitude_weights, METH_VARARGS, "Conservative latitude weights on a regular latitude axis." },
    { "conservative_longitude_weights", conservative_longitude_weights, METH_VARARGS, "Conservative longitude weights on a regular longitude axis." },
    { "conservative_regrid", conservative_regrid, METH_VARARGS, "NaN-aware conservative regridding on regular lon/lat grids." },
    { "nearest_regrid_apply", nearest_regrid_apply, METH_VARARGS, "Apply cached nearest-neighbor indices across 2D/3D/4D regular-grid fields." },
    { "nearest_neighbor_indices", nearest_neighbor_indices, METH_VARARGS, "Compute nearest-neighbor source indices for regular lon/lat grids." },
    { NULL, NULL, 0, NULL },
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "regrid",
    "Native nearest-neighbor regridding helpers for regular rectilinear grids.",
    -1,
    nearest_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyMODINIT_FUNC
PyInit_regrid(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
