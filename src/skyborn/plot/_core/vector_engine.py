"""Pure vector-engine helpers extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np

from .geometry import _point_within_grid_data
from .thinning import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
)


def _finite_plot_field_values(field, field_name):
    """Return finite values from a style field or fail fast with a clear error."""
    finite_values = np.asarray(field, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        raise ValueError(
            f"{field_name} field must contain at least one finite value after masking"
        )
    return finite_values


def _resolve_artist_coordinate_context(axes, transform):
    """Choose the final artist transform for the rendered curly-vector geometry."""
    artist_transform = transform
    artist_inverse_transform = None
    bake_display_geometry = bool(
        hasattr(axes, "projection")
        and transform is not None
        and transform is not axes.transData
    )
    if not bake_display_geometry:
        return artist_transform, artist_inverse_transform, False

    try:
        artist_inverse_transform = axes.transData.inverted()
    except Exception:
        return transform, None, False

    return axes.transData, artist_inverse_transform, True


def _default_ncl_max_length_px(axes_bbox, density):
    density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
    density_xy = np.maximum(density_xy, 0.1)
    nx = max(30.0 * density_xy[0], 1.0)
    ny = max(30.0 * density_xy[1], 1.0)
    sx = float(axes_bbox.width) / nx
    sy = float(axes_bbox.height) / ny
    return max(np.sqrt((sx * sx + sy * sy) / 2.0), 1.0)


def _resolve_ncl_reference_length_px(
    min_mag,
    max_mag,
    ref_mag,
    requested_ref_length_px,
    min_frac_length,
    default_max_length_px,
):
    min_mag = float(min_mag)
    max_mag = float(max_mag)
    ref_mag = float(ref_mag)
    requested_ref_length_px = float(requested_ref_length_px)
    min_frac_length = float(np.clip(min_frac_length, 0.0, 1.0))
    default_max_length_px = max(float(default_max_length_px), 1.0)

    if requested_ref_length_px > 0.0:
        return requested_ref_length_px
    if max_mag - min_mag <= 1e-12:
        return default_max_length_px
    if ref_mag > 0.0 and min_frac_length > 0.0:
        ratio = np.clip((ref_mag - min_mag) / max(max_mag - min_mag, 1e-12), 0.0, 1.0)
        denominator = 1.0 - min_frac_length + min_frac_length * ratio
        if denominator > 1e-12:
            return default_max_length_px * ratio / denominator
    if ref_mag > 0.0:
        return default_max_length_px * ref_mag / max(max_mag, 1e-12)
    return default_max_length_px


def _resolve_ncl_length_scale(
    min_mag,
    max_mag,
    ref_mag,
    requested_ref_length_px,
    min_frac_length,
    default_max_length_px,
):
    uvmn = float(min_mag)
    uvmx = max(float(max_mag), 1e-12)
    rvrm = float(ref_mag)
    rvrl = float(requested_ref_length_px)
    vfr = float(np.clip(min_frac_length, 0.0, 1.0))
    rdmx = max(float(default_max_length_px), 1.0)
    rdmn = 0.0
    vrl = rdmx
    vfl = rdmn
    iav = 0

    if uvmx - uvmn <= 1e-12:
        if rvrl > 0.0 and rvrm > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
        elif rvrm > 0.0:
            vrl = rdmx * rvrm / uvmx
        elif rvrl > 0.0:
            rdmx = rvrl
            vrl = rdmx
        vfl = vrl
        rdmn = vfl
    elif rvrm <= 0.0:
        if rvrl > 0.0:
            rdmx = rvrl
        vrl = rdmx
        if vfr > 0.0:
            iav = 1
            vfl = vfr * rdmx
            rdmn = vfl
        else:
            rdmn = rdmx * (uvmn / uvmx)
            vfl = rdmn
    elif rvrm <= uvmn:
        iav = 1
        if rvrl > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
        elif vfr > 0.0:
            vrl = rdmx * vfr
        else:
            vrl = rdmx * rvrm / uvmx
        rdmn = vrl * uvmn / max(rvrm, 1e-12)
        vfl = rdmn
    elif vfr > 0.0:
        iav = 1
        if rvrl > 0.0:
            vrl = rvrl
            vfl = vfr * vrl
            rdmn = vfl
            rdmx = rdmn + (vrl - rdmn) * (uvmx - uvmn) / max(rvrm - uvmn, 1e-12)
        else:
            ratio = (rvrm - uvmn) / max(uvmx - uvmn, 1e-12)
            denominator = 1.0 - vfr + vfr * ratio
            if denominator > 1e-12:
                vrl = rdmx * ratio / denominator
            vfl = vfr * vrl
            rdmn = vfl
    else:
        if rvrl > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
            vfl = vrl * uvmn / max(rvrm, 1e-12)
            rdmn = vfl
        else:
            vrl = rdmx * rvrm / uvmx
            vfl = rdmx * uvmn / uvmx
            rdmn = vfl

    return {
        "min_mag": uvmn,
        "max_mag": uvmx,
        "ref_mag": rvrm,
        "ref_length_px": vrl,
        "min_length_px": rdmn,
        "max_length_px": rdmx,
        "adjust_min": bool(iav),
        "min_frac_length": vfr,
    }


def _curve_length_from_magnitude(magnitude_value, length_scale):
    magnitude_value = float(magnitude_value)
    min_mag = length_scale["min_mag"]
    max_mag = length_scale["max_mag"]
    max_length_px = length_scale["max_length_px"]
    min_length_px = length_scale["min_length_px"]

    if max_mag - min_mag <= 1e-12:
        return max_length_px

    if not length_scale["adjust_min"]:
        return np.clip(magnitude_value / max_mag, 0.0, 1.0) * max_length_px

    scale = np.clip(
        (magnitude_value - min_mag) / max(max_mag - min_mag, 1e-12), 0.0, 1.0
    )
    return min_length_px + (max_length_px - min_length_px) * scale


def _density_xy(density):
    density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
    return np.maximum(density_xy, 0.1)


def _density_scalar(density):
    if np.isscalar(density):
        return max(float(density), 0.1)
    return float(np.mean(_density_xy(density)))


def _select_ncl_centers(
    grid,
    magnitude,
    transform,
    axes,
    density,
    start_points,
    min_distance,
    display_sampler=None,
    ncl_preset=None,
    sample_grid_field_array=None,
    thin_ncl_mapped_candidates=None,
):
    candidates, candidate_magnitudes = _prepare_ncl_center_candidates(
        grid=grid,
        magnitude=magnitude,
        density=density,
        start_points=start_points,
        ncl_preset=ncl_preset,
        sample_grid_field_array=sample_grid_field_array,
    )
    if display_sampler is not None:
        display_points = display_sampler.sample_display_points(candidates)
    else:
        display_points = transform.transform(candidates)
    valid = _valid_ncl_center_candidates(
        grid=grid,
        candidates=candidates,
        candidate_magnitudes=candidate_magnitudes,
        display_points=display_points,
        start_points=start_points,
    )

    candidates = candidates[valid]
    display_points = display_points[valid]
    candidate_magnitudes = candidate_magnitudes[valid]
    if candidates.size == 0:
        return []

    mapped_points = _map_ncl_display_points_to_viewport(display_points, axes.bbox)
    spacing_frac = _resolve_ncl_min_distance_fraction(
        density=density,
        min_distance=min_distance,
        ncl_preset=ncl_preset,
    )
    selected_indices = thin_ncl_mapped_candidates(mapped_points, spacing_frac)
    return [(candidates[idx], candidate_magnitudes[idx]) for idx in selected_indices]


def _prepare_ncl_center_candidates(
    grid,
    magnitude,
    density,
    start_points,
    ncl_preset,
    sample_grid_field_array,
):
    if start_points is None:
        candidates = _default_ncl_box_center_candidates(
            grid,
            density=density,
            ncl_preset=ncl_preset,
        )
        candidate_magnitudes = sample_grid_field_array(grid, magnitude, candidates)
        return candidates, candidate_magnitudes

    candidates = np.asanyarray(start_points, dtype=float)
    if candidates.ndim != 2 or candidates.shape[1] != 2:
        raise ValueError("'start_points' must be an (N, 2) array")
    candidate_magnitudes = sample_grid_field_array(grid, magnitude, candidates)
    return candidates, candidate_magnitudes


def _default_ncl_candidate_shape(grid, density):
    density_xy = _density_xy(density)
    candidate_nx = min(grid.nx - 1, max(int(np.ceil(30.0 * density_xy[0])), 1))
    candidate_ny = min(grid.ny - 1, max(int(np.ceil(30.0 * density_xy[1])), 1))
    return candidate_nx, candidate_ny


def _default_ncl_box_center_candidates(grid, density=1, ncl_preset=None):
    del ncl_preset
    if grid.nx < 2 or grid.ny < 2:
        return np.empty((0, 2), dtype=float)

    candidate_nx, candidate_ny = _default_ncl_candidate_shape(grid, density)
    if candidate_nx == grid.nx - 1 and candidate_ny == grid.ny - 1:
        xs = grid.x_origin + (np.arange(grid.nx - 1) + 0.5) * grid.dx
        ys = grid.y_origin + (np.arange(grid.ny - 1) + 0.5) * grid.dy
    else:
        x_edges = np.linspace(
            grid.x_origin, grid.x_origin + grid.width, candidate_nx + 1
        )
        y_edges = np.linspace(
            grid.y_origin, grid.y_origin + grid.height, candidate_ny + 1
        )
        xs = 0.5 * (x_edges[:-1] + x_edges[1:])
        ys = 0.5 * (y_edges[:-1] + y_edges[1:])
    grid_x, grid_y = np.meshgrid(xs, ys, indexing="xy")
    return np.column_stack([grid_x.ravel(), grid_y.ravel()])


def _valid_ncl_center_candidates(
    grid, candidates, candidate_magnitudes, display_points, start_points
):
    valid = np.isfinite(display_points).all(axis=1) & np.isfinite(candidate_magnitudes)

    for idx, (xd, yd) in enumerate(candidates):
        if _point_within_grid_data(grid, np.array([xd, yd], dtype=float)):
            continue
        if start_points is not None:
            raise ValueError(f"Starting point ({xd}, {yd}) outside of data boundaries")
        valid[idx] = False

    return valid
