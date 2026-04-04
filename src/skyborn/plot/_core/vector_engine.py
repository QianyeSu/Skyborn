"""Pure vector-engine helpers extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np


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
