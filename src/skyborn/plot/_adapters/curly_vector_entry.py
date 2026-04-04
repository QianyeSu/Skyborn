"""Internal entry orchestration for the public curly-vector facade."""

from __future__ import annotations

from collections.abc import Hashable
from typing import Any


def _prepare_curly_vector_dataset_inputs_impl(
    ax,
    x,
    y,
    u,
    v,
    transform,
    regrid_shape,
    curvilinear_regrid_shape,
    target_extent,
    color,
    linewidth,
    density,
    *,
    is_cartopy_crs_like_fn,
    maybe_as_scalar_field_fn,
    is_curvilinear_grid_fn,
    normalize_regrid_shape_fn,
    default_curvilinear_regrid_shape_fn,
    regrid_curvilinear_vectors_fn,
    regrid_cartopy_vectors_fn,
    default_cartopy_target_extent_fn,
    extract_regular_grid_from_regridded_vectors_fn,
):
    src_crs = transform if is_cartopy_crs_like_fn(transform) else None
    color_field = maybe_as_scalar_field_fn(color, u.shape)
    linewidth_field = maybe_as_scalar_field_fn(linewidth, u.shape)
    if color_field is not None:
        color = color_field
    if linewidth_field is not None:
        linewidth = linewidth_field

    if is_curvilinear_grid_fn(x, y):
        scalar_fields = []
        scalar_keys = []
        if color_field is not None:
            scalar_fields.append(color_field)
            scalar_keys.append("color")
        if linewidth_field is not None:
            scalar_fields.append(linewidth_field)
            scalar_keys.append("linewidth")

        target_shape = (
            normalize_regrid_shape_fn(curvilinear_regrid_shape)
            if curvilinear_regrid_shape is not None
            else default_curvilinear_regrid_shape_fn(x, density)
        )
        regrid_result = regrid_curvilinear_vectors_fn(
            x,
            y,
            u,
            v,
            *scalar_fields,
            target_shape=target_shape,
        )
        x, y, u, v, *scalar_grids = regrid_result
        for key, scalar_grid in zip(scalar_keys, scalar_grids):
            if key == "color":
                color = scalar_grid
            elif key == "linewidth":
                linewidth = scalar_grid

    if regrid_shape is None:
        if src_crs is not None:
            transform = src_crs._as_mpl_transform(ax)
        return x, y, u, v, color, linewidth, transform

    if src_crs is None:
        raise ValueError("regrid_shape requires a Cartopy CRS passed via transform")
    if not hasattr(ax, "projection"):
        raise ValueError("regrid_shape requires a Cartopy GeoAxes with a projection")

    scalar_fields = []
    scalar_keys = []
    color_field = maybe_as_scalar_field_fn(color, u.shape)
    linewidth_field = maybe_as_scalar_field_fn(linewidth, u.shape)
    if color_field is not None:
        scalar_fields.append(color_field)
        scalar_keys.append("color")
    if linewidth_field is not None:
        scalar_fields.append(linewidth_field)
        scalar_keys.append("linewidth")

    regrid_result = regrid_cartopy_vectors_fn(
        src_crs=src_crs,
        target_proj=ax.projection,
        regrid_shape=regrid_shape,
        x=x,
        y=y,
        u=u,
        v=v,
        *scalar_fields,
        target_extent=default_cartopy_target_extent_fn(ax, target_extent),
    )

    x_grid, y_grid, u_grid, v_grid, *scalar_grids = regrid_result
    x, y, u, v, scalar_grids = extract_regular_grid_from_regridded_vectors_fn(
        x_grid,
        y_grid,
        u_grid,
        v_grid,
        scalar_grids,
    )

    for key, scalar_grid in zip(scalar_keys, scalar_grids):
        if key == "color":
            color = scalar_grid
        elif key == "linewidth":
            linewidth = scalar_grid

    return x, y, u, v, color, linewidth, None


def _curly_vector_from_dataset_impl(
    ds,
    x: Hashable,
    y: Hashable,
    u: Hashable,
    v: Hashable,
    ax=None,
    density: Any = 1,
    linewidth: Any = None,
    linewidths: Any = None,
    color: Any = None,
    c: Any = None,
    cmap: Any = None,
    norm: Any = None,
    vmin: float | None = None,
    vmax: float | None = None,
    alpha: float | None = None,
    facecolor: Any = None,
    facecolors: Any = None,
    edgecolor: Any = None,
    edgecolors: Any = None,
    rasterized: bool | None = None,
    arrowsize=1,
    arrowstyle="->",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
    pivot: str | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length=0.0,
    min_distance: float | None = None,
    ncl_preset: str | None = None,
    regrid_shape: Any = None,
    curvilinear_regrid_shape: Any = None,
    target_extent: Any = None,
    isel: Any = None,
    *,
    resolve_curly_style_aliases_fn,
    extract_curly_vector_dataset_source_fn,
    prepare_dataset_style_field_fn,
    gca_fn,
    prepare_curly_vector_dataset_inputs_fn,
    collect_named_kwargs_fn,
    array_curly_vector_kwarg_names,
    array_curly_vector_fn,
):
    color, linewidth, facecolor, edgecolor, vmin, vmax = resolve_curly_style_aliases_fn(
        color=color,
        c=c,
        linewidth=linewidth,
        linewidths=linewidths,
        facecolor=facecolor,
        facecolors=facecolors,
        edgecolor=edgecolor,
        edgecolors=edgecolors,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
    )

    x, y, u, v, source_metadata = extract_curly_vector_dataset_source_fn(
        ds,
        x,
        y,
        u,
        v,
        isel=isel,
        return_metadata=True,
    )
    color = prepare_dataset_style_field_fn(
        color,
        isel=isel,
        expected_shape=u.shape,
        vector_dims=source_metadata["vector_dims"],
        x_descending=bool(source_metadata["x_descending"]),
        y_descending=bool(source_metadata["y_descending"]),
        role="color",
    )
    linewidth = prepare_dataset_style_field_fn(
        linewidth,
        isel=isel,
        expected_shape=u.shape,
        vector_dims=source_metadata["vector_dims"],
        x_descending=bool(source_metadata["x_descending"]),
        y_descending=bool(source_metadata["y_descending"]),
        role="linewidth",
    )

    if ax is None:
        ax = gca_fn()
    x, y, u, v, color, linewidth, transform = prepare_curly_vector_dataset_inputs_fn(
        ax=ax,
        x=x,
        y=y,
        u=u,
        v=v,
        transform=transform,
        regrid_shape=regrid_shape,
        curvilinear_regrid_shape=curvilinear_regrid_shape,
        target_extent=target_extent,
        color=color,
        linewidth=linewidth,
        density=density,
    )

    scope = {
        "density": density,
        "linewidth": linewidth,
        "color": color,
        "vmin": vmin,
        "vmax": vmax,
        "cmap": cmap,
        "norm": norm,
        "alpha": alpha,
        "facecolor": facecolor,
        "edgecolor": edgecolor,
        "rasterized": rasterized,
        "arrowsize": arrowsize,
        "arrowstyle": arrowstyle,
        "transform": transform,
        "zorder": zorder,
        "start_points": start_points,
        "integration_direction": integration_direction,
        "grains": grains,
        "broken_streamlines": broken_streamlines,
        "anchor": anchor,
        "pivot": pivot,
        "ref_magnitude": ref_magnitude,
        "ref_length": ref_length,
        "min_frac_length": min_frac_length,
        "min_distance": min_distance,
        "ncl_preset": ncl_preset,
    }
    return array_curly_vector_fn(
        ax,
        x,
        y,
        u,
        v,
        **collect_named_kwargs_fn(scope, array_curly_vector_kwarg_names),
    )


def _curly_vector_from_arrays_impl(
    ax,
    x: Any,
    y: Any,
    u: Any,
    v: Any,
    density: Any = 1,
    linewidth: Any = None,
    linewidths: Any = None,
    color: Any = None,
    c: Any = None,
    cmap: Any = None,
    norm: Any = None,
    vmin: float | None = None,
    vmax: float | None = None,
    alpha: float | None = None,
    facecolor: Any = None,
    facecolors: Any = None,
    edgecolor: Any = None,
    edgecolors: Any = None,
    rasterized: bool | None = None,
    arrowsize=1,
    arrowstyle="->",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
    pivot: str | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length=0.0,
    min_distance: float | None = None,
    ncl_preset: str | None = None,
    regrid_shape: Any = None,
    curvilinear_regrid_shape: Any = None,
    target_extent: Any = None,
    *,
    resolve_curly_style_aliases_fn,
    asarray_fn,
    prepare_curly_vector_dataset_inputs_fn,
    collect_named_kwargs_fn,
    array_curly_vector_kwarg_names,
    array_curly_vector_fn,
):
    color, linewidth, facecolor, edgecolor, vmin, vmax = resolve_curly_style_aliases_fn(
        color=color,
        c=c,
        linewidth=linewidth,
        linewidths=linewidths,
        facecolor=facecolor,
        facecolors=facecolors,
        edgecolor=edgecolor,
        edgecolors=edgecolors,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
    )

    x, y, u, v, color, linewidth, transform = prepare_curly_vector_dataset_inputs_fn(
        ax=ax,
        x=asarray_fn(x),
        y=asarray_fn(y),
        u=asarray_fn(u),
        v=asarray_fn(v),
        transform=transform,
        regrid_shape=regrid_shape,
        curvilinear_regrid_shape=curvilinear_regrid_shape,
        target_extent=target_extent,
        color=color,
        linewidth=linewidth,
        density=density,
    )

    scope = {
        "density": density,
        "linewidth": linewidth,
        "color": color,
        "vmin": vmin,
        "vmax": vmax,
        "cmap": cmap,
        "norm": norm,
        "alpha": alpha,
        "facecolor": facecolor,
        "edgecolor": edgecolor,
        "rasterized": rasterized,
        "arrowsize": arrowsize,
        "arrowstyle": arrowstyle,
        "transform": transform,
        "zorder": zorder,
        "start_points": start_points,
        "integration_direction": integration_direction,
        "grains": grains,
        "broken_streamlines": broken_streamlines,
        "anchor": anchor,
        "pivot": pivot,
        "ref_magnitude": ref_magnitude,
        "ref_length": ref_length,
        "min_frac_length": min_frac_length,
        "min_distance": min_distance,
        "ncl_preset": ncl_preset,
    }
    return array_curly_vector_fn(
        ax,
        x,
        y,
        u,
        v,
        **collect_named_kwargs_fn(scope, array_curly_vector_kwarg_names),
    )
