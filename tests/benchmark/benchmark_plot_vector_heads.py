"""Benchmarks for Python-vs-native curly-vector head geometry paths.

This script measures accuracy, timing, and RSS behavior for both the open
``"->"`` and filled ``"-|>"`` head paths under:

1. the scalar Python helper/reference path
2. the native batched helper path
"""

from __future__ import annotations

import gc
import json
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import psutil
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from skyborn.plot import curly_vector
from skyborn.plot import vector as vector_module
from skyborn.plot._core.result import CurlyVectorPlotSet
from skyborn.plot._core.vector_engine import (
    _curly_vector_ncl_impl,
    _curve_length_from_magnitude,
    _resolve_curly_anchor,
    _resolve_ncl_length_scale,
    _resolve_ncl_reference_length_px,
)
from skyborn.plot._shared.style import _ncl_arrow_edge_size_px, _resolve_open_arrow_size
from tests.plot_reference_geometry import (
    filled_arrow_geometry_reference,
    open_arrow_geometry_reference,
)

DATA_PATH = r"M:\CESM2LE\model\model_1011_001\Month\model_1011_001_Month_1850.nc"


@dataclass
class CurveSample:
    """One traced curly-vector shaft plus the head size used for it."""

    curve: np.ndarray
    display_curve: np.ndarray
    head_length_px: float
    head_width_px: float


def _rss_mb() -> float:
    """Return current process RSS in MiB."""

    return float(psutil.Process().memory_info().rss) / (1024.0 * 1024.0)


def _collect_curve_samples() -> tuple[list[CurveSample], dict[str, float]]:
    """Build a representative list of real CESM open-arrow shaft curves."""

    ds = xr.open_dataset(DATA_PATH)
    try:
        u = (
            ds["u"]
            .sel(level=200.0)
            .mean("time")
            .isel(lat=slice(None, None, 6), lon=slice(None, None, 6))
        )
        v = (
            ds["v"]
            .sel(level=200.0)
            .mean("time")
            .isel(lat=slice(None, None, 6), lon=slice(None, None, 6))
        )

        fig = plt.figure(figsize=(10, 4))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
        ax.set_global()
        ax.coastlines(linewidth=0.4)
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)

        x_values, y_values, u_values, v_values, _, _ = (
            vector_module._normalize_regular_grid_orientation(
                u["lon"].values,
                u["lat"].values,
                u.values,
                v.values,
                color=None,
                linewidth=None,
            )
        )
        grid = vector_module.Grid(x_values, y_values, allow_non_uniform=False)
        u_values = np.asarray(
            np.ma.masked_invalid(u_values).filled(np.nan), dtype=float
        )
        v_values = np.asarray(
            np.ma.masked_invalid(v_values).filled(np.nan), dtype=float
        )
        magnitude = np.hypot(u_values, v_values)
        valid_magnitude = magnitude[np.isfinite(magnitude)]

        display_sampler = vector_module._prepare_ncl_display_sampler(grid, transform)
        native_trace_context = vector_module._prepare_ncl_native_trace_context(
            grid,
            u_values,
            v_values,
            ax.bbox,
            display_sampler,
        )
        selected_centers = vector_module._select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=transform,
            axes=ax,
            density=1.0,
            start_points=None,
            min_distance=None,
            display_sampler=display_sampler,
            ncl_preset=None,
        )

        min_mag = float(np.min(valid_magnitude))
        max_mag = float(np.max(valid_magnitude))
        default_max_length_px = max(ax.bbox.width * 0.06, 1.0)
        ref_length_px = _resolve_ncl_reference_length_px(
            min_mag=min_mag,
            max_mag=max_mag,
            ref_mag=0.0,
            requested_ref_length_px=0.0,
            min_frac_length=0.0,
            default_max_length_px=default_max_length_px,
        )
        length_scale = _resolve_ncl_length_scale(
            min_mag=min_mag,
            max_mag=max_mag,
            ref_mag=0.0,
            requested_ref_length_px=0.0,
            min_frac_length=0.0,
            default_max_length_px=default_max_length_px,
        )
        step_px = max(1.5, ax.bbox.width * 0.0045)
        arrow_min_edge_px = max(ax.bbox.width * 0.003 * 1.0, 1.2)
        arrow_max_edge_px = max(ax.bbox.width * 0.012 * 1.0, arrow_min_edge_px)
        anchor = _resolve_curly_anchor(None, "both")

        samples: list[CurveSample] = []
        for center, center_mag in selected_centers:
            target_length_px = _curve_length_from_magnitude(center_mag, length_scale)
            curve_result = vector_module._build_ncl_curve(
                start_point=np.asarray(center, dtype=float),
                total_length_px=target_length_px,
                anchor=anchor,
                grid=grid,
                u=u_values,
                v=v_values,
                transform=transform,
                step_px=step_px,
                speed_scale=max_mag,
                viewport=ax.bbox,
                display_sampler=display_sampler,
                native_trace_context=native_trace_context,
            )
            if curve_result is None:
                continue
            curve, display_curve = curve_result
            if len(curve) < 2:
                continue

            head_length_px, head_width_px = _resolve_open_arrow_size(
                _ncl_arrow_edge_size_px(
                    center_mag,
                    max_mag=max_mag,
                    min_edge_px=arrow_min_edge_px,
                    max_edge_px=arrow_max_edge_px,
                )
            )
            samples.append(
                CurveSample(
                    curve=np.asarray(curve, dtype=float),
                    display_curve=np.asarray(display_curve, dtype=float),
                    head_length_px=float(head_length_px),
                    head_width_px=float(head_width_px),
                )
            )
        plt.close(fig)
    finally:
        ds.close()

    summary = {
        "curve_count": float(len(samples)),
        "ref_length_px": float(ref_length_px),
        "step_px": float(step_px),
        "max_mag": float(max_mag),
    }
    return samples, summary


def _benchmark_time(fn, *, loops: int, warmup: int = 3) -> dict[str, float]:
    """Time a callable repeatedly and summarize per-loop timings in milliseconds."""

    for _ in range(warmup):
        fn()

    timings = []
    for _ in range(loops):
        start = time.perf_counter()
        fn()
        timings.append((time.perf_counter() - start) * 1000.0)
    return {
        "loops": float(loops),
        "min_ms": float(min(timings)),
        "median_ms": float(statistics.median(timings)),
        "mean_ms": float(statistics.mean(timings)),
        "max_ms": float(max(timings)),
    }


def _measure_rss_stability(
    fn, *, iterations: int, collect_each: bool, warmup: int = 3
) -> dict[str, float]:
    """Measure RSS drift and range during repeated execution."""

    rss_values = []
    gc.collect()
    for _ in range(max(int(warmup), 0)):
        fn()
        if collect_each:
            gc.collect()
    gc.collect()
    start_rss = _rss_mb()
    for _ in range(iterations):
        fn()
        if collect_each:
            gc.collect()
        rss_values.append(_rss_mb())
    end_rss = rss_values[-1] if rss_values else start_rss
    return {
        "iterations": float(iterations),
        "start_rss_mb": float(start_rss),
        "min_rss_mb": float(min(rss_values, default=start_rss)),
        "max_rss_mb": float(max(rss_values, default=start_rss)),
        "end_rss_mb": float(end_rss),
        "drift_mb": float(end_rss - start_rss),
        "range_mb": float(
            max(rss_values, default=start_rss) - min(rss_values, default=start_rss)
        ),
    }


def main() -> None:
    """Run the baseline benchmark suite and print a JSON report."""

    samples, sample_summary = _collect_curve_samples()
    if not samples:
        raise RuntimeError("failed to collect any open-arrow curve samples")

    helper_fig = plt.figure(figsize=(10, 4))
    helper_ax = helper_fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    helper_ax.set_global()
    helper_transform = helper_ax.transData
    helper_inverse_transform = helper_transform.inverted()

    max_display_abs_diff = 0.0
    max_segment_abs_diff = 0.0
    max_filled_abs_diff = 0.0
    segment_count = 0
    for sample in samples:
        helper_geom = vector_module._open_arrow_geometry(
            sample.curve,
            helper_transform,
            sample.head_length_px,
            sample.head_width_px,
            display_curve=sample.display_curve,
        )
        ref_geom = open_arrow_geometry_reference(
            sample.display_curve,
            sample.head_length_px,
            sample.head_width_px,
        )
        if helper_geom is None or ref_geom is None:
            raise RuntimeError(
                "unexpected None geometry during baseline accuracy check"
            )

        for key in (
            "tip_display",
            "base_center_display",
            "left_display",
            "right_display",
        ):
            diff = np.max(
                np.abs(np.asarray(helper_geom[key]) - np.asarray(ref_geom[key]))
            )
            max_display_abs_diff = max(max_display_abs_diff, float(diff))

        helper_segments = vector_module._build_open_arrow_segments(
            curve=sample.curve,
            grid=None,
            transform=helper_transform,
            head_length_px=sample.head_length_px,
            head_width_px=sample.head_width_px,
            display_curve=sample.display_curve,
            inverse_transform=helper_inverse_transform,
        )
        display_vertices = np.vstack(
            [
                ref_geom["left_display"],
                ref_geom["tip_display"],
                ref_geom["right_display"],
            ]
        )
        ref_vertices = helper_inverse_transform.transform(display_vertices)
        ref_segments = [
            np.vstack([ref_vertices[0], ref_vertices[1]]),
            np.vstack([ref_vertices[2], ref_vertices[1]]),
        ]
        if len(helper_segments) != 2:
            raise RuntimeError("unexpected helper segment count for open arrow")
        for helper_seg, ref_seg in zip(helper_segments, ref_segments, strict=True):
            segment_count += 1
            diff = np.max(np.abs(np.asarray(helper_seg) - np.asarray(ref_seg)))
            max_segment_abs_diff = max(max_segment_abs_diff, float(diff))

        filled_polygon = vector_module._build_arrow_polygon(
            sample.curve,
            None,
            helper_transform,
            sample.head_length_px,
            sample.head_width_px,
            "k",
            "k",
            1.0,
            1.0,
            1.0,
            display_curve=sample.display_curve,
            inverse_transform=helper_inverse_transform,
        )
        ref_filled = filled_arrow_geometry_reference(
            sample.display_curve,
            sample.head_length_px,
            sample.head_width_px,
        )
        if filled_polygon is None or ref_filled is None:
            raise RuntimeError(
                "unexpected None geometry during filled-arrow accuracy check"
            )
        polygon_vertices = np.asarray(filled_polygon.get_xy()[:3], dtype=float)
        ref_polygon = helper_inverse_transform.transform(
            np.vstack(
                [
                    ref_filled["tip_display"],
                    ref_filled["left_display"],
                    ref_filled["right_display"],
                ]
            )
        )
        diff = np.max(np.abs(polygon_vertices - ref_polygon))
        max_filled_abs_diff = max(max_filled_abs_diff, float(diff))

    def _run_open_arrow_geometry_batch() -> None:
        for sample in samples:
            vector_module._open_arrow_geometry(
                sample.curve,
                helper_transform,
                sample.head_length_px,
                sample.head_width_px,
                display_curve=sample.display_curve,
            )

    def _run_open_arrow_segments_python_batch() -> None:
        for sample in samples:
            vector_module._build_open_arrow_segments(
                curve=sample.curve,
                grid=None,
                transform=helper_transform,
                head_length_px=sample.head_length_px,
                head_width_px=sample.head_width_px,
                display_curve=sample.display_curve,
                inverse_transform=helper_inverse_transform,
            )

    def _run_open_arrow_segments_native_batch() -> None:
        vector_module._build_open_arrow_segments_batch(
            transform=helper_transform,
            display_curves=[sample.display_curve for sample in samples],
            head_lengths_px=np.asarray(
                [sample.head_length_px for sample in samples], dtype=float
            ),
            head_widths_px=np.asarray(
                [sample.head_width_px for sample in samples], dtype=float
            ),
            inverse_transform=helper_inverse_transform,
        )

    def _run_filled_arrow_polygons_python_batch() -> None:
        for sample in samples:
            vector_module._build_arrow_polygon(
                sample.curve,
                None,
                helper_transform,
                sample.head_length_px,
                sample.head_width_px,
                "k",
                "k",
                1.0,
                1.0,
                1.0,
                display_curve=sample.display_curve,
                inverse_transform=helper_inverse_transform,
            )

    def _run_filled_arrow_polygons_native_batch() -> None:
        vector_module._build_filled_arrow_polygons_batch(
            transform=helper_transform,
            display_curves=[sample.display_curve for sample in samples],
            head_lengths_px=np.asarray(
                [sample.head_length_px for sample in samples], dtype=float
            ),
            head_widths_px=np.asarray(
                [sample.head_width_px for sample in samples], dtype=float
            ),
            inverse_transform=helper_inverse_transform,
        )

    ds = xr.open_dataset(DATA_PATH)
    try:
        u = (
            ds["u"]
            .sel(level=200.0)
            .mean("time")
            .isel(lat=slice(None, None, 6), lon=slice(None, None, 6))
        )
        v = (
            ds["v"]
            .sel(level=200.0)
            .mean("time")
            .isel(lat=slice(None, None, 6), lon=slice(None, None, 6))
        )
        x_render, y_render, u_render, v_render, _, _ = (
            vector_module._normalize_regular_grid_orientation(
                u["lon"].values,
                u["lat"].values,
                u.values,
                v.values,
                color=None,
                linewidth=None,
            )
        )

        def _render_internal(arrowstyle: str, *, use_native_batch: bool):
            fig_local = plt.figure(figsize=(10, 4))
            ax_local = fig_local.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            ax_local.set_global()
            ax_local.coastlines(linewidth=0.4)
            transform = ccrs.PlateCarree()._as_mpl_transform(ax_local)
            result = _curly_vector_ncl_impl(
                axes=ax_local,
                x=x_render,
                y=y_render,
                u=u_render,
                v=v_render,
                density=1.0,
                color="k",
                arrowsize=1.0,
                arrowstyle=arrowstyle,
                transform=transform,
                grid_cls=vector_module.Grid,
                prepare_ncl_display_sampler_fn=vector_module._prepare_ncl_display_sampler,
                prepare_ncl_native_trace_context_fn=vector_module._prepare_ncl_native_trace_context,
                select_ncl_centers_fn=vector_module._select_ncl_centers,
                build_ncl_curve_fn=vector_module._build_ncl_curve,
                sample_grid_field_fn=vector_module._sample_grid_field,
                build_ncl_arrow_artists_fn=vector_module._build_ncl_arrow_artists,
                build_open_arrow_segments_batch_fn=(
                    vector_module._build_open_arrow_segments_batch
                    if use_native_batch and arrowstyle == "->"
                    else None
                ),
                build_filled_arrow_polygons_batch_fn=(
                    vector_module._build_filled_arrow_polygons_batch
                    if use_native_batch and arrowstyle == "-|>"
                    else None
                ),
                display_points_to_data_fn=vector_module._display_points_to_data,
                result_cls=CurlyVectorPlotSet,
            )
            return fig_local, result

        def _compare_render_geometry(arrowstyle: str) -> dict[str, float]:
            fig_py, result_py = _render_internal(arrowstyle, use_native_batch=False)
            fig_native, result_native = _render_internal(
                arrowstyle, use_native_batch=True
            )
            try:
                line_py = [
                    np.asarray(seg, dtype=float)
                    for seg in result_py.lines.get_segments()
                ]
                line_native = [
                    np.asarray(seg, dtype=float)
                    for seg in result_native.lines.get_segments()
                ]
                if len(line_py) != len(line_native):
                    raise RuntimeError(
                        "line segment shapes differ between Python and native renders"
                    )
                line_diff = 0.0
                for py_seg, native_seg in zip(line_py, line_native, strict=True):
                    if py_seg.shape != native_seg.shape:
                        raise RuntimeError(
                            "individual line segment shapes differ between Python and native renders"
                        )
                    if py_seg.size:
                        line_diff = max(
                            line_diff,
                            float(np.max(np.abs(py_seg - native_seg))),
                        )

                if arrowstyle == "-|>":
                    polys_py = [
                        np.asarray(patch.get_xy()[:3], dtype=float)
                        for patch in result_py.arrows
                    ]
                    polys_native = []
                    for artist in result_native.arrows:
                        if hasattr(artist, "get_paths"):
                            for path in artist.get_paths():
                                polys_native.append(
                                    np.asarray(path.vertices[:3], dtype=float)
                                )
                        else:
                            polys_native.append(
                                np.asarray(artist.get_xy()[:3], dtype=float)
                            )
                    if len(polys_py) != len(polys_native):
                        raise RuntimeError(
                            "filled-arrow polygon counts differ between Python and native renders"
                        )
                    if polys_py:
                        poly_diff = max(
                            float(np.max(np.abs(py_poly - native_poly)))
                            for py_poly, native_poly in zip(
                                polys_py, polys_native, strict=True
                            )
                        )
                    else:
                        poly_diff = 0.0
                else:
                    poly_diff = 0.0
                return {
                    "line_segment_count": float(len(line_py)),
                    "line_segment_max_abs_diff": line_diff,
                    "filled_polygon_count": float(len(result_native.arrows)),
                    "filled_polygon_max_abs_diff": poly_diff,
                }
            finally:
                plt.close(fig_py)
                plt.close(fig_native)

        def _run_full_render_open() -> None:
            fig_local = plt.figure(figsize=(10, 4))
            ax_local = fig_local.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            ax_local.set_global()
            ax_local.coastlines(linewidth=0.4)
            curly_vector(
                ax_local,
                u["lon"].values,
                u["lat"].values,
                u.values,
                v.values,
                transform=ccrs.PlateCarree(),
                density=1.0,
                color="k",
                arrowsize=1.0,
                arrowstyle="->",
            )
            plt.close(fig_local)

        def _run_full_render_filled() -> None:
            fig_local = plt.figure(figsize=(10, 4))
            ax_local = fig_local.add_subplot(1, 1, 1, projection=ccrs.Robinson())
            ax_local.set_global()
            ax_local.coastlines(linewidth=0.4)
            curly_vector(
                ax_local,
                u["lon"].values,
                u["lat"].values,
                u.values,
                v.values,
                transform=ccrs.PlateCarree(),
                density=1.0,
                color="k",
                arrowsize=1.0,
                arrowstyle="-|>",
            )
            plt.close(fig_local)

        def _run_full_render_open_python() -> None:
            fig_local, _ = _render_internal("->", use_native_batch=False)
            plt.close(fig_local)

        def _run_full_render_open_native() -> None:
            fig_local, _ = _render_internal("->", use_native_batch=True)
            plt.close(fig_local)

        def _run_full_render_filled_python() -> None:
            fig_local, _ = _render_internal("-|>", use_native_batch=False)
            plt.close(fig_local)

        def _run_full_render_filled_native() -> None:
            fig_local, _ = _render_internal("-|>", use_native_batch=True)
            plt.close(fig_local)

        timing_report = {
            "open_arrow_geometry_batch": _benchmark_time(
                _run_open_arrow_geometry_batch,
                loops=20,
            ),
            "open_arrow_segments_python_batch": _benchmark_time(
                _run_open_arrow_segments_python_batch,
                loops=20,
            ),
            "open_arrow_segments_native_batch": _benchmark_time(
                _run_open_arrow_segments_native_batch,
                loops=20,
            ),
            "filled_arrow_polygons_python_batch": _benchmark_time(
                _run_filled_arrow_polygons_python_batch,
                loops=20,
            ),
            "filled_arrow_polygons_native_batch": _benchmark_time(
                _run_filled_arrow_polygons_native_batch,
                loops=20,
            ),
            "full_curly_vector_render_open": _benchmark_time(
                _run_full_render_open,
                loops=12,
            ),
            "full_curly_vector_render_filled": _benchmark_time(
                _run_full_render_filled,
                loops=12,
            ),
            "full_curly_vector_render_open_python": _benchmark_time(
                _run_full_render_open_python,
                loops=12,
            ),
            "full_curly_vector_render_open_native": _benchmark_time(
                _run_full_render_open_native,
                loops=12,
            ),
            "full_curly_vector_render_filled_python": _benchmark_time(
                _run_full_render_filled_python,
                loops=12,
            ),
            "full_curly_vector_render_filled_native": _benchmark_time(
                _run_full_render_filled_native,
                loops=12,
            ),
        }
        memory_report = {
            "open_arrow_segments_python_batch": _measure_rss_stability(
                _run_open_arrow_segments_python_batch,
                iterations=40,
                collect_each=True,
            ),
            "open_arrow_segments_native_batch": _measure_rss_stability(
                _run_open_arrow_segments_native_batch,
                iterations=40,
                collect_each=True,
            ),
            "filled_arrow_polygons_python_batch": _measure_rss_stability(
                _run_filled_arrow_polygons_python_batch,
                iterations=40,
                collect_each=True,
            ),
            "filled_arrow_polygons_native_batch": _measure_rss_stability(
                _run_filled_arrow_polygons_native_batch,
                iterations=40,
                collect_each=True,
            ),
            "full_curly_vector_render_open": _measure_rss_stability(
                _run_full_render_open,
                iterations=18,
                collect_each=True,
            ),
            "full_curly_vector_render_filled": _measure_rss_stability(
                _run_full_render_filled,
                iterations=18,
                collect_each=True,
            ),
            "full_curly_vector_render_open_python": _measure_rss_stability(
                _run_full_render_open_python,
                iterations=18,
                collect_each=True,
            ),
            "full_curly_vector_render_open_native": _measure_rss_stability(
                _run_full_render_open_native,
                iterations=18,
                collect_each=True,
            ),
            "full_curly_vector_render_filled_python": _measure_rss_stability(
                _run_full_render_filled_python,
                iterations=18,
                collect_each=True,
            ),
            "full_curly_vector_render_filled_native": _measure_rss_stability(
                _run_full_render_filled_native,
                iterations=18,
                collect_each=True,
            ),
        }
        render_accuracy = {
            "open": _compare_render_geometry("->"),
            "filled": _compare_render_geometry("-|>"),
        }
    finally:
        ds.close()
        plt.close(helper_fig)

    report = {
        "benchmark_case": {
            "dataset": DATA_PATH,
            "field": "annual-mean CESM2LE 200 hPa u/v",
            "projection": "Robinson",
            "arrowstyle": "->",
            "stride_lat_lon": 6.0,
            "density": 1.0,
        },
        "sample_summary": sample_summary,
        "accuracy": {
            "curve_count": float(len(samples)),
            "segment_count": float(segment_count),
            "max_display_abs_diff": float(max_display_abs_diff),
            "max_segment_abs_diff": float(max_segment_abs_diff),
            "max_filled_abs_diff": float(max_filled_abs_diff),
        },
        "render_accuracy": render_accuracy,
        "timings_ms_per_loop": timing_report,
        "memory_rss_mb": memory_report,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
