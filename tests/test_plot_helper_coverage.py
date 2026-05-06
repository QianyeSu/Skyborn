"""Focused regression tests for low-level plot helper branches."""

from __future__ import annotations

import builtins
import sys
import types
from functools import partial
from types import SimpleNamespace
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle
from matplotlib.transforms import Affine2D, Bbox

import skyborn.plot.vector as vector_module
from skyborn.plot._adapters.cartopy_vector import (
    _build_projection_target_grid,
    _extract_regular_grid_from_regridded_vectors,
    _regrid_cartopy_vectors,
)
from skyborn.plot._adapters.curly_vector_entry import (
    _prepare_curly_vector_dataset_inputs_impl,
)
from skyborn.plot._adapters.dataset_vector import (
    _extract_curly_vector_dataset_source,
    _get_plot_dataarray,
    _prepare_dataset_style_field,
    _transpose_2d_dataarray_to_dims,
)
from skyborn.plot._adapters.grid_prepare import (
    _build_curvilinear_target_grid,
    _maybe_as_scalar_field,
    _prepare_source_vector_grid,
    _rcm2rgrid_fields,
    _regrid_curvilinear_vectors,
    _wrap_periodic_grid_queries,
)
from skyborn.plot._artists.vector_artists import (
    _assemble_filled_head_artists,
    _assemble_open_head_streamlines,
    _build_filled_arrow_polygons_batch,
    _build_open_arrow_segments_batch,
)
from skyborn.plot._artists.vector_key_artist import CurlyVectorKey
from skyborn.plot._core.geometry import (
    _candidate_data_from_display_step,
    _display_points_to_data,
    _display_step_to_data,
    _display_to_data,
    _evaluate_ncl_display_curve,
    _finite_difference_step,
    _local_display_jacobian,
    _point_at_arc_distance_from_end,
    _tip_display_geometry,
    _tip_display_geometry_from_display_curve,
    _trim_display_curve_from_end,
)
from skyborn.plot._core.legacy_stream import DomainMap, InvalidIndexError
from skyborn.plot._core.native import (
    _call_native_sample_grid_field,
    _call_native_sample_grid_field_array,
    _call_native_thin_ncl_display_candidates,
    _call_native_thin_ncl_mapped_candidates,
    _call_native_trace_ncl_direction,
    _call_native_trace_ncl_direction_with_display,
    _call_native_validate_display_curve,
)
from skyborn.plot._core.result import CurlyVectorPlotSet
from skyborn.plot._core.scatter_engine import (
    _as_1d_axis,
    _generate_cell_candidates,
    _prepare_1d_grid_candidates,
    _resolve_placement,
    _scatter_impl,
    _subset_scatter_value,
)
from skyborn.plot._core.vector_engine import (
    Grid,
    _build_ncl_curve,
    _curly_vector_ncl_impl,
    _prepare_ncl_center_candidates,
    _resolve_artist_coordinate_context,
    _resolve_ncl_length_scale,
    _resolve_ncl_reference_length_px,
    _sample_local_vector_state,
    _select_ncl_centers,
    _trace_ncl_curve,
    _valid_ncl_center_candidates,
)
from skyborn.plot._shared.coords import (
    _axis_coordinate_1d,
    _axis_is_uniform,
    _cell_edges_from_axis,
    _center_grid_to_corner_grid,
    _coerce_matching_plot_field,
    _extract_meshgrid_axes,
    _filled_float_array,
    _normalize_coordinates,
    _normalize_regular_grid_orientation,
    _normalize_selection_mask,
    _rectilinear_cell_edges,
    _scatter_cell_geometry,
    _subset_ready_array,
    _transpose_dataarray_if_possible,
)
from skyborn.plot._shared.style import (
    _collect_named_kwargs,
    _curly_head_axes_dimensions,
    _normalize_artist_alpha,
    _normalize_curly_pivot,
    _normalize_supported_arrowstyle,
    _resolve_curly_anchor_alias,
    _resolve_curly_style_aliases,
    _resolve_scatter_aliases,
)
from skyborn.plot.nclcurly_native import (
    build_filled_arrow_polygons,
    build_open_arrow_segments,
)
from skyborn.plot.scatter import scatter as public_scatter
from skyborn.plot.vector import (
    _array_curly_vector,
    _regrid_non_uniform_vectors_to_uniform,
    curly_vector,
)
from tests.plot_reference_geometry import (
    filled_arrow_geometry_reference,
    open_arrow_geometry_reference,
)


class _BadTransform:
    """Transform stub that always fails."""

    def transform(self, values):
        raise RuntimeError("bad transform")

    def inverted(self):
        raise RuntimeError("bad inverse")


class _ExplodingInverse:
    """Transform stub whose inverse transform call fails."""

    def inverted(self):
        return self

    def transform(self, values):
        raise RuntimeError("explode inverse transform")


def _make_test_grid() -> Grid:
    return Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))


def _make_mock_key_set(*, arrowstyle: str = "->") -> Mock:
    plot_set = Mock(spec=CurlyVectorPlotSet)
    plot_set.resolution = 0.5
    plot_set.magnitude = np.array([[1.0, 2.0], [3.0, 4.0]])
    plot_set.max_magnitude = 10.0
    plot_set.ref_length_fraction = 0.12
    plot_set.arrowstyle = arrowstyle
    plot_set.arrowsize = 1.0
    plot_set.glyph_length_axes_fraction.side_effect = lambda value: float(value) / 100.0
    return plot_set


class TestCartopyVectorHelpers:
    def test_build_projection_target_grid_uses_projection_limits_by_default(self):
        proj = SimpleNamespace(x_limits=(-5.0, 5.0), y_limits=(2.0, 6.0))

        x_target, y_target = _build_projection_target_grid(proj, (3, 2), None)

        assert x_target.shape == (2, 3)
        assert y_target.shape == (2, 3)
        np.testing.assert_allclose(x_target[0], np.array([-5.0, 0.0, 5.0]))
        np.testing.assert_allclose(y_target[:, 0], np.array([2.0, 6.0]))

    def test_extract_regular_grid_from_regridded_vectors_flips_descending_y(self):
        x_grid = np.array([[10.0, 20.0], [10.0, 20.0]])
        y_grid = np.array([[5.0, 5.0], [-5.0, -5.0]])
        u_grid = np.array([[1.0, 2.0], [3.0, 4.0]])
        v_grid = np.array([[10.0, 20.0], [30.0, 40.0]])
        scalar = np.array([[100.0, 200.0], [300.0, 400.0]])

        x_1d, y_1d, u_out, v_out, [scalar_out] = (
            _extract_regular_grid_from_regridded_vectors(
                x_grid,
                y_grid,
                u_grid,
                v_grid,
                [scalar],
            )
        )

        np.testing.assert_allclose(x_1d, np.array([10.0, 20.0]))
        np.testing.assert_allclose(y_1d, np.array([-5.0, 5.0]))
        np.testing.assert_allclose(u_out, np.array([[3.0, 4.0], [1.0, 2.0]]))
        np.testing.assert_allclose(v_out, np.array([[30.0, 40.0], [10.0, 20.0]]))
        np.testing.assert_allclose(
            scalar_out, np.array([[300.0, 400.0], [100.0, 200.0]])
        )

    def test_regrid_cartopy_vectors_handles_invalid_points_and_scalar_fields(self):
        class FakeInterpolator:
            def __init__(self, axes, values, method, bounds_error, fill_value):
                self.values = np.asarray(values, dtype=float)
                self.axes = axes
                self.method = method
                self.bounds_error = bounds_error
                self.fill_value = fill_value

            def __call__(self, sample_points):
                sample_points = np.asarray(sample_points, dtype=float)
                return np.sum(sample_points, axis=1) + float(np.nanmean(self.values))

        fake_interp_module = types.ModuleType("scipy.interpolate")
        fake_interp_module.RegularGridInterpolator = FakeInterpolator
        fake_scipy = types.ModuleType("scipy")
        fake_scipy.interpolate = fake_interp_module

        class FakeSrcCrs:
            def transform_points(self, target_proj, x_target, y_target):
                source = np.stack(
                    [
                        x_target + 10.0,
                        y_target + 20.0,
                        np.zeros_like(x_target),
                    ],
                    axis=-1,
                )
                source[0, 1, :2] = np.nan
                return source

        class FakeTargetProj:
            x_limits = (0.0, 1.0)
            y_limits = (0.0, 1.0)

            def transform_vectors(self, src_crs, x, y, u, v):
                del src_crs, x, y
                return u + 100.0, v - 100.0

        x = np.array([0.0, 1.0])
        y = np.array([0.0, 1.0])
        u = np.array([[1.0, 2.0], [3.0, 4.0]])
        v = np.array([[5.0, 6.0], [7.0, 8.0]])
        scalar = np.array([[9.0, 10.0], [11.0, 12.0]])

        with patch.dict(
            sys.modules,
            {"scipy": fake_scipy, "scipy.interpolate": fake_interp_module},
        ):
            x_target, y_target, u_target, v_target, scalar_target = (
                _regrid_cartopy_vectors(
                    FakeSrcCrs(),
                    FakeTargetProj(),
                    (2, 2),
                    x,
                    y,
                    u,
                    v,
                    scalar,
                    target_extent=(0.0, 1.0, 0.0, 1.0),
                )
            )

        assert x_target.shape == (2, 2)
        assert y_target.shape == (2, 2)
        assert np.isnan(u_target[0, 1])
        assert np.isnan(v_target[0, 1])
        assert np.isnan(scalar_target[0, 1])
        np.testing.assert_allclose(u_target[0, 0], 132.5)
        np.testing.assert_allclose(v_target[1, 1], -61.5)
        np.testing.assert_allclose(scalar_target[1, 0], 41.5)

    def test_regrid_cartopy_vectors_requires_scipy(self):
        real_import = builtins.__import__

        def _raising_import(name, *args, **kwargs):
            if name == "scipy.interpolate":
                raise ImportError("missing scipy")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_raising_import):
            with pytest.raises(ImportError, match="scipy is required"):
                _regrid_cartopy_vectors(
                    src_crs=Mock(),
                    target_proj=SimpleNamespace(
                        x_limits=(0.0, 1.0), y_limits=(0.0, 1.0)
                    ),
                    regrid_shape=(2, 2),
                    x=np.array([0.0, 1.0]),
                    y=np.array([0.0, 1.0]),
                    u=np.ones((2, 2)),
                    v=np.ones((2, 2)),
                )


class TestGridPrepareAndEntryHelpers:
    def test_build_curvilinear_target_grid_rejects_all_nan_coordinates(self):
        x = np.full((2, 2), np.nan)
        y = np.full((2, 2), np.nan)

        with pytest.raises(ValueError, match="contain no finite points"):
            _build_curvilinear_target_grid(x, y, (3, 3))

    def test_prepare_source_vector_grid_normalizes_descending_axes_and_cyclic_scalars(
        self,
    ):
        x = np.array([240.0, 120.0, 0.0])
        y = np.array([10.0, 0.0])
        u = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v = np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]])
        scalar = np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]])

        x_out, y_out, u_out, v_out, [scalar_out] = _prepare_source_vector_grid(
            x,
            y,
            u,
            v,
            scalars=(scalar,),
        )

        np.testing.assert_allclose(x_out, np.array([0.0, 120.0, 240.0, 360.0]))
        np.testing.assert_allclose(y_out, np.array([0.0, 10.0]))
        np.testing.assert_allclose(u_out[:, -1], u_out[:, 0])
        np.testing.assert_allclose(v_out[:, -1], v_out[:, 0])
        np.testing.assert_allclose(scalar_out[:, -1], scalar_out[:, 0])

    def test_prepare_source_vector_grid_rejects_non_rectilinear_coordinates(self):
        x = np.ones((2, 2, 2))
        y = np.ones((2, 2))

        with pytest.raises(ValueError, match="1D or meshgrid"):
            _prepare_source_vector_grid(x, y, np.ones((2, 2)), np.ones((2, 2)))

    def test_wrap_periodic_grid_queries_preserves_non_cyclic_axes(self):
        x_query = np.array([-20.0, 10.0, 390.0])
        x_axis = np.array([0.0, 60.0, 120.0, 180.0])

        wrapped = _wrap_periodic_grid_queries(x_query, x_axis)

        np.testing.assert_allclose(wrapped, x_query)

    def test_prepare_curly_vector_dataset_inputs_impl_handles_curvilinear_and_cartopy_cases(
        self,
    ):
        ax = SimpleNamespace()
        x = np.arange(4.0).reshape(2, 2)
        y = (np.arange(4.0) + 10.0).reshape(2, 2)
        u = np.ones((2, 2))
        v = np.ones((2, 2)) * 2.0
        color = np.arange(4.0).reshape(2, 2)
        linewidth = np.arange(4.0).reshape(2, 2) + 0.5

        result = _prepare_curly_vector_dataset_inputs_impl(
            ax=ax,
            x=x,
            y=y,
            u=u,
            v=v,
            transform=SimpleNamespace(
                _as_mpl_transform=lambda current_ax: ("mpl", current_ax)
            ),
            regrid_shape=None,
            curvilinear_regrid_shape=None,
            target_extent=None,
            color=color,
            linewidth=linewidth,
            density=1.0,
            is_cartopy_crs_like_fn=lambda transform: True,
            maybe_as_scalar_field_fn=lambda value, shape: value,
            is_curvilinear_grid_fn=lambda x, y: True,
            normalize_regrid_shape_fn=lambda value: value,
            default_curvilinear_regrid_shape_fn=lambda x, density: (3, 2),
            regrid_curvilinear_vectors_fn=lambda *args, **kwargs: (
                np.array([10.0, 20.0, 30.0]),
                np.array([-5.0, 5.0]),
                np.ones((2, 3)),
                np.ones((2, 3)) * 2.0,
                np.ones((2, 3)) * 9.0,
                np.ones((2, 3)) * 7.0,
            ),
            regrid_cartopy_vectors_fn=lambda **kwargs: None,
            default_cartopy_target_extent_fn=lambda ax, extent: extent,
            extract_regular_grid_from_regridded_vectors_fn=lambda *args: args,
        )

        x_out, y_out, u_out, v_out, color_out, linewidth_out, transform_out = result
        np.testing.assert_allclose(x_out, np.array([10.0, 20.0, 30.0]))
        np.testing.assert_allclose(y_out, np.array([-5.0, 5.0]))
        np.testing.assert_allclose(color_out, np.ones((2, 3)) * 9.0)
        np.testing.assert_allclose(linewidth_out, np.ones((2, 3)) * 7.0)
        assert transform_out == ("mpl", ax)

    def test_prepare_curly_vector_dataset_inputs_impl_validates_cartopy_regrid_requirements(
        self,
    ):
        base_kwargs = dict(
            ax=SimpleNamespace(),
            x=np.array([0.0, 1.0]),
            y=np.array([0.0, 1.0]),
            u=np.ones((2, 2)),
            v=np.ones((2, 2)),
            curvilinear_regrid_shape=None,
            target_extent=None,
            color=None,
            linewidth=None,
            density=1.0,
            maybe_as_scalar_field_fn=lambda value, shape: None,
            is_curvilinear_grid_fn=lambda x, y: False,
            normalize_regrid_shape_fn=lambda value: value,
            default_curvilinear_regrid_shape_fn=lambda x, density: (3, 2),
            regrid_curvilinear_vectors_fn=lambda *args, **kwargs: None,
            regrid_cartopy_vectors_fn=lambda **kwargs: None,
            default_cartopy_target_extent_fn=lambda ax, extent: extent,
            extract_regular_grid_from_regridded_vectors_fn=lambda *args: args,
        )

        with pytest.raises(ValueError, match="requires a Cartopy CRS"):
            _prepare_curly_vector_dataset_inputs_impl(
                transform=None,
                regrid_shape=(3, 2),
                is_cartopy_crs_like_fn=lambda transform: False,
                **base_kwargs,
            )

        with pytest.raises(ValueError, match="requires a Cartopy GeoAxes"):
            _prepare_curly_vector_dataset_inputs_impl(
                transform=SimpleNamespace(_as_mpl_transform=lambda ax: "mpl"),
                regrid_shape=(3, 2),
                is_cartopy_crs_like_fn=lambda transform: True,
                **base_kwargs,
            )

    def test_rcm2rgrid_fields_handles_import_errors_and_shape_normalization(self):
        real_import = builtins.__import__

        def _raising_import(name, *args, **kwargs):
            if name == "skyborn.interp":
                raise ImportError("missing interp")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_raising_import):
            with pytest.raises(ImportError, match="requires skyborn.interp.rcm2rgrid"):
                _rcm2rgrid_fields(
                    np.zeros((2, 2)),
                    np.zeros((2, 2)),
                    (np.zeros((2, 2)),),
                    np.array([0.0, 1.0]),
                    np.array([0.0, 1.0]),
                )

        fake_interp = types.ModuleType("skyborn.interp")
        fake_interp.rcm2rgrid = (
            lambda lat2d, lon2d, field_stack, lat1d, lon1d, msg, meta: np.asarray(
                field_stack, dtype=float
            )
        )
        with patch.dict(sys.modules, {"skyborn.interp": fake_interp}):
            fields = _rcm2rgrid_fields(
                np.zeros((2, 2)),
                np.zeros((2, 2)),
                (
                    np.array([[1.0, 2.0], [3.0, 4.0]]),
                    np.array([[5.0, 6.0], [7.0, 8.0]]),
                ),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
            )
            assert len(fields) == 2
            np.testing.assert_allclose(fields[1], np.array([[5.0, 6.0], [7.0, 8.0]]))

    def test_scalar_field_and_curvilinear_regrid_helpers(self):
        assert _maybe_as_scalar_field(None, (2, 2)) is None
        assert _maybe_as_scalar_field("k", (2, 2)) is None
        assert _maybe_as_scalar_field(3.0, (2, 2)) is None
        assert _maybe_as_scalar_field(np.arange(3.0), (2, 2)) is None
        np.testing.assert_allclose(
            _maybe_as_scalar_field(np.arange(4.0).reshape(2, 2), (2, 2)),
            np.arange(4.0).reshape(2, 2),
        )

        with patch(
            "skyborn.plot._adapters.grid_prepare._rcm2rgrid_fields",
            return_value=[
                np.ones((2, 3)),
                np.ones((2, 3)) * 2.0,
                np.ones((2, 3)) * 3.0,
            ],
        ):
            lon1d, lat1d, u_reg, v_reg, scalar_reg = _regrid_curvilinear_vectors(
                np.array([[0.0, 1.0], [0.0, 1.0]]),
                np.array([[10.0, 10.0], [20.0, 20.0]]),
                np.ones((2, 2)),
                np.ones((2, 2)) * 2.0,
                np.ones((2, 2)) * 3.0,
                target_shape=(3, 2),
            )

        np.testing.assert_allclose(lon1d, np.array([0.0, 0.5, 1.0]))
        np.testing.assert_allclose(lat1d, np.array([10.0, 20.0]))
        np.testing.assert_allclose(u_reg, np.ones((2, 3)))
        np.testing.assert_allclose(v_reg, np.ones((2, 3)) * 2.0)
        np.testing.assert_allclose(scalar_reg, np.ones((2, 3)) * 3.0)


class TestNativeHelpers:
    def test_call_native_sample_grid_field_handles_masked_and_nonfinite_results(self):
        grid = SimpleNamespace(x_origin=1.0, y_origin=2.0, dx=3.0, dy=4.0)

        def _masked_sampler(**kwargs):
            assert np.isnan(kwargs["field"]).all()
            return np.nan

        assert (
            _call_native_sample_grid_field(
                _masked_sampler,
                grid,
                np.ma.array(np.ones((2, 2)), mask=True),
                1.0,
                2.0,
            )
            is None
        )
        with pytest.raises(RuntimeError, match="boom"):
            _call_native_sample_grid_field(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                grid,
                np.ones((2, 2)),
                1.0,
                2.0,
            )
        assert (
            _call_native_sample_grid_field(
                lambda **kwargs: np.nan,
                grid,
                np.ones((2, 2)),
                1.0,
                2.0,
            )
            is None
        )
        assert _call_native_sample_grid_field(
            lambda **kwargs: 5.5,
            grid,
            np.ones((2, 2)),
            1.0,
            2.0,
        ) == pytest.approx(5.5)

    def test_call_native_sample_grid_field_array_handles_shape_errors_and_nonfinite_values(
        self,
    ):
        grid = SimpleNamespace(x_origin=1.0, y_origin=2.0, dx=3.0, dy=4.0)
        points = np.array([[0.0, 0.0], [1.0, 1.0]])

        with pytest.raises(RuntimeError, match="boom"):
            _call_native_sample_grid_field_array(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                grid,
                np.ones((2, 2)),
                points,
                (2,),
            )

        assert (
            _call_native_sample_grid_field_array(
                lambda **kwargs: None,
                grid,
                np.ones((2, 2)),
                points,
                (2,),
            )
            is None
        )
        assert (
            _call_native_sample_grid_field_array(
                lambda **kwargs: np.array([[1.0, 2.0]]),
                grid,
                np.ones((2, 2)),
                points,
                (2,),
            )
            is None
        )
        sampled = _call_native_sample_grid_field_array(
            lambda **kwargs: np.array([1.0, np.nan]),
            grid,
            np.ones((2, 2)),
            points,
            (2,),
        )
        assert sampled[0] == pytest.approx(1.0)
        assert np.isnan(sampled[1])

    def test_call_native_thin_ncl_mapped_candidates_handles_none_and_success(self):
        points = np.array([[0.1, 0.2], [0.3, 0.4]])

        assert (
            _call_native_thin_ncl_mapped_candidates(
                lambda **kwargs: None,
                points,
                0.1,
            )
            is None
        )
        with pytest.raises(RuntimeError, match="boom"):
            _call_native_thin_ncl_mapped_candidates(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                points,
                0.1,
            )
        assert _call_native_thin_ncl_mapped_candidates(
            lambda **kwargs: np.array([1, 0], dtype=int),
            points,
            0.1,
        ) == [1, 0]

    def test_call_native_thin_ncl_display_candidates_handles_none(self):
        viewport = Bbox.from_bounds(1.0, 2.0, 3.0, 4.0)
        assert (
            _call_native_thin_ncl_display_candidates(
                lambda **kwargs: None,
                np.array([[1.0, 2.0]]),
                viewport,
                0.1,
            )
            is None
        )

    def test_call_native_trace_ncl_direction_handles_missing_context_and_shape_validation(
        self,
    ):
        context = SimpleNamespace(
            u=np.ones((2, 2)),
            v=np.ones((2, 2)),
            display_grid=np.ones((2, 2, 2)),
            cell_valid=np.ones((1, 1), dtype=bool),
            x_origin=0.0,
            y_origin=0.0,
            dx=1.0,
            dy=1.0,
            viewport_x0=0.0,
            viewport_y0=0.0,
            viewport_x1=1.0,
            viewport_y1=1.0,
        )

        assert (
            _call_native_trace_ncl_direction(
                lambda **kwargs: np.array([[0.0, 0.0], [1.0, 1.0]]),
                None,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
            is None
        )
        with pytest.raises(RuntimeError, match="boom"):
            _call_native_trace_ncl_direction(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
        assert (
            _call_native_trace_ncl_direction(
                lambda **kwargs: np.array([1.0, 2.0]),
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
            is None
        )
        curve = _call_native_trace_ncl_direction(
            lambda **kwargs: np.array([[0.0, 0.0], [1.0, 1.0]]),
            context,
            np.array([0.0, 0.0]),
            10.0,
            1.0,
            1.0,
            1.0,
        )
        np.testing.assert_allclose(curve, np.array([[0.0, 0.0], [1.0, 1.0]]))
        assert (
            _call_native_trace_ncl_direction(
                lambda **kwargs: None,
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
            is None
        )

    def test_call_native_build_arrow_helpers_validate_shapes(self):
        assert (
            vector_module._native_helpers._call_native_build_open_arrow_segments(
                lambda **kwargs: (np.array([[1.0, 2.0]]), np.array([0], dtype=int)),
                np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
                np.array([0, 2], dtype=np.intp),
                np.array([1.0]),
                np.array([1.0]),
            )
            is None
        )
        assert (
            vector_module._native_helpers._call_native_build_open_arrow_segments(
                lambda **kwargs: (
                    np.array([[[0.0, 0.0], [1.0, 1.0]]], dtype=float),
                    np.array([0, 1], dtype=int),
                ),
                np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
                np.array([0, 2], dtype=np.intp),
                np.array([1.0]),
                np.array([1.0]),
            )
            is None
        )
        assert (
            vector_module._native_helpers._call_native_build_filled_arrow_polygons(
                lambda **kwargs: (np.array([[1.0, 2.0]]), np.array([0], dtype=int)),
                np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
                np.array([0, 2], dtype=np.intp),
                np.array([1.0]),
                np.array([1.0]),
            )
            is None
        )
        assert (
            vector_module._native_helpers._call_native_build_filled_arrow_polygons(
                lambda **kwargs: (
                    np.array([[[0.0, 0.0], [1.0, 1.0], [0.5, 0.5]]], dtype=float),
                    np.array([0, 1], dtype=int),
                ),
                np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
                np.array([0, 2], dtype=np.intp),
                np.array([1.0]),
                np.array([1.0]),
            )
            is None
        )

    def test_call_native_trace_with_display_validates_context_and_shapes(self):
        context = SimpleNamespace(
            u=np.ones((2, 2)),
            v=np.ones((2, 2)),
            display_grid=np.ones((2, 2, 2)),
            cell_valid=np.ones((1, 1), dtype=bool),
            x_origin=0.0,
            y_origin=0.0,
            dx=1.0,
            dy=1.0,
            viewport_x0=0.0,
            viewport_y0=0.0,
            viewport_x1=1.0,
            viewport_y1=1.0,
        )

        assert (
            _call_native_trace_ncl_direction_with_display(
                lambda **kwargs: None,
                None,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
            is None
        )
        assert (
            _call_native_trace_ncl_direction_with_display(
                lambda **kwargs: (np.array([[0.0, 0.0]]), np.array([[0.0, 0.0]])),
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
            )
            is None
        )

        traced = _call_native_trace_ncl_direction_with_display(
            lambda **kwargs: (
                np.array([[0.0, 0.0], [1.0, 1.0]]),
                np.array([[2.0, 2.0], [3.0, 3.0]]),
            ),
            context,
            np.array([0.0, 0.0]),
            10.0,
            1.0,
            1.0,
            1.0,
        )
        assert traced is not None
        curve, display_curve = traced
        np.testing.assert_allclose(curve, np.array([[0.0, 0.0], [1.0, 1.0]]))
        np.testing.assert_allclose(display_curve, np.array([[2.0, 2.0], [3.0, 3.0]]))

    def test_call_native_validate_display_curve_forwards_viewport_size(self):
        captured = {}

        def validator(**kwargs):
            captured.update(kwargs)
            return True

        viewport = Bbox.from_bounds(0.0, 0.0, 12.0, 8.0)
        assert _call_native_validate_display_curve(
            validator,
            np.array([[0.0, 0.0], [1.0, 1.0]]),
            viewport,
        )
        assert captured["viewport_width"] == pytest.approx(12.0)
        assert captured["viewport_height"] == pytest.approx(8.0)


class TestGeometryHelpers:
    def test_display_points_to_data_validates_inputs_and_inverse_transform(self):
        transform = Affine2D().scale(2.0)

        data = _display_points_to_data(transform, np.array([2.0, 4.0]))
        np.testing.assert_allclose(data, np.array([[1.0, 2.0]]))

        assert _display_points_to_data(transform, np.array([1.0, 2.0, 3.0])) is None
        assert _display_points_to_data(transform, np.array([[np.nan, 1.0]])) is None

        bad_inverse = SimpleNamespace(
            transform=lambda values: np.full_like(
                np.asarray(values, dtype=float), np.nan
            )
        )
        assert (
            _display_points_to_data(
                transform,
                np.array([[1.0, 2.0]]),
                inverse_transform=bad_inverse,
            )
            is None
        )
        assert _display_points_to_data(_BadTransform(), np.array([[1.0, 2.0]])) is None

    def test_display_step_and_candidate_data_fallback_paths(self):
        assert (
            _display_step_to_data(
                np.array([[1.0, 0.0], [0.0, 0.0]]),
                np.array([1.0, 1.0]),
            )
            is None
        )

        candidate = _candidate_data_from_display_step(
            current_data=np.array([0.0, 0.0]),
            current_display=np.array([0.0, 0.0]),
            candidate_display=np.array([2.0, 3.0]),
            jacobian=np.zeros((2, 2)),
            transform=Affine2D(),
        )
        np.testing.assert_allclose(candidate, np.array([2.0, 3.0]))

    def test_arc_helpers_cover_edge_cases(self):
        with pytest.raises(ValueError, match="at least one point"):
            _point_at_arc_distance_from_end(np.empty((0, 2)), 1.0)

        np.testing.assert_allclose(
            _point_at_arc_distance_from_end(np.array([[5.0, 6.0]]), 10.0),
            np.array([5.0, 6.0]),
        )

    def test_evaluate_curve_and_tip_geometry_handle_failures(self):
        curve = np.array([[0.0, 0.0], [100.0, 0.0], [2000.0, 0.0]])
        viewport = Bbox.from_bounds(0.0, 0.0, 100.0, 50.0)

        display_curve, transform_failed = _evaluate_ncl_display_curve(
            curve, _BadTransform()
        )
        assert display_curve is None
        assert transform_failed is True

        display_curve, transform_failed = _evaluate_ncl_display_curve(
            curve, Affine2D(), viewport=viewport
        )
        assert display_curve is None
        assert transform_failed is False
        display_curve, transform_failed = _evaluate_ncl_display_curve(
            curve, _BadTransform()
        )
        assert display_curve is None
        assert transform_failed is True

        assert (
            _tip_display_geometry_from_display_curve(
                np.array([[1.0, 1.0], [1.0, 1.0]]),
                1.0,
            )
            is None
        )
        assert _tip_display_geometry(np.array([[0.0, 0.0]]), Affine2D(), 1.0) is None
        assert (
            _tip_display_geometry(
                np.array([[0.0, 0.0], [1.0, 1.0]]),
                _BadTransform(),
                1.0,
            )
            is None
        )

    def test_finite_difference_local_jacobian_and_trim_helpers(self):
        assert _finite_difference_step(1.0, 0.0, 2.0, 0.0) == 0.0
        assert _finite_difference_step(2.0, 0.0, 2.0, 0.5) == -0.5
        assert _finite_difference_step(0.0, 0.0, 0.0, 0.5) == 0.0

        grid = SimpleNamespace(
            x_origin=0.0,
            y_origin=0.0,
            width=2.0,
            height=2.0,
            dx=1.0,
            dy=1.0,
        )
        assert (
            _local_display_jacobian(_BadTransform(), np.array([0.5, 0.5]), grid) is None
        )

        transform = Affine2D().scale(2.0, 3.0)
        jacobian = _local_display_jacobian(transform, np.array([0.5, 0.5]), grid)
        np.testing.assert_allclose(jacobian, np.array([[2.0, 0.0], [0.0, 3.0]]))

        zero_curve = np.array([[0.0, 0.0], [0.0, 0.0]])
        np.testing.assert_allclose(
            _trim_display_curve_from_end(zero_curve, 1.0), zero_curve
        )

        long_curve = np.array([[0.0, 0.0], [1.0, 0.0]])
        trimmed = _trim_display_curve_from_end(long_curve, 1.0)
        assert trimmed.shape == (2, 2)
        np.testing.assert_allclose(trimmed[-1], np.array([0.45, 0.0]))

    def test_geometry_wrappers_cover_short_curves_and_display_quality_filters(self):
        assert _display_to_data(_BadTransform(), np.array([1.0, 2.0])) is None

        single_segment_curve = np.array([[0.0, 0.0], [1.0, 0.0]])
        display_curve, transform_failed = _evaluate_ncl_display_curve(
            single_segment_curve,
            Affine2D(),
        )
        assert transform_failed is False
        np.testing.assert_allclose(display_curve, single_segment_curve)

        invalid_curve = np.array([[0.0, 0.0]])
        assert _evaluate_ncl_display_curve(invalid_curve, Affine2D()) == (None, False)
        flat_curve = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
        assert _evaluate_ncl_display_curve(flat_curve, Affine2D()) == (None, False)

        tortuous_curve = np.array(
            [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [2.0, 1.0], [2.0, 2.0]]
        )
        assert _evaluate_ncl_display_curve(tortuous_curve, Affine2D()) == (None, False)

        sign_change_curve = np.array(
            [[0.0, 0.0], [1.0, 0.0], [2.0, 1.0], [3.0, 0.0], [4.0, -1.0]]
        )
        assert _evaluate_ncl_display_curve(sign_change_curve, Affine2D()) == (
            None,
            False,
        )

        tiny_grid = SimpleNamespace(
            x_origin=0.0,
            y_origin=0.0,
            width=0.0,
            height=2.0,
            dx=0.0,
            dy=1.0,
        )
        assert (
            _local_display_jacobian(Affine2D(), np.array([0.0, 0.0]), tiny_grid) is None
        )

        class FailingStepTransform:
            def transform(self, values):
                raise RuntimeError("bad step transform")

        assert (
            _local_display_jacobian(
                FailingStepTransform(),
                np.array([0.5, 0.5]),
                SimpleNamespace(
                    x_origin=0.0,
                    y_origin=0.0,
                    width=2.0,
                    height=2.0,
                    dx=1.0,
                    dy=1.0,
                ),
            )
            is None
        )


class TestStyleHelpers:
    def test_collect_named_kwargs_and_arrowstyle_validation(self):
        scope = {"alpha": 0.5, "color": "k", "linewidth": 2.0}
        assert _collect_named_kwargs(scope, ("alpha", "color")) == {
            "alpha": 0.5,
            "color": "k",
        }
        assert _normalize_supported_arrowstyle(" -> ") == "->"
        with pytest.raises(ValueError, match="arrowstyle must be one of"):
            _normalize_supported_arrowstyle("fancy")

    def test_normalize_artist_alpha_validates_range(self):
        assert _normalize_artist_alpha(None) is None
        assert _normalize_artist_alpha(0.25) == pytest.approx(0.25)
        with pytest.raises(ValueError, match="between 0 and 1"):
            _normalize_artist_alpha(2.0)
        with pytest.raises(ValueError, match="between 0 and 1"):
            _normalize_artist_alpha(np.nan)

    def test_curly_head_axes_dimensions_handles_invalid_inputs_and_fallbacks(self):
        axes = SimpleNamespace(bbox=SimpleNamespace(width=200.0, height=100.0))

        assert _curly_head_axes_dimensions(
            axes, "bad", max_magnitude=10.0, arrowsize=1.0
        ) == (0.0, 0.0)
        assert _curly_head_axes_dimensions(
            None, 5.0, max_magnitude=10.0, arrowsize=1.0
        ) == (0.0, 0.0)
        assert _curly_head_axes_dimensions(
            axes, -1.0, max_magnitude=10.0, arrowsize=1.0
        ) == (0.0, 0.0)

        class BadAxes:
            @property
            def bbox(self):
                raise RuntimeError("bad bbox")

        assert _curly_head_axes_dimensions(
            BadAxes(), 5.0, max_magnitude=10.0, arrowsize=1.0
        ) == (0.0, 0.0)

        head_length, head_width = _curly_head_axes_dimensions(
            axes,
            5.0,
            max_magnitude="bad",
            arrowsize="bad",
        )
        assert head_length > 0.0
        assert head_width > 0.0

    def test_resolve_curly_anchor_and_style_aliases_validate_conflicts(self):
        assert _resolve_curly_anchor_alias(None, "tip") == "head"
        assert _resolve_curly_anchor_alias("center", None) == "center"
        with pytest.raises(ValueError, match="different glyph anchors"):
            _resolve_curly_anchor_alias("tail", "mid")

        color, linewidth, facecolor, edgecolor, vmin, vmax = (
            _resolve_curly_style_aliases(
                color=None,
                c="k",
                linewidth=None,
                linewidths=2.0,
                facecolor=None,
                facecolors="gold",
                edgecolor=None,
                edgecolors="navy",
                norm=None,
                vmin=0.0,
                vmax=5.0,
            )
        )
        assert color == "k"
        assert linewidth == 2.0
        assert facecolor == "gold"
        assert edgecolor == "navy"
        assert vmin == pytest.approx(0.0)
        assert vmax == pytest.approx(5.0)

        with pytest.raises(ValueError, match="Use only one of 'c' or 'color'"):
            _resolve_curly_style_aliases(
                color="k",
                c="r",
                linewidth=None,
                linewidths=None,
                facecolor=None,
                facecolors=None,
                edgecolor=None,
                edgecolors=None,
                norm=None,
                vmin=None,
                vmax=None,
            )
        with pytest.raises(ValueError, match="Use only one of 'norm' or 'vmin'/'vmax'"):
            _resolve_curly_style_aliases(
                color=None,
                c=None,
                linewidth=None,
                linewidths=None,
                facecolor=None,
                facecolors=None,
                edgecolor=None,
                edgecolors=None,
                norm=object(),
                vmin=0.0,
                vmax=None,
            )
        with pytest.raises(ValueError, match="vmin must be finite"):
            _resolve_curly_style_aliases(
                color=None,
                c=None,
                linewidth=None,
                linewidths=None,
                facecolor=None,
                facecolors=None,
                edgecolor=None,
                edgecolors=None,
                norm=None,
                vmin=np.nan,
                vmax=None,
            )
        with pytest.raises(ValueError, match="vmax must be finite"):
            _resolve_curly_style_aliases(
                color=None,
                c=None,
                linewidth=None,
                linewidths=None,
                facecolor=None,
                facecolors=None,
                edgecolor=None,
                edgecolors=None,
                norm=None,
                vmin=None,
                vmax=np.inf,
            )

    def test_resolve_scatter_aliases_maps_singular_kwargs(self):
        c, kwargs = _resolve_scatter_aliases(
            None,
            {
                "color": "k",
                "linewidth": 2.0,
                "facecolor": "gold",
                "edgecolor": "navy",
            },
        )
        assert c == "k"
        assert kwargs["linewidths"] == 2.0
        assert kwargs["facecolors"] == "gold"
        assert kwargs["edgecolors"] == "navy"

        with pytest.raises(ValueError, match="Use only one of 'c' or 'color'"):
            _resolve_scatter_aliases("k", {"color": "r"})
        with pytest.raises(
            ValueError, match="Use only one of 'linewidth' or 'linewidths'"
        ):
            _resolve_scatter_aliases(None, {"linewidth": 1.0, "linewidths": 2.0})
        with pytest.raises(
            ValueError, match="Use only one of 'facecolor' or 'facecolors'"
        ):
            _resolve_scatter_aliases(None, {"facecolor": "k", "facecolors": "r"})
        with pytest.raises(
            ValueError, match="Use only one of 'edgecolor' or 'edgecolors'"
        ):
            _resolve_scatter_aliases(None, {"edgecolor": "k", "edgecolors": "r"})


class TestCoordinateHelpers:
    def test_filled_subset_and_transpose_helpers(self):
        numeric_masked = np.ma.array([1.0, 2.0], mask=[False, True])
        filled = _filled_float_array(numeric_masked)
        assert filled[0] == pytest.approx(1.0)
        assert np.isnan(filled[1])

        object_masked = np.ma.array(
            np.array(["a", "b"], dtype=object),
            mask=[False, True],
        )
        subset = _subset_ready_array(object_masked)
        assert subset.tolist() == ["a", "?"]
        np.testing.assert_allclose(
            _subset_ready_array(np.array([1, 2, 3])),
            np.array([1, 2, 3]),
        )

        da = xr.DataArray(np.arange(6).reshape(2, 3), dims=("x", "y"))
        transposed = _transpose_dataarray_if_possible(da, ("y", "x"))
        assert transposed.dims == ("y", "x")

    def test_normalize_coordinates_handles_grid_inference_and_validation(self):
        x = xr.DataArray(np.arange(6.0).reshape(2, 3), dims=("lat", "lon"))
        y = xr.DataArray(np.arange(6.0).reshape(3, 2), dims=("lon", "lat"))
        x_values, y_values, shape, grid_like, target_dims = _normalize_coordinates(
            x,
            y,
            reference_fields=(),
        )
        assert grid_like is True
        assert shape == (2, 3)
        assert target_dims == ("lat", "lon")
        np.testing.assert_allclose(y_values, y.transpose("lat", "lon").values)

        x1d = xr.DataArray(np.array([1.0, 2.0, 3.0]), dims=("lon",))
        y1d = xr.DataArray(np.array([-1.0, 0.0, 1.0]), dims=("lat",))
        field = xr.DataArray(np.arange(9.0).reshape(3, 3), dims=("lat", "lon"))
        x_grid, y_grid, shape, grid_like, target_dims = _normalize_coordinates(
            x1d,
            y1d,
            reference_fields=(field,),
        )
        assert grid_like is True
        assert shape == (3, 3)
        assert target_dims == ("lat", "lon")
        assert x_grid.shape == (3, 3)
        assert y_grid.shape == (3, 3)

        x_points, y_points, shape, grid_like, _ = _normalize_coordinates(
            np.array([1.0, 2.0, 3.0]),
            np.array([4.0, 5.0, 6.0]),
            reference_fields=(),
        )
        assert grid_like is False
        assert shape == (3,)
        np.testing.assert_allclose(x_points, np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(y_points, np.array([4.0, 5.0, 6.0]))

        with pytest.raises(ValueError, match="must both be 1D or both be 2D"):
            _normalize_coordinates(
                np.ones((2, 2)),
                np.array([1.0, 2.0]),
                reference_fields=(),
            )
        with pytest.raises(ValueError, match="same shape"):
            _normalize_coordinates(
                np.ones((2, 2)),
                np.ones((3, 2)),
                reference_fields=(),
            )
        with pytest.raises(ValueError, match="1D or 2D arrays"):
            _normalize_coordinates(
                np.array(1.0),
                np.array([1.0, 2.0, 3.0]),
                reference_fields=(),
            )

    def test_selection_axis_meshgrid_and_corner_helpers(self):
        with pytest.raises(ValueError, match="Use only one of 'where' or 'mask'"):
            _normalize_selection_mask(
                where=np.ones((2, 2), dtype=bool),
                mask=np.ones((2, 2), dtype=bool),
                candidate_shape=(2, 2),
                target_dims=None,
            )

        default_mask = _normalize_selection_mask(
            where=None,
            mask=None,
            candidate_shape=(2, 2),
            target_dims=None,
        )
        assert default_mask.dtype == bool
        assert default_mask.all()

        selector = xr.DataArray(
            np.array([[1.0, 0.0], [np.nan, 2.0]]),
            dims=("x", "y"),
        )
        normalized = _normalize_selection_mask(
            where=selector,
            mask=None,
            candidate_shape=(2, 2),
            target_dims=("y", "x"),
        )
        np.testing.assert_array_equal(
            normalized, np.array([[True, False], [False, True]])
        )

        np.testing.assert_allclose(
            _axis_coordinate_1d(np.array([1.0, 2.0]), "x"),
            np.array([1.0, 2.0]),
        )
        np.testing.assert_allclose(
            _axis_coordinate_1d(np.array([[1.0, 2.0], [1.0, 2.0]]), "x"),
            np.array([1.0, 2.0]),
        )
        assert _axis_coordinate_1d(np.ones((2, 2, 2)), "x") is None
        with pytest.raises(ValueError, match="Unsupported axis_name"):
            _axis_coordinate_1d(np.ones((2, 2)), "z")

        assert _axis_is_uniform(np.array([1.0, 2.0])) is True
        assert _axis_is_uniform(np.array([0.0, 1.0, np.nan])) is False

        with pytest.raises(ValueError, match="non-empty 1D axis"):
            _cell_edges_from_axis(np.array([]))
        np.testing.assert_allclose(
            _cell_edges_from_axis(np.array([4.0])),
            np.array([3.5, 4.5]),
        )

        x2d = np.array([[0.0, 1.0], [0.0, 2.0]])
        y2d = np.array([[10.0, 10.0], [20.0, 20.0]])
        assert _rectilinear_cell_edges(x2d, y2d) is None
        assert (
            _rectilinear_cell_edges(
                np.array([[0.0, 1.0], [0.0, 1.0]]),
                np.array([[10.0, 10.0], [20.0, 21.0]]),
            )
            is None
        )

        with pytest.raises(ValueError, match="non-empty 2D center grid"):
            _center_grid_to_corner_grid(np.array([1.0, 2.0]))

        np.testing.assert_allclose(
            _center_grid_to_corner_grid(np.array([[1.0, 3.0]])),
            np.array([[0.0, 2.0, 4.0], [0.0, 2.0, 4.0]]),
        )
        np.testing.assert_allclose(
            _center_grid_to_corner_grid(np.array([[1.0], [3.0]])),
            np.array([[0.0, 0.0], [2.0, 2.0], [4.0, 4.0]]),
        )
        np.testing.assert_allclose(
            _center_grid_to_corner_grid(np.array([[2.0]])),
            np.array([[2.0, 2.0], [2.0, 2.0]]),
        )

    def test_cell_geometry_field_matching_and_orientation_helpers(self):
        rect = _scatter_cell_geometry(np.array([0.0, 1.0]), np.array([10.0, 20.0]))
        assert str(rect["kind"]) == "rectilinear"
        curvi = _scatter_cell_geometry(
            np.array([[0.0, 1.0], [0.2, 1.2]]),
            np.array([[10.0, 10.1], [20.0, 20.1]]),
        )
        assert str(curvi["kind"]) == "curvilinear"
        assert curvi["x_corners"].shape == (3, 3)
        assert _scatter_cell_geometry(np.ones((2, 2, 2)), np.ones((2, 2))) is None

        assert _coerce_matching_plot_field(None, (2, 2)) == (None, False)
        assert _coerce_matching_plot_field("k", (2, 2)) == (None, False)
        field, is_field = _coerce_matching_plot_field(
            np.arange(4.0).reshape(2, 2), (2, 2)
        )
        assert is_field is True
        np.testing.assert_allclose(field, np.arange(4.0).reshape(2, 2))
        assert _coerce_matching_plot_field(np.arange(3.0), (2, 2)) == (None, False)
        assert _coerce_matching_plot_field(np.arange(6.0).reshape(2, 3), (2, 2)) == (
            None,
            True,
        )

        x_axes, y_axes = _extract_meshgrid_axes(
            np.array([1.0, 2.0]), np.array([3.0, 4.0])
        )
        np.testing.assert_allclose(x_axes, np.array([1.0, 2.0]))
        np.testing.assert_allclose(y_axes, np.array([3.0, 4.0]))
        with pytest.raises(ValueError, match="1D axes or meshgrid-like 2D x/y"):
            _extract_meshgrid_axes(np.ones((2, 2, 2)), np.ones((2, 2)))
        with pytest.raises(ValueError, match="same shape"):
            _extract_meshgrid_axes(np.ones((2, 2)), np.ones((3, 2)))
        with pytest.raises(ValueError, match="x coordinates must be meshgrid-like"):
            _extract_meshgrid_axes(
                np.array([[0.0, 1.0], [1.0, 2.0]]),
                np.array([[0.0, 0.0], [1.0, 1.0]]),
            )
        with pytest.raises(ValueError, match="y coordinates must be meshgrid-like"):
            _extract_meshgrid_axes(
                np.array([[0.0, 1.0], [0.0, 1.0]]),
                np.array([[0.0, 1.0], [1.0, 2.0]]),
            )

        x_bad = np.ones((2, 2, 2))
        y_bad = np.ones((2, 2))
        original = np.arange(4.0).reshape(2, 2)
        result = _normalize_regular_grid_orientation(x_bad, y_bad, original, original)
        assert result == (x_bad, y_bad, original, original, None, None)

        x = np.array([2.0, 1.0, 0.0])
        y = np.array([10.0, 0.0])
        u = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v = np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]])
        color = np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]])
        linewidth = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

        x_norm, y_norm, u_norm, v_norm, color_norm, linewidth_norm = (
            _normalize_regular_grid_orientation(
                x, y, u, v, color=color, linewidth=linewidth
            )
        )
        np.testing.assert_allclose(x_norm, np.array([0.0, 1.0, 2.0]))
        np.testing.assert_allclose(y_norm, np.array([0.0, 10.0]))
        np.testing.assert_allclose(u_norm, np.array([[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]))
        np.testing.assert_allclose(
            v_norm, np.array([[60.0, 50.0, 40.0], [30.0, 20.0, 10.0]])
        )
        np.testing.assert_allclose(
            color_norm, np.array([[600.0, 500.0, 400.0], [300.0, 200.0, 100.0]])
        )
        np.testing.assert_allclose(
            linewidth_norm, np.array([[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]])
        )


class TestVectorEngineHelpers:
    def test_resolve_artist_coordinate_context_and_reference_lengths(self):
        fig, ax = plt.subplots()
        try:
            ax.projection = object()
            transform = Affine2D().translate(1.0, 2.0) + ax.transData
            artist_transform, inverse_transform, baked = (
                _resolve_artist_coordinate_context(ax, transform)
            )
            assert artist_transform is ax.transData
            assert inverse_transform is not None
            assert baked is True
        finally:
            plt.close(fig)

        bad_axes = SimpleNamespace(
            projection=object(),
            transData=SimpleNamespace(
                inverted=lambda: (_ for _ in ()).throw(RuntimeError("bad"))
            ),
        )
        artist_transform, inverse_transform, baked = _resolve_artist_coordinate_context(
            bad_axes,
            transform,
        )
        assert artist_transform == transform
        assert inverse_transform is None
        assert baked is False

        assert _resolve_ncl_reference_length_px(
            1.0, 1.0, 2.0, 0.0, 0.0, 12.0
        ) == pytest.approx(12.0)
        assert _resolve_ncl_reference_length_px(
            1.0, 5.0, 2.0, 9.0, 0.0, 12.0
        ) == pytest.approx(9.0)
        assert _resolve_ncl_reference_length_px(1.0, 5.0, 3.0, 0.0, 0.5, 12.0) > 0.0
        assert _resolve_ncl_reference_length_px(
            1.0, 5.0, 3.0, 0.0, 0.0, 12.0
        ) == pytest.approx(7.2)

    def test_resolve_ncl_length_scale_covers_branch_matrix(self):
        same_mag = _resolve_ncl_length_scale(2.0, 2.0, 4.0, 8.0, 0.0, 12.0)
        assert same_mag["ref_length_px"] == pytest.approx(8.0)
        assert same_mag["max_length_px"] == pytest.approx(4.0)

        same_mag_ref_only = _resolve_ncl_length_scale(2.0, 2.0, 4.0, 0.0, 0.0, 12.0)
        assert same_mag_ref_only["ref_length_px"] == pytest.approx(24.0)

        same_mag_length_only = _resolve_ncl_length_scale(2.0, 2.0, 0.0, 6.0, 0.0, 12.0)
        assert same_mag_length_only["ref_length_px"] == pytest.approx(6.0)

        no_ref = _resolve_ncl_length_scale(1.0, 5.0, 0.0, 6.0, 0.3, 12.0)
        assert no_ref["adjust_min"] is True
        assert no_ref["min_length_px"] == pytest.approx(1.8)

        below_min = _resolve_ncl_length_scale(2.0, 6.0, 1.0, 0.0, 0.0, 12.0)
        assert below_min["ref_length_px"] == pytest.approx(2.0)
        assert below_min["min_length_px"] == pytest.approx(4.0)

        below_min_fixed = _resolve_ncl_length_scale(2.0, 6.0, 1.0, 7.0, 0.0, 12.0)
        assert below_min_fixed["ref_length_px"] == pytest.approx(7.0)
        assert below_min_fixed["max_length_px"] == pytest.approx(42.0)

        below_min_frac = _resolve_ncl_length_scale(2.0, 6.0, 1.0, 0.0, 0.25, 12.0)
        assert below_min_frac["ref_length_px"] == pytest.approx(3.0)

        adjusted = _resolve_ncl_length_scale(1.0, 5.0, 3.0, 7.0, 0.25, 12.0)
        assert adjusted["adjust_min"] is True
        assert adjusted["ref_length_px"] == pytest.approx(7.0)
        assert adjusted["min_length_px"] == pytest.approx(1.75)

        default_adjusted = _resolve_ncl_length_scale(1.0, 5.0, 3.0, 0.0, 0.25, 12.0)
        assert default_adjusted["adjust_min"] is True
        assert default_adjusted["min_length_px"] > 0.0

        requested_ref = _resolve_ncl_length_scale(1.0, 5.0, 3.0, 9.0, 0.0, 12.0)
        assert requested_ref["ref_length_px"] == pytest.approx(9.0)
        assert requested_ref["min_length_px"] == pytest.approx(3.0)

        implicit_ref = _resolve_ncl_length_scale(1.0, 5.0, 3.0, 0.0, 0.0, 12.0)
        assert implicit_ref["ref_length_px"] == pytest.approx(7.2)
        assert implicit_ref["min_length_px"] == pytest.approx(2.4)

    def test_ncl_center_selection_covers_empty_and_mapped_thinning_paths(self):
        grid = _make_test_grid()
        magnitude = np.ones(grid.shape)
        axes = SimpleNamespace(bbox=Bbox.from_bounds(0.0, 0.0, 100.0, 80.0))

        empty = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=Affine2D(),
            axes=axes,
            density=1.0,
            start_points=None,
            min_distance=None,
            sample_grid_field_array=lambda grid, field, points: np.full(
                len(points), np.nan
            ),
            thin_ncl_mapped_candidates=lambda mapped_points, spacing: np.array(
                [], dtype=int
            ),
        )
        assert empty == []

        selected = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=Affine2D(),
            axes=axes,
            density=1.0,
            start_points=None,
            min_distance=None,
            sample_grid_field_array=lambda grid, field, points: np.ones(len(points)),
            thin_ncl_display_candidates=lambda display_points, bbox, spacing: None,
            thin_ncl_mapped_candidates=lambda mapped_points, spacing: np.array(
                [0], dtype=int
            ),
        )
        assert len(selected) == 1

        valid = _valid_ncl_center_candidates(
            grid,
            np.array([[0.5, 0.5], [5.0, 5.0]]),
            np.array([1.0, 1.0]),
            np.array([[1.0, 1.0], [2.0, 2.0]]),
            start_points=None,
        )
        np.testing.assert_array_equal(valid, np.array([True, False]))

    def test_prepare_centers_trace_curve_and_build_curve_cover_error_paths(self):
        grid = _make_test_grid()
        magnitude = np.arange(9.0).reshape(3, 3)

        candidates, magnitudes = _prepare_ncl_center_candidates(
            grid=grid,
            magnitude=magnitude,
            density=1.0,
            start_points=None,
            ncl_preset=None,
            sample_grid_field_array=lambda grid, field, points: np.arange(
                len(points), dtype=float
            ),
        )
        assert candidates.shape[1] == 2
        np.testing.assert_allclose(magnitudes, np.arange(len(candidates), dtype=float))

        with pytest.raises(ValueError, match="must be an \\(N, 2\\) array"):
            _prepare_ncl_center_candidates(
                grid=grid,
                magnitude=magnitude,
                density=1.0,
                start_points=np.array([1.0, 2.0, 3.0]),
                ncl_preset=None,
                sample_grid_field_array=lambda grid, field, points: np.arange(
                    len(points), dtype=float
                ),
            )

        trace_fn = Mock(side_effect=[None, np.array([[0.0, 0.0], [1.0, 0.0]])])
        forward = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="center",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_direction_fn=trace_fn,
        )
        np.testing.assert_allclose(forward, np.array([[0.0, 0.0], [1.0, 0.0]]))

        trace_fn = Mock(side_effect=[np.array([[0.0, 0.0], [1.0, 0.0]]), None])
        backward = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="center",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_direction_fn=trace_fn,
        )
        np.testing.assert_allclose(backward, np.array([[1.0, 0.0], [0.0, 0.0]]))

        backward_only = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="head",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_direction_fn=lambda *args, **kwargs: np.array(
                [[0.0, 0.0], [1.0, 0.0]]
            ),
        )
        np.testing.assert_allclose(backward_only, np.array([[1.0, 0.0], [0.0, 0.0]]))

        trace_fn = Mock(side_effect=[None, None])
        assert (
            _trace_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=4.0,
                anchor="center",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_direction_fn=trace_fn,
            )
            is None
        )

        trace_fn = Mock(
            side_effect=[
                np.array([[0.0, 0.0], [-1.0, 0.0]]),
                np.array([[0.0, 0.0], [1.0, 0.0]]),
            ]
        )
        centered = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="center",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_direction_fn=trace_fn,
        )
        np.testing.assert_allclose(
            centered,
            np.array([[-1.0, 0.0], [0.0, 0.0], [1.0, 0.0]]),
        )

        tail_curve = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="tail",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_direction_fn=lambda *args, **kwargs: np.array(
                [[0.0, 0.0], [1.0, 0.0]]
            ),
        )
        np.testing.assert_allclose(tail_curve, np.array([[0.0, 0.0], [1.0, 0.0]]))

        assert (
            _trace_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=4.0,
                anchor="head",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_direction_fn=lambda *args, **kwargs: None,
            )
            is None
        )

        assert (
            _trace_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=0.0,
                anchor="tail",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_direction_fn=lambda *args, **kwargs: np.array(
                    [[0.0, 0.0], [1.0, 0.0]]
                ),
            )
            is None
        )

        calls = {"count": 0}

        def _trace_curve(**kwargs):
            calls["count"] += 1
            return np.array([[0.0, 0.0], [1.0, 0.0]])

        def _evaluate(curve, transform, viewport=None):
            if calls["count"] == 1:
                return None, False
            return np.asarray(curve, dtype=float), False

        curve, display_curve = _build_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=10.0,
            anchor="tail",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_curve_fn=_trace_curve,
            evaluate_ncl_display_curve_fn=_evaluate,
        )
        np.testing.assert_allclose(curve, np.array([[0.0, 0.0], [1.0, 0.0]]))
        np.testing.assert_allclose(display_curve, np.array([[0.0, 0.0], [1.0, 0.0]]))

        calls = {"count": 0}

        def _trace_curve_short(**kwargs):
            calls["count"] += 1
            return np.array([[0.0, 0.0], [0.05, 0.0]])

        assert (
            _build_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=1.2,
                anchor="tail",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_curve_fn=_trace_curve_short,
                evaluate_ncl_display_curve_fn=lambda curve, transform, viewport=None: (
                    None,
                    False,
                ),
            )
            is None
        )

        curve, display_curve = _build_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=10.0,
            anchor="tail",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
            trace_ncl_curve_fn=lambda **kwargs: np.array([[0.0, 0.0], [1.0, 0.0]]),
            evaluate_ncl_display_curve_fn=lambda curve, transform, viewport=None: (
                None,
                True,
            ),
        )
        assert display_curve is None

        with pytest.raises(TypeError, match="unexpected"):
            _build_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=10.0,
                anchor="tail",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.0, 0.0], [1.0, 0.0]]),
                    np.array([[0.0, 0.0], [1.0, 0.0]]),
                ),
                evaluate_ncl_display_curve_fn=lambda *args, **kwargs: (
                    _ for _ in ()
                ).throw(TypeError("unexpected failure")),
            )

    def test_tip_display_geometry_uses_explicit_display_curve(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float)
        display_curve = np.array([[10.0, 10.0], [20.0, 10.0]], dtype=float)
        tip, direction = _tip_display_geometry(
            curve,
            Affine2D(),
            2.0,
            display_curve=display_curve,
        )
        np.testing.assert_allclose(tip, np.array([20.0, 10.0]))
        np.testing.assert_allclose(direction, np.array([1.0, 0.0]))

    def test_sample_local_vector_state_handles_sampler_transform_and_norm_failures(
        self,
    ):
        grid = _make_test_grid()
        point = np.array([0.5, 0.5])
        transform = Affine2D()

        assert (
            _sample_local_vector_state(
                grid,
                np.ones(grid.shape),
                np.ones(grid.shape),
                transform,
                point,
                sample_grid_field_fn=lambda grid, field, x, y: None,
                local_display_jacobian_fn=lambda transform, point, grid: np.eye(2),
            )
            is None
        )
        assert (
            _sample_local_vector_state(
                grid,
                np.ones(grid.shape),
                np.ones(grid.shape),
                transform,
                point,
                display_sampler=SimpleNamespace(sample=lambda point: None),
                sample_grid_field_fn=lambda grid, field, x, y: 1.0,
                local_display_jacobian_fn=lambda transform, point, grid: np.eye(2),
            )
            is None
        )
        assert (
            _sample_local_vector_state(
                grid,
                np.ones(grid.shape),
                np.ones(grid.shape),
                transform,
                point,
                sample_grid_field_fn=lambda grid, field, x, y: 1.0,
                local_display_jacobian_fn=lambda transform, point, grid: None,
            )
            is None
        )

        class BadPointTransform:
            def transform(self, values):
                return np.array([[np.nan, np.nan]])

        assert (
            _sample_local_vector_state(
                grid,
                np.ones(grid.shape),
                np.ones(grid.shape),
                BadPointTransform(),
                point,
                sample_grid_field_fn=lambda grid, field, x, y: 1.0,
                local_display_jacobian_fn=lambda transform, point, grid: np.eye(2),
            )
            is None
        )
        assert (
            _sample_local_vector_state(
                grid,
                np.ones(grid.shape),
                np.ones(grid.shape),
                transform,
                point,
                sample_grid_field_fn=lambda grid, field, x, y: 0.0,
                local_display_jacobian_fn=lambda transform, point, grid: np.eye(2),
            )
            is None
        )

        state = _sample_local_vector_state(
            grid,
            np.ones(grid.shape),
            np.ones(grid.shape),
            transform,
            point,
            display_sampler=SimpleNamespace(
                sample=lambda point: (np.array([1.0, 2.0]), np.eye(2))
            ),
            sample_grid_field_fn=lambda grid, field, x, y: 2.0,
            local_display_jacobian_fn=lambda transform, point, grid: np.eye(2),
        )
        assert state[0].tolist() == [1.0, 2.0]
        assert state[3] == pytest.approx(np.hypot(2.0, 2.0))

    def test_grid_validation_covers_remaining_error_paths(self):
        with pytest.raises(ValueError, match="'x' can have at maximum 2 dimensions"):
            Grid(np.ones((2, 2, 2)), np.array([0.0, 1.0]))

        with pytest.raises(ValueError, match="The columns of 'y' must be equal"):
            Grid(
                np.array([0.0, 1.0]),
                np.array([[0.0, 1.0], [0.0, 2.0]]),
            )

        with pytest.raises(ValueError, match="'y' can have at maximum 2 dimensions"):
            Grid(np.array([0.0, 1.0]), np.ones((2, 2, 2)))

        with pytest.raises(ValueError, match="'y' must be strictly increasing"):
            Grid(np.array([0.0, 1.0]), np.array([1.0, 0.0]))

        with pytest.raises(ValueError, match="'y' values must be equally spaced"):
            Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 0.7, 2.0]))

        grid = Grid(
            np.array([0.0, 1.0, 3.0]),
            np.array([0.0, 2.0, 5.0]),
            allow_non_uniform=True,
        )
        assert grid.dx == pytest.approx(1.5)
        assert grid.dy == pytest.approx(2.5)

    def test_curly_vector_ncl_impl_covers_empty_short_and_color_default_paths(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            empty = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.full((2, 2), np.nan),
                np.full((2, 2), np.nan),
                rasterized=True,
                result_cls=CurlyVectorPlotSet,
            )
            assert empty.lines.get_rasterized() is True
            assert len(empty.lines.get_segments()) == 0

            common_kwargs = dict(
                prepare_ncl_display_sampler_fn=lambda grid, transform: None,
                prepare_ncl_native_trace_context_fn=lambda **kwargs: None,
                select_ncl_centers_fn=lambda **kwargs: [(np.array([0.5, 0.5]), 1.0)],
                sample_grid_field_fn=lambda grid, field, x, y: None,
                build_open_arrow_segments_batch_fn=lambda **kwargs: (
                    np.empty((0, 2, 2), dtype=float),
                    np.empty(0, dtype=int),
                ),
                build_filled_arrow_polygons_batch_fn=lambda **kwargs: (
                    np.empty((0, 3, 2), dtype=float),
                    np.empty(0, dtype=int),
                ),
                display_points_to_data_fn=lambda *args, **kwargs: None,
                result_cls=CurlyVectorPlotSet,
            )

            short_curve = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.ones((2, 2)),
                np.ones((2, 2)),
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.5, 0.5]]),
                    np.array([[0.5, 0.5]]),
                ),
                **common_kwargs,
            )
            assert len(short_curve.lines.get_segments()) == 0

            multicolor = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.ones((2, 2)),
                np.ones((2, 2)),
                color=np.array([[1.0, 2.0], [3.0, 4.0]]),
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                ),
                **common_kwargs,
            )
            assert len(multicolor.lines.get_segments()) == 1
        finally:
            plt.close(fig)

    def test_curly_vector_ncl_impl_uses_batched_open_arrow_segments(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            open_head_segments = np.array(
                [
                    [[0.70, 0.70], [0.75, 0.75]],
                    [[0.80, 0.80], [0.75, 0.75]],
                ],
                dtype=float,
            )
            batch_builder = Mock(
                return_value=(open_head_segments, np.array([0, 0], dtype=int))
            )

            result = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.ones((2, 2)),
                np.ones((2, 2)),
                arrowstyle="->",
                prepare_ncl_display_sampler_fn=lambda grid, transform: None,
                prepare_ncl_native_trace_context_fn=lambda **kwargs: None,
                select_ncl_centers_fn=lambda **kwargs: [(np.array([0.5, 0.5]), 1.0)],
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                ),
                sample_grid_field_fn=lambda grid, field, x, y: None,
                build_open_arrow_segments_batch_fn=batch_builder,
                build_filled_arrow_polygons_batch_fn=lambda **kwargs: (
                    np.empty((0, 3, 2), dtype=float),
                    np.empty(0, dtype=int),
                ),
                display_points_to_data_fn=lambda *args, **kwargs: None,
                result_cls=CurlyVectorPlotSet,
            )

            assert len(result.lines.get_segments()) == 3
            batch_builder.assert_called_once()
            assert len(batch_builder.call_args.kwargs["display_curves"]) == 1
        finally:
            plt.close(fig)


class TestResultHelpers:
    def _make_plot_set(
        self,
        *,
        magnitude,
        integration_direction="both",
        length_scale=None,
    ):
        axes = SimpleNamespace(bbox=SimpleNamespace(width=200.0, height=100.0))
        return CurlyVectorPlotSet(
            lines=None,
            arrows=(),
            resolution=0.5,
            magnitude=np.asarray(magnitude, dtype=float),
            zorder=1.0,
            transform=None,
            axes=axes,
            linewidth=1.0,
            color="k",
            cmap=None,
            arrowsize=1.0,
            arrowstyle="->",
            start_points=None,
            integration_direction=integration_direction,
            grains=10,
            broken_streamlines=True,
            allow_non_uniform_grid=False,
            density=1.0,
            anchor="center",
            ncl_preset=None,
            length_scale=length_scale,
            rasterized=None,
        )

    def test_curly_vector_plot_set_glyph_length_handles_all_fallbacks(self):
        plot_set = self._make_plot_set(magnitude=[[1.0, 2.0], [3.0, 4.0]])
        assert plot_set.get_scale_factor() == pytest.approx(2.0)
        assert plot_set.scale_value(8.0) == pytest.approx(4.0)
        assert plot_set.unscale_value(4.0) == pytest.approx(8.0)
        assert plot_set.glyph_length_axes_fraction("bad") == pytest.approx(0.0)
        assert plot_set.glyph_length_axes_fraction(-1.0) == pytest.approx(0.0)
        assert plot_set.glyph_length_axes_fraction(2.0) == pytest.approx(0.25)

        scaled = self._make_plot_set(
            magnitude=[[1.0, 2.0], [3.0, 4.0]],
            length_scale={
                "min_mag": 0.0,
                "max_mag": 4.0,
                "min_length_px": 0.0,
                "max_length_px": 40.0,
                "adjust_min": False,
                "ref_mag": 4.0,
                "ref_length_px": 40.0,
                "min_frac_length": 0.0,
            },
        )
        assert scaled.glyph_length_axes_fraction(2.0) == pytest.approx(0.1)

        empty = self._make_plot_set(magnitude=[[np.nan]])
        assert empty.max_magnitude is None
        assert empty.glyph_length_axes_fraction(1.0) == pytest.approx(0.0)


class TestDatasetVectorCoverage:
    def test_dataset_vector_helpers_cover_remaining_validation_paths(self):
        ds = xr.Dataset({"u": (("lat", "lon"), np.ones((2, 3)))})
        dataarray = _get_plot_dataarray(ds, ds["u"], role="u")
        assert dataarray.identical(ds["u"])

        with pytest.raises(ValueError, match="2D array before plotting"):
            _transpose_2d_dataarray_to_dims(
                xr.DataArray(np.array([1.0, 2.0]), dims=("lon",)),
                ("lat", "lon"),
                role="v",
            )

        x = xr.DataArray(np.array([0.0, 1.0, 2.0]), dims=("lon",))
        y = xr.DataArray(np.array([10.0, 20.0]), dims=("lat",))
        u = xr.DataArray(np.arange(6.0).reshape(2, 3), dims=("lat", "lon"))
        v = xr.DataArray(np.arange(6.0, 12.0).reshape(2, 3), dims=("lat", "lon"))

        x_out, y_out, u_out, v_out = _extract_curly_vector_dataset_source(
            xr.Dataset(),
            x,
            y,
            u,
            v,
        )
        np.testing.assert_allclose(x_out, x.values)
        np.testing.assert_allclose(y_out, y.values)
        np.testing.assert_allclose(u_out, u.values)
        np.testing.assert_allclose(v_out, v.values)

        with pytest.raises(ValueError, match="share the same 2D shape"):
            _extract_curly_vector_dataset_source(
                xr.Dataset(),
                x,
                y,
                u,
                xr.DataArray(np.ones((3, 3)), dims=("lat", "lon")),
            )

        with pytest.raises(ValueError, match="does not match the rectilinear x/y grid"):
            _extract_curly_vector_dataset_source(
                xr.Dataset(),
                xr.DataArray(np.array([0.0, 1.0, 2.0, 3.0]), dims=("lon",)),
                y,
                u,
                v,
            )

        with pytest.raises(ValueError, match="same shape"):
            _extract_curly_vector_dataset_source(
                xr.Dataset(),
                xr.DataArray(np.arange(6.0).reshape(2, 3), dims=("lat", "lon")),
                xr.DataArray(np.arange(4.0).reshape(2, 2), dims=("lat", "lon")),
                u,
                v,
            )

        with pytest.raises(ValueError, match="must match the u/v shape"):
            _extract_curly_vector_dataset_source(
                xr.Dataset(),
                xr.DataArray(np.arange(4.0).reshape(2, 2), dims=("lat", "lon")),
                xr.DataArray(np.arange(4.0).reshape(2, 2), dims=("lat", "lon")),
                u,
                v,
            )

    def test_prepare_dataset_style_field_handles_1d_dataarray_and_shape_mismatch(self):
        from_dataarray = _prepare_dataset_style_field(
            xr.DataArray(np.array([1.0, 2.0, 3.0]), dims=("lon",)),
            isel=None,
            expected_shape=(3,),
            vector_dims=("lat", "lon"),
            x_descending=False,
            y_descending=False,
            role="color",
        )
        np.testing.assert_allclose(from_dataarray, np.array([1.0, 2.0, 3.0]))

        original = np.arange(4.0).reshape(2, 2)
        returned = _prepare_dataset_style_field(
            original,
            isel=None,
            expected_shape=(2, 3),
            vector_dims=("lat", "lon"),
            x_descending=False,
            y_descending=False,
            role="linewidth",
        )
        assert returned is original


class TestScatterCoverage:
    def test_scatter_grid_axis_and_display_filter_branches(self):
        assert _as_1d_axis(xr.DataArray(np.ones((2, 2)), dims=("y", "x"))) is None

        grid_state = _prepare_1d_grid_candidates(
            np.array([0.0, np.nan, 2.0]),
            np.array([10.0, 20.0]),
            reference_fields=(np.ones((2, 3), dtype=bool),),
            where=None,
            mask=None,
            resolved_placement="points",
        )
        assert grid_state is not None
        np.testing.assert_allclose(
            grid_state["source_points"],
            np.array(
                [
                    [0.0, 10.0],
                    [2.0, 10.0],
                    [0.0, 20.0],
                    [2.0, 20.0],
                ]
            ),
        )

        class _PartlyNonfiniteTransform(Affine2D):
            def transform(self, values):
                transformed = np.asarray(super().transform(values), dtype=float)
                if transformed.shape[0] > 1:
                    transformed[1, 0] = np.nan
                return transformed

        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            collection = _scatter_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                transform=_PartlyNonfiniteTransform(),
            )
            assert collection.get_offsets().shape[0] == 1
        finally:
            plt.close(fig)

    def test_public_scatter_argument_validation(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            with pytest.raises(
                TypeError, match="expects at least x and y positional arguments"
            ):
                public_scatter()

            with pytest.raises(TypeError, match="requires x and y arguments"):
                public_scatter(ax, np.array([0.0, 1.0]))

            with pytest.raises(
                TypeError, match="received too many positional arguments"
            ):
                public_scatter(
                    ax,
                    np.array([0.0, 1.0]),
                    np.array([0.0, 1.0]),
                    5.0,
                    "k",
                    "extra",
                )
        finally:
            plt.close(fig)

    def test_scatter_engine_helpers_cover_remaining_branches(self):
        assert (
            _resolve_placement(
                None,
                cell_geometry=None,
            )
            == "points"
        )

        with pytest.raises(ValueError, match="placement must be one of"):
            _resolve_placement(
                "invalid",
                cell_geometry={"kind": np.array("rectilinear")},
            )

        with pytest.raises(
            ValueError, match="placement='cells' requires gridded coordinates"
        ):
            _resolve_placement(
                "cells",
                cell_geometry=None,
            )

        candidate_shape = (2, 2)
        source_indices = np.array([0, 2], dtype=int)
        candidate_source_positions = np.array([1, 0, 1], dtype=int)
        retained_indices = np.array([2], dtype=int)

        flattened = _subset_scatter_value(
            np.array([0.0, 1.0, 2.0, 3.0]),
            candidate_shape=candidate_shape,
            source_indices=source_indices,
            candidate_source_positions=candidate_source_positions,
            retained_indices=retained_indices,
            target_dims=None,
        )
        np.testing.assert_allclose(flattened, np.array([2.0]))

        generated = _subset_scatter_value(
            np.array([10.0, 11.0, 12.0]),
            candidate_shape=candidate_shape,
            source_indices=source_indices,
            candidate_source_positions=candidate_source_positions,
            retained_indices=retained_indices,
            target_dims=None,
        )
        np.testing.assert_allclose(generated, np.array([12.0]))

        odd = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        assert (
            _subset_scatter_value(
                odd,
                candidate_shape=candidate_shape,
                source_indices=source_indices,
                candidate_source_positions=candidate_source_positions,
                retained_indices=retained_indices,
                target_dims=None,
            )
            is odd
        )

    def test_generate_cell_candidates_and_scatter_impl_cover_error_paths(self):
        bbox_axes = SimpleNamespace(bbox=Bbox.from_bounds(0.0, 0.0, 100.0, 100.0))
        rect_geometry = {
            "kind": np.array("rectilinear"),
            "x_edges": np.array([0.0, 1.0]),
            "y_edges": np.array([0.0, 1.0]),
        }

        empty_points, empty_sources = _generate_cell_candidates(
            ax=bbox_axes,
            transform=Affine2D(),
            source_points=np.empty((0, 2), dtype=float),
            source_indices=np.empty(0, dtype=int),
            candidate_shape=(1, 1),
            cell_geometry=rect_geometry,
            spacing_fraction=0.1,
        )
        assert empty_points.shape == (0, 2)
        assert empty_sources.shape == (0,)

        with pytest.raises(ValueError, match="Unsupported cell geometry kind"):
            _generate_cell_candidates(
                ax=bbox_axes,
                transform=Affine2D(),
                source_points=np.array([[0.5, 0.5]]),
                source_indices=np.array([0], dtype=int),
                candidate_shape=(1, 1),
                cell_geometry={"kind": np.array("weird")},
                spacing_fraction=0.1,
            )

        fallback_points, fallback_sources = _generate_cell_candidates(
            ax=bbox_axes,
            transform=_BadTransform(),
            source_points=np.array([[0.5, 0.5]]),
            source_indices=np.array([0], dtype=int),
            candidate_shape=(1, 1),
            cell_geometry=rect_geometry,
            spacing_fraction=0.1,
        )
        np.testing.assert_allclose(fallback_points, np.array([[0.5, 0.5]]))
        np.testing.assert_array_equal(fallback_sources, np.array([0]))

        with patch(
            "skyborn.plot._core.scatter_engine._map_ncl_display_points_to_viewport",
            return_value=np.full((4, 2), np.nan),
        ):
            nan_points, nan_sources = _generate_cell_candidates(
                ax=bbox_axes,
                transform=Affine2D(),
                source_points=np.array([[0.5, 0.5]]),
                source_indices=np.array([0], dtype=int),
                candidate_shape=(1, 1),
                cell_geometry=rect_geometry,
                spacing_fraction=0.1,
            )
        np.testing.assert_allclose(nan_points, np.array([[0.5, 0.5]]))
        np.testing.assert_array_equal(nan_sources, np.array([0]))

        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            empty_collection = _scatter_impl(
                ax, np.array([np.nan]), np.array([0.0]), zorder=7.0
            )
            assert empty_collection.get_zorder() == pytest.approx(7.0)

            collection = _scatter_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                zorder=5.0,
            )
            assert collection.get_zorder() == pytest.approx(5.0)
        finally:
            plt.close(fig)


class TestStyleAndGeometryCoverage:
    def test_geometry_style_and_thinning_helpers_cover_remaining_branches(self):
        curve = np.array([[0.0, 0.0], [1.0, 1.0]])

        class _NonfiniteTransform:
            def transform(self, values):
                del values
                return np.array([[0.0, 0.0], [np.nan, 1.0]])

        assert _tip_display_geometry(curve, _NonfiniteTransform(), 1.0) is None
        assert (
            _tip_display_geometry(
                curve,
                Affine2D(),
                1.0,
                display_curve=np.array([[0.0, 0.0], [np.nan, 1.0]]),
            )
            is None
        )

        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            head_length, head_width = _curly_head_axes_dimensions(
                ax,
                3.0,
                max_magnitude=0.0,
                arrowsize=1.0,
            )
            assert head_length > 0.0
            assert head_width > 0.0
        finally:
            plt.close(fig)

        with pytest.raises(ValueError, match="pivot must be one of"):
            _normalize_curly_pivot("corner")

        assert _resolve_curly_anchor_alias("center", "mid") == "center"


class TestLegacyStreamCoverage:
    def test_domain_map_rejects_out_of_grid_trajectory_updates(self):
        grid = SimpleNamespace(nx=3, ny=3, dx=1.0, dy=1.0)
        grid.within_grid = lambda xg, yg: 0.0 <= xg <= 2.0 and 0.0 <= yg <= 2.0
        mask = SimpleNamespace(
            nx=3,
            ny=3,
            _update_trajectory=Mock(),
        )

        domain = DomainMap(grid, mask)
        with pytest.raises(InvalidIndexError):
            domain.update_trajectory(3.0, 1.0)


class TestVectorArtistCoverage:
    def test_batched_open_arrow_helpers_cover_empty_and_failure_paths(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float)

        segments, sources = _build_open_arrow_segments_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_open_arrow_segments_batch_fn=None,
        )
        assert segments.shape == (0, 2, 2)
        assert sources.shape == (0,)

        segments, sources = _build_open_arrow_segments_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_open_arrow_segments_batch_fn=lambda *args: None,
        )
        assert segments.shape == (0, 2, 2)
        assert sources.shape == (0,)

        segments, sources = _build_open_arrow_segments_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_open_arrow_segments_batch_fn=lambda *args: (
                np.empty((0, 2, 2), dtype=float),
                np.empty(0, dtype=int),
            ),
        )
        assert segments.shape == (0, 2, 2)
        assert sources.shape == (0,)

        segments, sources = _build_open_arrow_segments_batch(
            transform=_ExplodingInverse(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_open_arrow_segments_batch_fn=lambda *args: (
                np.array([[[0.0, 0.0], [1.0, 1.0]]], dtype=float),
                np.array([0], dtype=int),
            ),
        )
        assert segments.shape == (0, 2, 2)
        assert sources.shape == (0,)

        segments, sources = _build_open_arrow_segments_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_open_arrow_segments_batch_fn=lambda *args: (
                np.array([[[np.nan, np.nan], [np.nan, np.nan]]], dtype=float),
                np.array([0], dtype=int),
            ),
        )
        assert segments.shape == (0, 2, 2)
        assert sources.shape == (0,)

    def test_batched_filled_arrow_helpers_cover_error_and_empty_paths(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float)

        with pytest.raises(ValueError, match="must match display_curves"):
            _build_filled_arrow_polygons_batch(
                transform=Affine2D(),
                display_curves=[curve],
                head_lengths_px=np.array([1.0, 2.0]),
                head_widths_px=np.array([1.0]),
                build_filled_arrow_polygons_batch_fn=lambda *args: None,
                display_points_to_data_fn=_display_points_to_data,
            )

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[],
            head_lengths_px=np.array([], dtype=float),
            head_widths_px=np.array([], dtype=float),
            build_filled_arrow_polygons_batch_fn=lambda *args: None,
            display_points_to_data_fn=_display_points_to_data,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=None,
            display_points_to_data_fn=_display_points_to_data,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=lambda *args: (
                np.array([[[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]]]),
                np.array([0], dtype=int),
            ),
            display_points_to_data_fn=_display_points_to_data,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=lambda *args: (
                np.array([[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]]]),
                np.array([0], dtype=int),
            ),
            display_points_to_data_fn=lambda *args, **kwargs: np.array(
                [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]],
                dtype=float,
            ),
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=lambda *args: None,
            display_points_to_data_fn=_display_points_to_data,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=lambda *args: (
                np.array([[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]]], dtype=float),
                np.array([0], dtype=int),
            ),
            display_points_to_data_fn=lambda *args, **kwargs: None,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

        polygons, sources = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve],
            head_lengths_px=np.array([1.0]),
            head_widths_px=np.array([1.0]),
            build_filled_arrow_polygons_batch_fn=lambda *args: (
                np.array([[[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]]]),
                np.array([0], dtype=int),
            ),
            display_points_to_data_fn=_display_points_to_data,
        )
        assert polygons.shape == (0, 3, 2)
        assert sources.shape == (0,)

    def test_assemble_filled_head_artists_covers_empty_and_shape_errors(self):
        with pytest.raises(ValueError, match="polygon_sources must match"):
            _assemble_filled_head_artists(
                shafts=[np.array([[0.0, 0.0], [1.0, 0.0]])],
                shaft_colors=["k"],
                shaft_linewidths=[1.0],
                filled_polygons=np.array([[[0.0, 0.0], [1.0, 0.0], [0.5, 0.5]]]),
                polygon_sources=np.array([0, 1], dtype=int),
                facecolors=["k"],
                edgecolors=["k"],
                use_multicolor_lines=True,
                use_linewidth_field=True,
            )

        payload = _assemble_filled_head_artists(
            shafts=[np.array([[0.0, 0.0], [1.0, 0.0]])],
            shaft_colors=["k"],
            shaft_linewidths=[1.0],
            filled_polygons=np.empty((0, 3, 2), dtype=float),
            polygon_sources=np.array([], dtype=int),
            facecolors=["k"],
            edgecolors=["k"],
            use_multicolor_lines=False,
            use_linewidth_field=False,
        )
        assert payload[0]
        assert payload[1] is None
        assert payload[2] is None
        assert payload[3] == []
        assert payload[4] == []
        assert payload[5] == []

    def test_assemble_open_head_streamlines_preserves_source_order(self):
        shafts = [
            np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
            np.array([[2.0, 0.0], [3.0, 0.0]], dtype=float),
        ]
        head_segments = np.array(
            [
                [[0.9, 0.1], [1.0, 0.0]],
                [[0.9, -0.1], [1.0, 0.0]],
                [[2.9, 0.1], [3.0, 0.0]],
                [[2.9, -0.1], [3.0, 0.0]],
            ],
            dtype=float,
        )
        source_positions = np.array([1, 0, 1, 0], dtype=int)

        streamlines, colors, widths = _assemble_open_head_streamlines(
            shafts=shafts,
            shaft_colors=["red", "blue"],
            shaft_linewidths=[1.5, 2.0],
            head_segments=head_segments,
            source_positions=source_positions,
            use_multicolor_lines=True,
            use_linewidth_field=True,
        )

        assert len(streamlines) == 6
        np.testing.assert_allclose(streamlines[0], shafts[0])
        np.testing.assert_allclose(streamlines[1], head_segments[1])
        np.testing.assert_allclose(streamlines[2], head_segments[3])
        np.testing.assert_allclose(streamlines[3], shafts[1])
        np.testing.assert_allclose(streamlines[4], head_segments[0])
        np.testing.assert_allclose(streamlines[5], head_segments[2])
        assert colors == ["red", "red", "red", "blue", "blue", "blue"]
        assert widths == [1.5, 1.5, 1.5, 2.0, 2.0, 2.0]

    def test_assemble_open_head_streamlines_rejects_bad_source_shape(self):
        with pytest.raises(ValueError, match="source_positions must match"):
            _assemble_open_head_streamlines(
                shafts=[np.array([[0.0, 0.0], [1.0, 0.0]])],
                shaft_colors=["k"],
                shaft_linewidths=[1.0],
                head_segments=np.array([[[0.0, 0.0], [1.0, 0.0]]], dtype=float),
                source_positions=np.array([0, 1], dtype=int),
                use_multicolor_lines=False,
                use_linewidth_field=False,
            )

    def test_assemble_filled_head_artists_builds_polygon_styles(self):
        shafts = [
            np.array([[0.0, 0.0], [1.0, 0.0]], dtype=float),
            np.array([[2.0, 0.0], [3.0, 0.0]], dtype=float),
        ]
        polygons = np.array(
            [
                [[1.0, 0.0], [0.9, 0.1], [0.9, -0.1]],
                [[3.0, 0.0], [2.9, 0.1], [2.9, -0.1]],
            ],
            dtype=float,
        )

        (
            streamlines,
            line_colors,
            line_widths,
            polygon_facecolors,
            polygon_edgecolors,
            polygon_linewidths,
        ) = _assemble_filled_head_artists(
            shafts=shafts,
            shaft_colors=["red", "blue"],
            shaft_linewidths=[1.0, 3.0],
            filled_polygons=polygons,
            polygon_sources=np.array([1, 0], dtype=int),
            facecolors=["pink", "cyan"],
            edgecolors=["darkred", "darkblue"],
            use_multicolor_lines=True,
            use_linewidth_field=True,
        )

        assert len(streamlines) == 2
        np.testing.assert_allclose(streamlines[0], shafts[0])
        np.testing.assert_allclose(streamlines[1], shafts[1])
        assert line_colors == ["red", "blue"]
        assert line_widths == [1.0, 3.0]
        assert polygon_facecolors == ["cyan", "pink"]
        assert polygon_edgecolors == ["darkblue", "darkred"]
        assert polygon_linewidths == [1.5, 0.5]

    def test_batched_open_arrow_helpers_match_scalar_geometry(self):
        curve_1 = np.array([[0.0, 0.0], [2.0, 0.0], [3.0, 0.0]], dtype=float)
        curve_2 = np.array([[1.0, 1.0], [2.0, 2.0]], dtype=float)
        head_lengths = np.array([1.0, 0.5], dtype=float)
        head_widths = np.array([2.0, 1.0], dtype=float)

        display_points = np.vstack([curve_1, curve_2])
        curve_offsets = np.array([0, len(curve_1), len(curve_1) + len(curve_2)])
        expected_1 = open_arrow_geometry_reference(
            curve_1,
            head_lengths[0],
            head_widths[0],
        )
        expected_2 = open_arrow_geometry_reference(
            curve_2,
            head_lengths[1],
            head_widths[1],
        )
        native_builder = partial(
            vector_module._native_helpers._call_native_build_open_arrow_segments,
            vector_module._build_open_arrow_segments_native,
        )

        segments, source_positions = _build_open_arrow_segments_batch(
            transform=Affine2D(),
            display_curves=[curve_1, curve_2],
            head_lengths_px=head_lengths,
            head_widths_px=head_widths,
            build_open_arrow_segments_batch_fn=native_builder,
        )

        assert segments.shape == (4, 2, 2)
        np.testing.assert_array_equal(source_positions, np.array([0, 0, 1, 1]))
        np.testing.assert_allclose(
            segments[0],
            np.vstack([expected_1["left_display"], expected_1["tip_display"]]),
        )
        np.testing.assert_allclose(
            segments[1],
            np.vstack([expected_1["right_display"], expected_1["tip_display"]]),
        )
        np.testing.assert_allclose(
            segments[2],
            np.vstack([expected_2["left_display"], expected_2["tip_display"]]),
        )
        np.testing.assert_allclose(
            segments[3],
            np.vstack([expected_2["right_display"], expected_2["tip_display"]]),
        )

        native_segments = build_open_arrow_segments(
            display_points=display_points,
            curve_offsets=curve_offsets,
            head_lengths_px=head_lengths,
            head_widths_px=head_widths,
        )
        native_display_segments, native_sources = native_segments
        np.testing.assert_array_equal(native_sources, np.array([0, 0, 1, 1]))
        np.testing.assert_allclose(native_display_segments, segments)

    def test_batched_filled_arrow_helpers_match_scalar_geometry(self):
        curve_1 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 1.0]], dtype=float)
        curve_2 = np.array([[2.0, 0.0], [3.0, 0.0]], dtype=float)
        head_lengths = np.array([0.8, 0.5], dtype=float)
        head_widths = np.array([1.2, 0.8], dtype=float)

        display_points = np.vstack([curve_1, curve_2])
        curve_offsets = np.array([0, len(curve_1), len(curve_1) + len(curve_2)])
        polygons, sources = build_filled_arrow_polygons(
            display_points=display_points,
            curve_offsets=curve_offsets,
            head_lengths_px=head_lengths,
            head_widths_px=head_widths,
        )

        assert polygons.shape == (2, 3, 2)
        np.testing.assert_array_equal(sources, np.array([0, 1]))

        expected_1 = filled_arrow_geometry_reference(
            curve_1,
            head_lengths[0],
            head_widths[0],
        )
        expected_2 = filled_arrow_geometry_reference(
            curve_2,
            head_lengths[1],
            head_widths[1],
        )
        np.testing.assert_allclose(
            polygons[0],
            np.vstack(
                [
                    expected_1["tip_display"],
                    expected_1["left_display"],
                    expected_1["right_display"],
                ]
            ),
        )
        np.testing.assert_allclose(
            polygons[1],
            np.vstack(
                [
                    expected_2["tip_display"],
                    expected_2["left_display"],
                    expected_2["right_display"],
                ]
            ),
        )
        native_builder = partial(
            vector_module._native_helpers._call_native_build_filled_arrow_polygons,
            vector_module._build_filled_arrow_polygons_native,
        )

        batched_polygons, source_positions = _build_filled_arrow_polygons_batch(
            transform=Affine2D(),
            display_curves=[curve_1, curve_2],
            head_lengths_px=head_lengths,
            head_widths_px=head_widths,
            inverse_transform=Affine2D().inverted(),
            build_filled_arrow_polygons_batch_fn=native_builder,
            display_points_to_data_fn=_display_points_to_data,
        )
        np.testing.assert_array_equal(source_positions, np.array([0, 1]))
        np.testing.assert_allclose(batched_polygons, polygons)

    def test_vector_engine_uses_batched_filled_head_path(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            polygon = np.array(
                [[[0.70, 0.75], [0.65, 0.70], [0.75, 0.70]]],
                dtype=float,
            )
            batch_builder = Mock(return_value=(polygon, np.array([0], dtype=int)))

            result = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.ones((2, 2)),
                np.ones((2, 2)),
                arrowstyle="-|>",
                prepare_ncl_display_sampler_fn=lambda grid, transform: None,
                prepare_ncl_native_trace_context_fn=lambda **kwargs: None,
                select_ncl_centers_fn=lambda **kwargs: [(np.array([0.5, 0.5]), 1.0)],
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                ),
                sample_grid_field_fn=lambda grid, field, x, y: None,
                build_open_arrow_segments_batch_fn=None,
                build_filled_arrow_polygons_batch_fn=batch_builder,
                display_points_to_data_fn=lambda *args, **kwargs: None,
                result_cls=CurlyVectorPlotSet,
            )

            assert len(result.lines.get_segments()) == 1
            assert len(result.arrows) == 1
            assert isinstance(result.arrows[0], PolyCollection)
            batch_builder.assert_called_once()
            assert len(batch_builder.call_args.kwargs["display_curves"]) == 1
        finally:
            plt.close(fig)

    def test_vector_engine_rejects_missing_batch_builders_and_bad_arrowstyle(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            common = dict(
                axes=ax,
                x=np.array([0.0, 1.0]),
                y=np.array([0.0, 1.0]),
                u=np.ones((2, 2)),
                v=np.ones((2, 2)),
                prepare_ncl_display_sampler_fn=lambda grid, transform: None,
                prepare_ncl_native_trace_context_fn=lambda **kwargs: None,
                select_ncl_centers_fn=lambda **kwargs: [(np.array([0.5, 0.5]), 1.0)],
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                ),
                sample_grid_field_fn=lambda grid, field, x, y: None,
                display_points_to_data_fn=lambda *args, **kwargs: None,
                result_cls=CurlyVectorPlotSet,
            )
            with pytest.raises(RuntimeError, match="Open-arrow rendering requires"):
                _curly_vector_ncl_impl(
                    arrowstyle="->",
                    build_open_arrow_segments_batch_fn=None,
                    build_filled_arrow_polygons_batch_fn=lambda **kwargs: None,
                    **common,
                )

            with pytest.raises(RuntimeError, match="Filled-arrow rendering requires"):
                _curly_vector_ncl_impl(
                    arrowstyle="-|>",
                    build_open_arrow_segments_batch_fn=lambda **kwargs: None,
                    build_filled_arrow_polygons_batch_fn=None,
                    **common,
                )

            with pytest.raises(ValueError, match="Unsupported arrowstyle"):
                _curly_vector_ncl_impl(
                    arrowstyle="x",
                    build_open_arrow_segments_batch_fn=lambda **kwargs: None,
                    build_filled_arrow_polygons_batch_fn=lambda **kwargs: None,
                    **common,
                )
        finally:
            plt.close(fig)

    def test_vector_engine_uses_display_curve_typeerror_fallback(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            curve = np.array([[0.25, 0.25], [0.75, 0.75]], dtype=float)
            calls = {"count": 0}

            def _trace_curve(**kwargs):
                return curve, curve

            def _evaluate(curve_arg, transform, viewport=None):
                calls["count"] += 1
                return np.asarray(curve_arg, dtype=float), False

            built = _build_ncl_curve(
                start_point=np.array([0.0, 0.0]),
                total_length_px=5.0,
                anchor="tail",
                grid=_make_test_grid(),
                u=np.ones((3, 3)),
                v=np.ones((3, 3)),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=Bbox.from_bounds(0.0, 0.0, 10.0, 10.0),
                trace_ncl_curve_fn=_trace_curve,
                evaluate_ncl_display_curve_fn=_evaluate,
            )
            assert built is not None
            assert calls["count"] == 1
        finally:
            plt.close(fig)

    def test_vector_engine_extends_filled_line_styles(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            polygon = np.array(
                [[[0.70, 0.75], [0.65, 0.70], [0.75, 0.70]]],
                dtype=float,
            )
            result = _curly_vector_ncl_impl(
                ax,
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
                np.ones((2, 2)),
                np.ones((2, 2)),
                arrowstyle="-|>",
                color=np.array([[1.0, 2.0], [3.0, 4.0]]),
                linewidth=np.array([[1.0, 2.0], [3.0, 4.0]]),
                prepare_ncl_display_sampler_fn=lambda grid, transform: None,
                prepare_ncl_native_trace_context_fn=lambda **kwargs: None,
                select_ncl_centers_fn=lambda **kwargs: [(np.array([0.5, 0.5]), 1.0)],
                build_ncl_curve_fn=lambda **kwargs: (
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                    np.array([[0.25, 0.25], [0.75, 0.75]]),
                ),
                sample_grid_field_fn=lambda grid, field, x, y: 2.5,
                build_open_arrow_segments_batch_fn=None,
                build_filled_arrow_polygons_batch_fn=lambda **kwargs: (
                    polygon,
                    np.array([0], dtype=int),
                ),
                display_points_to_data_fn=lambda *args, **kwargs: None,
                result_cls=CurlyVectorPlotSet,
            )
            assert len(result.lines.get_segments()) == 1
            assert len(result.arrows) == 1
        finally:
            plt.close(fig)

    def test_vector_engine_adds_patch_for_non_collection_arrow_artist(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            added_patches = []
            original_add_patch = ax.add_patch

            def _record_patch(artist):
                added_patches.append(artist)
                return original_add_patch(artist)

            ax.add_patch = _record_patch
            ax.add_patch(Circle((0.5, 0.5), 0.1))
            assert added_patches
        finally:
            plt.close(fig)


class TestVectorKeyCoverage:
    def test_vector_key_helpers_cover_remaining_layout_and_fallback_paths(self):
        fig, ax = plt.subplots(figsize=(6, 4))
        try:
            with pytest.raises(ValueError, match="finite axes coordinates"):
                CurlyVectorKey(
                    ax=ax,
                    curly_vector_set=_make_mock_key_set(),
                    U=2.0,
                    x=np.nan,
                    y=0.1,
                )

            legend_s = CurlyVectorKey(
                ax=ax,
                curly_vector_set=_make_mock_key_set(),
                U=2.5,
                description="Reference Vector",
                labelpos="S",
            )
            assert legend_s._format_reference_value() == "2.5 m/s"
            assert legend_s._resolve_text_blocks() == (
                "Reference Vector",
                "2.5 m/s",
                None,
            )

            with patch(
                "skyborn.plot._artists.vector_key_artist._curly_head_axes_dimensions",
                side_effect=RuntimeError("bad head geometry"),
            ):
                head_length, head_width = legend_s._calculate_head_geometry()
            assert head_length > 0.0
            assert head_width > 0.0

            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

            filled = CurlyVectorKey(
                ax=ax,
                curly_vector_set=_make_mock_key_set(arrowstyle="-|>"),
                U=4.0,
                labelpos="W",
                arrow_props={"arrowstyle": "-|>"},
            )
            filled._update_arrow_geometry(0.1, 0.3, 0.4, 0.06, 0.05)
            assert filled.head_left.get_visible() is False
            assert filled.head_right.get_visible() is False
            assert filled.head_fill.get_visible() is True

            filled._layout_vertical(renderer, None, "Bottom")
            assert filled.text.get_visible() is False
            assert filled.text2.get_visible() is True

            filled._layout_horizontal(renderer, "Side")
            assert filled.text.get_horizontalalignment() == "right"

            bad_fraction = CurlyVectorKey(
                ax=ax,
                curly_vector_set=_make_mock_key_set(),
                U=4.0,
                reference_speed=2.0,
                max_arrow_length=0.08,
            )
            bad_fraction.curly_vector_set.glyph_length_axes_fraction.side_effect = (
                ValueError("bad")
            )
            bad_fraction.curly_vector_set.max_magnitude = "bad"
            bad_fraction.curly_vector_set.ref_length_fraction = object()
            assert bad_fraction._calculate_arrow_length() == pytest.approx(0.16)
        finally:
            plt.close(fig)


class TestVectorWrapperCoverage:
    def test_vector_module_trace_wrappers_cover_display_branches(self, monkeypatch):
        grid = _make_test_grid()
        viewport = Bbox.from_bounds(0.0, 0.0, 10.0, 10.0)
        curve = np.array([[0.0, 0.0], [1.0, 0.0]])
        display_curve = np.array([[2.0, 2.0], [3.0, 2.0]])
        original_trace_with_display = vector_module._trace_ncl_direction_with_display

        wrapped = _trace_ncl_curve(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="tail",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=viewport,
            trace_ncl_direction_fn=lambda *args, **kwargs: curve,
        )
        np.testing.assert_allclose(wrapped, curve)

        assert (
            vector_module._trace_ncl_curve_with_display(
                start_point=np.array([0.0, 0.0]),
                total_length_px=0.0,
                anchor="tail",
                grid=grid,
                u=np.ones(grid.shape),
                v=np.ones(grid.shape),
                transform=Affine2D(),
                step_px=1.0,
                speed_scale=1.0,
                viewport=viewport,
            )
            is None
        )

        forward_pair = (curve, display_curve)
        trace_with_display = Mock(side_effect=[None, forward_pair])
        monkeypatch.setattr(
            vector_module,
            "_trace_ncl_direction_with_display",
            trace_with_display,
        )
        forward = vector_module._trace_ncl_curve_with_display(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="center",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=viewport,
        )
        assert forward is forward_pair

        backward_pair = (curve, display_curve)
        trace_with_display = Mock(side_effect=[backward_pair, None])
        monkeypatch.setattr(
            vector_module,
            "_trace_ncl_direction_with_display",
            trace_with_display,
        )
        backward = vector_module._trace_ncl_curve_with_display(
            start_point=np.array([0.0, 0.0]),
            total_length_px=4.0,
            anchor="center",
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            transform=Affine2D(),
            step_px=1.0,
            speed_scale=1.0,
            viewport=viewport,
        )
        assert backward is not None
        backward_curve, backward_display = backward
        np.testing.assert_allclose(backward_curve, curve[::-1])
        np.testing.assert_allclose(backward_display, display_curve[::-1])

        monkeypatch.setattr(
            vector_module,
            "_trace_ncl_direction_with_display",
            original_trace_with_display,
        )
        sampler = SimpleNamespace()
        context = object()

        def prepare_context(**kwargs):
            assert kwargs["display_sampler"] is sampler
            return context

        def trace_via_native(**kwargs):
            assert kwargs["native_trace_context"] is context
            return curve, display_curve

        monkeypatch.setattr(
            vector_module,
            "_prepare_ncl_native_trace_context",
            prepare_context,
        )
        monkeypatch.setattr(
            vector_module._native_helpers,
            "_call_native_trace_ncl_direction_with_display",
            lambda *args, **kwargs: trace_via_native(
                native_trace_context=(
                    kwargs.get("native_trace_context", args[1])
                    if len(args) > 1
                    else None
                )
            ),
        )
        traced = vector_module._trace_ncl_direction_with_display(
            start_point=np.array([0.0, 0.0]),
            max_length_px=4.0,
            direction_sign=1.0,
            grid=grid,
            u=np.ones(grid.shape),
            v=np.ones(grid.shape),
            step_px=1.0,
            speed_scale=1.0,
            viewport=viewport,
            display_sampler=sampler,
        )
        assert traced is not None

    def test_regrid_non_uniform_vectors_to_uniform_covers_import_and_validation_paths(
        self,
    ):
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([10.0, 20.0])
        u = np.ones((2, 3))
        v = np.ones((2, 3)) * 2.0

        real_import = builtins.__import__

        def _raising_import(name, *args, **kwargs):
            if name == "scipy.interpolate":
                raise ImportError("missing scipy")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_raising_import):
            with pytest.raises(
                ImportError, match="scipy is required for non-uniform grid support"
            ):
                _regrid_non_uniform_vectors_to_uniform(x, y, u, v)

        fake_interp_module = types.ModuleType("scipy.interpolate")
        fake_interp_module.RegularGridInterpolator = object
        fake_scipy = types.ModuleType("scipy")
        fake_scipy.interpolate = fake_interp_module

        with patch.dict(
            sys.modules,
            {"scipy": fake_scipy, "scipy.interpolate": fake_interp_module},
        ):
            with pytest.raises(
                ValueError, match="must match the non-uniform grid shape"
            ):
                _regrid_non_uniform_vectors_to_uniform(x, y, np.ones((2, 2)), v)

            with pytest.raises(ValueError, match="Non-uniform scalar style fields"):
                _regrid_non_uniform_vectors_to_uniform(x, y, u, v, np.ones((2, 2)))

            with pytest.raises(
                ValueError,
                match="x coordinates must be strictly monotonic after sorting",
            ):
                _regrid_non_uniform_vectors_to_uniform(
                    np.array([1.0, 1.0, 2.0]),
                    y,
                    u,
                    v,
                )

            with pytest.raises(
                ValueError,
                match="y coordinates must be strictly monotonic after sorting",
            ):
                _regrid_non_uniform_vectors_to_uniform(
                    x,
                    np.array([10.0, 10.0]),
                    u,
                    v,
                )

    def test_array_curly_vector_and_public_wrapper_cover_remaining_argument_paths(
        self,
    ):
        fig, ax = plt.subplots(figsize=(4, 3))
        try:
            x = np.array([0.0, 2.0, 5.0])
            y = np.array([1000.0, 500.0])
            u = np.ones((2, 3))
            v = np.ones((2, 3)) * 2.0

            with pytest.raises(ValueError, match="If 'color' is given"):
                _array_curly_vector(
                    ax,
                    x,
                    y,
                    u,
                    v,
                    allow_non_uniform_grid=True,
                    color=np.ones((2, 2)),
                )

            with pytest.raises(ValueError, match="If 'linewidth' is given"):
                _array_curly_vector(
                    ax,
                    x,
                    y,
                    u,
                    v,
                    allow_non_uniform_grid=True,
                    linewidth=np.ones((2, 2)),
                )

            mock_result = Mock(spec=CurlyVectorPlotSet)
            with patch(
                "skyborn.plot.vector._curly_vector_from_arrays",
                return_value=mock_result,
            ) as mock_arrays:
                result = curly_vector(ax, x, y, u, v, density=0.8)

            assert result is mock_result
            assert mock_arrays.call_args.args[0] is ax
        finally:
            plt.close(fig)
