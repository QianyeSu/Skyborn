"""Focused regression tests for low-level plot helper branches."""

from __future__ import annotations

import builtins
import sys
import types
from types import SimpleNamespace
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.transforms import Affine2D, Bbox

from skyborn.plot._adapters.cartopy_vector import (
    _build_projection_target_grid,
    _extract_regular_grid_from_regridded_vectors,
    _regrid_cartopy_vectors,
)
from skyborn.plot._adapters.curly_vector_entry import (
    _prepare_curly_vector_dataset_inputs_impl,
)
from skyborn.plot._adapters.grid_prepare import (
    _build_curvilinear_target_grid,
    _maybe_as_scalar_field,
    _prepare_source_vector_grid,
    _rcm2rgrid_2d,
    _rcm2rgrid_fields,
    _regrid_curvilinear_vectors,
    _wrap_periodic_grid_queries,
)
from skyborn.plot._artists.vector_artists import (
    _build_arrow_polygon,
    _build_ncl_arrow_artists,
    _build_open_arrow_segments,
    _open_arrow_geometry,
    _trim_curve_for_open_head,
    _uses_open_arrow_head,
)
from skyborn.plot._core.geometry import (
    _acceptable_ncl_display_curve,
    _candidate_data_from_display_step,
    _chaikin_smooth_display_curve,
    _curve_shape_is_acceptable,
    _curve_to_display,
    _display_points_to_data,
    _display_step_to_data,
    _display_to_data,
    _evaluate_ncl_display_curve,
    _finite_difference_step,
    _fit_single_bend_display_curve,
    _local_display_jacobian,
    _point_at_arc_distance_from_end,
    _point_at_arc_fraction,
    _postprocess_ncl_curve,
    _tip_display_geometry,
    _tip_display_geometry_from_display_curve,
    _trim_display_curve_from_end,
)
from skyborn.plot._core.native import (
    _try_native_sample_grid_field,
    _try_native_sample_grid_field_array,
    _try_native_thin_ncl_mapped_candidates,
    _try_native_trace_ncl_direction,
)
from skyborn.plot._core.result import CurlyVectorPlotSet
from skyborn.plot._core.sampling import (
    _sample_grid_field_array_python,
    _sample_grid_field_python,
)
from skyborn.plot._core.vector_engine import (
    Grid,
    _build_ncl_curve,
    _prepare_ncl_center_candidates,
    _resolve_artist_coordinate_context,
    _resolve_ncl_length_scale,
    _resolve_ncl_reference_length_px,
    _sample_local_vector_state,
    _trace_ncl_curve,
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
    _normalize_supported_arrowstyle,
    _resolve_curly_anchor_alias,
    _resolve_curly_style_aliases,
    _resolve_scatter_aliases,
)


class _StopSampling(Exception):
    """Marker exception for masked sampling tests."""


class _BadTransform:
    """Transform stub that always fails."""

    def transform(self, values):
        raise RuntimeError("bad transform")

    def inverted(self):
        raise RuntimeError("bad inverse")


def _make_test_grid() -> Grid:
    return Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))


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

    def test_rcm2rgrid_helpers_handle_import_errors_and_shape_normalization(self):
        real_import = builtins.__import__

        def _raising_import(name, *args, **kwargs):
            if name == "skyborn.interp":
                raise ImportError("missing interp")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_raising_import):
            with pytest.raises(ImportError, match="requires skyborn.interp.rcm2rgrid"):
                _rcm2rgrid_2d(
                    np.zeros((2, 2)),
                    np.zeros((2, 2)),
                    np.zeros((2, 2)),
                    np.array([0.0, 1.0]),
                    np.array([0.0, 1.0]),
                )
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
            regridded_2d = _rcm2rgrid_2d(
                np.zeros((2, 2)),
                np.zeros((2, 2)),
                np.array([[1.0, 2.0], [3.0, 4.0]]),
                np.array([0.0, 1.0]),
                np.array([0.0, 1.0]),
            )
            np.testing.assert_allclose(
                regridded_2d,
                np.array([[1.0, 2.0], [3.0, 4.0]]),
            )

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


class TestSamplingHelpers:
    def test_sample_grid_field_python_handles_out_of_bounds_and_masked_paths(self):
        grid = _make_test_grid()
        field = np.arange(9.0).reshape(3, 3)

        assert (
            _sample_grid_field_python(
                grid,
                field,
                -1.0,
                0.5,
                interpgrid_fn=lambda field, xi, yi: 0.0,
                terminate_trajectory_exc=_StopSampling,
            )
            is None
        )

        masked = np.ma.array(field, mask=False)

        assert (
            _sample_grid_field_python(
                grid,
                masked,
                0.5,
                0.5,
                interpgrid_fn=lambda field, xi, yi: (_ for _ in ()).throw(
                    _StopSampling()
                ),
                terminate_trajectory_exc=_StopSampling,
            )
            is None
        )
        assert (
            _sample_grid_field_python(
                grid,
                masked,
                0.5,
                0.5,
                interpgrid_fn=lambda field, xi, yi: np.ma.masked,
                terminate_trajectory_exc=_StopSampling,
            )
            is None
        )
        assert _sample_grid_field_python(
            grid,
            masked,
            0.5,
            0.5,
            interpgrid_fn=lambda field, xi, yi: 7.0,
            terminate_trajectory_exc=_StopSampling,
        ) == pytest.approx(7.0)

    def test_sample_grid_field_python_returns_none_for_nonfinite_bilinear_value(self):
        grid = _make_test_grid()
        field = np.arange(9.0).reshape(3, 3)
        field[1, 1] = np.nan

        assert (
            _sample_grid_field_python(
                grid,
                field,
                0.5,
                0.5,
                interpgrid_fn=lambda field, xi, yi: 0.0,
                terminate_trajectory_exc=_StopSampling,
            )
            is None
        )

    def test_sample_grid_field_array_python_handles_empty_1d_and_invalid_points(self):
        grid = _make_test_grid()
        field = np.arange(9.0).reshape(3, 3)

        empty = _sample_grid_field_array_python(
            grid,
            field,
            np.empty((0, 2)),
            interpgrid_fn=lambda field, xi, yi: np.array([], dtype=float),
        )
        assert empty.size == 0

        one_point = _sample_grid_field_array_python(
            grid,
            field,
            np.array([0.5, 0.5]),
            interpgrid_fn=lambda field, xi, yi: np.array([1.25], dtype=float),
        )
        np.testing.assert_allclose(one_point, np.array([1.25]))

        invalid = _sample_grid_field_array_python(
            grid,
            field,
            np.array([[-5.0, -5.0], [9.0, 9.0]]),
            interpgrid_fn=lambda field, xi, yi: np.array([], dtype=float),
        )
        assert np.isnan(invalid).all()

    def test_sample_grid_field_array_python_replaces_nonfinite_samples(self):
        grid = _make_test_grid()
        field = np.arange(9.0).reshape(3, 3)
        points = np.array([[0.5, 0.5], [1.0, 1.0]])

        sampled = _sample_grid_field_array_python(
            grid,
            field,
            points,
            interpgrid_fn=lambda field, xi, yi: np.ma.array([3.5, np.nan]),
        )

        assert sampled[0] == pytest.approx(3.5)
        assert np.isnan(sampled[1])


class TestNativeHelpers:
    def test_try_native_sample_grid_field_handles_error_and_nonfinite_results(self):
        grid = SimpleNamespace(x_origin=1.0, y_origin=2.0, dx=3.0, dy=4.0)
        on_error = Mock()

        assert (
            _try_native_sample_grid_field(
                None,
                grid,
                np.ones((2, 2)),
                1.0,
                2.0,
                on_error,
            )
            is None
        )
        assert (
            _try_native_sample_grid_field(
                lambda **kwargs: 1.0,
                grid,
                np.ma.array(np.ones((2, 2)), mask=True),
                1.0,
                2.0,
                on_error,
            )
            is None
        )
        assert (
            _try_native_sample_grid_field(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                grid,
                np.ones((2, 2)),
                1.0,
                2.0,
                on_error,
            )
            is None
        )
        on_error.assert_called_once()
        assert (
            _try_native_sample_grid_field(
                lambda **kwargs: np.nan,
                grid,
                np.ones((2, 2)),
                1.0,
                2.0,
                Mock(),
            )
            is None
        )
        assert _try_native_sample_grid_field(
            lambda **kwargs: 5.5,
            grid,
            np.ones((2, 2)),
            1.0,
            2.0,
            Mock(),
        ) == pytest.approx(5.5)

    def test_try_native_sample_grid_field_array_handles_shape_errors_and_nonfinite_values(
        self,
    ):
        grid = SimpleNamespace(x_origin=1.0, y_origin=2.0, dx=3.0, dy=4.0)
        on_error = Mock()
        points = np.array([[0.0, 0.0], [1.0, 1.0]])

        assert (
            _try_native_sample_grid_field_array(
                None, grid, np.ones((2, 2)), points, (2,), on_error
            )
            is None
        )
        assert (
            _try_native_sample_grid_field_array(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                grid,
                np.ones((2, 2)),
                points,
                (2,),
                on_error,
            )
            is None
        )
        on_error.assert_called_once()
        assert (
            _try_native_sample_grid_field_array(
                lambda **kwargs: None,
                grid,
                np.ones((2, 2)),
                points,
                (2,),
                Mock(),
            )
            is None
        )
        assert (
            _try_native_sample_grid_field_array(
                lambda **kwargs: np.array([[1.0, 2.0]]),
                grid,
                np.ones((2, 2)),
                points,
                (2,),
                Mock(),
            )
            is None
        )
        sampled = _try_native_sample_grid_field_array(
            lambda **kwargs: np.array([1.0, np.nan]),
            grid,
            np.ones((2, 2)),
            points,
            (2,),
            Mock(),
        )
        assert sampled[0] == pytest.approx(1.0)
        assert np.isnan(sampled[1])

    def test_try_native_thin_ncl_mapped_candidates_handles_error_and_success(self):
        on_error = Mock()
        points = np.array([[0.1, 0.2], [0.3, 0.4]])

        assert (
            _try_native_thin_ncl_mapped_candidates(None, points, 0.1, on_error) is None
        )
        assert (
            _try_native_thin_ncl_mapped_candidates(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                points,
                0.1,
                on_error,
            )
            is None
        )
        on_error.assert_called_once()
        assert (
            _try_native_thin_ncl_mapped_candidates(
                lambda **kwargs: None,
                points,
                0.1,
                Mock(),
            )
            is None
        )
        assert _try_native_thin_ncl_mapped_candidates(
            lambda **kwargs: np.array([1, 0], dtype=int),
            points,
            0.1,
            Mock(),
        ) == [1, 0]

    def test_try_native_trace_ncl_direction_handles_missing_context_error_and_shape_validation(
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
        on_error = Mock()

        assert (
            _try_native_trace_ncl_direction(
                None,
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
                on_error,
            )
            is None
        )
        assert (
            _try_native_trace_ncl_direction(
                lambda **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
                on_error,
            )
            is None
        )
        on_error.assert_called_once()
        assert (
            _try_native_trace_ncl_direction(
                lambda **kwargs: np.array([1.0, 2.0]),
                context,
                np.array([0.0, 0.0]),
                10.0,
                1.0,
                1.0,
                1.0,
                Mock(),
            )
            is None
        )
        curve = _try_native_trace_ncl_direction(
            lambda **kwargs: np.array([[0.0, 0.0], [1.0, 1.0]]),
            context,
            np.array([0.0, 0.0]),
            10.0,
            1.0,
            1.0,
            1.0,
            Mock(),
        )
        np.testing.assert_allclose(curve, np.array([[0.0, 0.0], [1.0, 1.0]]))


class TestGeometryHelpers:
    def test_curve_to_display_uses_sampler_then_falls_back_to_transform(self):
        curve = np.array([[0.0, 0.0], [1.0, 1.0]])
        transform = Affine2D().scale(2.0)

        sampler = SimpleNamespace(sample_display_points=lambda curve: curve + 3.0)
        np.testing.assert_allclose(
            _curve_to_display(curve, transform, sampler), curve + 3.0
        )

        sampler = SimpleNamespace(
            sample_display_points=lambda curve: (_ for _ in ()).throw(
                RuntimeError("bad")
            )
        )
        np.testing.assert_allclose(
            _curve_to_display(curve, transform, sampler), curve * 2.0
        )

        assert _curve_to_display(curve, _BadTransform(), None) is None

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

    def test_postprocess_curve_and_arc_helpers_cover_edge_cases(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])

        class NoInverseTransform:
            def transform(self, values):
                return np.asarray(values, dtype=float)

            def inverted(self):
                raise RuntimeError("no inverse")

        np.testing.assert_allclose(
            _postprocess_ncl_curve(curve, NoInverseTransform()), curve
        )

        with pytest.raises(ValueError, match="at least one point"):
            _point_at_arc_fraction(np.empty((0, 2)), 0.5)
        with pytest.raises(ValueError, match="at least one point"):
            _point_at_arc_distance_from_end(np.empty((0, 2)), 1.0)

        np.testing.assert_allclose(
            _point_at_arc_fraction(np.array([[5.0, 6.0]]), 0.5),
            np.array([5.0, 6.0]),
        )
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
        assert _curve_shape_is_acceptable(curve, _BadTransform()) is True

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
        assert _postprocess_ncl_curve(np.array([[0.0, 0.0]]), _BadTransform()) is None

        short_curve = np.array([[0.0, 0.0], [1.0, 1.0]])
        np.testing.assert_allclose(
            _chaikin_smooth_display_curve(short_curve),
            short_curve,
        )
        np.testing.assert_allclose(
            _fit_single_bend_display_curve(short_curve),
            short_curve,
        )
        degenerate_curve = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
        np.testing.assert_allclose(
            _fit_single_bend_display_curve(degenerate_curve),
            degenerate_curve,
        )

        flat_curve = np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]])
        np.testing.assert_allclose(
            _point_at_arc_fraction(flat_curve, 0.5), flat_curve[0]
        )

        single_segment_curve = np.array([[0.0, 0.0], [1.0, 0.0]])
        display_curve, transform_failed = _evaluate_ncl_display_curve(
            single_segment_curve,
            Affine2D(),
        )
        assert transform_failed is False
        np.testing.assert_allclose(display_curve, single_segment_curve)
        np.testing.assert_allclose(
            _acceptable_ncl_display_curve(single_segment_curve, Affine2D()),
            single_segment_curve,
        )

        invalid_curve = np.array([[0.0, 0.0]])
        assert _evaluate_ncl_display_curve(invalid_curve, Affine2D()) == (None, False)
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


class TestVectorArtistHelpers:
    def test_uses_open_arrow_head_normalizes_arrowstyle_text(self):
        assert _uses_open_arrow_head("->")
        assert _uses_open_arrow_head(" -> ")
        assert not _uses_open_arrow_head("-|>")

    def test_open_arrow_geometry_handles_failure_and_width_variants(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0]])
        tip_fn = lambda display_curve, backoff: (
            display_curve[-1],
            np.array([1.0, 0.0]),
        )

        assert (
            _open_arrow_geometry(
                curve[:1],
                Affine2D(),
                2.0,
                tip_display_geometry_from_display_curve_fn=tip_fn,
            )
            is None
        )
        assert (
            _open_arrow_geometry(
                curve,
                _BadTransform(),
                2.0,
                tip_display_geometry_from_display_curve_fn=tip_fn,
            )
            is None
        )
        assert (
            _open_arrow_geometry(
                curve,
                Affine2D(),
                2.0,
                display_curve=np.array([[0.0, 0.0], [np.nan, 1.0]]),
                tip_display_geometry_from_display_curve_fn=tip_fn,
            )
            is None
        )
        assert (
            _open_arrow_geometry(
                curve,
                Affine2D(),
                2.0,
                tip_display_geometry_from_display_curve_fn=lambda display_curve, backoff: None,
            )
            is None
        )

        geometry = _open_arrow_geometry(
            curve,
            Affine2D(),
            2.0,
            head_width_px=4.0,
            tip_display_geometry_from_display_curve_fn=tip_fn,
        )
        assert "left_display" in geometry
        assert "right_display" in geometry

    def test_trim_curve_for_open_head_handles_geometry_and_inverse_failures(self):
        curve = np.array([[0.0, 0.0], [3.0, 0.0]])

        np.testing.assert_allclose(
            _trim_curve_for_open_head(
                curve[:1],
                Affine2D(),
                2.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: None,
                trim_display_curve_from_end_fn=lambda display_curve, distance: display_curve,
            ),
            curve[:1],
        )

        np.testing.assert_allclose(
            _trim_curve_for_open_head(
                curve,
                Affine2D(),
                2.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: None,
                trim_display_curve_from_end_fn=lambda display_curve, distance: display_curve,
            ),
            curve,
        )

        np.testing.assert_allclose(
            _trim_curve_for_open_head(
                curve,
                Affine2D(),
                2.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: {
                    "display_curve": curve.copy(),
                    "base_center_display": np.array([2.0, 0.0]),
                },
                trim_display_curve_from_end_fn=lambda display_curve, distance: np.array(
                    [[0.0, 0.0]]
                ),
            ),
            curve,
        )

        bad_inverse = SimpleNamespace(
            inverted=lambda: SimpleNamespace(
                transform=lambda values: (_ for _ in ()).throw(RuntimeError("bad"))
            )
        )
        np.testing.assert_allclose(
            _trim_curve_for_open_head(
                curve,
                bad_inverse,
                2.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: {
                    "display_curve": curve.copy(),
                    "base_center_display": np.array([2.0, 0.0]),
                },
                trim_display_curve_from_end_fn=lambda display_curve, distance: display_curve,
            ),
            curve,
        )

    def test_build_open_arrow_segments_and_polygon_cover_success_and_failure_paths(
        self,
    ):
        curve = np.array([[0.0, 0.0], [3.0, 0.0]])
        geometry = {
            "display_curve": curve.copy(),
            "tip_display": np.array([3.0, 0.0]),
            "left_display": np.array([1.0, 1.0]),
            "right_display": np.array([1.0, -1.0]),
        }

        assert (
            _build_open_arrow_segments(
                curve,
                None,
                Affine2D(),
                2.0,
                4.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: None,
                display_points_to_data_fn=lambda transform, values, inverse_transform=None: values,
            )
            == []
        )
        assert (
            _build_open_arrow_segments(
                curve,
                None,
                Affine2D(),
                2.0,
                4.0,
                open_arrow_geometry_fn=lambda *args, **kwargs: geometry,
                display_points_to_data_fn=lambda transform, values, inverse_transform=None: None,
            )
            == []
        )

        segments = _build_open_arrow_segments(
            curve,
            None,
            Affine2D(),
            2.0,
            4.0,
            open_arrow_geometry_fn=lambda *args, **kwargs: geometry,
            display_points_to_data_fn=lambda transform, values, inverse_transform=None: np.asarray(
                values, dtype=float
            ),
        )
        assert len(segments) == 2

        assert (
            _build_arrow_polygon(
                curve,
                None,
                Affine2D(),
                2.0,
                4.0,
                facecolor="k",
                edgecolor="w",
                linewidth=2.0,
                alpha=0.5,
                zorder=3,
                tip_display_geometry_fn=lambda *args, **kwargs: None,
                display_points_to_data_fn=lambda *args, **kwargs: None,
            )
            is None
        )

        polygon = _build_arrow_polygon(
            curve,
            None,
            Affine2D(),
            2.0,
            4.0,
            facecolor="k",
            edgecolor="w",
            linewidth=0.2,
            alpha=0.5,
            zorder=3,
            tip_display_geometry_fn=lambda *args, **kwargs: (
                np.array([3.0, 0.0]),
                np.array([1.0, 0.0]),
            ),
            display_points_to_data_fn=lambda transform, values, inverse_transform=None: np.asarray(
                values, dtype=float
            ),
        )
        assert polygon.get_linewidth() == pytest.approx(0.5)

    def test_build_ncl_arrow_artists_dispatches_open_and_filled_heads(self):
        head_segments, polygon = _build_ncl_arrow_artists(
            curve=np.array([[0.0, 0.0], [1.0, 0.0]]),
            grid=None,
            transform=Affine2D(),
            arrowstyle="->",
            head_length_px=2.0,
            head_width_px=4.0,
            facecolor="k",
            edgecolor="w",
            linewidth=1.0,
            alpha=0.5,
            zorder=3,
            uses_open_arrow_head_fn=lambda arrowstyle: True,
            build_open_arrow_segments_fn=lambda **kwargs: ["segment"],
            build_arrow_polygon_fn=lambda **kwargs: "polygon",
        )
        assert head_segments == ["segment"]
        assert polygon is None

        head_segments, polygon = _build_ncl_arrow_artists(
            curve=np.array([[0.0, 0.0], [1.0, 0.0]]),
            grid=None,
            transform=Affine2D(),
            arrowstyle="-|>",
            head_length_px=2.0,
            head_width_px=4.0,
            facecolor="k",
            edgecolor="w",
            linewidth=1.0,
            alpha=0.5,
            zorder=3,
            uses_open_arrow_head_fn=lambda arrowstyle: False,
            build_open_arrow_segments_fn=lambda **kwargs: ["segment"],
            build_arrow_polygon_fn=lambda **kwargs: "polygon",
        )
        assert head_segments == []
        assert polygon == "polygon"


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
