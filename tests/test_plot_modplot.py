"""Tests for the low-level curly-vector rendering engine."""

from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.collections import LineCollection
from matplotlib.transforms import Bbox

import skyborn.plot.vector_plot as vector_plot_module
from skyborn.plot.vector_plot import (
    CurlyVectorPlotSet,
    DomainMap,
    Grid,
    InvalidIndexError,
    OutOfBounds,
    StreamMask,
    TerminateTrajectory,
    _apply_ncl_preset_defaults,
    _axis_coordinate_1d,
    _axis_is_uniform,
    _candidate_data_from_display_step,
    _curve_shape_is_acceptable,
    _default_ncl_box_center_candidates,
    _default_ncl_candidate_shape,
    _density_scalar,
    _density_xy,
    _display_step_to_data,
    _evaluate_ncl_display_curve,
    _finite_difference_step,
    _fit_single_bend_display_curve,
    _gen_starting_points,
    _infer_profile_ncl_ref_magnitude,
    _local_display_jacobian,
    _map_ncl_display_points_to_viewport,
    _ncl_step_length_px,
    _NCLDisplaySampler,
    _normalize_ncl_preset,
    _open_arrow_geometry,
    _point_at_arc_distance_from_end,
    _point_within_grid_data,
    _prepare_ncl_display_sampler,
    _prepare_ncl_native_trace_context,
    _resolve_default_ncl_preset,
    _resolve_ncl_min_distance_fraction,
    _sample_grid_field,
    _sample_grid_field_array,
    _select_ncl_centers,
    _thin_ncl_mapped_candidates,
    _tip_display_geometry_from_display_curve,
    _trim_curve_for_open_head,
    _trim_display_curve_from_end,
    curly_vector,
    interpgrid,
)


class TestCurlyVector:
    """Test the curly_vector function."""

    @pytest.fixture
    def sample_vector_field(self):
        """Create sample vector field data for testing."""
        # Create coordinate arrays
        x = np.linspace(-2, 2, 10)
        y = np.linspace(-2, 2, 8)

        # Create 2D grid
        X, Y = np.meshgrid(x, y)

        # Create simple circular flow
        u = -Y * 0.5
        v = X * 0.5

        return x, y, u, v

    def test_curly_vector_basic(self, sample_vector_field):
        """Test basic curly_vector functionality."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(ax, x, y, u, v)

        # Check return type
        assert isinstance(result, CurlyVectorPlotSet)

        # Check that result has required attributes
        assert hasattr(result, "lines")
        assert hasattr(result, "arrows")
        assert hasattr(result, "resolution")
        assert hasattr(result, "magnitude")
        assert hasattr(result, "axes")

        # Check types of components
        assert isinstance(result.lines, LineCollection)
        assert isinstance(result.arrows, tuple)
        assert result.render_mode == "ncl_curly"

        plt.close(fig)

    def test_curly_vector_adds_line_collection_without_autolim(
        self, sample_vector_field, monkeypatch
    ):
        """The NCL-like renderer should skip expensive collection autolim work."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        autolim_calls = []
        original_add_collection = ax.add_collection

        def _record_add_collection(collection, autolim=True):
            autolim_calls.append(bool(autolim))
            return original_add_collection(collection, autolim=autolim)

        monkeypatch.setattr(ax, "add_collection", _record_add_collection)

        curly_vector(ax, x, y, u, v, color="k")

        assert autolim_calls
        assert autolim_calls[-1] is False
        assert ax.get_xlim()[0] <= np.min(x)
        assert ax.get_xlim()[1] >= np.max(x)
        assert ax.get_ylim()[0] <= np.min(y)
        assert ax.get_ylim()[1] >= np.max(y)

        plt.close(fig)

    def test_curly_vector_with_parameters(self, sample_vector_field):
        """Test curly_vector with various parameters."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            density=2,
            linewidth=2.0,
            color="red",
            arrowsize=1.5,
            arrowstyle="->",
            zorder=5,
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert result.linewidth == 2.0
        assert result.color == "red"
        assert result.arrowsize == 1.5
        assert result.arrowstyle == "->"
        assert result.zorder == 5

        plt.close(fig)

    def test_curly_vector_default_ncl_curly_basic(self, sample_vector_field):
        """Test the default NCL-like curved glyph rendering path."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            density=1.2,
            min_frac_length=0.5,
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert isinstance(result.lines, LineCollection)
        assert isinstance(result.arrows, tuple)
        assert result.render_mode == "ncl_curly"
        assert result.density == 1.2
        assert result.anchor == "center"
        assert len(ax.collections) > 0

        plt.close(fig)

    def test_curly_vector_anchor_mapping(self, sample_vector_field):
        """Test automatic anchor mapping in the NCL-like branch."""
        x, y, u, v = sample_vector_field
        expected = {
            "forward": "tail",
            "backward": "head",
            "both": "center",
        }

        for direction, anchor in expected.items():
            fig, ax = plt.subplots(figsize=(6, 4))
            result = curly_vector(
                ax,
                x,
                y,
                u,
                v,
                integration_direction=direction,
            )
            assert result.anchor == anchor
            plt.close(fig)

    def test_curly_vector_open_arrowstyle_uses_open_heads(self, sample_vector_field):
        """Test that ``arrowstyle='->'`` uses open line heads."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            arrowstyle="->",
            density=0.8,
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert len(result.lines.get_segments()) > 0
        assert len(result.arrows) == 0

        plt.close(fig)

    def test_curly_vector_filled_arrowstyle_creates_arrow_paths(
        self, sample_vector_field
    ):
        """Filled arrow styles should return the actual arrow-head patches."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            arrowstyle="-|>",
            density=0.8,
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert len(result.lines.get_segments()) > 0
        assert len(result.arrows) > 0
        assert all(patch in ax.patches for patch in result.arrows)

        plt.close(fig)

    def test_curly_vector_all_nan_field_returns_empty_artists(
        self, sample_vector_field
    ):
        """All-NaN vector fields should return empty artists without crashing."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        result = curly_vector(
            ax,
            x,
            y,
            np.full_like(u, np.nan),
            np.full_like(v, np.nan),
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert len(result.lines.get_segments()) == 0
        assert len(result.arrows) == 0

        plt.close(fig)

    def test_curly_vector_all_nan_color_field_raises_clear_error(
        self, sample_vector_field
    ):
        """All-NaN color fields should fail fast with a clear error."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(
            ValueError,
            match="color field must contain at least one finite value after masking",
        ):
            curly_vector(
                ax,
                x,
                y,
                u,
                v,
                color=np.full_like(u, np.nan),
                cmap="viridis",
            )

        plt.close(fig)

    def test_curly_vector_all_nan_linewidth_field_raises_clear_error(
        self, sample_vector_field
    ):
        """All-NaN linewidth fields should fail fast with a clear error."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(
            ValueError,
            match="linewidth field must contain at least one finite value after masking",
        ):
            curly_vector(
                ax,
                x,
                y,
                u,
                v,
                linewidth=np.full_like(u, np.nan),
            )

        plt.close(fig)

    def test_curly_vector_partial_nan_linewidth_field_falls_back_to_finite_mean(
        self, sample_vector_field
    ):
        """NaN gaps in a linewidth field should not leak NaN into the artists."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        linewidth = np.full_like(u, np.nan, dtype=float)
        linewidth[::2, ::2] = 2.5

        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            linewidth=linewidth,
            density=0.8,
        )

        assert len(result.lines.get_segments()) > 0
        assert not np.isnan(
            np.asarray(result.lines.get_linewidths(), dtype=float)
        ).any()

        plt.close(fig)

    def test_curly_vector_no_longer_accepts_render_mode_argument(
        self, sample_vector_field
    ):
        """The legacy render-mode switch should be removed from the public API."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(TypeError):
            curly_vector(ax, x, y, u, v, render_mode="ncl_curly")

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_profile_preset_applies_conservative_defaults(
        self, mock_ncl_curly, sample_vector_field
    ):
        """The profile preset should opt into non-uniform support and milder scaling."""
        x, y, u, v = sample_vector_field
        v = v.copy()
        v[0, 0] = 12.0
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(ax, x, y, u, v, ncl_preset="profile")

        call_kwargs = mock_ncl_curly.call_args.kwargs
        expected_ref = np.nanpercentile(np.hypot(u, v), 97.0)

        assert call_kwargs["allow_non_uniform_grid"] is True
        assert call_kwargs["ref_length"] == pytest.approx(0.06)
        assert call_kwargs["min_distance"] is None
        assert call_kwargs["ref_magnitude"] == pytest.approx(expected_ref)
        assert call_kwargs["ncl_preset"] == "profile"

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_profile_preset_preserves_explicit_overrides(
        self, mock_ncl_curly, sample_vector_field
    ):
        """Explicit NCL length controls should win over the profile preset defaults."""
        x, y, u, v = sample_vector_field
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(
            ax,
            x,
            y,
            u,
            v,
            ncl_preset="profile",
            allow_non_uniform_grid=True,
            ref_magnitude=8.0,
            ref_length=0.04,
            min_distance=0.03,
        )

        call_kwargs = mock_ncl_curly.call_args.kwargs
        assert call_kwargs["allow_non_uniform_grid"] is True
        assert call_kwargs["ref_magnitude"] == 8.0
        assert call_kwargs["ref_length"] == 0.04
        assert call_kwargs["min_distance"] == 0.03

        plt.close(fig)

    def test_curly_vector_rejects_unknown_ncl_preset(self, sample_vector_field):
        """Unknown presets should fail fast with a clear error."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(ValueError):
            curly_vector(ax, x, y, u, v, ncl_preset="bad-preset")

        plt.close(fig)

    def test_profile_ref_magnitude_uses_robust_percentile(self):
        """Profile reference magnitude should ignore a single strong outlier."""
        u = np.ones((3, 4))
        v = np.ones((3, 4))
        v[0, 0] = 20.0

        ref = _infer_profile_ncl_ref_magnitude(u, v)

        assert ref == pytest.approx(np.nanpercentile(np.hypot(u, v), 97.0))

    def test_apply_profile_preset_defaults_only_when_missing(self):
        """The profile preset should fill only the missing NCL-like controls."""
        u = np.ones((3, 4))
        v = np.ones((3, 4))

        resolved = _apply_ncl_preset_defaults(
            ncl_preset="profile",
            allow_non_uniform_grid=False,
            ref_magnitude=None,
            ref_length=None,
            min_distance=None,
            u=u,
            v=v,
        )

        assert resolved[0] is True
        assert resolved[1] == pytest.approx(np.nanpercentile(np.hypot(u, v), 97.0))
        assert resolved[2] == pytest.approx(0.06)
        assert resolved[3] is None
        assert resolved[4] == "profile"

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_auto_profile_preset_for_non_uniform_coordinates(
        self, mock_ncl_curly, sample_vector_field
    ):
        """Non-uniform profile-like axes should not require an explicit preset."""
        x, _, u, v = sample_vector_field
        y = np.array([100.0, 150.0, 225.0, 325.0, 450.0, 600.0, 800.0, 1000.0])
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(ax, x, y, u, v)

        call_kwargs = mock_ncl_curly.call_args.kwargs
        expected_ref = np.nanpercentile(np.hypot(u, v), 97.0)

        assert call_kwargs["allow_non_uniform_grid"] is True
        assert call_kwargs["ncl_preset"] == "profile"
        assert call_kwargs["ref_length"] == pytest.approx(0.06)
        assert call_kwargs["min_distance"] is None
        assert call_kwargs["ref_magnitude"] == pytest.approx(expected_ref)

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_allow_non_uniform_grid_auto_uses_profile_defaults(
        self, mock_ncl_curly, sample_vector_field
    ):
        """Explicit non-uniform mode should opt into profile defaults by default."""
        x, y, u, v = sample_vector_field
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(ax, x, y, u, v, allow_non_uniform_grid=True)

        call_kwargs = mock_ncl_curly.call_args.kwargs
        expected_ref = np.nanpercentile(np.hypot(u, v), 97.0)

        assert call_kwargs["allow_non_uniform_grid"] is True
        assert call_kwargs["ncl_preset"] == "profile"
        assert call_kwargs["ref_length"] == pytest.approx(0.06)
        assert call_kwargs["min_distance"] is None
        assert call_kwargs["ref_magnitude"] == pytest.approx(expected_ref)

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_non_uniform_2d_meshgrid_preserves_original_resolution(
        self, mock_ncl_curly
    ):
        """2D meshgrid profile input should keep its native (ny, nx) resolution."""
        x_1d = np.linspace(-60.0, 60.0, 5)
        y_1d = np.array([1000.0, 700.0, 400.0])
        x_2d, y_2d = np.meshgrid(x_1d, y_1d, indexing="xy")
        u = np.arange(15.0).reshape(3, 5)
        v = -u
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(ax, x_2d, y_2d, u, v, allow_non_uniform_grid=True)

        call_args = mock_ncl_curly.call_args.args
        assert np.asarray(call_args[1]).shape == (5,)
        assert np.asarray(call_args[2]).shape == (3,)
        assert np.asarray(call_args[3]).shape == (3, 5)
        assert np.asarray(call_args[4]).shape == (3, 5)

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_non_uniform_masked_values_stay_nan_after_regridding(
        self, mock_ncl_curly
    ):
        """Masked profile data should not be converted into synthetic vectors."""
        x = np.array([-30.0, -10.0, 10.0, 30.0])
        y = np.array([1000.0, 700.0, 400.0])
        base = np.arange(12.0).reshape(3, 4)
        u = np.ma.array(base, mask=np.zeros((3, 4), dtype=bool))
        v = np.ma.array(-base, mask=np.zeros((3, 4), dtype=bool))
        u.mask[1, 2] = True
        v.mask[1, 2] = True
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(ax, x, y, u, v, allow_non_uniform_grid=True)

        u_grid = np.asarray(mock_ncl_curly.call_args.args[3], dtype=float)
        v_grid = np.asarray(mock_ncl_curly.call_args.args[4], dtype=float)
        assert np.isnan(u_grid[1, 2])
        assert np.isnan(v_grid[1, 2])

        plt.close(fig)

    @patch("skyborn.plot.vector_plot._curly_vector_ncl")
    def test_curly_vector_non_uniform_style_fields_follow_vector_regridding(
        self, mock_ncl_curly
    ):
        """Non-uniform style fields must be reordered/regridded with the vectors."""
        x = np.array([3.0, 2.0, 1.0, 0.0])
        y = np.array([2.0, 1.0, 0.0])
        u = np.arange(12.0).reshape(3, 4)
        v = -u
        color = 100.0 + u
        linewidth = 200.0 + u
        mock_ncl_curly.return_value = Mock(spec=CurlyVectorPlotSet)

        fig, ax = plt.subplots(figsize=(6, 4))

        curly_vector(
            ax,
            x,
            y,
            u,
            v,
            color=color,
            linewidth=linewidth,
            allow_non_uniform_grid=True,
        )

        call_args = mock_ncl_curly.call_args
        np.testing.assert_allclose(
            np.asarray(call_args.args[1], dtype=float), [0, 1, 2, 3]
        )
        np.testing.assert_allclose(
            np.asarray(call_args.args[2], dtype=float), [0, 1, 2]
        )
        np.testing.assert_allclose(
            np.asarray(call_args.args[3], dtype=float),
            u[::-1, ::-1],
        )
        np.testing.assert_allclose(
            np.asarray(call_args.args[4], dtype=float),
            v[::-1, ::-1],
        )
        np.testing.assert_allclose(
            np.asarray(call_args.kwargs["color"], dtype=float),
            color[::-1, ::-1],
        )
        np.testing.assert_allclose(
            np.asarray(call_args.kwargs["linewidth"], dtype=float),
            linewidth[::-1, ::-1],
        )

        plt.close(fig)

    def test_curly_vector_rejects_unsupported_arrowstyle(self, sample_vector_field):
        """Unsupported arrow styles should fail fast instead of silently degrading."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(6, 4))

        with pytest.raises(ValueError, match="arrowstyle must be one of"):
            curly_vector(ax, x, y, u, v, arrowstyle="fancy")

        plt.close(fig)

    def test_ncl_step_length_scales_with_local_speed(self):
        """The NCL-like stepper should shorten steps rapidly in weaker flow."""
        assert _ncl_step_length_px(5.0, 10.0, 10.0) == pytest.approx(5.0)
        assert _ncl_step_length_px(5.0, 5.0, 10.0) == pytest.approx(1.25)
        assert _ncl_step_length_px(5.0, 0.1, 10.0) == pytest.approx(0.35)

    def test_trace_ncl_direction_prefers_native_when_available(self, monkeypatch):
        """The wrapper should use the native tracer when it returns a curve."""
        grid = Grid(np.linspace(0.0, 2.0, 3), np.linspace(0.0, 2.0, 3))
        u = np.ones(grid.shape, dtype=float)
        v = np.zeros(grid.shape, dtype=float)
        display_grid = np.stack(
            np.meshgrid(
                np.linspace(0.0, 20.0, 3), np.linspace(0.0, 20.0, 3), indexing="xy"
            ),
            axis=2,
        )
        display_sampler = Mock()
        display_sampler.display_grid = display_grid
        display_sampler.cell_valid = np.ones((2, 2), dtype=bool)

        expected = np.array([[0.5, 0.5], [1.0, 0.5], [1.5, 0.5]], dtype=float)

        monkeypatch.setattr(
            vector_plot_module,
            "_trace_ncl_direction_native",
            lambda **kwargs: expected.copy(),
        )

        def _should_not_fallback(**kwargs):
            raise AssertionError(
                "native tracer should have short-circuited the Python fallback"
            )

        monkeypatch.setattr(
            vector_plot_module, "_trace_ncl_direction_python", _should_not_fallback
        )

        curve = vector_plot_module._trace_ncl_direction(
            start_point=np.array([0.5, 0.5]),
            max_length_px=10.0,
            direction_sign=1.0,
            grid=grid,
            u=u,
            v=v,
            transform=Mock(),
            step_px=2.0,
            speed_scale=1.0,
            viewport=Bbox.from_extents(0.0, 0.0, 20.0, 20.0),
            display_sampler=display_sampler,
        )

        np.testing.assert_allclose(curve, expected)

    def test_trace_ncl_direction_falls_back_when_native_returns_none(self, monkeypatch):
        """The wrapper should preserve the Python path when the native call does not yield a curve."""
        grid = Grid(np.linspace(0.0, 2.0, 3), np.linspace(0.0, 2.0, 3))
        u = np.ones(grid.shape, dtype=float)
        v = np.zeros(grid.shape, dtype=float)
        display_grid = np.stack(
            np.meshgrid(
                np.linspace(0.0, 20.0, 3), np.linspace(0.0, 20.0, 3), indexing="xy"
            ),
            axis=2,
        )
        display_sampler = Mock()
        display_sampler.display_grid = display_grid
        display_sampler.cell_valid = np.ones((2, 2), dtype=bool)

        expected = np.array([[0.5, 0.5], [1.0, 0.5]], dtype=float)
        monkeypatch.setattr(
            vector_plot_module, "_trace_ncl_direction_native", lambda **kwargs: None
        )
        monkeypatch.setattr(
            vector_plot_module,
            "_trace_ncl_direction_python",
            lambda **kwargs: expected.copy(),
        )

        curve = vector_plot_module._trace_ncl_direction(
            start_point=np.array([0.5, 0.5]),
            max_length_px=10.0,
            direction_sign=1.0,
            grid=grid,
            u=u,
            v=v,
            transform=Mock(),
            step_px=2.0,
            speed_scale=1.0,
            viewport=Bbox.from_extents(0.0, 0.0, 20.0, 20.0),
            display_sampler=display_sampler,
        )

        np.testing.assert_allclose(curve, expected)

    def test_sample_grid_field_prefers_native_when_available(self, monkeypatch):
        """Scalar sampling should use the native helper when it returns a value."""
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        field = np.arange(9.0).reshape(3, 3)

        monkeypatch.setattr(
            vector_plot_module, "_sample_grid_field_native", lambda **kwargs: 4.25
        )

        value = vector_plot_module._sample_grid_field(grid, field, 0.5, 0.5)

        assert value == pytest.approx(4.25)

    def test_sample_grid_field_array_prefers_native_when_available(self, monkeypatch):
        """Vectorized scalar sampling should use the native helper when available."""
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        field = np.arange(9.0).reshape(3, 3)
        points = np.array([[0.5, 0.5], [1.5, 1.5]])
        expected = np.array([1.25, 7.75], dtype=float)

        monkeypatch.setattr(
            vector_plot_module,
            "_sample_grid_field_array_native",
            lambda **kwargs: expected.copy(),
        )

        sampled = vector_plot_module._sample_grid_field_array(grid, field, points)

        np.testing.assert_allclose(sampled, expected)

    def test_sample_grid_field_array_falls_back_when_native_returns_none(
        self, monkeypatch
    ):
        """Vectorized scalar sampling should preserve the Python fallback path."""
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        field = np.array(
            [
                [0.0, 1.0, 2.0],
                [2.0, 3.0, 4.0],
                [4.0, 5.0, 6.0],
            ]
        )
        points = np.array([[0.5, 0.5], [1.5, 1.5], [3.0, 3.0]])

        monkeypatch.setattr(
            vector_plot_module, "_sample_grid_field_array_native", lambda **kwargs: None
        )

        sampled = vector_plot_module._sample_grid_field_array(grid, field, points)

        assert sampled[0] == pytest.approx(
            vector_plot_module._sample_grid_field(grid, field, 0.5, 0.5)
        )
        assert sampled[1] == pytest.approx(
            vector_plot_module._sample_grid_field(grid, field, 1.5, 1.5)
        )
        assert np.isnan(sampled[2])

    def test_thin_ncl_mapped_candidates_prefers_native_when_available(
        self, monkeypatch
    ):
        """Mapped thinning should use the native helper when it returns indices."""
        mapped_points = np.array(
            [
                [0.10, 0.10],
                [0.12, 0.11],
                [0.55, 0.55],
            ]
        )

        monkeypatch.setattr(
            vector_plot_module,
            "_thin_ncl_mapped_candidates_native",
            lambda **kwargs: np.array([1], dtype=int),
        )

        selected = vector_plot_module._thin_ncl_mapped_candidates(
            mapped_points, spacing_frac=0.05
        )

        assert selected == [1]

    def test_thin_ncl_mapped_candidates_falls_back_when_native_returns_none(
        self, monkeypatch
    ):
        """Mapped thinning should preserve the Python fallback when native declines."""
        mapped_points = np.array(
            [
                [0.10, 0.10],
                [0.12, 0.11],
                [0.55, 0.55],
                [0.57, 0.54],
            ]
        )

        monkeypatch.setattr(
            vector_plot_module,
            "_thin_ncl_mapped_candidates_native",
            lambda **kwargs: None,
        )

        selected = vector_plot_module._thin_ncl_mapped_candidates(
            mapped_points, spacing_frac=0.05
        )

        assert selected == [0, 2]

    def test_map_ncl_display_points_to_viewport_normalizes_bbox(self):
        """Mapped thinning coordinates should be normalized to the active viewport."""
        fig, ax = plt.subplots(figsize=(6, 4))
        viewport = ax.bbox
        display_points = np.array(
            [
                [viewport.x0, viewport.y0],
                [
                    viewport.x0 + viewport.width / 2.0,
                    viewport.y0 + viewport.height / 4.0,
                ],
                [viewport.x1, viewport.y1],
            ]
        )

        mapped = _map_ncl_display_points_to_viewport(display_points, viewport)

        assert mapped[0] == pytest.approx([0.0, 0.0])
        assert mapped[1] == pytest.approx([0.5, 0.25])
        assert mapped[2] == pytest.approx([1.0, 1.0])
        plt.close(fig)

    def test_thin_ncl_mapped_candidates_culls_later_neighbors_in_scan_order(self):
        """The STTHIN-style pass should keep the first candidate and cull later nearby ones."""
        mapped_points = np.array(
            [
                [0.10, 0.10],
                [0.12, 0.11],
                [0.55, 0.55],
                [0.57, 0.54],
            ]
        )

        selected = _thin_ncl_mapped_candidates(mapped_points, spacing_frac=0.05)

        assert selected == [0, 2]

    def test_default_ncl_box_center_candidates_use_cell_centers(self):
        """Default NCL-style starts should come from grid-box centers, not grid nodes."""
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]))

        candidates = _default_ncl_box_center_candidates(grid)

        assert candidates.shape == (2, 2)
        np.testing.assert_allclose(candidates, [[0.5, 0.5], [1.5, 0.5]])

    def test_default_ncl_candidate_shape_scales_with_density(self):
        """Density should cap the default candidate lattice before full-grid sampling."""
        grid = Grid(np.linspace(0.0, 10.0, 101), np.linspace(0.0, 10.0, 101))

        assert _default_ncl_candidate_shape(grid, 0.5) == (15, 15)
        assert _default_ncl_candidate_shape(grid, 2.0) == (60, 60)
        assert _default_ncl_candidate_shape(grid, 5.0) == (100, 100)

    def test_default_ncl_box_center_candidates_density_reduces_large_grid_work(self):
        """Low density should produce fewer default candidates on large grids."""
        grid = Grid(np.linspace(0.0, 10.0, 101), np.linspace(0.0, 10.0, 101))

        sparse = _default_ncl_box_center_candidates(grid, density=0.5)
        dense = _default_ncl_box_center_candidates(grid, density=2.0)

        assert sparse.shape == (225, 2)
        assert dense.shape == (3600, 2)
        assert dense.shape[0] > sparse.shape[0]

    def test_profile_spacing_fraction_still_tracks_density(self):
        """Profile preset should remain density-aware instead of pinning one fixed spacing."""
        low_density = _resolve_ncl_min_distance_fraction(
            0.5, None, ncl_preset="profile"
        )
        high_density = _resolve_ncl_min_distance_fraction(
            2.0, None, ncl_preset="profile"
        )

        assert low_density == pytest.approx(0.036)
        assert high_density == pytest.approx(0.009)
        assert low_density > high_density

    def test_select_ncl_centers_keeps_scan_order_after_mapped_thinning(self):
        """Mapped thinning should preserve scan order when suppressing close start points."""
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0.0, 4.0)
        ax.set_ylim(0.0, 1.0)

        x = np.linspace(0.0, 4.0, 5)
        y = np.array([0.0, 1.0])
        grid = Grid(x, y)
        magnitude = np.ones(grid.shape)
        start_points = np.array(
            [
                [0.0, 0.5],
                [0.1, 0.5],
                [2.0, 0.5],
            ]
        )

        selected = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=ax.transData,
            axes=ax,
            density=1.0,
            start_points=start_points,
            min_distance=0.08,
        )

        assert len(selected) == 2
        np.testing.assert_allclose(selected[0][0], [0.0, 0.5])
        np.testing.assert_allclose(selected[1][0], [2.0, 0.5])
        plt.close(fig)

    def test_select_ncl_centers_defaults_to_box_centers(self):
        """Without user start points, center selection should operate on box-center candidates."""
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 1.0)

        x = np.array([0.0, 1.0, 2.0])
        y = np.array([0.0, 1.0])
        grid = Grid(x, y)
        magnitude = np.ones(grid.shape)

        selected = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=ax.transData,
            axes=ax,
            density=1.0,
            start_points=None,
            min_distance=0.0,
        )

        assert len(selected) == 2
        np.testing.assert_allclose(selected[0][0], [0.5, 0.5])
        np.testing.assert_allclose(selected[1][0], [1.5, 0.5])
        plt.close(fig)

    def test_select_ncl_centers_profile_density_changes_candidate_count(self):
        """Profile center selection should respond to density when min_distance is implicit."""
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0.0, 10.0)
        ax.set_ylim(0.0, 10.0)

        x = np.linspace(0.0, 10.0, 101)
        y = np.linspace(0.0, 10.0, 101)
        grid = Grid(x, y)
        magnitude = np.ones(grid.shape)

        sparse = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=ax.transData,
            axes=ax,
            density=0.5,
            start_points=None,
            min_distance=None,
            ncl_preset="profile",
        )
        dense = _select_ncl_centers(
            grid=grid,
            magnitude=magnitude,
            transform=ax.transData,
            axes=ax,
            density=2.0,
            start_points=None,
            min_distance=None,
            ncl_preset="profile",
        )

        assert len(sparse) < len(dense)
        plt.close(fig)

    def test_prepare_ncl_display_sampler_matches_identity_transform(self):
        """The cached display sampler should reproduce the local display mapping."""
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 1.0)
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]))

        sampler = _prepare_ncl_display_sampler(grid, ax.transData)
        display_point, jacobian = sampler.sample(np.array([0.5, 0.5]))

        expected_display = ax.transData.transform(np.asarray([[0.5, 0.5]]))[0]
        expected_jacobian = _local_display_jacobian(
            ax.transData, np.array([0.5, 0.5]), grid
        )

        np.testing.assert_allclose(display_point, expected_display)
        np.testing.assert_allclose(jacobian, expected_jacobian)
        plt.close(fig)

    def test_prepare_ncl_display_sampler_marks_projection_seam_cells_invalid(self):
        """Projection seam cells should not be sampled as contiguous display quads."""
        ccrs = pytest.importorskip("cartopy.crs")

        fig = plt.figure(figsize=(8, 4))
        ax = plt.axes(projection=ccrs.Robinson())
        fig.canvas.draw()

        grid = Grid(np.array([179.0, 180.0, 181.0]), np.array([-13.0, -12.0, -11.0]))
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        sampler = _prepare_ncl_display_sampler(grid, transform)

        assert sampler is not None
        assert sampler.cell_valid[:, 0].all()
        assert not np.any(sampler.cell_valid[:, 1])
        plt.close(fig)

    def test_candidate_data_from_display_step_uses_local_jacobian_inverse(self):
        """Display-space step inversion should stay consistent with the local transform."""
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 1.0)
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0]))
        current_data = np.array([0.6, 0.4])
        current_display = ax.transData.transform(np.asarray([current_data]))[0]
        jacobian = _local_display_jacobian(ax.transData, current_data, grid)
        candidate_display = current_display + np.array([14.0, 9.0])

        candidate = _candidate_data_from_display_step(
            current_data=current_data,
            current_display=current_display,
            candidate_display=candidate_display,
            jacobian=jacobian,
            transform=ax.transData,
        )
        expected = ax.transData.inverted().transform(np.asarray([candidate_display]))[0]

        np.testing.assert_allclose(candidate, expected)
        plt.close(fig)

    def test_sample_grid_field_array_matches_scalar_sampling(self):
        """Vectorized scalar-field sampling should agree with the scalar helper."""
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        field = np.array(
            [
                [0.0, 1.0, 2.0],
                [2.0, 3.0, 4.0],
                [4.0, 5.0, 6.0],
            ]
        )
        points = np.array([[0.5, 0.5], [1.5, 1.5], [3.0, 3.0]])

        sampled = _sample_grid_field_array(grid, field, points)

        assert sampled[0] == pytest.approx(_sample_grid_field(grid, field, 0.5, 0.5))
        assert sampled[1] == pytest.approx(_sample_grid_field(grid, field, 1.5, 1.5))
        assert np.isnan(sampled[2])

    def test_point_at_arc_distance_from_end_tracks_back_along_polyline(self):
        """Arrow orientation should be measured along the final arc, not only the last segment."""
        curve = np.array(
            [
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 4.0],
            ]
        )

        point = _point_at_arc_distance_from_end(curve, 2.0)

        assert point == pytest.approx([3.0, 2.0])

    def test_trim_display_curve_from_end_removes_head_overlap(self):
        """Open arrowheads should be able to trim the shaft before the tip."""
        curve = np.array(
            [
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 4.0],
            ]
        )

        trimmed = _trim_display_curve_from_end(curve, 1.5)

        assert trimmed[-1] == pytest.approx([3.0, 2.5])
        assert len(trimmed) == 3

    def test_curve_shape_accepts_smooth_arc(self):
        """Smooth single-bend glyphs should pass the quality gate."""
        fig, ax = plt.subplots(figsize=(6, 4))
        curve = np.array(
            [
                [0.0, 0.0],
                [0.4, 0.1],
                [0.8, 0.35],
                [1.2, 0.7],
            ]
        )
        assert _curve_shape_is_acceptable(curve, ax.transData)
        plt.close(fig)

    def test_curve_shape_rejects_wormy_s_curve(self):
        """Short S-shaped glyphs should be rejected by the quality gate."""
        fig, ax = plt.subplots(figsize=(6, 4))
        curve = np.array(
            [
                [0.0, 0.0],
                [0.3, 0.25],
                [0.6, -0.2],
                [0.9, 0.22],
                [1.2, -0.05],
            ]
        )
        assert not _curve_shape_is_acceptable(curve, ax.transData)
        plt.close(fig)

    def test_evaluate_ncl_display_curve_rejects_projection_seam_jump(self):
        """A dateline seam jump should be rejected instead of drawn as a straight line."""
        ccrs = pytest.importorskip("cartopy.crs")

        fig = plt.figure(figsize=(8, 4))
        ax = plt.axes(projection=ccrs.Robinson())
        fig.canvas.draw()

        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        curve = np.array([[179.0, -12.0], [180.0, -12.0], [181.0, -12.0]])

        display_curve, transform_failed = _evaluate_ncl_display_curve(
            curve,
            transform,
            viewport=ax.bbox,
        )

        assert display_curve is None
        assert not transform_failed
        plt.close(fig)

    @pytest.mark.parametrize("central_longitude", [179.0, 90.0])
    def test_curly_vector_geoaxes_bakes_segments_into_axes_data(
        self, central_longitude
    ):
        """GeoAxes rendering should store final line geometry in projection data space."""
        ccrs = pytest.importorskip("cartopy.crs")

        fig = plt.figure(figsize=(8, 4))
        ax = plt.axes(projection=ccrs.Robinson(central_longitude=central_longitude))
        ax.set_global()
        fig.canvas.draw()

        lon = np.linspace(0.0, 360.0, 36, endpoint=False)
        lat = np.linspace(-72.0, 72.0, 11)
        lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")
        u = np.cos(np.deg2rad(lat2d))
        v = 0.35 * np.sin(np.deg2rad(lon2d))

        result = curly_vector(
            ax,
            lon,
            lat,
            u,
            v,
            transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            density=0.55,
            color="k",
            arrowstyle="->",
        )

        assert result.lines.get_transform() is ax.transData
        assert result.transform is ax.transData

        max_jump = 0.0
        for segment in result.lines.get_segments():
            segment = np.asarray(segment, dtype=float)
            if segment.ndim != 2 or len(segment) < 2:
                continue
            display = ax.transData.transform(segment)
            jumps = np.hypot(np.diff(display[:, 0]), np.diff(display[:, 1]))
            if jumps.size:
                max_jump = max(max_jump, float(np.nanmax(jumps)))

        assert max_jump < float(np.hypot(ax.bbox.width, ax.bbox.height) * 0.35)
        plt.close(fig)

    def test_curly_vector_geoaxes_bakes_filled_arrow_patches_into_axes_data(self):
        """Filled arrowheads on GeoAxes should also use projection data coordinates."""
        ccrs = pytest.importorskip("cartopy.crs")

        fig = plt.figure(figsize=(8, 4))
        ax = plt.axes(projection=ccrs.Robinson(central_longitude=179.0))
        ax.set_global()
        fig.canvas.draw()

        lon = np.linspace(0.0, 360.0, 24, endpoint=False)
        lat = np.linspace(-60.0, 60.0, 9)
        lon2d, lat2d = np.meshgrid(lon, lat, indexing="xy")
        u = np.cos(np.deg2rad(lat2d))
        v = 0.25 * np.sin(np.deg2rad(lon2d))

        curly_vector(
            ax,
            lon,
            lat,
            u,
            v,
            transform=ccrs.PlateCarree()._as_mpl_transform(ax),
            density=0.45,
            color="k",
            arrowstyle="-|>",
        )

        assert ax.patches
        assert all(patch.get_transform() is ax.transData for patch in ax.patches)
        plt.close(fig)

    def test_fit_single_bend_display_curve_removes_turn_sign_changes(self):
        """The single-bend fitter should remove S-shaped inflection changes."""
        curve = np.array(
            [
                [0.0, 0.0],
                [0.3, 0.25],
                [0.6, -0.2],
                [0.9, 0.22],
                [1.2, -0.05],
            ]
        )
        fitted = _fit_single_bend_display_curve(curve)
        segments = np.diff(fitted, axis=0)
        lengths = np.hypot(segments[:, 0], segments[:, 1])
        directions = segments[lengths > 1e-6] / lengths[lengths > 1e-6, None]
        cross = (
            directions[:-1, 0] * directions[1:, 1]
            - directions[:-1, 1] * directions[1:, 0]
        )
        turn_sign = np.sign(cross[np.abs(cross) > 1e-6])
        assert np.count_nonzero(turn_sign[1:] * turn_sign[:-1] < 0) == 0

    def test_curly_vector_integration_directions(self, sample_vector_field):
        """Test different integration directions."""
        x, y, u, v = sample_vector_field
        directions = ["forward", "backward", "both"]

        for direction in directions:
            fig, ax = plt.subplots(figsize=(6, 4))

            result = curly_vector(ax, x, y, u, v, integration_direction=direction)

            assert isinstance(result, CurlyVectorPlotSet)
            assert result.integration_direction == direction

            plt.close(fig)

    def test_curly_vector_with_start_points(self, sample_vector_field):
        """Test curly_vector with custom start points."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Define custom start points
        start_points = np.array([[-1, -1], [0, 0], [1, 1]])

        result = curly_vector(ax, x, y, u, v, start_points=start_points)

        assert isinstance(result, CurlyVectorPlotSet)
        np.testing.assert_array_equal(result.start_points, start_points)

        plt.close(fig)

    def test_curly_vector_variable_linewidth(self, sample_vector_field):
        """Test curly_vector with variable linewidth."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create variable linewidth array
        linewidth = np.sqrt(u**2 + v**2)  # Linewidth proportional to speed

        result = curly_vector(ax, x, y, u, v, linewidth=linewidth)

        assert isinstance(result, CurlyVectorPlotSet)
        np.testing.assert_array_equal(result.linewidth, linewidth)

        plt.close(fig)

    def test_curly_vector_variable_color(self, sample_vector_field):
        """Test curly_vector with variable color."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create variable color array
        color = np.sqrt(u**2 + v**2)  # Color proportional to speed

        result = curly_vector(ax, x, y, u, v, color=color, cmap="viridis")

        assert isinstance(result, CurlyVectorPlotSet)
        np.testing.assert_array_equal(result.color, color)
        assert result.cmap.name == "viridis"

        plt.close(fig)

    def test_curly_vector_invalid_array_shapes(self, sample_vector_field):
        """Test curly_vector with invalid array shapes."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create mismatched u array
        u_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            curly_vector(ax, x, y, u_wrong, v)

        plt.close(fig)

    def test_curly_vector_invalid_color_shape(self, sample_vector_field):
        """Test curly_vector with invalid color array shape."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create wrong-shaped color array
        color_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            curly_vector(ax, x, y, u, v, color=color_wrong)

        plt.close(fig)

    def test_curly_vector_invalid_linewidth_shape(self, sample_vector_field):
        """Test curly_vector with invalid linewidth array shape."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create wrong-shaped linewidth array
        linewidth_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            curly_vector(ax, x, y, u, v, linewidth=linewidth_wrong)

        plt.close(fig)

    def test_curly_vector_outside_boundary_start_points(self, sample_vector_field):
        """Test curly_vector with start points outside data boundaries."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create start points outside the domain
        start_points = np.array(
            [[-10, -10], [10, 10]]  # Outside domain  # Outside domain
        )

        with pytest.raises(
            ValueError, match="Starting point .* outside of data boundaries"
        ):
            curly_vector(ax, x, y, u, v, start_points=start_points)

        plt.close(fig)


class TestInternalHelperFunctions:
    """Focused unit tests for helper branches in the low-level renderer."""

    def test_normalize_ncl_preset_accepts_aliases_and_rejects_unknown_values(self):
        assert _normalize_ncl_preset(None) is None
        assert _normalize_ncl_preset("profile") == "profile"
        assert _normalize_ncl_preset("vertical_profile") == "profile"
        assert _normalize_ncl_preset("vertical-profile") == "profile"
        assert _normalize_ncl_preset("lat_pressure") == "profile"

        with pytest.raises(ValueError, match="Unsupported ncl_preset"):
            _normalize_ncl_preset("bad")

    def test_axis_coordinate_1d_extracts_expected_axis_and_validates_name(self):
        x2d, y2d = np.meshgrid(np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0]))

        np.testing.assert_allclose(
            _axis_coordinate_1d(np.array([1.0, 2.0]), "x"), [1.0, 2.0]
        )
        np.testing.assert_allclose(_axis_coordinate_1d(x2d, "x"), [1.0, 2.0, 3.0])
        np.testing.assert_allclose(_axis_coordinate_1d(y2d, "y"), [4.0, 5.0])
        assert _axis_coordinate_1d(np.ones((2, 2, 2)), "x") is None

        with pytest.raises(ValueError, match="Unsupported axis_name"):
            _axis_coordinate_1d(x2d, "z")

    def test_axis_is_uniform_and_default_preset_resolution(self):
        assert bool(_axis_is_uniform(np.array([0.0, 1.0, 2.0, 3.0]))) is True
        assert bool(_axis_is_uniform(np.array([0.0, 1.0, 2.2, 3.0]))) is False
        assert bool(_axis_is_uniform(np.array([0.0, np.nan, 2.0]))) is False

        assert _resolve_default_ncl_preset(
            x=np.array([0.0, 1.0, 2.0]),
            y=np.array([1000.0, 900.0, 700.0]),
            allow_non_uniform_grid=False,
            ncl_preset=None,
        ) == (True, "profile")
        assert _resolve_default_ncl_preset(
            x=np.array([0.0, 1.0]),
            y=np.array([0.0, 1.0]),
            allow_non_uniform_grid=False,
            ncl_preset="profile",
        ) == (True, "profile")

    def test_density_helpers_clip_small_values(self):
        np.testing.assert_allclose(_density_xy((0.0, 0.05)), np.array([0.1, 0.1]))
        assert _density_scalar(0.0) == pytest.approx(0.1)
        assert _density_scalar((0.2, 0.4)) == pytest.approx(0.3)

    def test_finite_difference_step_chooses_valid_direction_or_span_fraction(self):
        assert _finite_difference_step(5.0, 0.0, 10.0, 2.0) == pytest.approx(2.0)
        assert _finite_difference_step(9.5, 0.0, 10.0, 1.0) == pytest.approx(-1.0)
        assert _finite_difference_step(0.4, 0.0, 1.0, 0.8) == pytest.approx(0.25)
        assert _finite_difference_step(0.0, 0.0, 0.0, 1.0) == pytest.approx(0.0)

    def test_display_step_to_data_handles_regular_and_singular_jacobians(self):
        jacobian = np.array([[2.0, 0.0], [0.0, 4.0]])
        display_step = np.array([4.0, 8.0])

        np.testing.assert_allclose(
            _display_step_to_data(jacobian, display_step),
            np.array([2.0, 2.0]),
        )
        assert (
            _display_step_to_data(np.array([[1.0, 2.0], [2.0, 4.0]]), display_step)
            is None
        )

    def test_point_within_grid_data_uses_data_space_bounds(self):
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([10.0, 20.0, 30.0]))

        assert bool(_point_within_grid_data(grid, np.array([1.5, 25.0]))) is True
        assert bool(_point_within_grid_data(grid, np.array([-0.1, 25.0]))) is False
        assert bool(_point_within_grid_data(grid, np.array([1.5, 35.0]))) is False

    def test_point_at_arc_distance_from_end_rejects_empty_curve(self):
        with pytest.raises(ValueError, match="at least one point"):
            _point_at_arc_distance_from_end(np.empty((0, 2)), 1.0)

    def test_tip_display_geometry_from_display_curve_rejects_degenerate_tip(self):
        display_curve = np.array([[1.0, 1.0], [1.0, 1.0]])

        assert _tip_display_geometry_from_display_curve(display_curve, 2.0) is None

    def test_open_arrow_geometry_rejects_nonfinite_display_curve(self):
        curve = np.array([[0.0, 0.0], [1.0, 0.0]])

        assert (
            _open_arrow_geometry(
                curve,
                transform=Mock(),
                head_length_px=5.0,
                display_curve=np.array([[0.0, 0.0], [np.nan, 1.0]]),
            )
            is None
        )

    def test_trim_curve_for_open_head_returns_original_when_inverse_transform_fails(
        self, monkeypatch
    ):
        curve = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        transform = Mock()
        transform.inverted.return_value.transform.side_effect = RuntimeError("boom")
        monkeypatch.setattr(
            vector_plot_module,
            "_open_arrow_geometry",
            lambda *args, **kwargs: {
                "display_curve": np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]),
                "base_center_display": np.array([1.5, 0.0]),
            },
        )

        trimmed = _trim_curve_for_open_head(curve, transform, 0.4)

        np.testing.assert_allclose(trimmed, curve)

    def test_ncl_display_sampler_invalid_points_and_empty_batches(self):
        grid = Grid(np.array([0.0, 1.0]), np.array([0.0, 1.0]))
        sampler = _NCLDisplaySampler(
            grid,
            np.array(
                [
                    [[0.0, 0.0], [1.0, 0.0]],
                    [[0.0, 1.0], [1.0, 1.0]],
                ]
            ),
        )

        assert sampler.sample(np.array([5.0, 5.0])) is None
        display_points, jacobians, valid = sampler.sample_points(
            np.empty((0, 2)),
            include_jacobian=True,
        )
        assert display_points.shape == (0, 2)
        assert jacobians.shape == (0, 2, 2)
        assert valid.shape == (0,)

    def test_prepare_ncl_display_sampler_handles_transform_failures(self):
        grid = Grid(np.array([0.0, 1.0]), np.array([0.0, 1.0]))

        class BadTransform:
            def transform(self, points):
                raise RuntimeError("boom")

        class NaNTransform:
            def transform(self, points):
                return np.full((len(points), 2), np.nan)

        assert _prepare_ncl_display_sampler(grid, BadTransform()) is None
        assert _prepare_ncl_display_sampler(grid, NaNTransform()) is None

    def test_prepare_ncl_native_trace_context_requires_native_backend_and_sampler(
        self, monkeypatch
    ):
        grid = Grid(np.array([0.0, 1.0]), np.array([0.0, 1.0]))
        viewport = Bbox.from_extents(0.0, 0.0, 10.0, 10.0)

        monkeypatch.setattr(vector_plot_module, "_trace_ncl_direction_native", None)
        assert (
            _prepare_ncl_native_trace_context(
                grid,
                np.ones((2, 2)),
                np.ones((2, 2)),
                viewport,
                display_sampler=Mock(),
            )
            is None
        )
        assert (
            _prepare_ncl_native_trace_context(
                grid,
                np.ones((2, 2)),
                np.ones((2, 2)),
                viewport,
                display_sampler=None,
            )
            is None
        )

    def test_sample_grid_field_array_returns_empty_array_for_empty_input(self):
        grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        field = np.arange(9.0).reshape(3, 3)

        sampled = _sample_grid_field_array(grid, field, np.empty((0, 2)))

        assert sampled.shape == (0,)


class TestCurlyVectorPlotSet:
    """Test the CurlyVectorPlotSet class."""

    @pytest.fixture
    def mock_collections(self):
        """Create mock collections for testing."""
        lines = Mock(spec=LineCollection)
        arrows = (Mock(name="arrow_patch_1"), Mock(name="arrow_patch_2"))
        return lines, arrows

    def test_curly_vector_plot_set_creation(self, mock_collections):
        """Test CurlyVectorPlotSet creation."""
        lines, arrows = mock_collections
        magnitude = np.array([[1, 2], [3, 4]])

        quiver_set = CurlyVectorPlotSet(
            lines=lines,
            arrows=arrows,
            resolution=0.5,
            magnitude=magnitude,
            zorder=1,
            transform=None,
            axes=None,
            linewidth=1.0,
            color="blue",
            cmap=None,
            arrowsize=1.0,
            arrowstyle="-|>",
            start_points=None,
            integration_direction="both",
            grains=15,
            broken_streamlines=True,
        )

        # Check attributes
        assert quiver_set.lines == lines
        assert quiver_set.arrows == arrows
        assert quiver_set.resolution == 0.5
        np.testing.assert_array_equal(quiver_set.magnitude, magnitude)
        assert quiver_set.linewidth == 1.0
        assert quiver_set.color == "blue"
        assert quiver_set.integration_direction == "both"

        # Check derived attributes
        assert quiver_set.max_magnitude == 4.0  # max of magnitude array
        assert quiver_set.scale_factor == 2.0  # for 'both' direction

    def test_curly_vector_plot_set_scale_factor_forward(self, mock_collections):
        """Test scale factor for forward integration."""
        lines, arrows = mock_collections

        quiver_set = CurlyVectorPlotSet(
            lines,
            arrows,
            0.5,
            np.array([[1, 2]]),
            1,
            None,
            None,
            1.0,
            "blue",
            None,
            1.0,
            "-|>",
            None,
            "forward",
            15,
            True,
        )

        assert quiver_set.scale_factor == 1.0  # for 'forward' direction

    def test_curly_vector_plot_set_scale_factor_backward(self, mock_collections):
        """Test scale factor for backward integration."""
        lines, arrows = mock_collections

        quiver_set = CurlyVectorPlotSet(
            lines,
            arrows,
            0.5,
            np.array([[1, 2]]),
            1,
            None,
            None,
            1.0,
            "blue",
            None,
            1.0,
            "-|>",
            None,
            "backward",
            15,
            True,
        )

        assert quiver_set.scale_factor == 1.0  # for 'backward' direction

    def test_curly_vector_plot_set_scaling_methods(self, mock_collections):
        """Test scaling methods."""
        lines, arrows = mock_collections

        quiver_set = CurlyVectorPlotSet(
            lines,
            arrows,
            0.5,
            np.array([[1, 2]]),
            1,
            None,
            None,
            1.0,
            "blue",
            None,
            1.0,
            "-|>",
            None,
            "both",
            15,
            True,
        )

        # Test get_scale_factor
        assert quiver_set.get_scale_factor() == 2.0

        # Test scale_value (physical to display)
        assert quiver_set.scale_value(10.0) == 5.0

        # Test unscale_value (display to physical)
        assert quiver_set.unscale_value(5.0) == 10.0


class TestGrid:
    """Test the Grid class."""

    def test_grid_1d_arrays(self):
        """Test Grid with 1D coordinate arrays."""
        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 5, 6)

        grid = Grid(x, y)

        assert grid.nx == 11
        assert grid.ny == 6
        assert grid.dx == 1.0
        assert grid.dy == 1.0
        assert grid.x_origin == 0.0
        assert grid.y_origin == 0.0
        assert grid.width == 10.0
        assert grid.height == 5.0
        assert grid.shape == (6, 11)

    def test_grid_2d_arrays(self):
        """Test Grid with 2D coordinate arrays."""
        x_1d = np.linspace(0, 4, 5)
        y_1d = np.linspace(0, 3, 4)
        x_2d, y_2d = np.meshgrid(x_1d, y_1d)

        grid = Grid(x_2d, y_2d)

        assert grid.nx == 5
        assert grid.ny == 4
        assert grid.dx == 1.0
        assert grid.dy == 1.0
        assert grid.shape == (4, 5)

    def test_grid_within_grid(self):
        """Test within_grid method."""
        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 5, 6)
        grid = Grid(x, y)

        # Test points within grid
        assert grid.within_grid(0, 0)
        assert grid.within_grid(5, 2.5)
        assert grid.within_grid(10, 5)

        # Test points outside grid
        assert not grid.within_grid(-1, 0)
        assert not grid.within_grid(11, 0)
        assert not grid.within_grid(0, -1)
        assert not grid.within_grid(0, 6)

    def test_grid_non_uniform_spacing_error(self):
        """Test Grid with non-uniform spacing raises error."""
        x = np.array([0, 1, 3, 4])  # Non-uniform spacing
        y = np.linspace(0, 3, 4)

        with pytest.raises(ValueError, match="values must be equally spaced"):
            Grid(x, y)

    def test_grid_non_increasing_error(self):
        """Test Grid with non-increasing values raises error."""
        x = np.array([0, 1, 0.5, 2])  # Not strictly increasing
        y = np.linspace(0, 3, 4)

        with pytest.raises(ValueError, match="must be strictly increasing"):
            Grid(x, y)

    def test_grid_inconsistent_2d_arrays_error(self):
        """Test Grid with inconsistent 2D arrays raises error."""
        # Create inconsistent 2D arrays
        x = np.array([[0, 1, 2], [0, 1, 3]])  # Different rows
        y = np.array([[0, 0, 0], [1, 1, 1]])

        with pytest.raises(ValueError, match="rows of 'x' must be equal"):
            Grid(x, y)

    def test_grid_requires_at_least_two_points_per_axis(self):
        """A grid with fewer than 2 points per axis should fail clearly."""
        with pytest.raises(ValueError, match="at least 2 points"):
            Grid(np.array([0.0]), np.array([0.0, 1.0]))

        with pytest.raises(ValueError, match="at least 2 points"):
            Grid(np.array([0.0, 1.0]), np.array([0.0]))


class TestStreamMask:
    """Test the StreamMask class."""

    def test_stream_mask_creation_scalar_density(self):
        """Test StreamMask creation with scalar density."""
        mask = StreamMask(1)

        assert mask.nx == 30
        assert mask.ny == 30
        assert mask.shape == (30, 30)
        assert np.all(mask._mask == 0)

    def test_stream_mask_creation_tuple_density(self):
        """Test StreamMask creation with tuple density."""
        mask = StreamMask((2, 1.5))

        assert mask.nx == 60
        assert mask.ny == 45
        assert mask.shape == (45, 60)

    def test_stream_mask_negative_density_error(self):
        """Test StreamMask with negative density raises error."""
        with pytest.raises(ValueError, match="must be positive"):
            StreamMask(-1)

    def test_stream_mask_invalid_density_error(self):
        """Test StreamMask with invalid density raises error."""
        with pytest.raises(ValueError, match="must be a scalar or be of length 2"):
            StreamMask([1, 2, 3])

    def test_stream_mask_trajectory_tracking(self):
        """Test trajectory tracking functionality."""
        mask = StreamMask(1)

        # Start trajectory
        mask._start_trajectory(5, 5)
        assert mask[5, 5] == 1
        assert (5, 5) in mask._traj

        # Update trajectory
        mask._update_trajectory(6, 5)
        assert mask[5, 6] == 1
        assert (5, 6) in mask._traj

        # Undo trajectory
        mask._undo_trajectory()
        assert mask[5, 5] == 0
        assert mask[5, 6] == 0

    def test_stream_mask_invalid_index_error(self):
        """Test InvalidIndexError when updating occupied cell."""
        mask = StreamMask(1)

        # Start trajectory
        mask._start_trajectory(5, 5)
        mask._update_trajectory(6, 5, broken_streamlines=True)

        # Try to re-enter an already occupied cell from a new position
        with pytest.raises(InvalidIndexError):
            mask._update_trajectory(5, 5, broken_streamlines=True)

    def test_stream_mask_broken_streamlines_false(self):
        """Test behavior when broken_streamlines=False."""
        mask = StreamMask(1)

        # Start trajectory
        mask._start_trajectory(5, 5)

        # Try to update to already occupied cell with broken_streamlines=False
        # Should not raise error
        mask._update_trajectory(5, 5, broken_streamlines=False)


class TestDomainMap:
    """Test the DomainMap class."""

    @pytest.fixture
    def sample_grid_and_mask(self):
        """Create sample grid and mask for testing."""
        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 5, 6)
        grid = Grid(x, y)
        mask = StreamMask(1)
        return grid, mask

    def test_domain_map_creation(self, sample_grid_and_mask):
        """Test DomainMap creation."""
        grid, mask = sample_grid_and_mask
        dmap = DomainMap(grid, mask)

        assert dmap.grid == grid
        assert dmap.mask == mask

        # Check conversion factors
        assert dmap.x_grid2mask == (mask.nx - 1) / (grid.nx - 1)
        assert dmap.y_grid2mask == (mask.ny - 1) / (grid.ny - 1)

    def test_domain_map_coordinate_conversions(self, sample_grid_and_mask):
        """Test coordinate conversion methods."""
        grid, mask = sample_grid_and_mask
        dmap = DomainMap(grid, mask)

        # Test grid2mask conversion
        xm, ym = dmap.grid2mask(5.0, 2.5)
        assert isinstance(xm, int)
        assert isinstance(ym, int)

        # Test mask2grid conversion (should be approximately inverse)
        xg, yg = dmap.mask2grid(xm, ym)
        assert abs(xg - 5.0) < 1.0  # Allow for rounding
        assert abs(yg - 2.5) < 1.0

        # Test data2grid conversion
        xg, yg = dmap.data2grid(5.0, 2.5)
        assert xg == 5.0  # Grid spacing is 1.0
        assert yg == 2.5

        # Test grid2data conversion (should be inverse)
        xd, yd = dmap.grid2data(xg, yg)
        assert xd == 5.0
        assert yd == 2.5

    def test_domain_map_trajectory_methods(self, sample_grid_and_mask):
        """Test trajectory-related methods."""
        grid, mask = sample_grid_and_mask
        dmap = DomainMap(grid, mask)

        # Test start_trajectory
        dmap.start_trajectory(5.0, 2.5)

        # Test reset_start_point
        dmap.reset_start_point(6.0, 3.0)

        # Test update_trajectory
        dmap.update_trajectory(6.5, 3.5)

        # Test undo_trajectory
        dmap.undo_trajectory()


class TestUtilityFunctions:
    """Test utility functions."""

    def test_interpgrid_basic(self):
        """Test basic interpgrid functionality."""
        # Create simple 2D array
        a = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)

        # Test interpolation at grid points
        assert interpgrid(a, 0, 0) == 1
        assert interpgrid(a, 1, 1) == 5
        assert interpgrid(a, 2, 2) == 9

        # Test interpolation between grid points
        result = interpgrid(a, 0.5, 0.5)
        expected = (1 + 2 + 4 + 5) / 4  # Bilinear interpolation
        assert abs(result - expected) < 1e-10

    def test_interpgrid_boundary_conditions(self):
        """Test interpgrid at boundaries."""
        a = np.array([[1, 2], [3, 4]], dtype=float)

        # Test at maximum indices
        assert interpgrid(a, 1, 1) == 4

        # Test beyond boundaries (should clip)
        assert interpgrid(a, 2, 2) == 4

    def test_interpgrid_array_input_clips_to_domain(self):
        """Array inputs beyond the grid should clip instead of indexing out."""
        a = np.array([[1, 2], [3, 4]], dtype=float)
        xi = np.array([0, 2, 3])
        yi = np.array([0, 2, 3])

        result = interpgrid(a, xi, yi)

        np.testing.assert_allclose(result, np.array([1.0, 4.0, 4.0]))

    def test_interpgrid_array_input(self):
        """Test interpgrid with array inputs."""
        a = np.array([[1, 2], [3, 4]], dtype=float)
        xi = np.array([0, 1])
        yi = np.array([0, 1])

        result = interpgrid(a, xi, yi)
        expected = np.array([1, 4])
        np.testing.assert_array_equal(result, expected)

    def test_interpgrid_masked_array(self):
        """Test interpgrid with masked arrays."""
        a = np.ma.array([[1, 2], [3, 4]], mask=[[False, True], [False, False]])

        # This should raise TerminateTrajectory for masked values
        with pytest.raises(TerminateTrajectory):
            interpgrid(a, 1, 0)  # Access masked location

    def test_gen_starting_points_integer_grains(self):
        """Test _gen_starting_points with integer grains."""
        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 5, 6)

        points = _gen_starting_points(x, y, 3)

        assert points.shape == (9, 2)  # 3x3 grid

        # Check that points are within bounds
        assert np.all(points[:, 0] >= x.min())
        assert np.all(points[:, 0] <= x.max())
        assert np.all(points[:, 1] >= y.min())
        assert np.all(points[:, 1] <= y.max())

    def test_gen_starting_points_tuple_grains(self):
        """Test _gen_starting_points with tuple grains."""
        x = np.linspace(0, 10, 11)
        y = np.linspace(0, 5, 6)

        points = _gen_starting_points(x, y, (4, 2))

        assert points.shape == (8, 2)  # 4x2 grid


class TestExceptionClasses:
    """Test custom exception classes."""

    def test_invalid_index_error(self):
        """Test InvalidIndexError can be raised and caught."""
        with pytest.raises(InvalidIndexError):
            raise InvalidIndexError("Test error")

    def test_terminate_trajectory(self):
        """Test TerminateTrajectory can be raised and caught."""
        with pytest.raises(TerminateTrajectory):
            raise TerminateTrajectory("Test termination")

    def test_out_of_bounds(self):
        """Test OutOfBounds can be raised and caught."""
        with pytest.raises(OutOfBounds):
            raise OutOfBounds("Test out of bounds")

    def test_exception_inheritance(self):
        """Test that exceptions have correct inheritance."""
        assert issubclass(InvalidIndexError, Exception)
        assert issubclass(TerminateTrajectory, Exception)
        assert issubclass(OutOfBounds, IndexError)


class TestIntegration:
    """Integration tests for curly-vector rendering."""

    @pytest.fixture
    def complex_vector_field(self):
        """Create a more complex vector field for integration testing."""
        # Create higher resolution grid
        x = np.linspace(-3, 3, 25)
        y = np.linspace(-2, 2, 20)
        X, Y = np.meshgrid(x, y)

        # Create complex flow pattern (double gyre)
        A = 0.1
        omega = 0.5
        epsilon = 0.1

        a = epsilon * np.sin(omega * 1.0)  # time = 1.0
        b = 1 - 2 * epsilon * np.sin(omega * 1.0)
        f = a * X**2 + b * X

        u = -np.pi * A * np.sin(np.pi * f) * np.cos(np.pi * Y)
        v = np.pi * A * np.cos(np.pi * f) * np.sin(np.pi * Y) * (2 * a * X + b)

        return x, y, u, v

    def test_full_curly_vector_integration(self, complex_vector_field):
        """Test full curly_vector functionality with complex field."""
        x, y, u, v = complex_vector_field
        fig, ax = plt.subplots(figsize=(10, 8))

        # Test with multiple parameters
        result = curly_vector(
            ax,
            x,
            y,
            u,
            v,
            density=1.5,
            linewidth=1.5,
            color=np.sqrt(u**2 + v**2),  # Color by speed
            cmap="plasma",
            arrowsize=1.2,
            arrowstyle="->",
            integration_direction="both",
            grains=20,
            broken_streamlines=True,
        )

        # Verify result
        assert isinstance(result, CurlyVectorPlotSet)
        assert result.density == 1.5
        assert result.integration_direction == "both"
        assert result.grains == 20
        assert result.broken_streamlines == True

        # Check that collections were added to axes
        assert len(ax.collections) > 0
        assert result.lines in ax.collections
        assert len(result.lines.get_segments()) > 0
        assert len(result.arrows) == 0
        assert len(ax.patches) == 0

        # Set plot properties
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title("Complex Vector Field")
        ax.grid(True, alpha=0.3)

        plt.close(fig)

    def test_error_handling_integration(self):
        """Test error handling in integration scenarios."""
        # Create simple field
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 4)
        u = np.ones((4, 5))
        v = np.zeros((4, 5))

        fig, ax = plt.subplots()

        # Test with various problematic inputs

        # Wrong shape arrays
        with pytest.raises(ValueError):
            curly_vector(ax, x, y, np.ones((3, 5)), v)

        # Invalid integration direction
        with pytest.raises(ValueError):
            curly_vector(ax, x, y, u, v, integration_direction="invalid")

        plt.close(fig)

    def test_masked_data_handling(self):
        """Test handling of masked data."""
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 4)

        # Create masked arrays
        u = np.ma.array(np.ones((4, 5)), mask=np.zeros((4, 5), dtype=bool))
        v = np.ma.array(np.zeros((4, 5)), mask=np.zeros((4, 5), dtype=bool))

        # Mask some points
        u.mask[1:3, 1:3] = True
        v.mask[1:3, 1:3] = True

        fig, ax = plt.subplots()

        result = curly_vector(ax, x, y, u, v, density=1)

        # Should handle masked data gracefully
        assert isinstance(result, CurlyVectorPlotSet)

        plt.close(fig)
