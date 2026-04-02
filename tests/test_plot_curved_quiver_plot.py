"""Tests for the dataset-level curly-vector plotting modules."""

from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

from skyborn.plot.ncl_vector import curly_vector
from skyborn.plot.vector_key import CurlyVectorKey, curly_vector_key
from skyborn.plot.vector_plot import CurlyVectorPlotSet


class TestDatasetCurlyVector:
    """Test the curly_vector function."""

    @pytest.fixture
    def sample_data(self):
        """Create sample wind data for testing."""
        # Create 2D coordinate arrays
        x = np.linspace(-10, 10, 15)
        y = np.linspace(-5, 5, 10)

        # Create 2D grid
        X, Y = np.meshgrid(x, y)

        # Create sample wind field (circular pattern)
        u = -Y * 0.5  # u component
        v = X * 0.3  # v component

        # Create xarray Dataset
        ds = xr.Dataset(
            {
                "u": (["y", "x"], u),
                "v": (["y", "x"], v),
            },
            coords={"x": (["x"], x), "y": (["y"], y)},
        )

        return ds

    def test_curly_vector_basic(self, sample_data):
        """Test basic curly_vector functionality."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create dataset curly-vector plot
        result = curly_vector(sample_data, x="x", y="y", u="u", v="v", ax=ax, density=1)

        # Check return type
        assert isinstance(result, CurlyVectorPlotSet)

        # Check that result has required attributes
        assert hasattr(result, "lines")
        assert hasattr(result, "arrows")
        assert hasattr(result, "resolution")
        assert hasattr(result, "magnitude")

        plt.close(fig)

    def test_curly_vector_with_parameters(self, sample_data):
        """Test curly_vector with various parameters."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            density=2,
            linewidth=2.0,
            color="red",
            arrowsize=1.5,
            arrowstyle="->",
        )

        assert isinstance(result, CurlyVectorPlotSet)
        assert result.linewidth == 2.0
        assert result.color == "red"
        assert result.arrowsize == 1.5
        assert result.arrowstyle == "->"

        plt.close(fig)

    def test_curly_vector_no_axes(self, sample_data):
        """Test curly_vector without providing axes."""
        # This should use current axes
        plt.figure(figsize=(8, 6))

        result = curly_vector(sample_data, x="x", y="y", u="u", v="v", density=1)

        assert isinstance(result, CurlyVectorPlotSet)
        plt.close()

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_array_input_uses_current_axes(self, mock_curly_vector):
        """Array input should default to the current axes when ax is omitted."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.linspace(-10.0, 10.0, 6)
        y = np.linspace(-5.0, 5.0, 5)
        u = np.ones((5, 6))
        v = np.ones((5, 6)) * 2.0

        fig, ax = plt.subplots(figsize=(8, 6))
        result = curly_vector(x, y, u, v, density=0.8)

        call_args = mock_curly_vector.call_args
        assert call_args[0][0] is ax
        np.testing.assert_allclose(call_args[0][1], x)
        np.testing.assert_allclose(call_args[0][2], y)
        np.testing.assert_allclose(call_args[0][3], u)
        np.testing.assert_allclose(call_args[0][4], v)
        assert call_args[1]["density"] == 0.8
        assert result is mock_result

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_array_input_accepts_ax_keyword(self, mock_curly_vector):
        """Array input should also accept ax=... like pyplot.quiver()."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.linspace(100.0, 110.0, 4)
        y = np.linspace(20.0, 26.0, 3)
        u = np.ones((3, 4)) * 3.0
        v = np.ones((3, 4)) * -1.5

        fig, ax = plt.subplots(figsize=(8, 6))
        result = curly_vector(x, y, u, v, ax=ax, color="k")

        call_args = mock_curly_vector.call_args
        assert call_args[0][0] is ax
        assert call_args[1]["color"] == "k"
        assert result is mock_result

        plt.close(fig)

    def test_curly_vector_integration_directions(self, sample_data):
        """Test different integration directions."""
        directions = ["forward", "backward", "both"]

        for direction in directions:
            fig, ax = plt.subplots(figsize=(6, 4))

            result = curly_vector(
                sample_data,
                x="x",
                y="y",
                u="u",
                v="v",
                ax=ax,
                integration_direction=direction,
            )

            assert isinstance(result, CurlyVectorPlotSet)
            assert result.integration_direction == direction

            plt.close(fig)

    def test_curly_vector_with_start_points(self, sample_data):
        """Test curly_vector with custom start points."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Define custom start points
        start_points = np.array([[-5, -2], [0, 0], [5, 2]])

        result = curly_vector(
            sample_data, x="x", y="y", u="u", v="v", ax=ax, start_points=start_points
        )

        assert isinstance(result, CurlyVectorPlotSet)
        np.testing.assert_array_equal(result.start_points, start_points)

        plt.close(fig)

    def test_curly_vector_with_transform(self, sample_data):
        """Test curly_vector with transform parameter."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Skip this test due to complex mock transform requirements
        plt.close(fig)
        pytest.skip("Transform test requires complex matplotlib transform mocking")

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_converts_crs_like_transform(
        self, mock_curly_vector, sample_data
    ):
        """Test CRS-like objects are converted through _as_mpl_transform."""
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result
        converted_transform = Mock()

        class DummyCRS:
            def _as_mpl_transform(self, target_ax):
                assert target_ax is ax
                return converted_transform

        curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            transform=DummyCRS(),
        )

        assert mock_curly_vector.call_args[1]["transform"] is converted_transform
        plt.close(fig)

    def test_curly_vector_with_colormap(self, sample_data):
        """Test curly_vector with colormap settings."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create color data
        color_data = np.random.randn(10, 15)

        result = curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            color=color_data,
            cmap="viridis",
            norm=plt.Normalize(vmin=0, vmax=1),
        )

        assert isinstance(result, CurlyVectorPlotSet)
        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._regrid_cartopy_vectors")
    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_regrids_cartopy_vectors(
        self, mock_curly_vector, mock_regrid_vectors, sample_data
    ):
        """Test projection-aware vector regridding before calling curly_vector."""
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x_grid, y_grid = np.meshgrid(np.linspace(-2, 2, 4), np.linspace(-1, 1, 3))
        u_grid = np.full_like(x_grid, 5.0)
        v_grid = np.full_like(y_grid, 1.5)
        color_grid = np.full_like(x_grid, 7.0)
        linewidth_grid = np.full_like(x_grid, 0.8)
        mock_regrid_vectors.return_value = (
            x_grid,
            y_grid,
            u_grid,
            v_grid,
            color_grid,
            linewidth_grid,
        )

        class DummyCRS:
            def _as_mpl_transform(self, target_ax):
                return Mock()

        ax.projection = Mock()
        curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            transform=DummyCRS(),
            regrid_shape=(4, 3),
            color=sample_data["u"].data.copy(),
            linewidth=np.ones_like(sample_data["u"].data),
        )

        mock_regrid_vectors.assert_called_once()
        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args[0][1], x_grid[0, :])
        np.testing.assert_allclose(call_args[0][2], y_grid[:, 0])
        np.testing.assert_allclose(call_args[0][3], u_grid)
        np.testing.assert_allclose(call_args[0][4], v_grid)
        np.testing.assert_allclose(call_args[1]["color"], color_grid)
        np.testing.assert_allclose(call_args[1]["linewidth"], linewidth_grid)
        assert call_args[1]["transform"] is None

        plt.close(fig)

    def test_curly_vector_regrid_requires_cartopy_crs(self, sample_data):
        """Test that regrid_shape requires a CRS-like transform."""
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(ValueError, match="regrid_shape requires a Cartopy CRS"):
            curly_vector(
                sample_data,
                x="x",
                y="y",
                u="u",
                v="v",
                ax=ax,
                regrid_shape=20,
            )

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_supports_isel_for_multidimensional_inputs(
        self, mock_curly_vector
    ):
        """Test that curly_vector can select a 2D slice via isel."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.linspace(100.0, 110.0, 6)
        y = np.linspace(20.0, 28.0, 5)
        u = np.arange(2 * 3 * 5 * 6, dtype=float).reshape(2, 3, 5, 6)
        v = -u
        ds = xr.Dataset(
            {
                "u": (["time", "level", "y", "x"], u),
                "v": (["time", "level", "y", "x"], v),
            },
            coords={"x": x, "y": y},
        )

        fig, ax = plt.subplots(figsize=(8, 6))
        curly_vector(
            ds,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            isel={"time": 1, "level": 2},
        )

        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args[0][1], x)
        np.testing.assert_allclose(call_args[0][2], y)
        np.testing.assert_allclose(call_args[0][3], u[1, 2])
        np.testing.assert_allclose(call_args[0][4], v[1, 2])

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_sorts_descending_rectilinear_coordinates(
        self, mock_curly_vector
    ):
        """Descending 1D x/y coordinates should be normalized for plotting."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.array([110.0, 108.0, 106.0, 104.0])
        y = np.array([30.0, 25.0, 20.0])
        u = np.arange(12, dtype=float).reshape(3, 4)
        v = -u
        ds = xr.Dataset(
            {"u": (["y", "x"], u), "v": (["y", "x"], v)},
            coords={"x": x, "y": y},
        )

        fig, ax = plt.subplots(figsize=(8, 6))
        curly_vector(ds, x="x", y="y", u="u", v="v", ax=ax)

        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args[0][1], x[::-1])
        np.testing.assert_allclose(call_args[0][2], y[::-1])
        np.testing.assert_allclose(call_args[0][3], u[::-1, ::-1])
        np.testing.assert_allclose(call_args[0][4], v[::-1, ::-1])

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._regrid_curvilinear_vectors")
    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_regularizes_curvilinear_coordinates(
        self, mock_curly_vector, mock_regrid_curvilinear
    ):
        """Test that 2D curvilinear coordinates are regularized before plotting."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.array(
            [
                [100.0, 102.0, 104.0, 106.0],
                [100.4, 102.4, 104.4, 106.4],
                [100.8, 102.8, 104.8, 106.8],
            ]
        )
        y = np.array(
            [
                [20.0, 20.2, 20.4, 20.6],
                [22.0, 22.2, 22.4, 22.6],
                [24.0, 24.2, 24.4, 24.6],
            ]
        )
        u = np.ones_like(x)
        v = np.ones_like(x) * 2.0
        ds = xr.Dataset(
            {
                "x2d": (["y", "x"], x),
                "y2d": (["y", "x"], y),
                "u": (["y", "x"], u),
                "v": (["y", "x"], v),
            }
        )

        lon1d = np.linspace(100.0, 107.0, 5)
        lat1d = np.linspace(20.0, 24.5, 4)
        u_reg = np.full((4, 5), 3.0)
        v_reg = np.full((4, 5), 1.0)
        mock_regrid_curvilinear.return_value = (lon1d, lat1d, u_reg, v_reg)

        fig, ax = plt.subplots(figsize=(8, 6))
        curly_vector(
            ds,
            x="x2d",
            y="y2d",
            u="u",
            v="v",
            ax=ax,
            curvilinear_regrid_shape=(5, 4),
        )

        mock_regrid_curvilinear.assert_called_once()
        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args[0][1], lon1d)
        np.testing.assert_allclose(call_args[0][2], lat1d)
        np.testing.assert_allclose(call_args[0][3], u_reg)
        np.testing.assert_allclose(call_args[0][4], v_reg)

        plt.close(fig)

    def test_curly_vector_rejects_unaligned_vector_components(self):
        """Vector components must already live on the same physical grid."""
        x = np.linspace(100.0, 110.0, 6)
        y = np.linspace(20.0, 28.0, 5)
        ds = xr.Dataset(
            {
                "u": (["y", "x_u"], np.ones((5, 7))),
                "v": (["y_v", "x"], np.ones((6, 6))),
            },
            coords={"x": x, "y": y},
        )

        fig, ax = plt.subplots(figsize=(8, 6))
        with pytest.raises(ValueError, match="same physical grid"):
            curly_vector(ds, x="x", y="y", u="u", v="v", ax=ax)
        plt.close(fig)

    def test_curly_vector_with_zorder(self, sample_data):
        """Test curly_vector with zorder parameter."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(sample_data, x="x", y="y", u="u", v="v", ax=ax, zorder=10)

        assert isinstance(result, CurlyVectorPlotSet)
        plt.close(fig)

    def test_curly_vector_grains_parameter(self, sample_data):
        """Test curly_vector with different grains settings."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            grains=20,
            broken_streamlines=False,
        )

        assert isinstance(result, CurlyVectorPlotSet)
        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_calls_curly_vector(self, mock_curly_vector, sample_data):
        """Test that curly_vector properly calls curly_vector."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector(
            sample_data, x="x", y="y", u="u", v="v", ax=ax, density=2, linewidth=1.5
        )

        # Check that curly_vector was called
        mock_curly_vector.assert_called_once()

        # Check some of the arguments passed to curly_vector
        call_args = mock_curly_vector.call_args
        assert call_args[0][0] == ax  # First argument should be axes
        assert call_args[1]["density"] == 2
        assert call_args[1]["linewidth"] == 1.5
        assert "render_mode" not in call_args[1]

        assert result == mock_result
        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_forwards_ncl_options(self, mock_curly_vector, sample_data):
        """Test forwarding of the NCL-like rendering options."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))

        curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            anchor="tail",
            ref_magnitude=12.0,
            ref_length=0.05,
            min_frac_length=0.6,
            min_distance=0.02,
            ncl_preset="profile",
        )

        call_args = mock_curly_vector.call_args
        assert "render_mode" not in call_args[1]
        assert call_args[1]["anchor"] == "tail"
        assert call_args[1]["ref_magnitude"] == 12.0
        assert call_args[1]["ref_length"] == 0.05
        assert call_args[1]["min_frac_length"] == 0.6
        assert call_args[1]["min_distance"] == 0.02
        assert call_args[1]["ncl_preset"] == "profile"

        plt.close(fig)

    def test_curly_vector_no_longer_accepts_render_mode_argument(self, sample_data):
        """The wrapper should expose only the NCL-like path."""
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(TypeError):
            curly_vector(
                sample_data,
                x="x",
                y="y",
                u="u",
                v="v",
                ax=ax,
                render_mode="ncl_curly",
            )

        plt.close(fig)


class TestCurlyVectorKey:
    """Focused tests for the NCL-like reference annotation."""

    @pytest.fixture
    def mock_curly_vector_set(self):
        mock_set = Mock(spec=CurlyVectorPlotSet)
        mock_set.resolution = 0.5
        mock_set.magnitude = np.array([[0.5, 1.0], [1.5, 2.0]])
        mock_set.max_magnitude = 20.0
        mock_set.ref_length_fraction = 0.12
        mock_set.arrowstyle = "->"
        mock_set.arrowsize = 1.0
        mock_set.glyph_length_axes_fraction.side_effect = (
            lambda value: float(value) / 100.0
        )
        mock_set.glyph_head_axes_dimensions.return_value = (0.018, 0.026)
        return mock_set

    def test_curly_vector_legend_default_ncl_layout(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax, curly_vector_set=mock_curly_vector_set, U=10.0, units="m/s"
        )
        fig.canvas.draw()

        assert legend in ax.artists
        assert legend.labelpos == "N"
        assert legend.text.get_text() == "10 m/s"
        assert legend.text2.get_text() == "Reference Vector"
        assert legend.patch.get_width() > 0.0
        assert legend.patch.get_height() > 0.0
        assert len(legend.arrow.get_xdata()) == 2
        assert legend.head_left.get_visible() is True
        assert legend.head_right.get_visible() is True

        plt.close(fig)

    def test_curly_vector_legend_uses_plot_glyph_length(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax, curly_vector_set=mock_curly_vector_set, U=12.0, units="m/s"
        )

        assert legend._calculate_arrow_length() == pytest.approx(0.12)

        plt.close(fig)

    def test_curly_vector_legend_supports_horizontal_layout(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=8.0,
            units="m/s",
            labelpos="E",
        )
        fig.canvas.draw()

        assert legend.text.get_text() == "8 m/s\nReference Vector"
        assert legend.text2.get_visible() is False
        assert legend.patch.get_width() > legend.patch.get_height()

        plt.close(fig)

    def test_curly_vector_legend_can_hide_frame(self, mock_curly_vector_set):
        """frameon=False should suppress the annotation box."""
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=8.0,
            units="m/s",
            frameon=False,
        )
        fig.canvas.draw()

        assert legend.patch.get_visible() is False

        plt.close(fig)

    def test_curly_vector_legend_can_hide_description(self, mock_curly_vector_set):
        """show_description=False should suppress the default Reference Vector text."""
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=8.0,
            units="m/s",
            show_description=False,
        )
        fig.canvas.draw()

        assert legend.text.get_text() == "8 m/s"
        assert legend.text2.get_visible() is False
        assert legend.text2.get_text() == ""

        plt.close(fig)

    def test_curly_vector_legend_explicit_label_and_description(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=5.0,
            units="m/s",
            label="20 m/s",
            description="Jet Core",
        )
        fig.canvas.draw()

        assert legend.text.get_text() == "20 m/s"
        assert legend.text2.get_text() == "Jet Core"

        plt.close(fig)

    def test_curly_vector_legend_invalid_location(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(ValueError):
            CurlyVectorKey(
                ax=ax,
                curly_vector_set=mock_curly_vector_set,
                U=2.0,
                units="m/s",
                loc="invalid_location",
            )

        plt.close(fig)


class TestCurlyVectorKeyHelper:
    """Test the curly_vector_key convenience function."""

    @pytest.fixture
    def mock_curly_vector_set(self):
        mock_set = Mock(spec=CurlyVectorPlotSet)
        mock_set.resolution = 0.5
        mock_set.magnitude = np.array([[0.5, 1.0], [1.5, 2.0]])
        mock_set.max_magnitude = 20.0
        mock_set.ref_length_fraction = 0.1
        mock_set.arrowstyle = "->"
        mock_set.arrowsize = 1.0
        mock_set.glyph_length_axes_fraction.return_value = 0.1
        mock_set.glyph_head_axes_dimensions.return_value = (0.016, 0.024)
        return mock_set

    def test_curly_vector_key_basic(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(
            ax=ax, curly_vector_set=mock_curly_vector_set, U=2.0, units="m/s"
        )

        assert isinstance(result, CurlyVectorKey)
        assert result.U == 2.0
        assert result.units == "m/s"
        assert result.loc == "lower right"
        assert result.labelpos == "N"

        plt.close(fig)

    def test_curly_vector_key_uses_current_axes_when_ax_is_omitted(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(mock_curly_vector_set, U=2.0, units="m/s")

        assert isinstance(result, CurlyVectorKey)
        assert result.ax is ax

        plt.close(fig)

    def test_curly_vector_key_accepts_positional_plotset_and_ax_keyword(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(mock_curly_vector_set, 3.0, ax=ax, labelpos="E")

        assert isinstance(result, CurlyVectorKey)
        assert result.ax is ax
        assert result.U == 3.0
        assert result.labelpos == "E"

        plt.close(fig)

    def test_curly_vector_key_kwargs_passed(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=3.0,
            arrow_props={"color": "green", "linewidth": 3},
            margin=0.05,
            description="Reference Vector",
        )

        assert isinstance(result, CurlyVectorKey)
        assert result.arrow_props["color"] == "green"
        assert result.arrow_props["linewidth"] == 3
        assert result.margin == 0.05

        plt.close(fig)

    def test_curly_vector_key_frame_and_description_options(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(
            mock_curly_vector_set,
            U=4.0,
            ax=ax,
            frameon=False,
            show_description=False,
        )
        fig.canvas.draw()

        assert result.patch.get_visible() is False
        assert result.text.get_text() == "4 m/s"
        assert result.text2.get_visible() is False

        plt.close(fig)


class TestDatasetCurlyVectorIntegration:
    """Integration tests for dataset curly-vector plotting."""

    @pytest.fixture
    def wind_data(self):
        """Create realistic wind data for integration testing."""
        # Create coordinate arrays
        lon = np.linspace(-10, 10, 20)
        lat = np.linspace(-5, 5, 15)

        # Create 2D grid
        LON, LAT = np.meshgrid(lon, lat)

        # Create realistic wind pattern (westerly jet with some curvature)
        u = 10 * (1 - (LAT / 5) ** 2) * np.exp(-0.1 * LON**2)  # Westerly component
        v = 2 * np.sin(LON * np.pi / 10) * np.exp(-0.1 * LAT**2)  # Meridional component

        # Add some noise
        np.random.seed(42)
        u += np.random.normal(0, 0.5, u.shape)
        v += np.random.normal(0, 0.3, v.shape)

        # Create xarray Dataset
        ds = xr.Dataset(
            {
                "u": (["lat", "lon"], u),
                "v": (["lat", "lon"], v),
            },
            coords={"lon": (["lon"], lon), "lat": (["lat"], lat)},
        )

        return ds

    def test_curly_vector_with_legend_integration(self, wind_data):
        """Test full integration of dataset curly vectors with the key."""
        fig, ax = plt.subplots(figsize=(10, 8))

        # Create dataset curly-vector plot
        quiver_set = curly_vector(
            wind_data,
            x="lon",
            y="lat",
            u="u",
            v="v",
            ax=ax,
            density=1.5,
            color="blue",
            linewidth=1.5,
            arrowsize=1.2,
        )

        # Add legend
        legend = curly_vector_key(
            ax=ax,
            curly_vector_set=quiver_set,
            U=10.0,
            units="m/s",
            loc="upper right",
        )

        # Verify everything was created properly
        assert isinstance(quiver_set, CurlyVectorPlotSet)
        assert isinstance(legend, CurlyVectorKey)

        assert len(ax.collections) > 0  # Should have LineCollection
        assert legend in ax.artists

        # Set some basic plot properties
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Dataset Curly Vector Plot with Key")
        ax.grid(True, alpha=0.3)

        plt.close(fig)

    def test_multiple_legends_different_positions(self, wind_data):
        """Test multiple legends in different positions."""
        fig, ax = plt.subplots(figsize=(12, 8))

        # Create dataset curly-vector plot
        quiver_set = curly_vector(
            wind_data, x="lon", y="lat", u="u", v="v", ax=ax, density=1
        )

        # Add multiple legends
        positions = ["lower left", "lower right", "upper left"]
        speeds = [5.0, 10.0, 15.0]

        legends = []
        for pos, speed in zip(positions, speeds):
            legend = curly_vector_key(
                ax=ax,
                curly_vector_set=quiver_set,
                U=speed,
                units="m/s",
                loc=pos,
                width=0.15,
                height=0.10,
            )
            legends.append(legend)

        # All legends should be created
        assert len(legends) == 3
        for legend in legends:
            assert isinstance(legend, CurlyVectorKey)

        plt.close(fig)


class TestDatasetCurlyVectorErrorHandling:
    """Test error handling in dataset curly-vector plotting."""

    def test_curly_vector_with_invalid_data(self):
        """Test curly_vector with invalid data shapes."""
        u_data = np.random.randn(10, 15)

        x = np.linspace(-10, 10, 15)
        y = np.linspace(-5, 5, 10)

        # Test with missing variables should raise error
        ds = xr.Dataset(
            {"u": (["y", "x"], u_data)}, coords={"x": x, "y": y}  # Missing v variable
        )

        fig, ax = plt.subplots()

        # Should handle gracefully or raise appropriate error
        with pytest.raises((ValueError, KeyError)):
            curly_vector(ds, x="x", y="y", u="u", v="v", ax=ax)

        plt.close(fig)

    def test_legend_with_invalid_location(self):
        """Test legend creation with invalid location."""
        fig, ax = plt.subplots()
        mock_set = Mock(spec=CurlyVectorPlotSet)
        mock_set.resolution = 0.5
        mock_set.magnitude = np.array([[1.0, 2.0]])

        # Should raise error for invalid location
        with pytest.raises(ValueError):
            CurlyVectorKey(
                ax=ax,
                curly_vector_set=mock_set,
                U=2.0,
                units="m/s",
                loc="invalid_location",
            )

        plt.close(fig)
