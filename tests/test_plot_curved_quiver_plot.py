"""Tests for the dataset-level curly-vector plotting modules."""

import warnings
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

import skyborn.plot.ncl_vector as ncl_vector_module
from skyborn.plot.ncl_vector import (
    _append_cyclic_column,
    _apply_dataset_isel,
    _build_curvilinear_target_grid,
    _default_cartopy_target_extent,
    _default_curvilinear_regrid_shape,
    _extract_curly_vector_dataset_source,
    _extract_regular_grid_from_regridded_vectors,
    _get_plot_dataarray,
    _grid_spans_full_longitude,
    _is_curvilinear_grid,
    _maybe_as_scalar_field,
    _normalize_density_pair,
    _normalize_regrid_shape,
    _prepare_source_vector_grid,
    _wrap_periodic_grid_queries,
    curly_vector,
)
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

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_array_input_forwards_quiver_style_kwargs(
        self, mock_curly_vector
    ):
        """The public wrapper should pass through supported quiver-like style kwargs."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.linspace(100.0, 110.0, 4)
        y = np.linspace(20.0, 26.0, 3)
        u = np.ones((3, 4)) * 3.0
        v = np.ones((3, 4)) * -1.5

        fig, ax = plt.subplots(figsize=(8, 6))
        result = curly_vector(
            x,
            y,
            u,
            v,
            ax=ax,
            alpha=0.4,
            facecolor="gold",
            edgecolor="firebrick",
            pivot="mid",
            rasterized=True,
        )

        call_args = mock_curly_vector.call_args
        assert call_args[0][0] is ax
        assert call_args[1]["alpha"] == pytest.approx(0.4)
        assert call_args[1]["facecolor"] == "gold"
        assert call_args[1]["edgecolor"] == "firebrick"
        assert call_args[1]["pivot"] == "mid"
        assert call_args[1]["rasterized"] is True
        assert result is mock_result

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_array_input_forwards_quiver_style_aliases(
        self, mock_curly_vector
    ):
        """The public wrapper should resolve Matplotlib-style aliases."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        x = np.linspace(100.0, 110.0, 4)
        y = np.linspace(20.0, 26.0, 3)
        u = np.ones((3, 4)) * 3.0
        v = np.ones((3, 4)) * -1.5
        color_field = np.arange(12, dtype=float).reshape(3, 4)

        fig, ax = plt.subplots(figsize=(8, 6))
        result = curly_vector(
            x,
            y,
            u,
            v,
            ax=ax,
            c=color_field,
            linewidths=1.75,
            facecolors="gold",
            edgecolors="navy",
            vmin=0.0,
            vmax=10.0,
        )

        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args[1]["color"], color_field)
        assert call_args[1]["linewidth"] == pytest.approx(1.75)
        assert call_args[1]["facecolor"] == "gold"
        assert call_args[1]["edgecolor"] == "navy"
        assert call_args[1]["vmin"] == pytest.approx(0.0)
        assert call_args[1]["vmax"] == pytest.approx(10.0)
        assert result is mock_result

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_dataset_input_resolves_quiver_style_aliases(
        self, mock_curly_vector, sample_data
    ):
        """Dataset wrapper should apply aliases before style-field preparation."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))
        color_field = sample_data["u"] ** 2 + sample_data["v"] ** 2
        linewidth_field = xr.full_like(sample_data["u"], 1.25)

        result = curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            c=color_field,
            linewidths=linewidth_field,
            facecolors="gold",
            edgecolors="navy",
            vmin=0.0,
            vmax=50.0,
        )

        call_args = mock_curly_vector.call_args
        assert call_args[0][0] is ax
        np.testing.assert_allclose(call_args[0][1], sample_data["x"].values)
        np.testing.assert_allclose(call_args[0][2], sample_data["y"].values)
        np.testing.assert_allclose(call_args[1]["color"], color_field.values)
        np.testing.assert_allclose(call_args[1]["linewidth"], linewidth_field.values)
        assert call_args[1]["facecolor"] == "gold"
        assert call_args[1]["edgecolor"] == "navy"
        assert call_args[1]["vmin"] == pytest.approx(0.0)
        assert call_args[1]["vmax"] == pytest.approx(50.0)
        assert result is mock_result

        plt.close(fig)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_wrapper_default_arrowstyle_matches_core(
        self, mock_curly_vector, sample_data
    ):
        """Dataset wrapper should use the same default arrowstyle as the core API."""
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        curly_vector(sample_data, x="x", y="y", u="u", v="v", ax=ax)

        assert mock_curly_vector.call_args.kwargs["arrowstyle"] == "->"

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

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_coerces_arraylike_color_and_linewidth(
        self, mock_curly_vector, sample_data
    ):
        """Array-like color/linewidth inputs should be normalized to ndarrays."""
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))
        color = sample_data["u"]
        linewidth = xr.DataArray(
            np.full(sample_data["u"].shape, 0.75),
            dims=sample_data["u"].dims,
            coords=sample_data["u"].coords,
        )

        curly_vector(
            sample_data,
            x="x",
            y="y",
            u="u",
            v="v",
            ax=ax,
            color=color,
            linewidth=linewidth,
        )

        call_args = mock_curly_vector.call_args
        assert isinstance(call_args[1]["color"], np.ndarray)
        assert isinstance(call_args[1]["linewidth"], np.ndarray)
        np.testing.assert_allclose(call_args[1]["color"], color.data)
        np.testing.assert_allclose(call_args[1]["linewidth"], linewidth.data)

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
            color=sample_data["u"],
            linewidth=xr.DataArray(
                np.ones_like(sample_data["u"].data),
                dims=sample_data["u"].dims,
                coords=sample_data["u"].coords,
            ),
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
        ncl_vector_module._ISSUED_PLOT_WARNINGS.clear()

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
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
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
        messages = [
            str(item.message)
            for item in caught
            if "Ignoring isel indexers" in str(item.message)
        ]
        assert len(messages) == 1

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


class TestDatasetCurlyVectorHelpers:
    """Unit tests for dataset/cartopy helper branches."""

    def test_apply_dataset_isel_warns_once_for_unused_indexers(self):
        da = xr.DataArray(
            np.arange(2 * 3 * 4).reshape(2, 3, 4),
            dims=("time", "level", "x"),
        )
        ncl_vector_module._ISSUED_PLOT_WARNINGS.clear()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            selected = _apply_dataset_isel(da, {"time": 1, "lat": 0})
            repeated = _apply_dataset_isel(da, {"time": 1, "lat": 0})

        assert selected.dims == ("level", "x")
        np.testing.assert_array_equal(selected.data, da.isel(time=1).data)
        np.testing.assert_array_equal(repeated.data, da.isel(time=1).data)
        messages = [
            str(item.message)
            for item in caught
            if "Ignoring isel indexers" in str(item.message)
        ]
        assert len(messages) == 1

    def test_get_plot_dataarray_rejects_scalar_after_selection(self):
        ds = xr.Dataset({"x": (["time"], np.array([1.0]))})

        with pytest.raises(
            ValueError, match="x must resolve to at least one dimension"
        ):
            _get_plot_dataarray(ds, "x", isel={"time": 0}, role="x")

    def test_get_plot_dataarray_rejects_more_than_2d_after_selection(self):
        data = np.arange(2 * 3 * 4 * 5, dtype=float).reshape(2, 3, 4, 5)
        ds = xr.Dataset({"u": (["time", "level", "y", "x"], data)})

        with pytest.raises(
            ValueError,
            match=r"u must resolve to 1D or 2D data after selection; remaining dims",
        ):
            _get_plot_dataarray(ds, "u", isel={"time": 0}, role="u")

    def test_extract_curly_vector_dataset_source_rejects_non_2d_vectors(self):
        ds = xr.Dataset(
            {
                "x": (["x"], np.linspace(0.0, 2.0, 3)),
                "y": (["y"], np.linspace(0.0, 1.0, 2)),
                "u": (["x"], np.linspace(1.0, 3.0, 3)),
                "v": (["y", "x"], np.ones((2, 3))),
            }
        )

        with pytest.raises(ValueError, match="u and v must each resolve to 2D arrays"):
            _extract_curly_vector_dataset_source(ds, "x", "y", "u", "v")

    def test_extract_curly_vector_dataset_source_rejects_mixed_coordinate_dims(self):
        x = xr.DataArray(np.linspace(0.0, 2.0, 3), dims=("x",))
        y = xr.DataArray(np.ones((2, 3)), dims=("y", "x"))
        ds = xr.Dataset(
            {
                "x": x,
                "y2d": y,
                "u": (["y", "x"], np.ones((2, 3))),
                "v": (["y", "x"], np.ones((2, 3))),
            }
        )

        with pytest.raises(ValueError, match="x and y must both be 1D or both be 2D"):
            _extract_curly_vector_dataset_source(ds, "x", "y2d", "u", "v")

    def test_default_cartopy_target_extent_uses_axes_when_available(self):
        ax = Mock()
        ax.projection = object()
        ax.get_extent.return_value = (-10.0, 10.0, -5.0, 5.0)

        assert _default_cartopy_target_extent(ax, None) == (-10.0, 10.0, -5.0, 5.0)
        assert _default_cartopy_target_extent(ax, (1.0, 2.0, 3.0, 4.0)) == (
            1.0,
            2.0,
            3.0,
            4.0,
        )

    def test_default_cartopy_target_extent_handles_axes_errors(self):
        ax = Mock()
        ax.projection = object()
        ax.get_extent.side_effect = RuntimeError("boom")

        assert _default_cartopy_target_extent(ax, None) is None
        assert _default_cartopy_target_extent(object(), None) is None

    def test_normalize_regrid_shape_and_density_pair_helpers(self):
        assert _normalize_regrid_shape(1) == (2, 2)
        assert _normalize_regrid_shape((7.9, 3.1)) == (7, 3)
        assert _normalize_density_pair(0.8) == (0.8, 0.8)
        assert _normalize_density_pair((0.6, 1.4)) == (0.6, 1.4)

        with pytest.raises(ValueError, match="density must be a scalar"):
            _normalize_density_pair((1.0, 2.0, 3.0))

    def test_is_curvilinear_grid_distinguishes_rectilinear_and_warped_meshes(self):
        x1d = np.linspace(100.0, 106.0, 4)
        y1d = np.linspace(20.0, 24.0, 3)
        x_rect, y_rect = np.meshgrid(x1d, y1d)

        assert bool(_is_curvilinear_grid(x_rect, y_rect)) is False

        x_warped = x_rect + np.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [0.2, 0.2, 0.2, 0.2],
                [0.4, 0.4, 0.4, 0.4],
            ]
        )
        assert bool(_is_curvilinear_grid(x_warped, y_rect)) is True

        with pytest.raises(ValueError, match="matching shapes"):
            _is_curvilinear_grid(np.ones((2, 3)), np.ones((3, 2)))

    def test_default_curvilinear_regrid_shape_respects_density_and_source_size(self):
        x = np.zeros((12, 20))

        assert _default_curvilinear_regrid_shape(x, density=0.4) == (20, 12)
        assert _default_curvilinear_regrid_shape(x, density=(2.0, 3.0)) == (20, 12)

    def test_build_curvilinear_target_grid_unwraps_longitude_and_rejects_all_nan(self):
        lon = np.array([[170.0, 175.0, -179.0], [170.0, 175.0, -179.0]])
        lat = np.array([[10.0, 10.0, 10.0], [20.0, 20.0, 20.0]])

        lon1d, lat1d = _build_curvilinear_target_grid(lon, lat, (4, 3))

        assert lon1d[0] == pytest.approx(170.0)
        assert lon1d[-1] == pytest.approx(181.0)
        assert lat1d[0] == pytest.approx(10.0)
        assert lat1d[-1] == pytest.approx(20.0)

        with pytest.raises(ValueError, match="contain no finite points"):
            _build_curvilinear_target_grid(
                np.full((2, 2), np.nan),
                np.full((2, 2), np.nan),
                (3, 3),
            )

    def test_maybe_as_scalar_field_accepts_only_matching_arrays(self):
        expected_shape = (2, 3)

        assert _maybe_as_scalar_field(None, expected_shape) is None
        assert _maybe_as_scalar_field("k", expected_shape) is None
        assert _maybe_as_scalar_field(2.0, expected_shape) is None
        assert _maybe_as_scalar_field(np.ones((3, 2)), expected_shape) is None
        np.testing.assert_allclose(
            _maybe_as_scalar_field([[1, 2, 3], [4, 5, 6]], expected_shape),
            np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        )

    def test_prepare_source_vector_grid_sorts_descending_coordinates(self):
        x = np.array([30.0, 20.0, 10.0])
        y = np.array([50.0, 40.0])
        u = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v = -u
        scalar = u + 100.0

        x_out, y_out, u_out, v_out, scalars_out = _prepare_source_vector_grid(
            x, y, u, v, scalars=(scalar,)
        )

        np.testing.assert_allclose(x_out, np.array([10.0, 20.0, 30.0]))
        np.testing.assert_allclose(y_out, np.array([40.0, 50.0]))
        np.testing.assert_allclose(u_out, np.array([[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]))
        np.testing.assert_allclose(v_out, -u_out)
        np.testing.assert_allclose(scalars_out[0], u_out + 100.0)

    def test_prepare_source_vector_grid_appends_cyclic_column_for_global_longitude(
        self,
    ):
        x = np.array([0.0, 90.0, 180.0, 270.0])
        y = np.array([-30.0, 30.0])
        u = np.arange(8.0).reshape(2, 4)
        v = -u
        scalar = np.full((2, 4), 7.0)

        x_out, y_out, u_out, v_out, scalars_out = _prepare_source_vector_grid(
            x, y, u, v, scalars=(scalar,)
        )

        np.testing.assert_allclose(x_out, np.array([0.0, 90.0, 180.0, 270.0, 360.0]))
        np.testing.assert_allclose(y_out, y)
        np.testing.assert_allclose(u_out[:, -1], u[:, 0])
        np.testing.assert_allclose(v_out[:, -1], v[:, 0])
        np.testing.assert_allclose(scalars_out[0][:, -1], scalar[:, 0])

    def test_prepare_source_vector_grid_rejects_non_rectilinear_axes(self):
        with pytest.raises(ValueError, match="requires 1D or meshgrid"):
            _prepare_source_vector_grid(
                np.ones((2, 2, 2)),
                np.ones((2, 2)),
                np.ones((2, 2)),
                np.ones((2, 2)),
            )

    def test_periodic_grid_helpers_wrap_queries_and_append_cyclic_data(self):
        x = np.array([0.0, 120.0, 240.0])
        field = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

        assert bool(_grid_spans_full_longitude(x)) is True
        x_cyclic, field_cyclic = _append_cyclic_column(x, field)
        np.testing.assert_allclose(x_cyclic, np.array([0.0, 120.0, 240.0, 360.0]))
        np.testing.assert_allclose(field_cyclic[:, -1], field[:, 0])
        np.testing.assert_allclose(
            _wrap_periodic_grid_queries(np.array([-10.0, 361.0]), x),
            np.array([-10.0, 1.0]),
        )
        np.testing.assert_allclose(
            _wrap_periodic_grid_queries(np.array([-10.0, 361.0]), np.array([0.0, 1.0])),
            np.array([-10.0, 361.0]),
        )

    def test_periodic_grid_helpers_wrap_queries_after_cyclic_append(self):
        x = np.array([0.0, 120.0, 240.0])
        field = np.array([[1.0, 2.0, 3.0]])

        x_cyclic, _ = _append_cyclic_column(x, field)

        np.testing.assert_allclose(
            _wrap_periodic_grid_queries(np.array([-10.0, 361.0]), x_cyclic),
            np.array([350.0, 1.0]),
        )

    def test_extract_curly_vector_dataset_source_transposes_2d_coords_to_match_uv_dims(
        self,
    ):
        x_values = np.array([[10.0, 20.0], [30.0, 40.0]])
        y_values = np.array([[50.0, 60.0], [70.0, 80.0]])
        u_values = np.array([[1.0, 2.0], [3.0, 4.0]])
        v_values = -u_values

        ds = xr.Dataset(
            {
                "u": (("y", "x"), u_values),
                "v": (("y", "x"), v_values),
                "lon": (("x", "y"), x_values),
                "lat": (("x", "y"), y_values),
            }
        )

        x_out, y_out, u_out, v_out = _extract_curly_vector_dataset_source(
            ds, "lon", "lat", "u", "v"
        )

        np.testing.assert_allclose(x_out, x_values.transpose())
        np.testing.assert_allclose(y_out, y_values.transpose())
        np.testing.assert_allclose(u_out, u_values)
        np.testing.assert_allclose(v_out, v_values)

    @patch("skyborn.plot.ncl_vector._array_curly_vector")
    def test_curly_vector_dataset_flips_external_style_fields_with_descending_axes(
        self, mock_curly_vector
    ):
        u_values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v_values = -u_values
        color = xr.DataArray(100.0 + u_values, dims=("y", "x"))
        linewidth = xr.DataArray((200.0 + u_values).transpose(), dims=("x", "y"))
        ds = xr.Dataset(
            {
                "u": (("y", "x"), u_values),
                "v": (("y", "x"), v_values),
            },
            coords={
                "x": ("x", np.array([2.0, 1.0, 0.0])),
                "y": ("y", np.array([10.0, 0.0])),
            },
        )
        mock_result = Mock(spec=CurlyVectorPlotSet)
        mock_curly_vector.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))

        curly_vector(
            ds, x="x", y="y", u="u", v="v", ax=ax, color=color, linewidth=linewidth
        )

        call_args = mock_curly_vector.call_args
        np.testing.assert_allclose(call_args.args[1], np.array([0.0, 1.0, 2.0]))
        np.testing.assert_allclose(call_args.args[2], np.array([0.0, 10.0]))
        expected = u_values[::-1, ::-1]
        np.testing.assert_allclose(
            np.asarray(call_args.kwargs["color"], dtype=float), 100.0 + expected
        )
        np.testing.assert_allclose(
            np.asarray(call_args.kwargs["linewidth"], dtype=float), 200.0 + expected
        )
        assert call_args.args[0] is ax
        assert mock_result is not None

        plt.close(fig)

    def test_extract_curly_vector_dataset_source_rejects_unmatched_2d_coord_dims(
        self,
    ):
        ds = xr.Dataset(
            {
                "u": (("y", "x"), np.ones((2, 3))),
                "v": (("y", "x"), np.ones((2, 3))),
                "lon": (("row", "col"), np.ones((2, 3))),
                "lat": (("row", "col"), np.ones((2, 3))),
            }
        )

        with pytest.raises(ValueError, match="dims"):
            _extract_curly_vector_dataset_source(ds, "lon", "lat", "u", "v")

    def test_extract_regular_grid_from_regridded_vectors_sorts_descending_axes(self):
        x_grid = np.array([[30.0, 20.0, 10.0], [30.0, 20.0, 10.0]])
        y_grid = np.array([[60.0, 60.0, 60.0], [50.0, 50.0, 50.0]])
        u_grid = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        v_grid = -u_grid
        scalar_grid = u_grid + 10.0

        x_1d, y_1d, u_out, v_out, scalar_grids = (
            _extract_regular_grid_from_regridded_vectors(
                x_grid,
                y_grid,
                u_grid,
                v_grid,
                [scalar_grid],
            )
        )

        np.testing.assert_allclose(x_1d, np.array([10.0, 20.0, 30.0]))
        np.testing.assert_allclose(y_1d, np.array([50.0, 60.0]))
        np.testing.assert_allclose(u_out, np.array([[6.0, 5.0, 4.0], [3.0, 2.0, 1.0]]))
        np.testing.assert_allclose(v_out, -u_out)
        np.testing.assert_allclose(scalar_grids[0], u_out + 10.0)

    def test_curly_vector_dispatch_rejects_missing_and_unsupported_arguments(self):
        with pytest.raises(TypeError, match="expects either"):
            curly_vector()

        with pytest.raises(TypeError, match="Unsupported arguments"):
            curly_vector(np.array([1.0]), np.array([2.0]), np.array([3.0]))


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

    def test_curly_vector_legend_open_arrow_shaft_reaches_tip(
        self, mock_curly_vector_set
    ):
        """Open key arrows should render as one connected glyph."""
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=8.0,
            units="m/s",
            labelpos="E",
        )
        fig.canvas.draw()

        arrow_tip = float(legend.head_left.get_xdata()[-1])
        shaft_tip = float(legend.arrow.get_xdata()[-1])
        assert shaft_tip == pytest.approx(arrow_tip)

        plt.close(fig)

    def test_curly_vector_legend_supports_explicit_xy_axes_position(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=8.0,
            x=0.8,
            y=1.05,
            labelpos="E",
        )
        fig.canvas.draw()

        assert legend._calculate_position() == pytest.approx((0.8, 1.05))
        assert legend.patch.get_x() == pytest.approx(0.8)
        assert legend.patch.get_y() == pytest.approx(1.05)

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


class TestCurlyVectorKeyInternals:
    """Focused unit tests for reference-key helper branches."""

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
        mock_set.glyph_head_axes_dimensions.return_value = (0.02, 0.03)
        return mock_set

    def test_curly_vector_key_rejects_invalid_labelpos(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(ValueError, match="labelpos must be one of"):
            CurlyVectorKey(
                ax=ax, curly_vector_set=mock_curly_vector_set, U=2.0, labelpos="Q"
            )

        plt.close(fig)

    def test_curly_vector_key_requires_x_and_y_together(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(ValueError, match="x and y must both be provided together"):
            CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=2.0, x=0.8)

        plt.close(fig)

    @pytest.mark.parametrize("bad_u", [np.nan, np.inf, -1.0, 0.0])
    def test_curly_vector_key_rejects_non_positive_or_non_finite_u(
        self, mock_curly_vector_set, bad_u
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(
            ValueError, match="U must be a finite positive reference magnitude"
        ):
            CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=bad_u)

        plt.close(fig)

    def test_curly_vector_key_rejects_unsupported_arrowstyle(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(ValueError, match="arrowstyle must be one of"):
            CurlyVectorKey(
                ax=ax,
                curly_vector_set=mock_curly_vector_set,
                U=4.0,
                arrow_props={"arrowstyle": "fancy"},
            )

        plt.close(fig)

    def test_curly_vector_key_center_label_hides_description(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=6.0,
            label="6 m/s",
            description="Reference Vector",
            center_label=True,
        )
        fig.canvas.draw()

        assert legend.text.get_text() == "6 m/s"
        assert legend.text2.get_visible() is False

        plt.close(fig)

    def test_calculate_head_geometry_falls_back_when_plotset_returns_invalid_values(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_curly_vector_set.glyph_head_axes_dimensions.return_value = (np.nan, np.nan)

        legend = CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=10.0)
        head_length, head_width = legend._calculate_head_geometry()

        assert 0.0 < head_length <= legend.arrow_length * 0.5
        assert head_width > 0.0

        plt.close(fig)

    def test_calculate_arrow_length_uses_ref_length_fraction_fallback(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_curly_vector_set.glyph_length_axes_fraction.side_effect = ValueError("bad")

        legend = CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=10.0)

        assert legend._calculate_arrow_length() == pytest.approx(0.06)

        plt.close(fig)

    def test_calculate_arrow_length_uses_reference_speed_fallback(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        mock_curly_vector_set.glyph_length_axes_fraction.side_effect = TypeError("bad")
        mock_curly_vector_set.max_magnitude = None
        mock_curly_vector_set.ref_length_fraction = None

        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=4.0,
            reference_speed=2.0,
            max_arrow_length=0.08,
        )

        assert legend._calculate_arrow_length() == pytest.approx(0.16)

        plt.close(fig)

    def test_calculate_arrow_length_handles_unexpected_exception(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=4.0,
            max_arrow_length=0.08,
        )
        legend.curly_vector_set = object()
        legend.reference_speed = object()

        assert legend._calculate_arrow_length() == pytest.approx(0.048)

        plt.close(fig)

    def test_calculate_position_supports_all_corner_locations(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        legend = CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=4.0)

        legend.loc = "lower left"
        assert legend._calculate_position() == pytest.approx(
            (legend.margin, legend.margin)
        )

        legend.loc = "upper left"
        assert legend._calculate_position() == pytest.approx(
            (legend.margin, 1.0 - legend.height - legend.margin)
        )

        legend.loc = "upper right"
        assert legend._calculate_position() == pytest.approx(
            (1.0 - legend.width - legend.margin, 1.0 - legend.height - legend.margin)
        )

        plt.close(fig)

    def test_setup_prop_helpers_merge_user_overrides(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))
        legend = CurlyVectorKey(
            ax=ax,
            curly_vector_set=mock_curly_vector_set,
            U=4.0,
            arrow_props={
                "color": "green",
                "linewidth": 4.0,
                "alpha": 0.3,
                "facecolor": "yellow",
                "edgecolor": "red",
            },
        )

        merged = legend._setup_props({"alpha": 0.5}, {"linewidth": 1.0, "alpha": 1.0})
        assert merged == {"linewidth": 1.0, "alpha": 0.5}
        assert legend._setup_line_props()["color"] == "green"
        assert legend._setup_line_props()["linewidth"] == 4.0
        assert legend._setup_head_fill_props()["facecolor"] == "yellow"
        assert legend._setup_head_fill_props()["edgecolor"] == "red"

        plt.close(fig)

    def test_set_figure_propagates_to_child_artists(self, mock_curly_vector_set):
        fig1, ax1 = plt.subplots(figsize=(8, 6))
        legend = CurlyVectorKey(ax=ax1, curly_vector_set=mock_curly_vector_set, U=4.0)
        legend.set_figure(fig1)

        assert legend.figure is fig1
        assert legend.patch.figure is fig1
        assert legend.arrow.figure is fig1
        assert legend.text.figure is fig1

        plt.close(fig1)

    def test_draw_returns_early_when_invisible(
        self, mock_curly_vector_set, monkeypatch
    ):
        fig, ax = plt.subplots(figsize=(8, 6))
        legend = CurlyVectorKey(ax=ax, curly_vector_set=mock_curly_vector_set, U=4.0)
        legend.set_visible(False)
        monkeypatch.setattr(
            legend,
            "_update_layout",
            lambda renderer: (_ for _ in ()).throw(AssertionError("should not run")),
        )

        legend.draw(fig.canvas.get_renderer())

        plt.close(fig)

    def test_curly_vector_key_rejects_missing_args_and_extra_positionals(
        self, mock_curly_vector_set
    ):
        fig, ax = plt.subplots(figsize=(8, 6))

        with pytest.raises(TypeError, match="expects either"):
            curly_vector_key()

        with pytest.raises(TypeError, match="missing required curly_vector_set"):
            curly_vector_key(ax)

        with pytest.raises(TypeError, match="too many positional arguments"):
            curly_vector_key(ax, mock_curly_vector_set, 2.0, "m/s", "extra")

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

    def test_curly_vector_key_accepts_explicit_xy_kwargs(self, mock_curly_vector_set):
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curly_vector_key(
            mock_curly_vector_set,
            U=4.0,
            ax=ax,
            x=0.8,
            y=1.05,
            labelpos="E",
        )
        fig.canvas.draw()

        assert result.x == pytest.approx(0.8)
        assert result.y == pytest.approx(1.05)
        assert result.patch.get_x() == pytest.approx(0.8)
        assert result.patch.get_y() == pytest.approx(1.05)

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
