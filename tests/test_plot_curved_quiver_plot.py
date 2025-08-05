"""
Tests for skyborn.plot.curved_quiver_plot module.

This module tests the curved quiver plotting functionality,
including the curved_quiver function and CurvedQuiverLegend class.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib
from unittest.mock import Mock, patch, MagicMock
from matplotlib.patches import FancyArrowPatch
from matplotlib.text import Text

from skyborn.plot.curved_quiver_plot import (
    curved_quiver,
    add_curved_quiverkey,
    CurvedQuiverLegend,
)
from skyborn.plot.modplot import CurvedQuiverplotSet


class TestCurvedQuiver:
    """Test the curved_quiver function."""

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

    def test_curved_quiver_basic(self, sample_data):
        """Test basic curved_quiver functionality."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create curved quiver plot
        result = curved_quiver(
            sample_data, x="x", y="y", u="u", v="v", ax=ax, density=1
        )

        # Check return type
        assert isinstance(result, CurvedQuiverplotSet)

        # Check that result has required attributes
        assert hasattr(result, "lines")
        assert hasattr(result, "arrows")
        assert hasattr(result, "resolution")
        assert hasattr(result, "magnitude")

        plt.close(fig)

    def test_curved_quiver_with_parameters(self, sample_data):
        """Test curved_quiver with various parameters."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = curved_quiver(
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

        assert isinstance(result, CurvedQuiverplotSet)
        assert result.linewidth == 2.0
        assert result.color == "red"
        assert result.arrowsize == 1.5
        assert result.arrowstyle == "->"

        plt.close(fig)

    def test_curved_quiver_no_axes(self, sample_data):
        """Test curved_quiver without providing axes."""
        # This should use current axes
        plt.figure(figsize=(8, 6))

        result = curved_quiver(sample_data, x="x", y="y", u="u", v="v", density=1)

        assert isinstance(result, CurvedQuiverplotSet)
        plt.close()

    def test_curved_quiver_integration_directions(self, sample_data):
        """Test different integration directions."""
        directions = ["forward", "backward", "both"]

        for direction in directions:
            fig, ax = plt.subplots(figsize=(6, 4))

            result = curved_quiver(
                sample_data,
                x="x",
                y="y",
                u="u",
                v="v",
                ax=ax,
                integration_direction=direction,
            )

            assert isinstance(result, CurvedQuiverplotSet)
            assert result.integration_direction == direction

            plt.close(fig)

    def test_curved_quiver_with_start_points(self, sample_data):
        """Test curved_quiver with custom start points."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Define custom start points
        start_points = np.array([[-5, -2], [0, 0], [5, 2]])

        result = curved_quiver(
            sample_data, x="x", y="y", u="u", v="v", ax=ax, start_points=start_points
        )

        assert isinstance(result, CurvedQuiverplotSet)
        np.testing.assert_array_equal(result.start_points, start_points)

        plt.close(fig)

    @patch("skyborn.plot.curved_quiver_plot.velovect")
    def test_curved_quiver_calls_velovect(self, mock_velovect, sample_data):
        """Test that curved_quiver properly calls velovect."""
        mock_result = Mock(spec=CurvedQuiverplotSet)
        mock_velovect.return_value = mock_result

        fig, ax = plt.subplots(figsize=(8, 6))

        result = curved_quiver(
            sample_data, x="x", y="y", u="u", v="v", ax=ax, density=2, linewidth=1.5
        )

        # Check that velovect was called
        mock_velovect.assert_called_once()

        # Check some of the arguments passed to velovect
        call_args = mock_velovect.call_args
        assert call_args[0][0] == ax  # First argument should be axes
        assert call_args[1]["density"] == 2
        assert call_args[1]["linewidth"] == 1.5

        assert result == mock_result
        plt.close(fig)


class TestCurvedQuiverLegend:
    """Test the CurvedQuiverLegend class."""

    @pytest.fixture
    def mock_curved_quiver_set(self):
        """Create a mock CurvedQuiverplotSet for testing."""
        mock_set = Mock(spec=CurvedQuiverplotSet)
        mock_set.resolution = 0.5
        mock_set.magnitude = np.array([[0.5, 1.0], [1.5, 2.0]])
        return mock_set

    @pytest.fixture
    def sample_axes(self):
        """Create sample matplotlib axes."""
        fig, ax = plt.subplots(figsize=(8, 6))
        return ax, fig

    def test_curved_quiver_legend_basic(self, sample_axes, mock_curved_quiver_set):
        """Test basic CurvedQuiverLegend creation."""
        ax, fig = sample_axes

        legend = CurvedQuiverLegend(
            ax=ax, curved_quiver_set=mock_curved_quiver_set, U=2.0, units="m/s"
        )

        # Check basic attributes
        assert legend.ax == ax
        assert legend.curved_quiver_set == mock_curved_quiver_set
        assert legend.U == 2.0
        assert legend.units == "m/s"
        assert legend.show_units == True

        # Check that components were created
        assert hasattr(legend, "patch")
        assert hasattr(legend, "arrow")
        assert hasattr(legend, "text")

        plt.close(fig)

    def test_curved_quiver_legend_positions(self, sample_axes, mock_curved_quiver_set):
        """Test different legend positions."""
        ax, fig = sample_axes
        positions = ["lower left", "lower right", "upper left", "upper right"]

        for pos in positions:
            legend = CurvedQuiverLegend(
                ax=ax,
                curved_quiver_set=mock_curved_quiver_set,
                U=2.0,
                units="m/s",
                loc=pos,
            )

            assert legend.loc == pos
            # Position should be calculated
            assert hasattr(legend, "x")
            assert hasattr(legend, "y")

        plt.close(fig)

    def test_curved_quiver_legend_label_positions(
        self, sample_axes, mock_curved_quiver_set
    ):
        """Test different label positions."""
        ax, fig = sample_axes
        label_positions = ["N", "S", "E", "W"]

        for labelpos in label_positions:
            legend = CurvedQuiverLegend(
                ax=ax,
                curved_quiver_set=mock_curved_quiver_set,
                U=2.0,
                units="m/s",
                labelpos=labelpos,
            )

            assert legend.labelpos == labelpos

        plt.close(fig)

    def test_curved_quiver_legend_no_units(self, sample_axes, mock_curved_quiver_set):
        """Test legend with no units (auto-centering)."""
        ax, fig = sample_axes

        legend = CurvedQuiverLegend(
            ax=ax,
            curved_quiver_set=mock_curved_quiver_set,
            U=5.0,
            units="",  # Empty units
        )

        assert legend.units == ""
        assert legend.show_units == False
        assert legend.center_label == True  # Should auto-center
        assert legend.text_content == "5.0"  # Should not include units

        plt.close(fig)

    def test_curved_quiver_legend_with_units(self, sample_axes, mock_curved_quiver_set):
        """Test legend with units."""
        ax, fig = sample_axes

        legend = CurvedQuiverLegend(
            ax=ax, curved_quiver_set=mock_curved_quiver_set, U=3.5, units="m/s"
        )

        assert legend.units == "m/s"
        assert legend.show_units == True
        assert legend.center_label == False

        plt.close(fig)

    def test_curved_quiver_legend_custom_properties(
        self, sample_axes, mock_curved_quiver_set
    ):
        """Test legend with custom properties."""
        ax, fig = sample_axes

        arrow_props = {"color": "blue", "linewidth": 2}
        patch_props = {"facecolor": "yellow", "alpha": 0.8}
        text_props = {"fontsize": 12, "color": "red"}

        legend = CurvedQuiverLegend(
            ax=ax,
            curved_quiver_set=mock_curved_quiver_set,
            U=2.0,
            units="m/s",
            arrow_props=arrow_props,
            patch_props=patch_props,
            text_props=text_props,
        )

        # Check that properties were applied
        assert legend.arrow_props["color"] == "blue"
        assert legend.arrow_props["linewidth"] == 2
        assert legend.patch_props["facecolor"] == "yellow"
        assert legend.patch_props["alpha"] == 0.8
        assert legend.text_props["fontsize"] == 12
        assert legend.text_props["color"] == "red"

        plt.close(fig)

    def test_curved_quiver_legend_arrow_length_calculation(
        self, sample_axes, mock_curved_quiver_set
    ):
        """Test arrow length calculation."""
        ax, fig = sample_axes

        # Test with different wind speeds
        wind_speeds = [1.0, 2.0, 5.0, 10.0]

        for U in wind_speeds:
            legend = CurvedQuiverLegend(
                ax=ax,
                curved_quiver_set=mock_curved_quiver_set,
                U=U,
                reference_speed=2.0,
                max_arrow_length=0.08,
            )

            # Arrow length should be proportional to wind speed
            expected_scale = U / 2.0  # reference_speed = 2.0
            assert hasattr(legend, "arrow_length")
            assert legend.arrow_length > 0

        plt.close(fig)

    def test_curved_quiver_legend_size_adjustment(
        self, sample_axes, mock_curved_quiver_set
    ):
        """Test that legend box size adjusts for content."""
        ax, fig = sample_axes

        # Create legend with long text
        legend = CurvedQuiverLegend(
            ax=ax,
            curved_quiver_set=mock_curved_quiver_set,
            U=15.5,
            units="very_long_unit_string",
            width=0.1,  # Initial small width
            height=0.05,  # Initial small height
        )

        # Box should have been resized to fit content
        assert legend.width >= 0.1  # Should be at least the initial width
        assert legend.height >= 0.05  # Should be at least the initial height

        plt.close(fig)

    def test_curved_quiver_legend_invalid_position(
        self, sample_axes, mock_curved_quiver_set
    ):
        """Test invalid legend position raises error."""
        ax, fig = sample_axes

        with pytest.raises(ValueError, match="loc must be one of"):
            CurvedQuiverLegend(
                ax=ax,
                curved_quiver_set=mock_curved_quiver_set,
                U=2.0,
                units="m/s",
                loc="invalid_position",
            )

        plt.close(fig)


class TestAddCurvedQuiverkey:
    """Test the add_curved_quiverkey convenience function."""

    @pytest.fixture
    def mock_curved_quiver_set(self):
        """Create a mock CurvedQuiverplotSet for testing."""
        mock_set = Mock(spec=CurvedQuiverplotSet)
        mock_set.resolution = 0.5
        mock_set.magnitude = np.array([[0.5, 1.0], [1.5, 2.0]])
        return mock_set

    def test_add_curved_quiverkey_basic(self, mock_curved_quiver_set):
        """Test basic add_curved_quiverkey functionality."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = add_curved_quiverkey(
            ax=ax, curved_quiver_set=mock_curved_quiver_set, U=2.0, units="m/s"
        )

        # Should return a CurvedQuiverLegend instance
        assert isinstance(result, CurvedQuiverLegend)
        assert result.U == 2.0
        assert result.units == "m/s"
        assert result.loc == "lower right"  # default
        assert result.labelpos == "E"  # default

        plt.close(fig)

    def test_add_curved_quiverkey_with_parameters(self, mock_curved_quiver_set):
        """Test add_curved_quiverkey with custom parameters."""
        fig, ax = plt.subplots(figsize=(8, 6))

        result = add_curved_quiverkey(
            ax=ax,
            curved_quiver_set=mock_curved_quiver_set,
            U=5.0,
            units="knots",
            loc="upper left",
            labelpos="W",
            width=0.2,
            height=0.1,
        )

        assert isinstance(result, CurvedQuiverLegend)
        assert result.U == 5.0
        assert result.units == "knots"
        assert result.loc == "upper left"
        assert result.labelpos == "W"
        assert result.width == 0.2
        assert result.height == 0.1

        plt.close(fig)

    def test_add_curved_quiverkey_kwargs_passed(self, mock_curved_quiver_set):
        """Test that additional kwargs are passed to CurvedQuiverLegend."""
        fig, ax = plt.subplots(figsize=(8, 6))

        custom_arrow_props = {"color": "green", "linewidth": 3}

        result = add_curved_quiverkey(
            ax=ax,
            curved_quiver_set=mock_curved_quiver_set,
            U=3.0,
            arrow_props=custom_arrow_props,
            margin=0.05,
        )

        assert isinstance(result, CurvedQuiverLegend)
        assert result.arrow_props["color"] == "green"
        assert result.arrow_props["linewidth"] == 3
        assert result.margin == 0.05

        plt.close(fig)


class TestCurvedQuiverIntegration:
    """Integration tests for curved quiver plotting."""

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

    def test_curved_quiver_with_legend_integration(self, wind_data):
        """Test full integration of curved quiver with legend."""
        fig, ax = plt.subplots(figsize=(10, 8))

        # Create curved quiver plot
        quiver_set = curved_quiver(
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
        legend = add_curved_quiverkey(
            ax=ax,
            curved_quiver_set=quiver_set,
            U=10.0,
            units="m/s",
            loc="upper right",
            labelpos="E",
        )

        # Verify everything was created properly
        assert isinstance(quiver_set, CurvedQuiverplotSet)
        assert isinstance(legend, CurvedQuiverLegend)

        # Check that axes has the expected artists
        # Lines and patches should be added to the axes
        assert len(ax.collections) > 0  # Should have LineCollection
        # Should have arrow patches and legend patch
        assert len(ax.patches) > 0

        # Set some basic plot properties
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Curved Quiver Plot with Legend")
        ax.grid(True, alpha=0.3)

        plt.close(fig)

    def test_multiple_legends_different_positions(self, wind_data):
        """Test multiple legends in different positions."""
        fig, ax = plt.subplots(figsize=(12, 8))

        # Create curved quiver plot
        quiver_set = curved_quiver(
            wind_data, x="lon", y="lat", u="u", v="v", ax=ax, density=1
        )

        # Add multiple legends
        positions = ["lower left", "lower right", "upper left"]
        speeds = [5.0, 10.0, 15.0]

        legends = []
        for pos, speed in zip(positions, speeds):
            legend = add_curved_quiverkey(
                ax=ax,
                curved_quiver_set=quiver_set,
                U=speed,
                units="m/s",
                loc=pos,
                width=0.15,
                height=0.08,
            )
            legends.append(legend)

        # All legends should be created
        assert len(legends) == 3
        for legend in legends:
            assert isinstance(legend, CurvedQuiverLegend)

        plt.close(fig)
