"""
Tests for skyborn.plot.plotting module.

This module tests the plotting functionality in skyborn.plot.plotting,
including matplotlib integration and image comparison tests.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from skyborn.plot.plotting import add_equal_axes, createFigure


class TestCreateFigure:
    """Test the createFigure function."""

    def test_createFigure_basic(self):
        """Test basic createFigure functionality."""
        fig = createFigure()

        # Check default parameters
        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 8
        assert fig.dpi == 300

        plt.close(fig)

    def test_createFigure_custom_figsize(self):
        """Test createFigure with custom figsize."""
        fig = createFigure(figsize=(8, 6))

        assert fig.get_figwidth() == 8
        assert fig.get_figheight() == 6
        assert fig.dpi == 300

        plt.close(fig)

    def test_createFigure_custom_dpi(self):
        """Test createFigure with custom DPI."""
        fig = createFigure(dpi=150)

        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 8
        assert fig.dpi == 150

        plt.close(fig)

    def test_createFigure_with_subplotAdj(self):
        """Test createFigure with subplot adjustments."""
        subplot_params = {
            "left": 0.1,
            "right": 0.9,
            "top": 0.9,
            "bottom": 0.1,
            "hspace": 0.2,
            "wspace": 0.2,
        }

        fig = createFigure(subplotAdj=subplot_params)

        # Check that figure was created successfully
        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 8

        plt.close(fig)

    def test_createFigure_kwargs(self):
        """Test createFigure with additional keyword arguments."""
        fig = createFigure(facecolor="lightblue", edgecolor="red")

        # Check that figure was created successfully
        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 8

        plt.close(fig)

    def test_createFigure_all_parameters(self):
        """Test createFigure with all parameters specified."""
        subplot_params = {"left": 0.2, "right": 0.8}

        fig = createFigure(
            figsize=(10, 5),
            dpi=200,
            subplotAdj=subplot_params,
            facecolor="white",
            edgecolor="black",
        )

        assert fig.get_figwidth() == 10
        assert fig.get_figheight() == 5
        assert fig.dpi == 200

        plt.close(fig)


class TestAddEqualAxes:
    """Test the add_equal_axes function."""

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_add_equal_axes_right(self):
        """Test adding axes to the right."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create some sample data
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        ax.plot(x, y, "b-", linewidth=2)
        ax.set_xlabel("X axis")
        ax.set_ylabel("Y axis")
        ax.set_title("Main Plot")
        ax.grid(True, alpha=0.3)

        # Add new axes to the right
        ax_new = add_equal_axes(ax, loc="right", pad=0.05, width=0.1)

        # Add content to new axes
        ax_new.plot(y, x, "r-", linewidth=2)
        ax_new.set_xlabel("Y values")
        ax_new.set_ylabel("X values")
        ax_new.grid(True, alpha=0.3)

        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_add_equal_axes_left(self):
        """Test adding axes to the left."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create some sample data
        x = np.linspace(0, 2 * np.pi, 100)
        y = np.cos(x)
        ax.plot(x, y, "g-", linewidth=2)
        ax.set_xlabel("X axis")
        ax.set_ylabel("Cosine values")
        ax.set_title("Main Plot")
        ax.grid(True, alpha=0.3)

        # Add new axes to the left
        ax_new = add_equal_axes(ax, loc="left", pad=0.05, width=0.15)

        # Add histogram to new axes
        ax_new.hist(y, bins=20, orientation="horizontal", alpha=0.7, color="orange")
        ax_new.set_xlabel("Frequency")
        ax_new.set_ylabel("Cosine values")

        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_add_equal_axes_top(self):
        """Test adding axes to the top."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create 2D data
        x = np.linspace(-2, 2, 50)
        y = np.linspace(-2, 2, 50)
        X, Y = np.meshgrid(x, y)
        Z = np.exp(-(X**2 + Y**2))

        # Main contour plot
        cs = ax.contourf(X, Y, Z, levels=20, cmap="viridis")
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
        ax.set_title("2D Gaussian")

        # Add colorbar-like axes at top
        ax_new = add_equal_axes(ax, loc="top", pad=0.05, width=0.05)

        # Show mean values along x-axis
        mean_z = np.mean(Z, axis=0)
        ax_new.plot(x, mean_z, "k-", linewidth=2)
        ax_new.set_xlabel("X coordinate")
        ax_new.set_ylabel("Mean Z")
        ax_new.grid(True, alpha=0.3)

        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_add_equal_axes_bottom(self):
        """Test adding axes to the bottom."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create scatter plot
        np.random.seed(42)
        x = np.random.randn(200)
        y = 2 * x + np.random.randn(200) * 0.5

        ax.scatter(x, y, alpha=0.6, c="blue", s=30)
        ax.set_xlabel("X values")
        ax.set_ylabel("Y values")
        ax.set_title("Scatter Plot with Correlation")
        ax.grid(True, alpha=0.3)

        # Add marginal histogram at bottom
        ax_new = add_equal_axes(ax, loc="bottom", pad=0.05, width=0.1)

        ax_new.hist(x, bins=30, alpha=0.7, color="red", edgecolor="black")
        ax_new.set_xlabel("X values")
        ax_new.set_ylabel("Frequency")
        ax_new.grid(True, alpha=0.3)

        return fig

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_add_equal_axes_multiple(self):
        """Test adding multiple axes around a central plot."""
        fig, ax = plt.subplots(figsize=(10, 8))

        # Central heatmap
        np.random.seed(123)
        data = np.random.randn(20, 20)
        im = ax.imshow(data, cmap="RdBu_r", aspect="equal")
        ax.set_title("Central Heatmap")
        ax.set_xlabel("X index")
        ax.set_ylabel("Y index")

        # Add axes on all sides
        ax_right = add_equal_axes(ax, loc="right", pad=0.02, width=0.08)
        ax_top = add_equal_axes(ax, loc="top", pad=0.02, width=0.08)

        # Right: row means
        row_means = np.mean(data, axis=1)
        ax_right.plot(row_means, np.arange(len(row_means)), "k-", linewidth=2)
        ax_right.set_xlabel("Row mean")
        ax_right.set_ylim(ax.get_ylim())
        ax_right.grid(True, alpha=0.3)

        # Top: column means
        col_means = np.mean(data, axis=0)
        ax_top.plot(np.arange(len(col_means)), col_means, "k-", linewidth=2)
        ax_top.set_ylabel("Col mean")
        ax_top.set_xlim(ax.get_xlim())
        ax_top.grid(True, alpha=0.3)

        return fig

    def test_add_equal_axes_array_input(self):
        """Test add_equal_axes with array of axes as input."""
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        # Plot on both axes
        x = np.linspace(0, 10, 100)
        axes[0].plot(x, np.sin(x))
        axes[1].plot(x, np.cos(x))

        # Add axes to the right of both subplots
        ax_new = add_equal_axes(axes, loc="right", pad=0.05, width=0.1)

        # Check that new axes was created
        assert ax_new is not None
        assert hasattr(ax_new, "plot")  # Should be a valid axes object

        plt.close(fig)

    def test_add_equal_axes_invalid_location(self):
        """Test add_equal_axes with invalid location parameter."""
        fig, ax = plt.subplots()

        with pytest.raises(ValueError):
            add_equal_axes(ax, loc="invalid_location", pad=0.05, width=0.1)

        plt.close(fig)

    def test_add_equal_axes_parameters(self):
        """Test add_equal_axes parameter validation."""
        fig, ax = plt.subplots()

        # Test with various parameter values
        ax_new1 = add_equal_axes(ax, loc="right", pad=0.1, width=0.2)
        ax_new2 = add_equal_axes(ax, loc="left", pad=0.01, width=0.05)

        # Check that axes were created successfully
        assert ax_new1 is not None
        assert ax_new2 is not None

        plt.close(fig)


class TestPlottingIntegration:
    """Integration tests for plotting module."""

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_plotting_with_climate_data(self, sample_2d_field):
        """Test plotting functions with realistic climate data."""
        fig, ax = plt.subplots(figsize=(10, 6))

        # Create a temperature map
        temp = sample_2d_field
        lat = temp.lat
        lon = temp.lon

        # Plot temperature field
        im = ax.contourf(lon, lat, temp, levels=20, cmap="RdYlBu_r")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Surface Temperature")

        # Add colorbar using add_equal_axes
        cax = add_equal_axes(ax, loc="right", pad=0.02, width=0.03)
        fig.colorbar(im, cax=cax, label="Temperature (K)")

        # Add marginal plots
        ax_top = add_equal_axes(ax, loc="top", pad=0.02, width=0.08)
        ax_right = add_equal_axes(
            ax, loc="right", pad=0.15, width=0.08
        )  # Further right

        # Zonal mean (average over longitude)
        zonal_mean = temp.mean(dim="lon")
        ax_right.plot(zonal_mean, lat, "k-", linewidth=2)
        ax_right.set_xlabel("Zonal mean (K)")
        ax_right.set_ylim(ax.get_ylim())
        ax_right.grid(True, alpha=0.3)

        # Meridional mean (average over latitude)
        meridional_mean = temp.mean(dim="lat")
        ax_top.plot(lon, meridional_mean, "k-", linewidth=2)
        ax_top.set_ylabel("Meridional mean (K)")
        ax_top.set_xlim(ax.get_xlim())
        ax_top.grid(True, alpha=0.3)

        return fig

    def test_plotting_error_handling(self):
        """Test error handling in plotting functions."""
        fig, ax = plt.subplots()

        # Test with invalid parameters
        with pytest.raises(ValueError):
            add_equal_axes(ax, loc="nowhere", pad=0.05, width=0.1)

        plt.close(fig)

    @pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
    def test_complex_layout(self, sample_climate_data):
        """Test complex multi-panel layout using add_equal_axes."""
        temp = sample_climate_data["temperature"]
        precip = sample_climate_data["precipitation"]

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle("Climate Data Analysis")

        # Temperature maps for different months
        temp_jan = temp.isel(time=0)  # January
        temp_jul = temp.isel(time=6)  # July

        im1 = axes[0, 0].contourf(
            temp.lon, temp.lat, temp_jan, levels=20, cmap="RdYlBu_r"
        )
        axes[0, 0].set_title("Temperature - January")
        axes[0, 0].set_xlabel("Longitude")
        axes[0, 0].set_ylabel("Latitude")

        im2 = axes[0, 1].contourf(
            temp.lon, temp.lat, temp_jul, levels=20, cmap="RdYlBu_r"
        )
        axes[0, 1].set_title("Temperature - July")
        axes[0, 1].set_xlabel("Longitude")

        # Precipitation maps
        precip_jan = precip.isel(time=0)
        precip_jul = precip.isel(time=6)

        im3 = axes[1, 0].contourf(
            precip.lon, precip.lat, precip_jan, levels=20, cmap="Blues"
        )
        axes[1, 0].set_title("Precipitation - January")
        axes[1, 0].set_xlabel("Longitude")
        axes[1, 0].set_ylabel("Latitude")

        im4 = axes[1, 1].contourf(
            precip.lon, precip.lat, precip_jul, levels=20, cmap="Blues"
        )
        axes[1, 1].set_title("Precipitation - July")
        axes[1, 1].set_xlabel("Longitude")

        # Add colorbars using add_equal_axes
        for i, (ax, im, label) in enumerate(
            [
                (axes[0, 0], im1, "Temperature (K)"),
                (axes[0, 1], im2, "Temperature (K)"),
                (axes[1, 0], im3, "Precipitation (mm/day)"),
                (axes[1, 1], im4, "Precipitation (mm/day)"),
            ]
        ):
            cax = add_equal_axes(ax, loc="right", pad=0.02, width=0.02)
            fig.colorbar(im, cax=cax, label=label if i % 2 == 0 else "")

        plt.tight_layout()
        return fig


if __name__ == "__main__":
    # Quick test runner
    pytest.main([__file__, "-v"])
