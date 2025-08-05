"""
Tests for skyborn.plot.modplot module.

This module tests the modular plotting utilities and streamline functionality,
including the velovect function and supporting classes.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from unittest.mock import Mock, patch
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import FancyArrowPatch

from skyborn.plot.modplot import (
    velovect,
    CurvedQuiverplotSet,
    Grid,
    StreamMask,
    DomainMap,
    interpgrid,
    _gen_starting_points,
    _get_integrator,
    InvalidIndexError,
    TerminateTrajectory,
    OutOfBounds,
)


class TestVelovect:
    """Test the velovect function."""

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

    def test_velovect_basic(self, sample_vector_field):
        """Test basic velovect functionality."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        result = velovect(ax, x, y, u, v)

        # Check return type
        assert isinstance(result, CurvedQuiverplotSet)

        # Check that result has required attributes
        assert hasattr(result, "lines")
        assert hasattr(result, "arrows")
        assert hasattr(result, "resolution")
        assert hasattr(result, "magnitude")
        assert hasattr(result, "axes")

        # Check types of components
        assert isinstance(result.lines, LineCollection)
        assert isinstance(result.arrows, PatchCollection)

        plt.close(fig)

    def test_velovect_with_parameters(self, sample_vector_field):
        """Test velovect with various parameters."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        result = velovect(
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

        assert isinstance(result, CurvedQuiverplotSet)
        assert result.linewidth == 2.0
        assert result.color == "red"
        assert result.arrowsize == 1.5
        assert result.arrowstyle == "->"
        assert result.zorder == 5

        plt.close(fig)

    def test_velovect_integration_directions(self, sample_vector_field):
        """Test different integration directions."""
        x, y, u, v = sample_vector_field
        directions = ["forward", "backward", "both"]

        for direction in directions:
            fig, ax = plt.subplots(figsize=(6, 4))

            result = velovect(ax, x, y, u, v, integration_direction=direction)

            assert isinstance(result, CurvedQuiverplotSet)
            assert result.integration_direction == direction

            plt.close(fig)

    def test_velovect_with_start_points(self, sample_vector_field):
        """Test velovect with custom start points."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Define custom start points
        start_points = np.array([[-1, -1], [0, 0], [1, 1]])

        result = velovect(ax, x, y, u, v, start_points=start_points)

        assert isinstance(result, CurvedQuiverplotSet)
        np.testing.assert_array_equal(result.start_points, start_points)

        plt.close(fig)

    def test_velovect_variable_linewidth(self, sample_vector_field):
        """Test velovect with variable linewidth."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create variable linewidth array
        linewidth = np.sqrt(u**2 + v**2)  # Linewidth proportional to speed

        result = velovect(ax, x, y, u, v, linewidth=linewidth)

        assert isinstance(result, CurvedQuiverplotSet)
        np.testing.assert_array_equal(result.linewidth, linewidth)

        plt.close(fig)

    def test_velovect_variable_color(self, sample_vector_field):
        """Test velovect with variable color."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create variable color array
        color = np.sqrt(u**2 + v**2)  # Color proportional to speed

        result = velovect(ax, x, y, u, v, color=color, cmap="viridis")

        assert isinstance(result, CurvedQuiverplotSet)
        np.testing.assert_array_equal(result.color, color)
        assert result.cmap.name == "viridis"

        plt.close(fig)

    def test_velovect_invalid_array_shapes(self, sample_vector_field):
        """Test velovect with invalid array shapes."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create mismatched u array
        u_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            velovect(ax, x, y, u_wrong, v)

        plt.close(fig)

    def test_velovect_invalid_color_shape(self, sample_vector_field):
        """Test velovect with invalid color array shape."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create wrong-shaped color array
        color_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            velovect(ax, x, y, u, v, color=color_wrong)

        plt.close(fig)

    def test_velovect_invalid_linewidth_shape(self, sample_vector_field):
        """Test velovect with invalid linewidth array shape."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create wrong-shaped linewidth array
        linewidth_wrong = np.ones((5, 5))  # Wrong shape

        with pytest.raises(ValueError, match="must match the shape"):
            velovect(ax, x, y, u, v, linewidth=linewidth_wrong)

        plt.close(fig)

    def test_velovect_outside_boundary_start_points(self, sample_vector_field):
        """Test velovect with start points outside data boundaries."""
        x, y, u, v = sample_vector_field
        fig, ax = plt.subplots(figsize=(8, 6))

        # Create start points outside the domain
        start_points = np.array(
            [[-10, -10], [10, 10]]  # Outside domain  # Outside domain
        )

        with pytest.raises(
            ValueError, match="Starting point .* outside of data boundaries"
        ):
            velovect(ax, x, y, u, v, start_points=start_points)

        plt.close(fig)


class TestCurvedQuiverplotSet:
    """Test the CurvedQuiverplotSet class."""

    @pytest.fixture
    def mock_collections(self):
        """Create mock collections for testing."""
        lines = Mock(spec=LineCollection)
        arrows = Mock(spec=PatchCollection)
        return lines, arrows

    def test_curved_quiverplot_set_creation(self, mock_collections):
        """Test CurvedQuiverplotSet creation."""
        lines, arrows = mock_collections
        magnitude = np.array([[1, 2], [3, 4]])

        quiver_set = CurvedQuiverplotSet(
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

    def test_curved_quiverplot_set_scale_factor_forward(self, mock_collections):
        """Test scale factor for forward integration."""
        lines, arrows = mock_collections

        quiver_set = CurvedQuiverplotSet(
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

    def test_curved_quiverplot_set_scale_factor_backward(self, mock_collections):
        """Test scale factor for backward integration."""
        lines, arrows = mock_collections

        quiver_set = CurvedQuiverplotSet(
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

    def test_curved_quiverplot_set_scaling_methods(self, mock_collections):
        """Test scaling methods."""
        lines, arrows = mock_collections

        quiver_set = CurvedQuiverplotSet(
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

        # Try to update to already occupied cell
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
    """Integration tests for modplot functionality."""

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

    def test_full_velovect_integration(self, complex_vector_field):
        """Test full velovect functionality with complex field."""
        x, y, u, v = complex_vector_field
        fig, ax = plt.subplots(figsize=(10, 8))

        # Test with multiple parameters
        result = velovect(
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
        assert isinstance(result, CurvedQuiverplotSet)
        assert result.density == 1.5
        assert result.integration_direction == "both"
        assert result.grains == 20
        assert result.broken_streamlines == True

        # Check that collections were added to axes
        assert len(ax.collections) > 0
        assert len(ax.patches) > 0

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
            velovect(ax, x, y, np.ones((3, 5)), v)

        # Invalid integration direction
        with pytest.raises(ValueError):
            velovect(ax, x, y, u, v, integration_direction="invalid")

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

        result = velovect(ax, x, y, u, v, density=1)

        # Should handle masked data gracefully
        assert isinstance(result, CurvedQuiverplotSet)

        plt.close(fig)
