"""
Tests for skyborn.gridfill module.

These tests verify the functionality of the gridfill module for filling
missing values in gridded data using Poisson equation solvers.
"""

import numpy as np
import numpy.ma as ma
import pytest

from skyborn.gridfill import fill


class TestGridfill:
    """Test suite for gridfill functionality."""

    def test_basic_fill(self):
        """Test basic filling functionality with a simple case."""
        # Create a simple 5x5 grid with one missing value in the center
        data = np.arange(25).reshape(5, 5).astype(float)
        mask = np.zeros_like(data, dtype=bool)
        mask[2, 2] = True  # Center point missing

        masked_data = ma.array(data, mask=mask)

        # Fill the missing value
        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

        # Check convergence
        assert converged[0], "Algorithm should converge for simple case"

        # Check that the filled value is reasonable (should be close to 12)
        filled_value = filled[2, 2]
        assert 10 < filled_value < 14, f"Filled value {filled_value} seems unreasonable"

        # Check that non-missing values are preserved
        assert np.allclose(
            filled[mask == False], data[mask == False]
        ), "Non-missing values should be preserved"

    def test_rectangular_gap(self):
        """Test filling a rectangular gap."""
        # Create a 10x10 grid with a 3x3 gap in the center
        data = np.random.rand(10, 10)
        original_data = data.copy()

        mask = np.zeros_like(data, dtype=bool)
        mask[4:7, 4:7] = True  # 3x3 gap in center

        masked_data = ma.array(data, mask=mask)

        # Fill the gap
        filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)

        # Check convergence
        assert converged[0], "Algorithm should converge"

        # Check that boundary values are preserved
        boundary_mask = ~mask
        assert np.allclose(
            filled[boundary_mask], original_data[boundary_mask]
        ), "Boundary values should be preserved"

        # Check that filled values are not NaN or infinite
        assert np.all(np.isfinite(filled[mask])), "Filled values should be finite"

    def test_zonal_initialization(self):
        """Test zonal mean initialization."""
        # Create data with strong zonal pattern
        lon = np.linspace(0, 360, 20)
        lat = np.linspace(-90, 90, 10)
        LON, LAT = np.meshgrid(lon, lat)

        # Create a sinusoidal pattern in latitude
        data = np.sin(np.radians(LAT)) + 0.1 * np.random.rand(*LAT.shape)

        # Create a gap in the middle
        mask = np.zeros_like(data, dtype=bool)
        mask[4:6, 8:12] = True

        masked_data = ma.array(data, mask=mask)

        # Fill with and without zonal initialization
        filled_zero, _ = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, initzonal=False, verbose=False
        )
        filled_zonal, _ = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, initzonal=True, verbose=False
        )

        # Both should be finite and different
        assert np.all(np.isfinite(filled_zero[mask]))
        assert np.all(np.isfinite(filled_zonal[mask]))

        # For strongly zonal data, zonal initialization often gives better results
        # but we just test that both work
        assert not np.allclose(
            filled_zero[mask], filled_zonal[mask]
        ), "Different initialization should give different results"

    def test_cyclic_boundary(self):
        """Test cyclic boundary conditions."""
        # Create a 8x16 grid (like a simplified global grid)
        data = np.random.rand(8, 16)

        # Make the data cyclic by ensuring left and right edges are similar
        data[:, 0] = data[:, -1]

        # Create a gap near the edge to test cyclicity
        mask = np.zeros_like(data, dtype=bool)
        mask[3:5, 0:2] = True  # Gap at left edge
        mask[3:5, -2:] = True  # Gap at right edge

        masked_data = ma.array(data, mask=mask)

        # Fill with cyclic boundary
        filled_cyclic, conv_cyc = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, cyclic=True, verbose=False
        )

        # Fill without cyclic boundary
        filled_non_cyclic, conv_non = fill(
            masked_data, xdim=1, ydim=0, eps=1e-4, cyclic=False, verbose=False
        )

        # Both should converge
        assert conv_cyc[0] and conv_non[0], "Both methods should converge"

        # Results should be different
        assert not np.allclose(
            filled_cyclic[mask], filled_non_cyclic[mask]
        ), "Cyclic and non-cyclic should give different results"

    def test_convergence_parameters(self):
        """Test different convergence parameters."""
        # Create test data
        data = np.random.rand(15, 20)
        mask = np.zeros_like(data, dtype=bool)
        mask[7:10, 9:12] = True

        masked_data = ma.array(data, mask=mask)

        # Test with different eps values
        filled_loose, conv_loose = fill(masked_data, xdim=1, ydim=0, eps=1e-2)
        filled_tight, conv_tight = fill(masked_data, xdim=1, ydim=0, eps=1e-6)

        # Both should converge but with different precision
        assert conv_loose[0] and conv_tight[0], "Both should converge"

        # Tighter tolerance should be closer to looser tolerance
        # (this is a weak test, but ensures numerical stability)
        assert np.allclose(
            filled_loose[~mask], filled_tight[~mask]
        ), "Boundary values should be identical regardless of tolerance"

    def test_error_conditions(self):
        """Test error handling."""
        # Test with non-masked array (should raise TypeError)
        regular_array = np.random.rand(5, 5)

        with pytest.raises(TypeError):
            fill(regular_array, xdim=1, ydim=0, eps=1e-4)

        # Test with invalid dimensions
        data = np.random.rand(5, 5)
        mask = np.zeros_like(data, dtype=bool)
        mask[2, 2] = True
        masked_data = ma.array(data, mask=mask)

        with pytest.raises(ValueError):
            fill(masked_data, xdim=5, ydim=0, eps=1e-4)  # Invalid xdim

        with pytest.raises(ValueError):
            fill(masked_data, xdim=1, ydim=5, eps=1e-4)  # Invalid ydim

    def test_different_relaxation_parameters(self):
        """Test different relaxation constants."""
        # Create test data
        data = np.random.rand(10, 15)
        mask = np.zeros_like(data, dtype=bool)
        mask[4:6, 6:9] = True

        masked_data = ma.array(data, mask=mask)

        # Test with different relaxation parameters
        relax_values = [0.5, 0.6, 0.8]
        results = []

        for relax in relax_values:
            filled, converged = fill(
                masked_data, xdim=1, ydim=0, eps=1e-4, relax=relax, verbose=False
            )
            assert converged[0], f"Should converge with relax={relax}"
            results.append(filled)

        # All results should preserve boundary values
        for result in results:
            assert np.allclose(
                result[~mask], data[~mask]
            ), "Boundary values should be preserved"

    def test_multidimensional_data(self):
        """Test with 3D data (time series of 2D grids)."""
        # Create 3D data: (time, lat, lon)
        nt, nlat, nlon = 5, 8, 12
        data = np.random.rand(nt, nlat, nlon)

        # Create the same mask for all time steps
        mask = np.zeros_like(data, dtype=bool)
        mask[:, 3:5, 5:8] = True  # Same spatial gap for all times

        masked_data = ma.array(data, mask=mask)

        # Fill (should handle each time slice independently)
        filled, converged = fill(masked_data, xdim=2, ydim=1, eps=1e-4)

        # Should converge for all time slices
        assert np.all(converged), "Should converge for all time slices"

        # Check that boundary values are preserved for all times
        assert np.allclose(
            filled[~mask], data[~mask]
        ), "Boundary values should be preserved across all time steps"

        # Check that filled values are finite
        assert np.all(np.isfinite(filled[mask])), "All filled values should be finite"


# Test fill_cube function if iris is available
try:
    import iris
    import iris.tests
    from skyborn.gridfill import fill_cube

    class TestGridfillCube:
        """Test suite for fill_cube functionality with iris cubes."""

        def test_fill_cube_basic(self):
            """Test basic cube filling functionality."""
            # Create a simple test cube
            data = np.random.rand(5, 10, 15)  # (time, lat, lon)

            # Create missing values
            data = ma.array(data)
            data.mask = False
            data.mask[2, 2:4, 6:9] = True  # Create a gap

            # Create coordinates
            time_coord = iris.coords.DimCoord(
                np.arange(5), standard_name="time", units="days since 2000-01-01"
            )
            lat_coord = iris.coords.DimCoord(
                np.linspace(-45, 45, 10), standard_name="latitude", units="degrees"
            )
            lon_coord = iris.coords.DimCoord(
                np.linspace(0, 360, 15), standard_name="longitude", units="degrees"
            )
            lon_coord.circular = True  # Make longitude cyclic

            # Create cube
            cube = iris.cube.Cube(
                data,
                dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
            )

            # Fill the cube
            filled_cube = fill_cube(cube, eps=1e-4, verbose=False)

            # Check that we got a cube back
            assert isinstance(filled_cube, iris.cube.Cube)

            # Check that filled data has no missing values where we had gaps
            original_mask = cube.data.mask
            filled_data = filled_cube.data

            # Where we had missing data, we should now have finite values
            if hasattr(filled_data, "mask"):
                assert not np.any(
                    filled_data.mask[original_mask]
                ), "Previously missing values should now be filled"

            # Check that original data is preserved where not missing
            if hasattr(filled_data, "mask"):
                good_original = ~original_mask
                assert np.allclose(
                    filled_data.data[good_original], cube.data.data[good_original]
                ), "Original data should be preserved"
            else:
                good_original = ~original_mask
                assert np.allclose(
                    filled_data[good_original], cube.data.data[good_original]
                ), "Original data should be preserved"

        def test_fill_cube_full_output(self):
            """Test cube filling with full output."""
            # Create test data with a gap that might not converge easily
            data = np.random.rand(3, 6, 8)
            data = ma.array(data)
            data.mask = False
            data.mask[1, 2:4, 3:6] = True  # Create a gap

            # Create a minimal cube
            cube = iris.cube.Cube(data)

            # Add required coordinates for auto-detection
            lat_coord = iris.coords.DimCoord(
                np.linspace(-30, 30, 6), standard_name="latitude", units="degrees"
            )
            lon_coord = iris.coords.DimCoord(
                np.linspace(0, 360, 8), standard_name="longitude", units="degrees"
            )

            cube.add_dim_coord(lat_coord, 1)
            cube.add_dim_coord(lon_coord, 2)

            # Fill with full output
            filled_cube, not_converged = fill_cube(
                cube, eps=1e-6, itermax=50, full_output=True, verbose=False
            )

            # Check return types
            assert isinstance(filled_cube, iris.cube.Cube)
            assert isinstance(not_converged, np.ndarray)
            assert not_converged.dtype == bool

except ImportError:
    # Skip iris tests if iris is not available
    pass


if __name__ == "__main__":
    # Run basic tests
    test_suite = TestGridfill()

    print("Running gridfill tests...")
    test_suite.test_basic_fill()
    print("✓ Basic fill test passed")

    test_suite.test_rectangular_gap()
    print("✓ Rectangular gap test passed")

    test_suite.test_zonal_initialization()
    print("✓ Zonal initialization test passed")

    test_suite.test_cyclic_boundary()
    print("✓ Cyclic boundary test passed")

    test_suite.test_convergence_parameters()
    print("✓ Convergence parameters test passed")

    test_suite.test_different_relaxation_parameters()
    print("✓ Relaxation parameters test passed")

    test_suite.test_multidimensional_data()
    print("✓ Multidimensional data test passed")

    print("All gridfill tests passed!")
