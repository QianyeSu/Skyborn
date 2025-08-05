"""
Comprehensive tests for skyborn.interp module.

This module tests all interpolation functionality including regridding,
hybrid-sigma interpolation, and multidimensional spatial interpolation.
"""

import pytest
import numpy as np
import xarray as xr
from numpy.testing import assert_array_almost_equal, assert_array_equal


class TestInterpModuleImports:
    """Test that all interpolation module components can be imported."""

    def test_regridding_imports(self):
        """Test regridding imports."""
        from skyborn.interp import (
            Grid,
            Regridder,
            NearestRegridder,
            BilinearRegridder,
            ConservativeRegridder,
            regrid_dataset,
            nearest_neighbor_indices,
        )

        # Basic sanity checks
        assert Grid is not None
        assert Regridder is not None
        assert NearestRegridder is not None
        assert BilinearRegridder is not None
        assert ConservativeRegridder is not None
        assert regrid_dataset is not None
        assert nearest_neighbor_indices is not None

    def test_interpolation_imports(self):
        """Test interpolation function imports."""
        from skyborn.interp import (
            interp_hybrid_to_pressure,
            interp_sigma_to_hybrid,
            interp_multidim,
        )

        # Basic sanity checks
        assert interp_hybrid_to_pressure is not None
        assert interp_sigma_to_hybrid is not None
        assert interp_multidim is not None

    def test_module_all_exports(self):
        """Test that __all__ exports are correct."""
        import skyborn.interp as interp_module

        expected_exports = [
            "Grid",
            "Regridder",
            "NearestRegridder",
            "BilinearRegridder",
            "ConservativeRegridder",
            "nearest_neighbor_indices",
            "regrid_dataset",
            "interp_hybrid_to_pressure",
            "interp_sigma_to_hybrid",
            "interp_multidim",
        ]

        # Check that all expected exports are available
        for export in expected_exports:
            assert hasattr(interp_module, export), f"Missing export: {export}"


class TestInterpModuleIntegration:
    """Integration tests for the interpolation module."""

    @pytest.fixture
    def sample_climate_dataset(self):
        """Create a realistic climate dataset for testing."""
        # Create CMIP-like data structure
        time = np.arange(12)  # 12 months
        lev = np.arange(10)  # 10 vertical levels
        lat = np.linspace(-60, 60, 25)  # 25 latitude points
        lon = np.linspace(0, 360, 36, endpoint=False)  # 36 longitude points

        # Create realistic temperature data
        temp_data = np.random.uniform(200, 300, (12, 10, 25, 36))

        # Create surface pressure with realistic pattern
        psfc_pattern = 101325 - 1000 * np.sin(np.deg2rad(lat)) ** 2
        psfc_data = np.tile(psfc_pattern[None, :, None], (12, 1, 36))

        # Create hybrid coordinates
        hya = np.linspace(0, 80000, 10)  # Pa
        hyb = np.linspace(1.0, 0.0, 10)  # dimensionless

        dataset = xr.Dataset(
            {
                "temperature": (["time", "lev", "lat", "lon"], temp_data),
                "surface_pressure": (["time", "lat", "lon"], psfc_data),
            },
            coords={
                "time": time,
                "lev": lev,
                "lat": lat,
                "lon": lon,
                "hya": (["lev"], hya),
                "hyb": (["lev"], hyb),
            },
            attrs={
                "title": "Test CMIP dataset",
                "institution": "Test Institution",
            },
        )

        return dataset

    def test_regrid_then_pressure_interpolation(self, sample_climate_dataset):
        """Test workflow: regrid spatially then interpolate vertically."""
        from skyborn.interp import Grid, regrid_dataset, interp_hybrid_to_pressure

        dataset = sample_climate_dataset

        # Step 1: Regrid spatially to coarser resolution
        target_lat = np.linspace(-45, 45, 10)
        target_lon = np.linspace(0, 360, 18, endpoint=False)
        target_grid = Grid.from_degrees(target_lon, target_lat)

        regridded = regrid_dataset(dataset, target_grid, method="bilinear")

        # Verify spatial regridding
        assert regridded["temperature"].shape == (
            12,
            10,
            10,
            18,
        )  # (time, lev, lat, lon)
        assert regridded["surface_pressure"].shape == (12, 10, 18)  # (time, lat, lon)

        # Step 2: Interpolate to pressure levels (single time step)
        single_time = regridded.isel(time=0)

        target_pressures = np.array([100000, 85000, 70000, 50000, 30000])  # Pa

        temp_on_pressure = interp_hybrid_to_pressure(
            single_time["temperature"],
            single_time["surface_pressure"],
            single_time["hya"],
            single_time["hyb"],
            target_pressures,
            method="log",
        )

        # Verify pressure interpolation
        assert temp_on_pressure.shape == (5, 10, 18)  # (pressure, lat, lon)

        # Check that results are reasonable
        assert np.all(temp_on_pressure[np.isfinite(temp_on_pressure)] > 150)
        assert np.all(temp_on_pressure[np.isfinite(temp_on_pressure)] < 350)

    def test_multidim_interpolation_workflow(self, sample_climate_dataset):
        """Test multidimensional interpolation workflow."""
        from skyborn.interp import interp_multidim

        dataset = sample_climate_dataset

        # Take surface pressure at one time step
        surface_data = dataset["surface_pressure"].isel(time=0)

        # Interpolate to irregular grid (e.g., station locations)
        station_lats = np.array([-45.2, -23.8, 0.0, 25.5, 55.7])
        station_lons = np.array([170.1, 135.6, 0.0, 78.3, 12.4])

        interpolated = interp_multidim(
            surface_data, station_lats, station_lons, extrap=True
        )

        # Verify results
        assert interpolated.shape == (5, 5)  # (station_lat, station_lon)
        assert isinstance(interpolated, xr.DataArray)

        # Should have reasonable pressure values
        finite_values = interpolated.values[np.isfinite(interpolated.values)]
        assert np.all(finite_values > 80000)  # Reasonable min pressure
        assert np.all(finite_values < 110000)  # Reasonable max pressure

    def test_combined_spatial_vertical_interpolation(self, sample_climate_dataset):
        """Test combined spatial and vertical interpolation."""
        from skyborn.interp import Grid, regrid_dataset, interp_hybrid_to_pressure

        dataset = sample_climate_dataset

        # Define target grids
        target_lat = np.linspace(-30, 30, 7)
        target_lon = np.linspace(0, 360, 12, endpoint=False)
        target_grid = Grid.from_degrees(target_lon, target_lat)
        target_pressures = np.array([100000, 70000, 50000, 30000, 10000])

        # Process multiple time steps
        results = []

        for t in range(0, dataset.sizes["time"], 3):  # Every 3rd time step
            # Spatial regridding
            single_time = dataset.isel(time=t)
            regridded = regrid_dataset(single_time, target_grid, method="conservative")

            # Vertical interpolation
            temp_on_pressure = interp_hybrid_to_pressure(
                regridded["temperature"],
                regridded["surface_pressure"],
                regridded["hya"],
                regridded["hyb"],
                target_pressures,
                method="linear",
            )

            results.append(temp_on_pressure)

        # Verify we processed multiple time steps
        assert len(results) == 4  # 0, 3, 6, 9

        # Check that all results have correct shape
        for result in results:
            assert result.shape == (5, 7, 12)  # (pressure, lat, lon)

        # Check that results are different (temporal variation)
        assert not np.allclose(results[0], results[-1], equal_nan=True)


class TestInterpModuleEdgeCases:
    """Test edge cases and error handling for interpolation module."""

    def test_dimension_order_preservation(self):
        """Test that dimension ordering is preserved across operations."""
        from skyborn.interp import Grid, regrid_dataset

        # Create data with specific dimension order
        data = np.random.randn(5, 8, 12)  # (time, lat, lon)

        dataset = xr.Dataset(
            {"var": (["time", "lat", "lon"], data)},
            coords={
                "time": np.arange(5),
                "lat": np.linspace(-40, 40, 8),
                "lon": np.linspace(0, 360, 12, endpoint=False),
            },
        )

        # Target grid
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 6, endpoint=False), np.linspace(-20, 20, 4)
        )

        # Regrid
        result = regrid_dataset(dataset, target_grid, method="nearest")

        # Check dimension order is preserved
        assert list(result["var"].dims) == ["time", "lat", "lon"]
        assert result["var"].shape == (5, 4, 6)

    def test_coordinate_attributes_preservation(self):
        """Test that coordinate attributes are preserved."""
        from skyborn.interp import Grid, regrid_dataset

        # Create dataset with coordinate attributes
        lat_vals = np.linspace(-45, 45, 10)
        lon_vals = np.linspace(0, 360, 15, endpoint=False)

        dataset = xr.Dataset(
            {"temperature": (["lat", "lon"], np.random.randn(10, 15))},
            coords={
                "lat": xr.DataArray(
                    lat_vals,
                    dims="lat",
                    attrs={"units": "degrees_north", "long_name": "Latitude"},
                ),
                "lon": xr.DataArray(
                    lon_vals,
                    dims="lon",
                    attrs={"units": "degrees_east", "long_name": "Longitude"},
                ),
            },
            attrs={"title": "Test dataset"},
        )

        # Regrid
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 8, endpoint=False), np.linspace(-30, 30, 6)
        )

        result = regrid_dataset(dataset, target_grid, method="bilinear")

        # Dataset attributes should be preserved
        assert result.attrs["title"] == "Test dataset"

        # Note: Coordinate attributes may be updated due to regridding

    def test_mixed_spatial_nonspatial_variables(self):
        """Test handling of datasets with mixed spatial/non-spatial variables."""
        from skyborn.interp import Grid, regrid_dataset

        lat_vals = np.linspace(-30, 30, 6)
        lon_vals = np.linspace(0, 360, 8, endpoint=False)
        time_vals = np.arange(3)

        dataset = xr.Dataset(
            {
                "temperature": (["time", "lat", "lon"], np.random.randn(3, 6, 8)),
                "pressure": (["lat", "lon"], np.random.randn(6, 8)),
                # No spatial dims
                "time_bnds": (["time", "bnds"], np.random.randn(3, 2)),
                # No spatial dims
                "global_mean": (["time"], np.random.randn(3)),
            },
            coords={"time": time_vals, "lat": lat_vals, "lon": lon_vals},
        )

        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-15, 15, 3)
        )

        result = regrid_dataset(dataset, target_grid, method="nearest")

        # Spatial variables should be regridded
        assert result["temperature"].shape == (3, 3, 4)
        assert result["pressure"].shape == (3, 4)

        # Non-spatial variables should remain unchanged
        assert result["time_bnds"].shape == (3, 2)
        assert result["global_mean"].shape == (3,)

        # Values of non-spatial variables should be identical
        assert_array_equal(result["time_bnds"].values, dataset["time_bnds"].values)
        assert_array_equal(result["global_mean"].values, dataset["global_mean"].values)

    def test_empty_dataset_handling(self):
        """Test handling of empty datasets."""
        from skyborn.interp import Grid, regrid_dataset

        # Empty dataset
        empty_dataset = xr.Dataset(
            coords={
                "lat": np.linspace(-45, 45, 5),
                "lon": np.linspace(0, 360, 8, endpoint=False),
            }
        )

        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 4, endpoint=False), np.linspace(-30, 30, 3)
        )

        # Should handle empty dataset gracefully
        result = regrid_dataset(empty_dataset, target_grid, method="nearest")

        # Should have updated coordinates but no data variables
        assert len(result.data_vars) == 0
        assert len(result.coords["lat"]) == 3
        assert len(result.coords["lon"]) == 4


class TestInterpModulePerformanceFeatures:
    """Test performance-related features of interpolation module."""

    def test_regridder_caching_behavior(self):
        """Test that regridders cache expensive computations."""
        from skyborn.interp import ConservativeRegridder, Grid

        source_grid = Grid.from_degrees(
            np.linspace(0, 360, 20, endpoint=False), np.linspace(-60, 60, 15)
        )
        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 10, endpoint=False), np.linspace(-30, 30, 8)
        )

        regridder = ConservativeRegridder(source_grid, target_grid)

        # Initially, weights should not be computed
        assert regridder._lon_weights is None
        assert regridder._lat_weights is None

        # First data array regridding should compute weights
        data1 = np.random.randn(*source_grid.shape)
        result1 = regridder.regrid_array(data1)

        # Now weights should be cached
        assert regridder._lon_weights is not None
        assert regridder._lat_weights is not None

        # Store references to cached weights
        cached_lon_weights = regridder._lon_weights
        cached_lat_weights = regridder._lat_weights

        # Second regridding should reuse cached weights
        data2 = np.random.randn(*source_grid.shape)
        result2 = regridder.regrid_array(data2)

        # Should be the same objects (not recomputed)
        assert regridder._lon_weights is cached_lon_weights
        assert regridder._lat_weights is cached_lat_weights

        # Results should have correct shapes
        assert result1.shape == target_grid.shape
        assert result2.shape == target_grid.shape

    def test_dataset_regridding_efficiency(self):
        """Test that dataset regridding is efficient for multiple variables."""
        from skyborn.interp import Grid, regrid_dataset

        # Create dataset with multiple variables
        lat_vals = np.linspace(-45, 45, 20)
        lon_vals = np.linspace(0, 360, 30, endpoint=False)

        dataset = xr.Dataset(
            {
                "temperature": (["lat", "lon"], np.random.randn(20, 30)),
                "humidity": (["lat", "lon"], np.random.randn(20, 30)),
                "pressure": (["lat", "lon"], np.random.randn(20, 30)),
            },
            coords={"lat": lat_vals, "lon": lon_vals},
        )

        target_grid = Grid.from_degrees(
            np.linspace(0, 360, 15, endpoint=False), np.linspace(-30, 30, 10)
        )

        # Regrid all variables at once
        result = regrid_dataset(dataset, target_grid, method="conservative")

        # All variables should be regridded to same target shape
        for var_name in ["temperature", "humidity", "pressure"]:
            assert result[var_name].shape == (10, 15)

        # Should preserve variable attributes if they exist
        if "temperature" in dataset and hasattr(dataset["temperature"], "attrs"):
            for attr_name in dataset["temperature"].attrs:
                if attr_name in result["temperature"].attrs:
                    assert (
                        result["temperature"].attrs[attr_name]
                        == dataset["temperature"].attrs[attr_name]
                    )


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
