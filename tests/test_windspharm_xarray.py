"""
Tests for skyborn.windspharm.xarray module.

This module tests the xarray-specific functionality of VectorWind,
including coordinate preservation, metadata handling, and xarray DataArray operations.
"""

import numpy as np
import pytest
import warnings

try:
    import xarray as xr

    XARRAY_AVAILABLE = True
except ImportError:
    XARRAY_AVAILABLE = False
    xr = None

if XARRAY_AVAILABLE:
    from skyborn.windspharm.xarray import VectorWind


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayVectorWindInitialization:
    """Test VectorWind xarray initialization."""

    def test_xarray_vectorwind_basic(self):
        """Test basic xarray VectorWind initialization."""
        # Create xarray DataArrays with coordinates
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(19, 36)
        v_data = np.random.randn(19, 36)

        u = xr.DataArray(
            u_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "long_name": "zonal wind"},
        )
        v = xr.DataArray(
            v_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "long_name": "meridional wind"},
        )

        vw = VectorWind(u, v)
        assert vw is not None
        assert hasattr(vw, "_api")

    def test_xarray_vectorwind_3d(self):
        """Test VectorWind with 3D xarray data."""
        time = np.arange(5)
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)

        u_data = np.random.randn(5, 19, 36)
        v_data = np.random.randn(5, 19, 36)

        u = xr.DataArray(
            u_data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )
        v = xr.DataArray(
            v_data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )

        vw = VectorWind(u, v)
        assert vw is not None

    def test_xarray_vectorwind_custom_params(self):
        """Test VectorWind with custom parameters."""
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(19, 36)
        v_data = np.random.randn(19, 36)

        u = xr.DataArray(u_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})
        v = xr.DataArray(v_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})

        vw = VectorWind(u, v, rsphere=6.37e6, legfunc="computed")
        assert vw is not None


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayVectorWindOperations:
    """Test VectorWind xarray operations."""

    @pytest.fixture
    def sample_xarray_data(self):
        """Create sample xarray wind data."""
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(19, 36) * 10
        v_data = np.random.randn(19, 36) * 10

        u = xr.DataArray(
            u_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "long_name": "zonal wind"},
        )
        v = xr.DataArray(
            v_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "long_name": "meridional wind"},
        )
        return u, v

    def test_u_v_accessors(self, sample_xarray_data):
        """Test u and v component accessors."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        u_out = vw.u()
        v_out = vw.v()

        assert isinstance(u_out, xr.DataArray)
        assert isinstance(v_out, xr.DataArray)
        assert "standard_name" in u_out.attrs
        assert "standard_name" in v_out.attrs
        assert u_out.attrs["standard_name"] == "eastward_wind"
        assert v_out.attrs["standard_name"] == "northward_wind"

    def test_magnitude(self, sample_xarray_data):
        """Test wind magnitude calculation."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        magnitude = vw.magnitude()

        assert isinstance(magnitude, xr.DataArray)
        assert magnitude.shape == u.shape
        assert "standard_name" in magnitude.attrs
        assert magnitude.attrs["standard_name"] == "wind_speed"
        assert magnitude.attrs["units"] == "m s**-1"

    def test_vrtdiv(self, sample_xarray_data):
        """Test combined vorticity and divergence calculation."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        vrt, div = vw.vrtdiv()

        assert isinstance(vrt, xr.DataArray)
        assert isinstance(div, xr.DataArray)
        assert vrt.shape == u.shape
        assert div.shape == u.shape
        assert "standard_name" in vrt.attrs
        assert "standard_name" in div.attrs

    def test_vorticity_xarray(self, sample_xarray_data):
        """Test vorticity calculation with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        vorticity = vw.vorticity()

        # Check it's an xarray DataArray
        assert isinstance(vorticity, xr.DataArray)
        # Check coordinates are preserved
        assert "lat" in vorticity.coords
        assert "lon" in vorticity.coords
        # Check shape
        assert vorticity.shape == u.shape

    def test_divergence_xarray(self, sample_xarray_data):
        """Test divergence calculation with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        divergence = vw.divergence()

        assert isinstance(divergence, xr.DataArray)
        assert "lat" in divergence.coords
        assert "lon" in divergence.coords
        assert divergence.shape == u.shape

    def test_rossby_wave_source_xarray(self, sample_xarray_data):
        """Test Rossby wave source calculation with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        # Test basic RWS calculation
        rws = vw.rossbywavesource()

        assert isinstance(rws, xr.DataArray)
        assert "lat" in rws.coords
        assert "lon" in rws.coords
        assert rws.shape == u.shape
        assert np.all(np.isfinite(rws.values))

        # Check units and metadata
        assert "units" in rws.attrs
        assert rws.attrs["units"] == "s**-2"
        assert "standard_name" in rws.attrs
        assert rws.attrs["standard_name"] == "rossby_wave_source"

        # Test with truncation (use smaller truncation for test grid)
        rws_t15 = vw.rossbywavesource(truncation=15)
        assert isinstance(rws_t15, xr.DataArray)
        assert rws_t15.shape == u.shape
        assert np.all(np.isfinite(rws_t15.values))

        # Test with custom omega (significantly different from default 7.292e-5)
        rws_custom = vw.rossbywavesource(omega=1.0e-4)
        assert isinstance(rws_custom, xr.DataArray)
        assert rws_custom.shape == u.shape
        assert np.all(np.isfinite(rws_custom.values))

        # Basic functionality test passed - omega parameter is correctly accepted

    def test_streamfunction_xarray(self, sample_xarray_data):
        """Test streamfunction calculation with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        sf = vw.streamfunction()

        assert isinstance(sf, xr.DataArray)
        assert "lat" in sf.coords
        assert "lon" in sf.coords
        assert sf.shape == u.shape

    def test_velocity_potential_xarray(self, sample_xarray_data):
        """Test velocity potential calculation with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        vp = vw.velocitypotential()

        assert isinstance(vp, xr.DataArray)
        assert "lat" in vp.coords
        assert "lon" in vp.coords
        assert vp.shape == u.shape

    def test_sfvp_xarray(self, sample_xarray_data):
        """Test combined streamfunction and velocity potential."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        sf, vp = vw.sfvp()

        assert isinstance(sf, xr.DataArray)
        assert isinstance(vp, xr.DataArray)
        assert sf.shape == u.shape
        assert vp.shape == u.shape

    def test_helmholtz_xarray(self, sample_xarray_data):
        """Test Helmholtz decomposition with xarray."""
        u, v = sample_xarray_data
        vw = VectorWind(u, v)

        u_chi, v_chi, u_psi, v_psi = vw.helmholtz()

        for component in [u_chi, v_chi, u_psi, v_psi]:
            assert isinstance(component, xr.DataArray)
            assert component.shape == u.shape
            assert "lat" in component.coords
            assert "lon" in component.coords


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayMetadataPreservation:
    """Test metadata and coordinate preservation."""

    def test_coordinate_preservation(self):
        """Test that coordinates are preserved in outputs."""
        # Create data with specific coordinate attributes
        lat = xr.DataArray(
            np.linspace(-90, 90, 19),
            dims=["lat"],
            attrs={"units": "degrees_north", "long_name": "latitude"},
        )
        lon = xr.DataArray(
            np.linspace(0, 357.5, 36),
            dims=["lon"],
            attrs={"units": "degrees_east", "long_name": "longitude"},
        )

        u_data = np.random.randn(19, 36)
        v_data = np.random.randn(19, 36)

        u = xr.DataArray(
            u_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "standard_name": "eastward_wind"},
        )
        v = xr.DataArray(
            v_data,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={"units": "m/s", "standard_name": "northward_wind"},
        )

        vw = VectorWind(u, v)
        vorticity = vw.vorticity()

        # Check coordinate attributes are preserved
        assert vorticity.coords["lat"].attrs["units"] == "degrees_north"
        assert vorticity.coords["lon"].attrs["units"] == "degrees_east"

    def test_3d_coordinate_preservation(self):
        """Test coordinate preservation with 3D data."""
        time = xr.DataArray(
            np.arange(3), dims=["time"], attrs={"units": "days since 2020-01-01"}
        )
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)

        u_data = np.random.randn(3, 19, 36)
        v_data = np.random.randn(3, 19, 36)

        u = xr.DataArray(
            u_data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )
        v = xr.DataArray(
            v_data,
            dims=["time", "lat", "lon"],
            coords={"time": time, "lat": lat, "lon": lon},
        )

        vw = VectorWind(u, v)
        vorticity = vw.vorticity()

        # Check all coordinates preserved
        assert "time" in vorticity.coords
        assert "lat" in vorticity.coords
        assert "lon" in vorticity.coords
        assert vorticity.shape == (3, 19, 36)


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayValidation:
    """Test validation and error handling for xarray inputs."""

    def test_mismatched_coordinates(self):
        """Test error with mismatched coordinates."""
        lat1 = np.linspace(-90, 90, 19)
        lat2 = np.linspace(-85, 85, 19)  # Different range
        lon = np.linspace(0, 357.5, 36)

        u = xr.DataArray(
            np.random.randn(19, 36),
            dims=["lat", "lon"],
            coords={"lat": lat1, "lon": lon},
        )
        v = xr.DataArray(
            np.random.randn(19, 36),
            dims=["lat", "lon"],
            coords={"lat": lat2, "lon": lon},
        )

        # Should raise an error for mismatched coordinates
        with pytest.raises(ValueError, match="identical coordinate values"):
            VectorWind(u, v)

    def test_different_dimensions(self):
        """Test error with different dimension names."""
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)

        u = xr.DataArray(
            np.random.randn(19, 36),
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
        )
        v = xr.DataArray(
            np.random.randn(19, 36),
            dims=["y", "x"],  # Different dimension names
            coords={"y": lat, "x": lon},
        )

        # Should raise an error for different dimensions
        with pytest.raises(ValueError, match="identical dimensions"):
            VectorWind(u, v)


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayAdvancedOperations:
    """Test advanced xarray operations."""

    @pytest.fixture
    def complex_xarray_data(self):
        """Create complex xarray data with multiple dimensions."""
        time = np.arange(4)
        level = np.array([1000, 850, 500, 200])
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)

        u_data = np.random.randn(4, 4, 19, 36) * 20
        v_data = np.random.randn(4, 4, 19, 36) * 20

        u = xr.DataArray(
            u_data,
            dims=["time", "level", "lat", "lon"],
            coords={"time": time, "level": level, "lat": lat, "lon": lon},
            attrs={"units": "m/s"},
        )
        v = xr.DataArray(
            v_data,
            dims=["time", "level", "lat", "lon"],
            coords={"time": time, "level": level, "lat": lat, "lon": lon},
            attrs={"units": "m/s"},
        )
        return u, v

    def test_4d_operations(self, complex_xarray_data):
        """Test operations with 4D xarray data."""
        u, v = complex_xarray_data
        vw = VectorWind(u, v)

        vorticity = vw.vorticity()

        assert isinstance(vorticity, xr.DataArray)
        assert vorticity.shape == u.shape
        assert "time" in vorticity.coords
        assert "level" in vorticity.coords
        assert "lat" in vorticity.coords
        assert "lon" in vorticity.coords

    def test_truncated_outputs(self, complex_xarray_data):
        """Test operations with truncation parameter."""
        u, v = complex_xarray_data
        vw = VectorWind(u, v)

        # Test with valid truncation (must be <= nlat-1 = 18)
        vorticity_t15 = vw.vorticity(truncation=15)
        vorticity_full = vw.vorticity()

        assert isinstance(vorticity_t15, xr.DataArray)
        assert vorticity_t15.shape == vorticity_full.shape

        # Test that invalid truncation raises error
        with pytest.raises(Exception):  # Could be ValidationError or similar
            vw.vorticity(truncation=21)  # Too high for 19 latitudes

    def test_gradient_operations(self, complex_xarray_data):
        """Test gradient operations."""
        u, v = complex_xarray_data
        vw = VectorWind(u, v)

        # Get a scalar field first
        vorticity = vw.vorticity()

        # Test gradients
        grad_u, grad_v = vw.gradient(vorticity)

        assert isinstance(grad_u, xr.DataArray)
        assert isinstance(grad_v, xr.DataArray)
        assert grad_u.shape == u.shape
        assert grad_v.shape == v.shape


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayErrorHandling:
    """Test error handling specific to xarray interface."""

    def test_invalid_input_types(self):
        """Test error handling for invalid input types."""
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)

        # Valid xarray
        u = xr.DataArray(
            np.random.randn(19, 36),
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
        )

        # Invalid - numpy array instead of xarray
        v = np.random.randn(19, 36)

        with pytest.raises(TypeError, match="must be xarray.DataArray"):
            VectorWind(u, v)

    def test_missing_coordinates(self):
        """Test handling of data without proper coordinates."""
        u_data = np.random.randn(19, 36)
        v_data = np.random.randn(19, 36)

        # Create DataArrays without proper lat/lon coordinates
        u = xr.DataArray(u_data, dims=["y", "x"])
        v = xr.DataArray(v_data, dims=["y", "x"])

        # Should handle gracefully or provide meaningful error
        try:
            vw = VectorWind(u, v)
            # If it works, basic operations should also work
            result = vw.vorticity()
            assert isinstance(result, xr.DataArray)
        except Exception as e:
            # Should be a meaningful error message
            assert len(str(e)) > 0

    def test_xarray_fallback_import(self):
        """Test xarray import fallback behavior."""
        # This is harder to test directly due to import caching
        # We can at least verify the module has the proper error handling
        import sys

        original_modules = sys.modules.copy()

        try:
            # Test that the module handles imports properly
            from skyborn.windspharm import xarray as xa_module

            # If we got here, xarray is available
            assert hasattr(xa_module, "VectorWind")
        except ImportError:
            # If xarray isn't available, that's expected
            pass
        finally:
            # Restore modules
            sys.modules.clear()
            sys.modules.update(original_modules)


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarrayDifferentGridTypes:
    """Test different grid types and coordinate systems."""

    def test_gaussian_grid_inference(self):
        """Test automatic detection of Gaussian grids."""
        # Use a proper Gaussian grid from the library
        try:
            from skyborn.spharm import gaussian_lats_wts

            lat_gauss, _ = gaussian_lats_wts(18)
        except ImportError:
            # Skip if gaussian_lats_wts not available
            pytest.skip("gaussian_lats_wts not available")

        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(18, 36)
        v_data = np.random.randn(18, 36)

        u = xr.DataArray(
            u_data, dims=["lat", "lon"], coords={"lat": lat_gauss, "lon": lon}
        )
        v = xr.DataArray(
            v_data, dims=["lat", "lon"], coords={"lat": lat_gauss, "lon": lon}
        )

        vw = VectorWind(u, v)
        vorticity = vw.vorticity()
        assert isinstance(vorticity, xr.DataArray)

    def test_south_to_north_latitude_ordering(self):
        """Test automatic reversal of south-to-north latitude ordering."""
        # Create south-to-north latitudes (wrong order)
        lat = np.linspace(-90, 90, 19)  # South to north
        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(19, 36)
        v_data = np.random.randn(19, 36)

        u = xr.DataArray(u_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})
        v = xr.DataArray(v_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})

        vw = VectorWind(u, v)
        # Should automatically reverse and work
        vorticity = vw.vorticity()
        assert isinstance(vorticity, xr.DataArray)
        # Output should have proper north-to-south ordering
        assert vorticity.coords["lat"].values[0] > vorticity.coords["lat"].values[-1]


@pytest.mark.skipif(not XARRAY_AVAILABLE, reason="xarray not available")
class TestXarraySpecialOperations:
    """Test special operations and edge cases."""

    @pytest.fixture
    def special_data(self):
        """Create data for special tests."""
        lat = np.linspace(-90, 90, 19)
        lon = np.linspace(0, 357.5, 36)
        u_data = np.random.randn(19, 36) * 10
        v_data = np.random.randn(19, 36) * 10

        u = xr.DataArray(u_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})
        v = xr.DataArray(v_data, dims=["lat", "lon"], coords={"lat": lat, "lon": lon})
        return u, v

    def test_irrotationalcomponent(self, special_data):
        """Test irrotational component calculation."""
        u, v = special_data
        vw = VectorWind(u, v)

        u_irrot, v_irrot = vw.irrotationalcomponent()

        assert isinstance(u_irrot, xr.DataArray)
        assert isinstance(v_irrot, xr.DataArray)
        assert u_irrot.shape == u.shape
        assert v_irrot.shape == u.shape

    def test_nondivergentcomponent(self, special_data):
        """Test nondivergent component calculation."""
        u, v = special_data
        vw = VectorWind(u, v)

        u_nondiv, v_nondiv = vw.nondivergentcomponent()

        assert isinstance(u_nondiv, xr.DataArray)
        assert isinstance(v_nondiv, xr.DataArray)
        assert u_nondiv.shape == u.shape
        assert v_nondiv.shape == u.shape

    def test_planetaryvorticity(self, special_data):
        """Test planetary vorticity calculation."""
        u, v = special_data
        vw = VectorWind(u, v)

        f = vw.planetaryvorticity()

        assert isinstance(f, xr.DataArray)
        assert f.shape == u.shape
        # Check that it has proper attributes
        assert "standard_name" in f.attrs

    def test_absolutevorticity(self, special_data):
        """Test absolute vorticity calculation."""
        u, v = special_data
        vw = VectorWind(u, v)

        abs_vort = vw.absolutevorticity()

        assert isinstance(abs_vort, xr.DataArray)
        assert abs_vort.shape == u.shape
        # Check that it has proper attributes
        assert "standard_name" in abs_vort.attrs
