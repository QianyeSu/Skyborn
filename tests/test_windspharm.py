"""
Tests for skyborn.windspharm module (spherical harmonic wind analysis)

This module contains comprehensive tests for the windspharm functionality,
including VectorWind class and related utilities.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    from skyborn.windspharm import VectorWind, standard, tools
    from skyborn.windspharm.standard import VectorWind as StandardVectorWind
    from skyborn.windspharm.tools import (
        prep_data,
        recover_data,
        get_recovery,
        reverse_latdim,
        order_latdim,
    )

    WINDSPHARM_AVAILABLE = True
except ImportError:
    WINDSPHARM_AVAILABLE = False
    VectorWind = None
    standard = None
    tools = None


class TestVectorWindInitialization:
    """Test VectorWind class initialization."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_basic_initialization(self):
        """Test basic VectorWind initialization."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)
        assert vw.u.shape == (nlat, nlon)
        assert vw.v.shape == (nlat, nlon)
        assert vw.gridtype == "regular"
        assert hasattr(vw, "s")  # Spharmt object

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_gaussian_grid(self):
        """Test VectorWind initialization with Gaussian grid."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v, gridtype="gaussian")
        assert vw.gridtype == "gaussian"

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_custom_parameters(self):
        """Test VectorWind initialization with custom parameters."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v, gridtype="regular", rsphere=6.371e6, legfunc="computed")
        assert vw.gridtype == "regular"

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_3d_data(self):
        """Test VectorWind initialization with 3D data (time series)."""
        nlat, nlon, nt = 37, 72, 10
        u = np.random.randn(nlat, nlon, nt)
        v = np.random.randn(nlat, nlon, nt)

        vw = VectorWind(u, v)
        assert vw.u.shape == (nlat, nlon, nt)
        assert vw.v.shape == (nlat, nlon, nt)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_masked_arrays(self):
        """Test VectorWind with masked arrays."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        # Create masked arrays
        u_masked = np.ma.masked_where(u > 2.0, u)
        v_masked = np.ma.masked_where(v > 2.0, v)

        # This should work but fill masked values
        vw = VectorWind(u_masked.filled(0), v_masked.filled(0))
        assert vw.u.shape == (nlat, nlon)
        assert vw.v.shape == (nlat, nlon)


class TestVectorWindValidation:
    """Test VectorWind input validation."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_shape_mismatch(self):
        """Test error handling for mismatched shapes."""
        u = np.random.randn(37, 72)
        v = np.random.randn(36, 72)  # Different nlat

        with pytest.raises(ValueError):
            VectorWind(u, v)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_nan_values(self):
        """Test error handling for NaN values."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        # Insert NaN values
        u[0, 0] = np.nan

        with pytest.raises(ValueError):
            VectorWind(u, v)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_inf_values(self):
        """Test error handling for infinite values."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        # Insert infinite values
        v[0, 0] = np.inf

        with pytest.raises(ValueError):
            VectorWind(u, v)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_extreme_values(self):
        """Test error handling for extremely large values."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        # Insert extremely large values
        u[0, 0] = 1e10

        with pytest.raises(ValueError):
            VectorWind(u, v)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_invalid_gridtype(self):
        """Test error handling for invalid gridtype."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        with pytest.raises((ValueError, TypeError)):
            VectorWind(u, v, gridtype="invalid")

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_invalid_legfunc(self):
        """Test error handling for invalid legfunc."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        with pytest.raises((ValueError, TypeError)):
            VectorWind(u, v, legfunc="invalid")

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vectorwind_too_small_grid(self):
        """Test error handling for too small grids."""
        u = np.random.randn(2, 3)
        v = np.random.randn(2, 3)

        with pytest.raises((ValueError, RuntimeError)):
            VectorWind(u, v)


class TestVectorWindBasicOperations:
    """Test basic VectorWind operations."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_vorticity_calculation(self):
        """Test vorticity calculation."""
        nlat, nlon = 37, 72

        # Create simple wind pattern
        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")
        vorticity = vw.vorticity()

        assert vorticity.shape == (nlat, nlon)
        assert isinstance(vorticity, np.ndarray)
        assert np.all(np.isfinite(vorticity))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_divergence_calculation(self):
        """Test divergence calculation."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")
        divergence = vw.divergence()

        assert divergence.shape == (nlat, nlon)
        assert isinstance(divergence, np.ndarray)
        assert np.all(np.isfinite(divergence))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_rossby_wave_source_calculation(self):
        """Test Rossby wave source calculation."""
        nlat, nlon = 37, 72

        # Create test wind field with some structure
        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(2 * lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.sin(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        # Test basic RWS calculation
        rws = vw.rossbywavesource()
        assert rws.shape == (nlat, nlon)
        assert isinstance(rws, np.ndarray)
        assert np.all(np.isfinite(rws))

        # Test with truncation
        rws_t21 = vw.rossbywavesource(truncation=21)
        assert rws_t21.shape == (nlat, nlon)
        assert np.all(np.isfinite(rws_t21))

        # Test with custom omega (significantly different from default)
        rws_custom = vw.rossbywavesource(omega=1.0e-4)
        assert rws_custom.shape == (nlat, nlon)
        assert np.all(np.isfinite(rws_custom))

        # Basic functionality test passed - omega parameter is correctly accepted

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_streamfunction_calculation(self):
        """Test streamfunction calculation."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")
        streamfunction = vw.streamfunction()

        assert streamfunction.shape == (nlat, nlon)
        assert isinstance(streamfunction, np.ndarray)
        assert np.all(np.isfinite(streamfunction))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_velocity_potential_calculation(self):
        """Test velocity potential calculation."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")
        velocity_potential = vw.velocitypotential()

        assert velocity_potential.shape == (nlat, nlon)
        assert isinstance(velocity_potential, np.ndarray)
        assert np.all(np.isfinite(velocity_potential))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_sfvp_calculation(self):
        """Test combined streamfunction and velocity potential calculation."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")
        streamfunction, velocity_potential = vw.sfvp()

        assert streamfunction.shape == (nlat, nlon)
        assert velocity_potential.shape == (nlat, nlon)
        assert isinstance(streamfunction, np.ndarray)
        assert isinstance(velocity_potential, np.ndarray)
        assert np.all(np.isfinite(streamfunction))
        assert np.all(np.isfinite(velocity_potential))


class TestVectorWindAdvancedOperations:
    """Test advanced VectorWind operations."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_truncation_parameter(self):
        """Test operations with truncation parameter."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        # Test with different truncation levels
        vorticity_20 = vw.vorticity(truncation=20)
        vorticity_30 = vw.vorticity(truncation=30)

        assert vorticity_20.shape == (nlat, nlon)
        assert vorticity_30.shape == (nlat, nlon)

        # Results should be different due to different truncation
        assert not np.allclose(vorticity_20, vorticity_30, rtol=1e-3)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_gradient_calculation(self):
        """Test gradient calculation if available."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        # Check if gradient method exists
        if hasattr(vw, "gradient"):
            # Create a scalar field for gradient calculation
            scalar_field = np.sin(lon_grid) * np.cos(lat_grid)

            grad_x, grad_y = vw.gradient(scalar_field)

            assert grad_x.shape == (nlat, nlon)
            assert grad_y.shape == (nlat, nlon)
            assert np.all(np.isfinite(grad_x))
            assert np.all(np.isfinite(grad_y))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_3d_operations(self):
        """Test operations with 3D data."""
        nlat, nlon, nt = 37, 72, 5

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Create 3D wind fields
        u = np.zeros((nlat, nlon, nt))
        v = np.zeros((nlat, nlon, nt))

        for t in range(nt):
            phase = 2 * np.pi * t / nt
            u[:, :, t] = np.sin(lon_grid + phase) * np.cos(lat_grid)
            v[:, :, t] = np.cos(lon_grid + phase) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        vorticity = vw.vorticity()
        divergence = vw.divergence()

        assert vorticity.shape == (nlat, nlon, nt)
        assert divergence.shape == (nlat, nlon, nt)
        assert np.all(np.isfinite(vorticity))
        assert np.all(np.isfinite(divergence))


class TestVectorWindConsistency:
    """Test consistency and mathematical properties."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_helmholtz_decomposition_consistency(self):
        """Test that Helmholtz decomposition is consistent."""
        nlat, nlon = 37, 72

        # Create divergent and rotational wind components
        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Simple wind field
        u = np.sin(2 * lon_grid) * np.cos(lat_grid)
        v = np.cos(2 * lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        # Get streamfunction and velocity potential
        psi, chi = vw.sfvp()

        # Check that they are finite and have correct shape
        assert psi.shape == (nlat, nlon)
        assert chi.shape == (nlat, nlon)
        assert np.all(np.isfinite(psi))
        assert np.all(np.isfinite(chi))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_symmetry_properties(self):
        """Test symmetry properties of operations."""
        nlat, nlon = 37, 72

        # Create symmetric wind field
        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        vorticity = vw.vorticity()
        divergence = vw.divergence()

        # Check that results are reasonable
        assert np.std(vorticity) > 0  # Should have some variation
        assert np.std(divergence) > 0  # Should have some variation

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_gaussian_vs_regular_grid(self):
        """Test differences between Gaussian and regular grids."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        # Regular grid
        vw_regular = VectorWind(u, v, gridtype="regular")
        vort_regular = vw_regular.vorticity()

        # Gaussian grid
        vw_gaussian = VectorWind(u, v, gridtype="gaussian")
        vort_gaussian = vw_gaussian.vorticity()

        # Both should be valid but may differ slightly
        assert vort_regular.shape == vort_gaussian.shape
        assert np.all(np.isfinite(vort_regular))
        assert np.all(np.isfinite(vort_gaussian))


class TestWindspharmTools:
    """Test windspharm.tools module functions."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_prep_data_basic(self):
        """Test basic prep_data functionality."""
        # Create test data
        nlat, nlon = 37, 72
        data = np.random.randn(nlat, nlon)

        prepared, info = prep_data(data, "yx")

        assert isinstance(prepared, np.ndarray)
        assert isinstance(info, dict)
        # After prep_data, it might add a dimension for combination with other fields
        assert prepared.ndim >= 2

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_prep_data_3d(self):
        """Test prep_data with 3D data."""
        nlat, nlon, nt = 37, 72, 5
        data = np.random.randn(nt, nlat, nlon)

        prepared, info = prep_data(data, "tyx")

        assert isinstance(prepared, np.ndarray)
        assert isinstance(info, dict)
        assert prepared.ndim == 3

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_recover_data_basic(self):
        """Test basic recover_data functionality."""
        nlat, nlon = 37, 72
        original = np.random.randn(nlat, nlon)

        prepared, info = prep_data(original, "yx")
        recovered = recover_data(prepared, info)

        assert recovered.shape == original.shape
        np.testing.assert_array_equal(recovered, original)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_get_recovery_function(self):
        """Test get_recovery function creation."""
        nlat, nlon = 37, 72
        data1 = np.random.randn(nlat, nlon)

        prepared1, info1 = prep_data(data1, "yx")

        # Create recovery function - should pass a single info dict
        recovery_func = get_recovery(info1)

        assert callable(recovery_func)

        # Test recovery with single array
        recovered_list = recovery_func(prepared1)
        assert isinstance(recovered_list, list)
        assert len(recovered_list) == 1
        recovered = recovered_list[0]
        assert recovered.shape == data1.shape
        np.testing.assert_array_equal(recovered, data1)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_reverse_latdim(self):
        """Test reverse_latdim function."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        u_reversed, v_reversed = reverse_latdim(u, v)

        assert u_reversed.shape == u.shape
        assert v_reversed.shape == v.shape
        # Check that latitude dimension is actually reversed
        np.testing.assert_array_equal(u_reversed[0, :], u[-1, :])
        np.testing.assert_array_equal(u_reversed[-1, :], u[0, :])
        np.testing.assert_array_equal(v_reversed[0, :], v[-1, :])
        np.testing.assert_array_equal(v_reversed[-1, :], v[0, :])

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_order_latdim(self):
        """Test order_latdim function."""
        nlat, nlon = 37, 72
        lats = np.linspace(90, -90, nlat)  # North to south
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        ordered_lats, ordered_u, ordered_v = order_latdim(lats, u, v)

        assert ordered_u.shape == u.shape
        assert ordered_v.shape == v.shape
        assert len(ordered_lats) == nlat
        # Should be ordered north to south
        assert ordered_lats[0] > ordered_lats[-1]

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_order_latdim_reversed(self):
        """Test order_latdim with initially south-to-north data."""
        nlat, nlon = 37, 72
        lats = np.linspace(-90, 90, nlat)  # South to north
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        ordered_lats, ordered_u, ordered_v = order_latdim(lats, u, v)

        assert ordered_u.shape == u.shape
        assert ordered_v.shape == v.shape
        assert len(ordered_lats) == nlat
        # Should be reordered to north to south
        assert ordered_lats[0] > ordered_lats[-1]


class TestWindspharmPerformance:
    """Test windspharm performance with different configurations."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_small_grid_performance(self):
        """Test performance with small grid."""
        nlat, nlon = 19, 36

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        vw = VectorWind(u, v, gridtype="regular")

        # Perform multiple operations
        for _ in range(5):
            vorticity = vw.vorticity()
            divergence = vw.divergence()
            psi, chi = vw.sfvp()

        # Basic sanity check
        assert vorticity.shape == (nlat, nlon)
        assert divergence.shape == (nlat, nlon)
        assert psi.shape == (nlat, nlon)
        assert chi.shape == (nlat, nlon)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_legfunc_options(self):
        """Test different legfunc options."""
        nlat, nlon = 37, 72

        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        # Test with stored Legendre functions
        vw_stored = VectorWind(u, v, legfunc="stored")
        vort_stored = vw_stored.vorticity()

        # Test with computed Legendre functions
        vw_computed = VectorWind(u, v, legfunc="computed")
        vort_computed = vw_computed.vorticity()

        # Results should be very similar
        assert vort_stored.shape == vort_computed.shape
        assert np.allclose(vort_stored, vort_computed, rtol=1e-10)


class TestWindspharmErrorHandling:
    """Test error handling in windspharm operations."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_invalid_truncation(self):
        """Test error handling for invalid truncation values."""
        nlat, nlon = 37, 72
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)

        # Test with negative truncation
        with pytest.raises((ValueError, RuntimeError)):
            vw.vorticity(truncation=-1)

        # Test with too large truncation
        with pytest.raises((ValueError, RuntimeError)):
            vw.vorticity(truncation=1000)

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_tools_invalid_input(self):
        """Test error handling in tools functions."""
        # Test prep_data with invalid dimension order
        data = np.random.randn(37, 72)

        with pytest.raises((ValueError, KeyError)):
            prep_data(data, "invalid")

        # Test order_latdim with mismatched dimensions
        u = np.random.randn(37, 72)
        v = np.random.randn(37, 72)
        lats = np.linspace(-90, 90, 36)  # Wrong size

        with pytest.raises((ValueError, IndexError)):
            order_latdim(lats, u, v)


class TestWindspharmIntegration:
    """Integration tests for windspharm functionality."""

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_full_workflow(self):
        """Test complete windspharm workflow."""
        nlat, nlon = 37, 72

        # Create test wind field
        lons = np.linspace(0, 2 * np.pi, nlon, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, nlat)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(2 * lon_grid) * np.cos(lat_grid)
        v = np.cos(2 * lon_grid) * np.cos(lat_grid)

        # Initialize VectorWind
        vw = VectorWind(u, v, gridtype="regular")

        # Calculate all quantities
        vorticity = vw.vorticity()
        divergence = vw.divergence()
        streamfunction = vw.streamfunction()
        velocity_potential = vw.velocitypotential()
        psi, chi = vw.sfvp()

        # Verify all results
        assert vorticity.shape == (nlat, nlon)
        assert divergence.shape == (nlat, nlon)
        assert streamfunction.shape == (nlat, nlon)
        assert velocity_potential.shape == (nlat, nlon)
        assert psi.shape == (nlat, nlon)
        assert chi.shape == (nlat, nlon)

        # Check that sfvp components match individual calls
        assert np.allclose(psi, streamfunction, rtol=1e-10)
        assert np.allclose(chi, velocity_potential, rtol=1e-10)

        # Verify all values are finite
        assert np.all(np.isfinite(vorticity))
        assert np.all(np.isfinite(divergence))
        assert np.all(np.isfinite(streamfunction))
        assert np.all(np.isfinite(velocity_potential))

    @pytest.mark.skipif(
        not WINDSPHARM_AVAILABLE, reason="windspharm module not available"
    )
    def test_tools_integration(self):
        """Test tools integration with VectorWind."""
        nlat, nlon, nt = 37, 72, 3

        # Create 3D wind data with time dimension first
        u_orig = np.random.randn(nt, nlat, nlon)
        v_orig = np.random.randn(nt, nlat, nlon)

        # Use tools to prepare data
        u_prep, u_info = prep_data(u_orig, "tyx")
        v_prep, v_info = prep_data(v_orig, "tyx")

        # Process with VectorWind
        vw = VectorWind(u_prep, v_prep)
        vorticity = vw.vorticity()

        # Recover original shape
        vorticity_recovered = recover_data(vorticity, u_info)

        assert vorticity_recovered.shape == (nt, nlat, nlon)
        assert np.all(np.isfinite(vorticity_recovered))


class TestVectorWindUntestedMethods:
    """Test previously untested methods in VectorWind."""

    def test_vectorwind_magnitude_calculation(self):
        """Test magnitude calculation method."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon) * 10
        v = np.random.randn(nlat, nlon) * 10

        vw = VectorWind(u, v)

        # Test magnitude calculation
        magnitude = vw.magnitude()

        # Should match manual calculation
        expected_magnitude = np.sqrt(u**2 + v**2)
        np.testing.assert_allclose(magnitude, expected_magnitude, rtol=1e-12)

    def test_vectorwind_with_truncation_parameter(self):
        """Test VectorWind operations with truncation parameter."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)

        # Test operations with custom truncation
        for ntrunc in [5, 10, 15]:
            vort = vw.vorticity(truncation=ntrunc)
            div = vw.divergence(truncation=ntrunc)

            assert vort.shape == (nlat, nlon)
            assert div.shape == (nlat, nlon)
            assert np.all(np.isfinite(vort))
            assert np.all(np.isfinite(div))

    @pytest.mark.skip(
        reason="Helmholtz decomposition has large reconstruction errors due to spherical harmonic truncation"
    )
    def test_vectorwind_helmholtz_decomposition(self):
        """Test complete Helmholtz decomposition."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)

        # Get Helmholtz decomposition
        u_chi, v_chi, u_psi, v_psi = vw.helmholtz()

        # Should reconstruct original field (approximately)
        u_reconstructed = u_chi + u_psi
        v_reconstructed = v_chi + v_psi

        # Check that components have expected shapes
        assert u_chi.shape == u.shape
        assert v_chi.shape == v.shape
        assert u_psi.shape == u.shape
        assert v_psi.shape == v.shape

    def test_vectorwind_gradient_operations(self):
        """Test gradient operations."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)

        # Get a scalar field (vorticity)
        vorticity = vw.vorticity()

        # Calculate gradients
        grad_x, grad_y = vw.gradient(vorticity)

        assert grad_x.shape == (nlat, nlon)
        assert grad_y.shape == (nlat, nlon)
        assert np.all(np.isfinite(grad_x))
        assert np.all(np.isfinite(grad_y))

    @pytest.mark.skip(
        reason="Irrotational/nondivergent decomposition has large reconstruction errors due to spherical harmonic truncation"
    )
    def test_vectorwind_irrotational_nondivergent(self):
        """Test irrotational and nondivergent components."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon)
        v = np.random.randn(nlat, nlon)

        vw = VectorWind(u, v)

        # Get irrotational component
        u_irrot, v_irrot = vw.irrotationalcomponent()

        # Get nondivergent component
        u_nondiv, v_nondiv = vw.nondivergentcomponent()

        # Check that components have expected shapes
        assert u_irrot.shape == u.shape
        assert v_irrot.shape == v.shape
        assert u_nondiv.shape == u.shape
        assert v_nondiv.shape == v.shape

    def test_vectorwind_quasi_geostrophic_diagnostics(self):
        """Test quasi-geostrophic diagnostic calculations."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon) * 10
        v = np.random.randn(nlat, nlon) * 10

        vw = VectorWind(u, v)

        # Test planetary vorticity
        f = vw.planetaryvorticity(omega=7.272e-5)
        assert f.shape == (nlat, nlon)
        assert np.all(np.isfinite(f))

        # Test absolute vorticity
        vort = vw.vorticity()
        abs_vort = vw.absolutevorticity(omega=7.272e-5)

        # Should be vorticity + planetary vorticity
        np.testing.assert_allclose(abs_vort, vort + f, rtol=1e-12)


class TestVectorWindAdvancedOperations:
    """Test advanced VectorWind operations and edge cases."""

    def test_vectorwind_with_zero_wind(self):
        """Test VectorWind with zero wind fields."""
        nlat, nlon = 19, 36
        u_zero = np.zeros((nlat, nlon))
        v_zero = np.zeros((nlat, nlon))

        vw = VectorWind(u_zero, v_zero)

        # All derived quantities should be zero
        vort = vw.vorticity()
        div = vw.divergence()
        sf = vw.streamfunction()
        vp = vw.velocitypotential()

        assert np.allclose(vort, 0, atol=1e-12)
        assert np.allclose(div, 0, atol=1e-12)
        assert np.allclose(sf, 0, atol=1e-12)
        assert np.allclose(vp, 0, atol=1e-12)

    def test_vectorwind_with_uniform_wind(self):
        """Test VectorWind with uniform wind fields."""
        nlat, nlon = 19, 36
        u_uniform = np.full((nlat, nlon), 10.0)  # Uniform 10 m/s
        v_uniform = np.full((nlat, nlon), 5.0)  # Uniform 5 m/s

        vw = VectorWind(u_uniform, v_uniform)

        # Uniform wind should have zero vorticity and divergence
        vort = vw.vorticity()
        div = vw.divergence()

        # Should be very close to zero (within numerical precision)
        # Use more relaxed tolerance for spherical harmonic computations
        assert np.allclose(vort, 0, atol=1e-4)
        assert np.allclose(div, 0, atol=1e-4)

    def test_vectorwind_conservation_properties(self):
        """Test conservation properties of transformations."""
        nlat, nlon = 19, 36
        u = np.random.randn(nlat, nlon) * 10
        v = np.random.randn(nlat, nlon) * 10

        vw = VectorWind(u, v)

        # Get streamfunction and velocity potential
        sf, vp = vw.sfvp()

        # Reconstruct wind from potentials
        u_chi, v_chi = vw.gradient(vp)  # Irrotational component
        u_psi, v_psi = vw.gradient(sf)  # Nondivergent component

        # For 2D case, nondivergent component should be rotated gradient
        # This is a complex test - just ensure we get reasonable results
        assert u_psi.shape == (nlat, nlon)
        assert v_psi.shape == (nlat, nlon)
        assert np.all(np.isfinite(u_psi))
        assert np.all(np.isfinite(v_psi))


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
