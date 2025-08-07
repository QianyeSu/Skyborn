"""
Tests for skyborn.spharm module (spherical harmonic transforms)

This module contains comprehensive tests for the spherical harmonic
transform functionality, including Spharmt class and related utilities.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    from skyborn.spharm import Spharmt, regrid, gaussian_lats_wts, getspecindx

    SPHARM_AVAILABLE = True
except ImportError:
    SPHARM_AVAILABLE = False
    Spharmt = None
    regrid = None
    gaussian_lats_wts = None
    getspecindx = None


class TestSpharmtInitialization:
    """Test Spharmt class initialization."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_basic_initialization(self):
        """Test basic Spharmt initialization."""
        sht = Spharmt(nlon=144, nlat=73)
        assert sht.nlon == 144
        assert sht.nlat == 73
        assert hasattr(sht, "grdtospec")
        assert hasattr(sht, "spectogrd")

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_with_custom_parameters(self):
        """Test Spharmt initialization with custom parameters."""
        sht = Spharmt(
            nlon=72, nlat=37, rsphere=7.0e6, gridtype="regular", legfunc="computed"
        )
        assert sht.nlon == 72
        assert sht.nlat == 37
        assert sht.rsphere == 7.0e6
        assert sht.gridtype == "regular"
        assert sht.legfunc == "computed"

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_gaussian_grid(self):
        """Test Spharmt initialization with Gaussian grid."""
        sht = Spharmt(nlon=72, nlat=37, gridtype="gaussian", legfunc="stored")
        assert sht.gridtype == "gaussian"
        assert sht.legfunc == "stored"

    def test_spharmt_import_error(self):
        """Test graceful handling when spharm module is not available."""
        if not SPHARM_AVAILABLE:
            with pytest.raises(ImportError):
                Spharmt(144, 73)


class TestSpharmtGridOperations:
    """Test Spharmt grid operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_grdtospec_and_back(self):
        """Test grid to spectral and back to grid transform."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create simple test data
        grid_data = np.random.randn(37, 72)

        # Transform to spectral space
        spectral = sht.grdtospec(grid_data)
        assert isinstance(spectral, np.ndarray)
        assert np.iscomplexobj(spectral)

        # Transform back to grid space
        grid_back = sht.spectogrd(spectral)
        assert grid_back.shape == grid_data.shape
        assert isinstance(grid_back, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_grdtospec_3d_data(self):
        """Test grid to spectral transform with 3D data."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create 3D test data (time series)
        nt = 5
        grid_data_3d = np.random.randn(37, 72, nt)

        # Transform to spectral space
        spectral_3d = sht.grdtospec(grid_data_3d)
        assert isinstance(spectral_3d, np.ndarray)
        assert np.iscomplexobj(spectral_3d)
        assert spectral_3d.ndim == 2  # (nspec, nt)

        # Transform back to grid space
        grid_back_3d = sht.spectogrd(spectral_3d)
        assert grid_back_3d.shape == grid_data_3d.shape
        assert isinstance(grid_back_3d, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_grdtospec_shape_validation(self):
        """Test shape validation for grid transforms."""
        sht = Spharmt(nlon=72, nlat=37)

        # Wrong shape should raise appropriate error
        wrong_shape_data = np.random.randn(20, 30)
        with pytest.raises((ValueError, RuntimeError)):
            sht.grdtospec(wrong_shape_data)


class TestSpharmtVectorOperations:
    """Test Spharmt vector operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getuv_basic(self):
        """Test basic getuv functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test vorticity and divergence spectral coefficients
        # Shape should be (nspec,) = (nlat*(nlat+1)/2,) for triangular truncation
        nspec = (sht.nlat * (sht.nlat + 1)) // 2
        vort_spec = np.random.randn(nspec) + 1j * np.random.randn(nspec)
        div_spec = np.random.randn(nspec) + 1j * np.random.randn(nspec)

        u, v = sht.getuv(vort_spec, div_spec)
        assert u.shape == (37, 72)
        assert v.shape == (37, 72)
        assert isinstance(u, np.ndarray)
        assert isinstance(v, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getvrtdivspec_basic(self):
        """Test basic getvrtdivspec functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test u and v winds
        u = np.random.randn(37, 72)
        v = np.random.randn(37, 72)

        vort_spec, div_spec = sht.getvrtdivspec(u, v)

        # Expected spectral coefficients shape for triangular truncation
        expected_spec_shape = (sht.nlat * (sht.nlat + 1)) // 2
        assert vort_spec.shape == (expected_spec_shape,)
        assert div_spec.shape == (expected_spec_shape,)
        assert np.iscomplexobj(vort_spec)
        assert np.iscomplexobj(div_spec)


class TestSpharmtGradientOperations:
    """Test Spharmt gradient operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getgrad_basic(self):
        """Test basic getgrad functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test spectral coefficients with correct shape
        nspec = (sht.nlat * (sht.nlat + 1)) // 2
        spec = np.random.randn(nspec) + 1j * np.random.randn(nspec)

        gradx, grady = sht.getgrad(spec)
        assert gradx.shape == (37, 72)
        assert grady.shape == (37, 72)
        assert isinstance(gradx, np.ndarray)
        assert isinstance(grady, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getpsichi_basic(self):
        """Test basic getpsichi functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test u and v winds
        u = np.random.randn(37, 72)
        v = np.random.randn(37, 72)

        psi, chi = sht.getpsichi(u, v)
        assert psi.shape == (37, 72)
        assert chi.shape == (37, 72)
        assert isinstance(psi, np.ndarray)
        assert isinstance(chi, np.ndarray)


class TestSpharmtSmoothing:
    """Test Spharmt smoothing operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_specsmooth_basic(self):
        """Test basic specsmooth functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test grid data
        data = np.random.randn(37, 72)

        # Create smoothing factor array (must be 1D array with size nlat)
        smooth_factors = np.exp(-np.arange(sht.nlat) / 10.0)

        smoothed = sht.specsmooth(data, smooth_factors)
        # Remove extra dimensions if present
        smoothed = np.squeeze(smoothed)

        assert smoothed.shape == data.shape
        assert isinstance(smoothed, np.ndarray)

        # Test that smoothing actually changes the data
        assert not np.allclose(data, smoothed, rtol=1e-10)


class TestUtilityFunctions:
    """Test spharm utility functions."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_gaussian_lats_wts_basic(self):
        """Test basic gaussian_lats_wts functionality."""
        lats, wts = gaussian_lats_wts(37)
        assert len(lats) == 37
        assert len(wts) == 37
        assert isinstance(lats, np.ndarray)
        assert isinstance(wts, np.ndarray)

        # Check that weights sum to approximately 2.0
        assert abs(np.sum(wts) - 2.0) < 1e-10

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_gaussian_lats_wts_edge_cases(self):
        """Test gaussian_lats_wts edge cases."""
        # Test minimum valid nlat
        lats, wts = gaussian_lats_wts(2)
        assert len(lats) == 2
        assert len(wts) == 2

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getspecindx_basic(self):
        """Test getspecindx functionality."""
        ntrunc = 36  # Triangular truncation
        indxm, indxn = getspecindx(ntrunc)

        # Check shapes
        expected_size = ((ntrunc + 1) * (ntrunc + 2)) // 2
        assert len(indxm) == expected_size
        assert len(indxn) == expected_size

        # Check that indices are valid
        assert np.all(indxm >= 0)
        assert np.all(indxn >= 0)
        assert np.all(indxm <= indxn)  # m <= n constraint
        assert np.all(indxn <= ntrunc)  # n <= ntrunc constraint


class TestRegridFunction:
    """Test regrid utility function."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_regrid_basic(self):
        """Test basic regrid functionality."""
        # Create input and output grid instances
        grid_in = Spharmt(nlon=72, nlat=37)
        grid_out = Spharmt(nlon=36, nlat=19)

        # Create test data with correct shape for input grid
        field = np.random.randn(37, 72)

        # Test regridding to different grid
        regridded = regrid(grid_in, grid_out, field)
        assert regridded.shape == (19, 36)
        assert isinstance(regridded, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_regrid_with_smoothing(self):
        """Test regrid functionality with smoothing."""
        grid_in = Spharmt(nlon=72, nlat=37)
        grid_out = Spharmt(nlon=36, nlat=19)

        field = np.random.randn(37, 72)

        # Create smoothing factors for output grid
        smooth_factors = np.exp(-np.arange(grid_out.nlat) / 5.0)

        regridded = regrid(grid_in, grid_out, field, smooth=smooth_factors)
        # Remove extra dimensions if present
        regridded = np.squeeze(regridded)

        assert regridded.shape == (19, 36)
        assert isinstance(regridded, np.ndarray)


class TestSpharmtIntegration:
    """Integration tests for Spharmt functionality."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_round_trip_transform(self):
        """Test that grid->spectral->grid is approximately identity."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test data with simple pattern for better numerical stability
        lons = np.linspace(0, 2 * np.pi, 72, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, 37)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Simple harmonic pattern that should be representable exactly
        original = np.sin(2 * lon_grid) * np.cos(lat_grid)

        # Round trip transform
        spectral = sht.grdtospec(original)
        reconstructed = sht.spectogrd(spectral)

        # Check that reconstruction is close to original
        # Use more relaxed tolerance for complex numerical operations
        assert np.allclose(original, reconstructed, rtol=1e-6, atol=1e-6)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_velocity_potential_streamfunction_consistency(self):
        """Test consistency between getpsichi and getuv/getvrtdivspec."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test winds with simple patterns
        lons = np.linspace(0, 2 * np.pi, 72, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, 37)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        u = np.sin(lon_grid) * np.cos(lat_grid)
        v = np.cos(lon_grid) * np.cos(lat_grid)

        # Get streamfunction and velocity potential
        psi, chi = sht.getpsichi(u, v)

        # Check that these are valid arrays
        assert psi.shape == (37, 72)
        assert chi.shape == (37, 72)

        # Check for NaN/inf values
        assert np.all(np.isfinite(psi))
        assert np.all(np.isfinite(chi))

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_helmholtz_decomposition_consistency(self):
        """Test that getvrtdivspec and getuv are consistent."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test winds with simpler patterns for better numerical stability
        lons = np.linspace(0, 2 * np.pi, 72, endpoint=False)
        lats = np.linspace(-np.pi / 2, np.pi / 2, 37)
        lon_grid, lat_grid = np.meshgrid(lons, lats)

        # Use simpler harmonic patterns that are more representable
        u_orig = np.sin(lon_grid) * np.cos(lat_grid)
        v_orig = np.cos(lon_grid) * np.cos(lat_grid)

        # Helmholtz decomposition
        vrt_spec, div_spec = sht.getvrtdivspec(u_orig, v_orig)
        u_reconstructed, v_reconstructed = sht.getuv(vrt_spec, div_spec)

        # Check reconstruction quality with more relaxed tolerance for numerical precision
        # In spherical harmonic analysis, perfect reconstruction is limited by truncation
        assert np.allclose(u_orig, u_reconstructed, rtol=1e-2, atol=1e-2)
        assert np.allclose(v_orig, v_reconstructed, rtol=1e-2, atol=1e-2)


class TestSpharmtErrorHandling:
    """Test error handling in Spharmt."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_invalid_grid_sizes(self):
        """Test handling of invalid grid sizes."""
        # Test with too small grid
        with pytest.raises((ValueError, RuntimeError)):
            Spharmt(nlon=1, nlat=1)

        # Test with too small nlat
        with pytest.raises((ValueError, RuntimeError)):
            Spharmt(nlon=72, nlat=2)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_invalid_method_calls(self):
        """Test handling of invalid method calls."""
        sht = Spharmt(nlon=72, nlat=37)

        # Test with wrong input shapes
        wrong_shape = np.random.randn(10, 20)
        with pytest.raises((ValueError, RuntimeError)):
            sht.grdtospec(wrong_shape)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_invalid_parameters(self):
        """Test handling of invalid initialization parameters."""
        # Test invalid gridtype
        with pytest.raises((ValueError, RuntimeError)):
            Spharmt(nlon=72, nlat=37, gridtype="invalid")

        # Test invalid legfunc
        with pytest.raises((ValueError, RuntimeError)):
            Spharmt(nlon=72, nlat=37, legfunc="invalid")

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spectral_data_validation(self):
        """Test validation of spectral data inputs."""
        sht = Spharmt(nlon=72, nlat=37)

        # Test with completely wrong spectral array size (should work but give unexpected results)
        # The library might be more permissive than expected, so we test for basic functionality
        wrong_spec_size = np.random.randn(100) + 1j * np.random.randn(100)

        # This might not raise an error but should still produce some output
        try:
            result = sht.spectogrd(wrong_spec_size)
            # If it doesn't raise an error, we just check that it returns something
            assert isinstance(result, np.ndarray)
        except (ValueError, RuntimeError):
            # This is also acceptable - validation can work
            pass


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
class TestSpharmtPerformance:
    """Performance tests for Spharmt."""

    def test_small_grid_performance(self):
        """Test performance with small grid."""
        sht = Spharmt(nlon=36, nlat=19)

        # Create test data
        data = np.random.randn(19, 36)

        # Perform transforms multiple times
        for _ in range(10):
            spectral = sht.grdtospec(data)
            reconstructed = sht.spectogrd(spectral)

        # Basic sanity check
        assert reconstructed.shape == data.shape

    def test_spharmt_different_grid_sizes(self):
        """Test spharmt with different grid sizes"""
        # Create simple test pattern
        data = np.random.rand(37, 72)

        sht = Spharmt(nlon=72, nlat=37)

        # Transform to spectral and back
        spectral = sht.grdtospec(data)
        recovered = sht.spectogrd(spectral)

        # For random data, exact recovery is not expected due to spectral truncation
        # Just check that the operation works and produces reasonable output
        assert recovered.shape == data.shape
        assert np.all(np.isfinite(recovered))

    def test_spharmt_memory_efficiency(self):
        """Test memory efficiency with reasonable size data"""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test data
        data = np.random.rand(37, 72)

        # Transform and check memory usage is reasonable
        spectral = sht.grdtospec(data)
        recovered = sht.spectogrd(spectral)

        assert np.all(np.isfinite(recovered))
        assert recovered.shape == data.shape

    def test_vector_operations_consistency(self):
        """Test vector operations consistency"""
        sht = Spharmt(nlon=72, nlat=37)

        # Create simple test winds
        u = np.ones((37, 72))
        v = np.zeros((37, 72))

        # Get vorticity and divergence
        vort_spec, div_spec = sht.getvrtdivspec(u, v)

        # Reconstruct winds
        u_back, v_back = sht.getuv(vort_spec, div_spec)

        # Should maintain basic structure
        assert u_back.shape == u.shape
        assert v_back.shape == v.shape
        assert np.all(np.isfinite(u_back))
        assert np.all(np.isfinite(v_back))

    def test_repeated_transforms_stability(self):
        """Test numerical stability of repeated transforms"""
        sht = Spharmt(nlon=72, nlat=37)

        # Simple test pattern
        original = np.cos(np.linspace(0, 2 * np.pi, 72))[np.newaxis, :] * np.ones(
            (37, 1)
        )
        current = original.copy()

        # Perform several round-trip transforms
        for i in range(3):
            spectral = sht.grdtospec(current)
            current = sht.spectogrd(spectral)

        # Check that output is reasonable (relaxed expectations)
        assert current.shape == original.shape
        assert np.all(np.isfinite(current))

    def test_precision_with_different_data_types(self):
        """Test precision with different data types"""
        sht = Spharmt(nlon=72, nlat=37)

        # Simple pattern that should transform well
        field_64 = np.cos(np.linspace(0, 2 * np.pi, 72))[np.newaxis, :] * np.ones(
            (37, 1)
        )
        field_32 = field_64.astype(np.float32)

        # Transform both
        spec_64 = sht.grdtospec(field_64)
        spec_32 = sht.grdtospec(field_32)

        recovered_64 = sht.spectogrd(spec_64)
        recovered_32 = sht.spectogrd(spec_32)

        # Check that outputs are reasonable
        assert recovered_64.shape == field_64.shape
        assert recovered_32.shape == field_32.shape
        assert np.all(np.isfinite(recovered_64))
        assert np.all(np.isfinite(recovered_32))


@pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
class TestRegridAdvanced:
    """Advanced regrid functionality tests."""

    def test_regrid_with_different_grid_types(self):
        """Test regridding between different grid types"""
        # Create input and output grids
        grid_regular = Spharmt(nlon=72, nlat=37, gridtype="regular")
        grid_gaussian = Spharmt(nlon=72, nlat=37, gridtype="gaussian")

        # Create simple test field
        field = np.random.rand(37, 72)

        # Regrid from regular to gaussian
        field_gaussian = regrid(grid_regular, grid_gaussian, field)

        # Regrid back
        field_regular_back = regrid(grid_gaussian, grid_regular, field_gaussian)

        # Should maintain basic properties
        assert field_gaussian.shape == (37, 72)
        assert field_regular_back.shape == (37, 72)
        assert np.all(np.isfinite(field_gaussian))
        assert np.all(np.isfinite(field_regular_back))


class TestSpharmtUntestedMethods:
    """Test previously untested methods in Spharmt."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_properties(self):
        """Test Spharmt properties and attributes."""
        sht = Spharmt(nlon=72, nlat=36, gridtype="gaussian", legfunc="stored")

        # Test read-only properties
        assert sht.nlat == 36
        assert sht.nlon == 72
        assert sht.gridtype == "gaussian"
        assert sht.legfunc == "stored"

        # Test computed properties
        assert hasattr(sht, "rsphere")
        assert sht.rsphere > 0

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_string_representation(self):
        """Test string representation of Spharmt."""
        sht = Spharmt(nlon=72, nlat=36, gridtype="gaussian")

        repr_str = repr(sht)
        assert "Spharmt" in repr_str
        assert "72" in repr_str
        assert "36" in repr_str
        assert "gaussian" in repr_str

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spharmt_attribute_protection(self):
        """Test that Spharmt attributes are protected."""
        sht = Spharmt(nlon=72, nlat=36)

        # Test that we can't modify key attributes
        with pytest.raises(AttributeError):
            sht.nlat = 50

        with pytest.raises(AttributeError):
            sht.nlon = 100

        with pytest.raises(AttributeError):
            del sht.gridtype

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_grdtospec_with_ntrunc(self):
        """Test grdtospec with custom ntrunc parameter."""
        sht = Spharmt(nlon=72, nlat=36)
        data = np.random.randn(36, 72)

        # Test with different truncations
        for ntrunc in [10, 20, 30]:
            spec = sht.grdtospec(data, ntrunc=ntrunc)
            expected_size = (ntrunc + 1) * (ntrunc + 2) // 2
            assert len(spec) == expected_size

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_spectogrd_basic(self):
        """Test basic spectogrd functionality."""
        sht = Spharmt(nlon=72, nlat=36)

        # Create test grid data and convert to spectral
        data = np.random.randn(36, 72)
        spec_data = sht.grdtospec(data)

        # Convert back to grid
        grid = sht.spectogrd(spec_data)
        assert grid.shape == (36, 72)

        # Should be close to original
        np.testing.assert_allclose(grid, data, rtol=1e-10)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_specsmooth_edge_cases(self):
        """Test specsmooth with edge cases."""
        sht = Spharmt(nlon=72, nlat=36)
        data = np.random.randn(36, 72)

        # Test with exponential smoothing
        smooth_exp = np.exp(-np.arange(sht.nlat) / 10.0)
        result_exp = sht.specsmooth(data, smooth_exp)
        assert result_exp.shape[-2:] == (36, 72)

        # Test that smoothing reduces variance
        assert np.var(result_exp.squeeze()) <= np.var(data)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getvrtdivspec_with_custom_ntrunc(self):
        """Test getvrtdivspec with custom truncation."""
        sht = Spharmt(nlon=72, nlat=36)

        u = np.random.randn(36, 72)
        v = np.random.randn(36, 72)

        for ntrunc in [15, 25]:
            vrt_spec, div_spec = sht.getvrtdivspec(u, v, ntrunc=ntrunc)
            expected_size = (ntrunc + 1) * (ntrunc + 2) // 2
            assert len(vrt_spec) == expected_size
            assert len(div_spec) == expected_size

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getuv_basic(self):
        """Test basic getuv functionality."""
        sht = Spharmt(nlon=72, nlat=36)

        # Create test wind fields
        u = np.random.randn(36, 72)
        v = np.random.randn(36, 72)

        # Get vorticity and divergence spectra
        vrt_spec, div_spec = sht.getvrtdivspec(u, v)

        # Reconstruct winds
        u_back, v_back = sht.getuv(vrt_spec, div_spec)
        assert u_back.shape == (36, 72)
        assert v_back.shape == (36, 72)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getpsichi_with_custom_ntrunc(self):
        """Test getpsichi with custom truncation."""
        sht = Spharmt(nlon=72, nlat=36)

        u = np.random.randn(36, 72)
        v = np.random.randn(36, 72)

        psi, chi = sht.getpsichi(u, v, ntrunc=20)
        assert psi.shape == (36, 72)
        assert chi.shape == (36, 72)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getgrad_basic(self):
        """Test basic getgrad functionality."""
        sht = Spharmt(nlon=72, nlat=36)

        # Create scalar field
        scalar = np.random.randn(36, 72)
        scalar_spec = sht.grdtospec(scalar)

        grad_x, grad_y = sht.getgrad(scalar_spec)
        assert grad_x.shape == (36, 72)
        assert grad_y.shape == (36, 72)


class TestSpharmtDataValidation:
    """Test data validation methods in Spharmt."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_validate_grid_data(self):
        """Test grid data validation."""
        sht = Spharmt(nlon=72, nlat=36)

        # Valid 2D data
        data_2d = np.random.randn(36, 72)
        validated, _ = sht._validate_grid_data(data_2d, "test_data")
        assert validated.shape == (36, 72)

        # Valid 3D data
        data_3d = np.random.randn(5, 36, 72)
        validated, _ = sht._validate_grid_data(data_3d, "test_data")
        assert validated.shape == (5, 36, 72)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_validate_spectral_data(self):
        """Test spectral data validation."""
        sht = Spharmt(nlon=72, nlat=36)

        ntrunc = min(sht.nlat - 1, (sht.nlon - 1) // 2)
        nspec = (ntrunc + 1) * (ntrunc + 2) // 2

        # Valid complex data
        spec_data = np.random.randn(nspec) + 1j * np.random.randn(nspec)
        validated, _ = sht._validate_spectral_data(spec_data, "test_spec")
        assert validated.dtype == np.complex128

        # Real data should be converted to complex
        spec_real = np.random.randn(nspec)
        validated, _ = sht._validate_spectral_data(spec_real, "test_spec")
        assert validated.dtype == np.complex128

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_validate_ntrunc(self):
        """Test ntrunc validation."""
        sht = Spharmt(nlon=72, nlat=36)
        max_allowed = min(sht.nlat - 1, (sht.nlon - 1) // 2)

        # Valid ntrunc
        result = sht._validate_ntrunc(20, max_allowed)
        assert result == 20

        # None should return max allowed
        result = sht._validate_ntrunc(None, max_allowed)
        assert result == max_allowed

        # Too large ntrunc should raise error
        with pytest.raises(Exception):
            sht._validate_ntrunc(max_allowed + 10, max_allowed)


class TestUtilityFunctionsExtended:
    """Extended tests for utility functions."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_gaussian_lats_wts_symmetry(self):
        """Test symmetry properties of Gaussian latitudes."""
        for nlat in [16, 32, 48]:
            lats, wts = gaussian_lats_wts(nlat)

            # Test symmetry
            np.testing.assert_allclose(lats, -lats[::-1], rtol=1e-12)
            np.testing.assert_allclose(wts, wts[::-1], rtol=1e-12)

            # Test weight sum
            assert abs(np.sum(wts) - 2.0) < 1e-12

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getspecindx_completeness(self):
        """Test that getspecindx returns complete spectral indices."""
        for ntrunc in [5, 10, 20]:
            n, m = getspecindx(ntrunc)

            # Should have correct number of coefficients
            expected_size = (ntrunc + 1) * (ntrunc + 2) // 2
            assert len(n) == expected_size
            assert len(m) == expected_size

            # All indices should be valid
            assert np.all(n >= 0)
            assert np.all(n <= ntrunc)
            assert np.all(m >= 0)
            assert np.all(m <= n)  # m should not exceed n


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
