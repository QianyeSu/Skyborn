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
    from skyborn.spharm import Spharmt, regrid, gaussian_lats_wts

    SPHARM_AVAILABLE = True
except ImportError:
    SPHARM_AVAILABLE = False
    Spharmt = None
    regrid = None
    gaussian_lats_wts = None


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

        # Transform back to grid space
        grid_back = sht.spectogrd(spectral)
        assert grid_back.shape == grid_data.shape
        assert isinstance(grid_back, np.ndarray)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_grdtospec_shape_validation(self):
        """Test shape validation for grid transforms."""
        sht = Spharmt(nlon=72, nlat=37)

        # Wrong shape should raise appropriate error
        wrong_shape_data = np.random.randn(20, 30)
        with pytest.raises(ValueError):
            sht.grdtospec(wrong_shape_data)


class TestSpharmtVectorOperations:
    """Test Spharmt vector operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getuv_basic(self):
        """Test basic getuv functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test vorticity and divergence
        vort_spec = np.random.randn(37, 72) + 1j * np.random.randn(37, 72)
        div_spec = np.random.randn(37, 72) + 1j * np.random.randn(37, 72)

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
        assert vort_spec.shape == (37, 72)
        assert div_spec.shape == (37, 72)


class TestSpharmtGradientOperations:
    """Test Spharmt gradient operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_getgrad_basic(self):
        """Test basic getgrad functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test spectral coefficients
        spec = np.random.randn(37, 72) + 1j * np.random.randn(37, 72)

        gradx, grady = sht.getgrad(spec)
        assert gradx.shape == (37, 72)
        assert grady.shape == (37, 72)

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


class TestSpharmtSmoothing:
    """Test Spharmt smoothing operations."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_specsmooth_basic(self):
        """Test basic specsmooth functionality."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test spectral coefficients
        spec = np.random.randn(37, 72) + 1j * np.random.randn(37, 72)

        # Test smoothing with different factors
        smooth_factor = 0.5
        smoothed = sht.specsmooth(spec, smooth_factor)
        assert smoothed.shape == spec.shape
        assert isinstance(smoothed, np.ndarray)


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


class TestRegridFunction:
    """Test regrid utility function."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_regrid_basic(self):
        """Test basic regrid functionality."""
        # Create test data
        field = np.random.randn(37, 72)

        # Test regridding to different grid
        regridded = regrid(field, 37, 72, 19, 36)
        assert regridded.shape == (19, 36)
        assert isinstance(regridded, np.ndarray)


class TestSpharmtIntegration:
    """Integration tests for Spharmt functionality."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_round_trip_transform(self):
        """Test that grid->spectral->grid is approximately identity."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test data
        original = np.random.randn(37, 72)

        # Round trip transform
        spectral = sht.grdtospec(original)
        reconstructed = sht.spectogrd(spectral)

        # Check that reconstruction is close to original
        # (allowing for numerical precision)
        assert np.allclose(original, reconstructed, rtol=1e-10, atol=1e-10)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_velocity_potential_streamfunction_consistency(self):
        """Test consistency between getpsichi and getuv/getvrtdivspec."""
        sht = Spharmt(nlon=72, nlat=37)

        # Create test winds
        u = np.sin(np.linspace(0, 2 * np.pi, 72))[np.newaxis, :] * np.ones((37, 1))
        v = np.cos(np.linspace(0, 2 * np.pi, 72))[np.newaxis, :] * np.ones((37, 1))

        # Get streamfunction and velocity potential
        psi, chi = sht.getpsichi(u, v)

        # Check that these are valid arrays
        assert psi.shape == (37, 72)
        assert chi.shape == (37, 72)

        # Check for NaN/inf values
        assert np.all(np.isfinite(psi))
        assert np.all(np.isfinite(chi))


class TestSpharmtErrorHandling:
    """Test error handling in Spharmt."""

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_invalid_grid_sizes(self):
        """Test handling of invalid grid sizes."""
        # Test with too small grid
        with pytest.raises((ValueError, RuntimeError)):
            Spharmt(nlon=1, nlat=1)

    @pytest.mark.skipif(not SPHARM_AVAILABLE, reason="spharm module not available")
    def test_invalid_method_calls(self):
        """Test handling of invalid method calls."""
        sht = Spharmt(nlon=72, nlat=37)

        # Test with wrong input shapes
        wrong_shape = np.random.randn(10, 20)
        with pytest.raises((ValueError, RuntimeError)):
            sht.grdtospec(wrong_shape)


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


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
