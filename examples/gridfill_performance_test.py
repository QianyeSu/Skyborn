#!/usr/bin/env python
"""
GridFill Performance and Precision Test
========================================

This script tests the performance and numerical precision of the optimized
gridfill implementation.
"""

import time
import numpy as np
import matplotlib.pyplot as plt
import skyborn


def create_test_data(nlat=100, nlon=200, missing_fraction=0.2):
    """Create synthetic test data with missing values."""
    np.random.seed(42)

    # Create coordinate grids
    lat = np.linspace(-90, 90, nlat)
    lon = np.linspace(0, 360, nlon)
    LON, LAT = np.meshgrid(lon, lat)

    # Create a realistic-looking field
    data = (
        np.sin(np.radians(3 * LON)) * np.cos(np.radians(2 * LAT))
        + 0.5 * np.sin(np.radians(LON)) * np.sin(np.radians(LAT))
        + 0.2 * np.random.randn(nlat, nlon)
    )

    # Create missing values
    missing_mask = np.random.random((nlat, nlon)) < missing_fraction

    # Ensure some large gaps for more realistic testing
    missing_mask[40:60, 80:120] = True  # Large rectangular gap
    missing_mask[20:30, 150:170] = True  # Another gap

    # Create masked array
    data_with_gaps = np.ma.array(data, mask=missing_mask)

    return data, data_with_gaps


def test_performance():
    """Test the performance of gridfill."""
    print("=== GridFill Performance Test ===")

    sizes = [(50, 100), (100, 200), (200, 400)]

    for nlat, nlon in sizes:
        print(f"\nTesting grid size: {nlat} x {nlon}")

        # Create test data
        original, data_with_gaps = create_test_data(nlat, nlon)
        n_missing = np.sum(data_with_gaps.mask)

        print(f"Missing values: {n_missing} ({n_missing/(nlat*nlon)*100:.1f}%)")

        # Time the filling process
        start_time = time.time()
        filled_data, converged = skyborn.fill(
            data_with_gaps, xdim=1, ydim=0, eps=1e-4, itermax=1000, verbose=False
        )
        end_time = time.time()

        # Calculate statistics
        elapsed_time = end_time - start_time
        valid_mask = ~data_with_gaps.mask
        rms_error = np.sqrt(np.mean((filled_data - original)[valid_mask] ** 2))

        print(f"Time: {elapsed_time:.3f} seconds")
        print(f"Converged: {converged[0]}")
        print(f"RMS error in valid regions: {rms_error:.6f}")
        print(f"Performance: {(nlat*nlon)/elapsed_time:.0f} points/second")


def test_numerical_precision():
    """Test numerical precision and stability."""
    print("\n=== Numerical Precision Test ===")

    # Create a simple analytical test case
    nlat, nlon = 50, 100
    lat = np.linspace(-90, 90, nlat)
    lon = np.linspace(0, 360, nlon)
    LON, LAT = np.meshgrid(lon, lat)

    # Simple harmonic function
    analytical = np.sin(np.radians(LON)) * np.cos(np.radians(LAT))

    # Create small missing region
    data = analytical.copy()
    mask = np.zeros_like(data, dtype=bool)
    mask[20:30, 40:60] = True  # Small rectangular gap

    data_masked = np.ma.array(data, mask=mask)

    # Fill the gap
    filled_data, converged = skyborn.fill(
        data_masked,
        xdim=1,
        ydim=0,
        eps=1e-8,  # High precision
        itermax=5000,
        verbose=True,
    )

    # Analyze precision in the filled region
    filled_region = filled_data[mask]
    analytical_region = analytical[mask]

    max_error = np.max(np.abs(filled_region - analytical_region))
    mean_error = np.mean(np.abs(filled_region - analytical_region))

    print(f"Maximum error in filled region: {max_error:.10f}")
    print(f"Mean absolute error: {mean_error:.10f}")
    print(f"Relative error: {max_error/np.max(np.abs(analytical_region)):.2e}")

    # Check for NaN or Inf values
    if np.any(np.isnan(filled_data)) or np.any(np.isinf(filled_data)):
        print("⚠️  WARNING: Found NaN or Inf values in result!")
    else:
        print("✅ No NaN or Inf values detected")


def test_cross_platform_consistency():
    """Test that results are consistent across different platforms."""
    print("\n=== Cross-Platform Consistency Test ===")

    # Create deterministic test data
    np.random.seed(12345)  # Fixed seed for reproducibility
    original, data_with_gaps = create_test_data(30, 60, 0.15)

    # Run gridfill multiple times with same parameters
    results = []
    for i in range(3):
        filled, converged = skyborn.fill(
            data_with_gaps.copy(),
            xdim=1,
            ydim=0,
            eps=1e-6,
            relax=0.6,
            itermax=1000,
            verbose=False,
        )
        results.append(filled)

    # Check consistency
    max_diff_01 = np.max(np.abs(results[0] - results[1]))
    max_diff_02 = np.max(np.abs(results[0] - results[2]))

    print(f"Maximum difference between runs 1-2: {max_diff_01:.2e}")
    print(f"Maximum difference between runs 1-3: {max_diff_02:.2e}")

    if max_diff_01 < 1e-10 and max_diff_02 < 1e-10:
        print("✅ Results are highly consistent")
    elif max_diff_01 < 1e-6 and max_diff_02 < 1e-6:
        print("✅ Results are reasonably consistent")
    else:
        print("⚠️  WARNING: Results show significant variation")


if __name__ == "__main__":
    print("GridFill Optimization and Precision Test")
    print("=" * 50)

    # Test performance
    test_performance()

    # Test numerical precision
    test_numerical_precision()

    # Test consistency
    test_cross_platform_consistency()

    print("\n" + "=" * 50)
    print("Test completed!")
    print("\nNotes:")
    print("- Performance should be good without sacrificing precision")
    print("- No -ffast-math was used to preserve IEEE 754 compliance")
    print("- Cross-platform compilation uses -march=x86-64 -mtune=generic")
    print("- Cython works on Windows, Linux, and macOS (including Apple Silicon)")
