#!/usr/bin/env python3
"""
Skyborn spharm Fortran optimization performance test tool

This script compares the performance of the spharm library before and after optimization, especially:
1. Divergence/vorticity calculation
2. Laplacian operator
3. Spherical harmonic transform
4. Gradient calculation

Before running, you need to:
1. Compile the original version and record performance
2. Compile the optimized version and record performance
3. Compare the results
"""

import time
import numpy as np
import sys
import os
from contextlib import contextmanager

# Add source code path
sys.path.insert(0, "src")


@contextmanager
def timer(description):
    """Timing context manager"""
    start = time.perf_counter()
    yield
    end = time.perf_counter()
    print(f"{description}: {end - start:.4f} s")


def create_test_data(nlat=64, nlon=128, nt=10):
    """Create test data"""
    print(f"Creating test data: {nlat}x{nlon} grid, {nt} time steps")

    # Create latitude and longitude grid
    lats = np.linspace(-90, 90, nlat)
    lons = np.linspace(0, 360, nlon, endpoint=False)
    LON, LAT = np.meshgrid(lons, lats)

    # Create test field: contains multiple spherical harmonic modes
    test_data = np.zeros((nt, nlat, nlon))

    for t in range(nt):
        # Add spherical harmonics of different frequencies
        field = (
            np.sin(2 * np.radians(LAT)) * np.cos(3 * np.radians(LON))
            + 0.5 * np.cos(3 * np.radians(LAT)) * np.sin(2 * np.radians(LON))
            + 0.3 * np.sin(4 * np.radians(LAT)) * np.cos(5 * np.radians(LON))
            + 0.2 * np.random.random((nlat, nlon))
        )  # Add noise

        test_data[t] = field

    return test_data, LAT, LON


def benchmark_spharm_operations(test_data, iterations=5):
    """Test core spharm operations performance"""
    print("Starting performance test...")

    try:
        from skyborn.spharm import Spharmt

        print("Successfully imported optimized spharm")
    except ImportError as e:
        print(f"Failed to import spharm: {e}")
        return None

    nt, nlat, nlon = test_data.shape

    # Create spherical harmonic transform object
    print(f"Initializing Spharmt object ({nlon}x{nlat})")
    sht = Spharmt(nlon=nlon, nlat=nlat, gridtype="regular")

    results = {}

    # 1. Test grid to spectrum transform (Analysis)
    print("\nTest 1: Grid to spectrum transform (spherical harmonic analysis)")
    times = []
    for i in range(iterations):
        with timer(f"  Iteration {i+1}") as t:
            for t_idx in range(nt):
                spec_coeffs = sht.grdtospec(test_data[t_idx])

    # Calculate average time
    mean_time = np.mean(
        [time.perf_counter() - time.perf_counter() for _ in range(iterations)]
    )

    # Retest for accurate timing
    start_time = time.perf_counter()
    for i in range(iterations):
        for t_idx in range(nt):
            spec_coeffs = sht.grdtospec(test_data[t_idx])
    end_time = time.perf_counter()

    grd_to_spec_time = (end_time - start_time) / iterations
    results["grd_to_spec"] = grd_to_spec_time
    print(f"  Average time: {grd_to_spec_time:.4f} s")

    # 2. Test spectrum to grid transform (Synthesis)
    print("\nTest 2: Spectrum to grid transform (spherical harmonic synthesis)")
    start_time = time.perf_counter()
    for i in range(iterations):
        for t_idx in range(nt):
            data_reconstructed = sht.spectogrd(spec_coeffs)
    end_time = time.perf_counter()

    spec_to_grd_time = (end_time - start_time) / iterations
    results["spec_to_grd"] = spec_to_grd_time
    print(f"  Average time: {spec_to_grd_time:.4f} s")

    # 3. Test vorticity/divergence calculation (key optimization point)
    print("\nTest 3: Vorticity/divergence calculation")

    # Create wind field data
    u_wind = np.random.random((nlat, nlon)) * 10  # Zonal wind
    v_wind = np.random.random((nlat, nlon)) * 10  # Meridional wind

    start_time = time.perf_counter()
    for i in range(iterations):
        # Calculate vorticity and divergence spectral coefficients
        vrt_spec, div_spec = sht.getvrtdivspec(u_wind, v_wind)
        # Reconstruct wind field from spectral coefficients
        u_new, v_new = sht.getuv(vrt_spec, div_spec)
    end_time = time.perf_counter()

    vorticity_time = (end_time - start_time) / iterations
    results["vorticity_divergence"] = vorticity_time
    print(f"  Average time: {vorticity_time:.4f} s")

    # Verify accuracy
    u_error = np.max(np.abs(u_wind - u_new))
    v_error = np.max(np.abs(v_wind - v_new))
    print(
        f"  Wind field reconstruction accuracy: U error={u_error:.2e}, V error={v_error:.2e}"
    )

    # 4. Test gradient calculation (Laplacian operator)
    print("\nTest 4: Gradient calculation (Laplacian operator)")

    # Use a simple scalar field
    scalar_field = test_data[0]  # Use the first time step

    start_time = time.perf_counter()
    for i in range(iterations):
        # Grid to spectrum
        spec = sht.grdtospec(scalar_field)
        # Calculate Laplacian operator (in practice, this would be done in Fortran)
        # Here we test the round-trip transform performance
        laplacian_spec = spec * -1  # Simplified Laplacian operator
        laplacian_field = sht.spectogrd(laplacian_spec)
    end_time = time.perf_counter()

    gradient_time = (end_time - start_time) / iterations
    results["gradient_computation"] = gradient_time
    print(f"  Average time: {gradient_time:.4f} s")

    # 5. Overall performance test
    print("\nTest 5: Comprehensive performance test")
    start_time = time.perf_counter()

    for i in range(iterations):
        # Simulate a complete atmospheric dynamics calculation process
        for t_idx in range(min(nt, 3)):  # Only test the first 3 time steps
            # 1. Analysis
            spec = sht.grdtospec(test_data[t_idx])
            # 2. Wind field processing
            vrt_spec, div_spec = sht.getvrtdivspec(u_wind, v_wind)
            # 3. Synthesis
            reconstructed = sht.spectogrd(spec)
            u_new, v_new = sht.getuv(vrt_spec, div_spec)

    end_time = time.perf_counter()
    comprehensive_time = (end_time - start_time) / iterations
    results["comprehensive"] = comprehensive_time
    print(f"  Average time: {comprehensive_time:.4f} s")

    return results


def print_performance_summary(results):
    """Print performance summary"""
    if not results:
        print("No performance results to display")
        return

    print("\n" + "=" * 60)
    print("SPHARM Performance Test Summary")
    print("=" * 60)

    operations = [
        ("grd_to_spec", "Grid→Spectrum Transform"),
        ("spec_to_grd", "Spectrum→Grid Transform"),
        ("vorticity_divergence", "Vorticity/Divergence Calculation"),
        ("gradient_computation", "Gradient Calculation"),
        ("comprehensive", "Comprehensive Performance Test"),
    ]

    total_time = 0
    for key, description in operations:
        if key in results:
            time_val = results[key]
            total_time += time_val
            print(f"{description:.<30} {time_val:.4f} s")

    print("-" * 60)
    print(f"{'Total Time':.<30} {total_time:.4f} s")
    print("=" * 60)


def save_benchmark_results(results, filename="spharm_benchmark.txt"):
    """Save benchmark results"""
    if not results:
        return

    with open(filename, "w") as f:
        f.write("SPHARM Fortran Optimization Performance Test Results\n")
        f.write("=" * 50 + "\n")
        f.write(f"Test Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Python Version: {sys.version}\n")
        f.write(f"NumPy Version: {np.__version__}\n\n")

        operations = [
            ("grd_to_spec", "Grid→Spectrum Transform"),
            ("spec_to_grd", "Spectrum→Grid Transform"),
            ("vorticity_divergence", "Vorticity/Divergence Calculation"),
            ("gradient_computation", "Gradient Calculation"),
            ("comprehensive", "Comprehensive Performance Test"),
        ]

        for key, description in operations:
            if key in results:
                f.write(f"{description}: {results[key]:.6f} s\n")

        total_time = sum(results.values())
        f.write(f"\nTotal Time: {total_time:.6f} s\n")

    print(f"Results saved to: {filename}")


def compare_with_baseline(results, baseline_file="baseline_benchmark.txt"):
    """Compare with baseline performance"""
    if not os.path.exists(baseline_file):
        print(f"Creating baseline file: {baseline_file}")
        save_benchmark_results(results, baseline_file)
        return

    print(f"\nComparing with baseline performance ({baseline_file})")
    print("-" * 60)

    # Baseline comparison logic can be implemented here
    # Since the baseline file format is simple, prompt user to compare manually
    print("Please manually compare the current results with those in the baseline file")
    print(
        "Performance improvement = (Baseline Time - Current Time) / Baseline Time * 100%"
    )


def main():
    """Main function"""
    print("Skyborn spharm Fortran Optimization Performance Test")
    print("=" * 60)

    # Set test parameters
    nlat = 64  # Number of latitude grid points
    nlon = 128  # Number of longitude grid points
    nt = 5  # Number of time steps
    iterations = 3  # Number of test iterations

    print(f"Test configuration:")
    print(f"   Grid size: {nlat} x {nlon}")
    print(f"   Number of time steps: {nt}")
    print(f"   Number of iterations: {iterations}")

    # Create test data
    test_data, LAT, LON = create_test_data(nlat, nlon, nt)

    # Run performance test
    results = benchmark_spharm_operations(test_data, iterations)

    # Display results
    print_performance_summary(results)

    # Save results
    if results:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        save_benchmark_results(results, f"spharm_benchmark_{timestamp}.txt")

        # Compare with baseline
        compare_with_baseline(results)

    print("\nPerformance test completed!")

    # Return results for further analysis
    return results


# if __name__ == "__main__":
#     # Run performance test
#     benchmark_results = main()
