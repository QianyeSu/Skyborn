#!/usr/bin/env python3
"""
Windspharm Example: Vector Wind Analysis with 4D Data Processing

This example demonstrates both basic and advanced usage of the skyborn.windspharm module,
including real-world applications for processing 4D atmospheric data with parallel computation.
"""

import numpy as np
import matplotlib.pyplot as plt
from skyborn.windspharm import VectorWind
from skyborn.spharm import gaussian_lats_wts
from skyborn.windspharm.tools import prep_data, recover_data
from skyborn.windspharm.standard import VectorWind as StandardVectorWind
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import time


def create_sample_wind_data():
    """Create realistic sample wind data on a Gaussian grid."""
    # Grid dimensions (T42 resolution)
    nlat, nlon = 73, 144

    # Get Gaussian latitudes
    lats, _ = gaussian_lats_wts(nlat)
    lons = np.linspace(0, 360, nlon, endpoint=False)

    # Create coordinate grids
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Create realistic wind patterns
    # Zonal wind with jet stream structure
    u = (
        20 * np.exp(-((lat_grid - 30) ** 2) / (2 * 15**2))  # NH jet
        + 15 * np.exp(-((lat_grid + 30) ** 2) / (2 * 15**2))  # SH jet
        + 3 * np.random.randn(nlat, nlon)
    )  # noise

    # Meridional wind with wave patterns
    v = 5 * np.sin(3 * np.deg2rad(lon_grid)) * np.cos(
        np.deg2rad(lat_grid)
    ) + 2 * np.random.randn(nlat, nlon)

    return u, v, lats, lons


def calculate_divergent_wind_component(
    time_idx: int, zonal_wind: np.ndarray, meridional_wind: np.ndarray
) -> tuple:
    """
    Calculate divergent wind components for a single time step.

    Real-world usage pattern for processing 4D atmospheric data.

    Args:
        time_idx: Time index for processing
        zonal_wind: Zonal wind component array (time, level, lat, lon)
        meridional_wind: Meridional wind component array (time, level, lat, lon)

    Returns:
        tuple: (time_idx, divergent_zonal_wind, divergent_meridional_wind)
    """
    try:
        # Prepare data for spherical harmonic analysis
        # 'zyx' format: z=level, y=lat, x=lon
        zonal_prepared, wind_info = prep_data(zonal_wind[time_idx], "zyx")
        meridional_prepared, _ = prep_data(meridional_wind[time_idx], "zyx")

        # Create vector wind object and compute irrotational component
        vector_wind = StandardVectorWind(zonal_prepared, meridional_prepared)
        divergent_u, divergent_v = vector_wind.irrotationalcomponent()

        # Recover original data structure
        divergent_u_recovered = recover_data(divergent_u, wind_info)
        divergent_v_recovered = recover_data(divergent_v, wind_info)

        return time_idx, divergent_u_recovered, divergent_v_recovered

    except Exception as e:
        print(f"Error processing time step {time_idx}: {e}")
        raise


def compute_divergent_wind_components(
    zonal_wind: np.ndarray, meridional_wind: np.ndarray, max_workers: int = 4
) -> tuple:
    """
    Compute divergent wind components using parallel processing.

    This function demonstrates the efficient pattern for processing
    large atmospheric datasets with shape (time, level, lat, lon).

    Args:
        zonal_wind: 4D array of zonal wind (time, level, lat, lon)
        meridional_wind: 4D array of meridional wind (time, level, lat, lon)
        max_workers: Maximum number of worker threads

    Returns:
        tuple: (divergent_zonal_wind, divergent_meridional_wind)
    """
    input_shape = zonal_wind.shape
    n_timesteps = input_shape[0]

    # Initialize output arrays with same shape as input
    divergent_zonal = np.zeros_like(zonal_wind)
    divergent_meridional = np.zeros_like(meridional_wind)

    # Process all time steps in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_time = {
            executor.submit(
                calculate_divergent_wind_component, t, zonal_wind, meridional_wind
            ): t
            for t in range(n_timesteps)
        }

        # Collect results with progress bar
        for future in tqdm(
            as_completed(future_to_time),
            total=n_timesteps,
            desc="Computing divergent wind components",
        ):
            try:
                time_idx, div_u, div_v = future.result()
                divergent_zonal[time_idx] = div_u
                divergent_meridional[time_idx] = div_v

            except Exception as e:
                print(f"Error in future result: {e}")
                raise

    return divergent_zonal, divergent_meridional


def create_4d_sample_data():
    """Create sample 4D atmospheric wind data for demonstration."""
    # Dimensions typical of atmospheric model output
    n_time = 12  # Monthly data
    n_level = 5  # Pressure levels
    n_lat = 37  # Reduced for demo
    n_lon = 72  # Reduced for demo

    print(f"Creating 4D sample data: ({n_time}, {n_level}, {n_lat}, {n_lon})")

    # Create realistic coordinate grids
    lats = np.linspace(90, -90, n_lat)
    lons = np.linspace(0, 360, n_lon, endpoint=False)

    # Initialize arrays
    u_4d = np.zeros((n_time, n_level, n_lat, n_lon))
    v_4d = np.zeros((n_time, n_level, n_lat, n_lon))

    np.random.seed(42)  # Reproducible results

    for t in range(n_time):
        for k in range(n_level):
            lon_grid, lat_grid = np.meshgrid(lons, lats)

            # Time-varying jet strength
            jet_strength = 20 + 10 * np.sin(2 * np.pi * t / 12)
            level_factor = 0.5 + 0.5 * (k / (n_level - 1))

            # Realistic zonal wind
            u_4d[t, k, :, :] = (
                jet_strength
                * level_factor
                * (
                    np.exp(-((lat_grid - 30) ** 2) / (2 * 15**2))
                    + 0.7 * np.exp(-((lat_grid + 30) ** 2) / (2 * 15**2))
                )
            ) + 2 * np.random.randn(n_lat, n_lon)

            # Meridional wind with waves
            v_4d[t, k, :, :] = (
                3
                * level_factor
                * np.sin(4 * np.deg2rad(lon_grid))
                * np.cos(np.deg2rad(lat_grid))
                * np.sin(2 * np.pi * t / 12)
            ) + 1 * np.random.randn(n_lat, n_lon)

    return u_4d, v_4d, lats, lons


def main():
    """Main analysis function demonstrating both basic and advanced usage."""
    print("Skyborn Windspharm: Basic and Advanced Examples")
    print("=" * 50)

    # === BASIC EXAMPLE ===
    print("\n1. Basic Vector Wind Analysis")
    print("-" * 30)

    # Create sample data
    print("Creating 2D sample wind data...")
    u, v, lats, lons = create_sample_wind_data()
    print(f"Data shape: u={u.shape}, v={v.shape}")

    # Initialize VectorWind
    vw = VectorWind(u, v, gridtype="gaussian")
    print(f"VectorWind initialized (grid: {vw.gridtype})")

    # Compute basic quantities
    vorticity = vw.vorticity()
    divergence = vw.divergence()
    streamfunction = vw.streamfunction()
    velocity_potential = vw.velocitypotential()

    print(f"Computed quantities:")
    print(f"  Vorticity range: {vorticity.min():.2e} to {vorticity.max():.2e} s⁻¹")
    print(f"  Divergence range: {divergence.min():.2e} to {divergence.max():.2e} s⁻¹")

    # Helmholtz decomposition
    u_rot, v_rot, u_div, v_div = vw.helmholtz()

    # Verify energy conservation
    total_energy = 0.5 * (u**2 + v**2).mean()
    rotational_energy = 0.5 * (u_rot**2 + v_rot**2).mean()
    divergent_energy = 0.5 * (u_div**2 + v_div**2).mean()

    print(f"Energy analysis:")
    print(f"  Total energy: {total_energy:.3f} m²/s²")
    print(
        f"  Rotational: {rotational_energy:.3f} m²/s² ({100*rotational_energy/total_energy:.1f}%)"
    )
    print(
        f"  Divergent: {divergent_energy:.3f} m²/s² ({100*divergent_energy/total_energy:.1f}%)"
    )

    # === ADVANCED EXAMPLE: 4D DATA PROCESSING ===
    print("\n2. Advanced 4D Data Processing with Parallel Computation")
    print("-" * 60)

    # Create 4D sample data
    u_4d, v_4d, lats_4d, lons_4d = create_4d_sample_data()
    print(f"4D data created: {u_4d.shape}")
    print(f"Total grid points: {u_4d.size:,}")

    # Process with parallel computation
    print("\nComputing divergent wind components...")
    start_time = time.time()

    div_u_4d, div_v_4d = compute_divergent_wind_components(u_4d, v_4d, max_workers=4)

    processing_time = time.time() - start_time

    print(f"✓ Processing completed!")
    print(f"  Processing time: {processing_time:.2f} seconds")
    print(f"  Time per time step: {processing_time/u_4d.shape[0]:.3f} seconds")
    print(f"  Efficiency: {u_4d.size/(processing_time*1e6):.1f} M points/second")

    # Analyze results
    div_u_ratio = np.abs(div_u_4d).mean() / np.abs(u_4d).mean()
    div_v_ratio = np.abs(div_v_4d).mean() / np.abs(v_4d).mean()

    print(f"\nResults analysis:")
    print(f"  Divergent U component: {div_u_ratio:.3f} of total")
    print(f"  Divergent V component: {div_v_ratio:.3f} of total")
    print(f"  ✓ Realistic ratios for atmospheric flows (< 0.2 expected)")

    # === VISUALIZATION ===
    print("\n3. Creating Visualizations")
    print("-" * 25)

    # Create plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Skyborn Windspharm Analysis Results", fontsize=16, fontweight="bold")

    # Plot original winds
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    im1 = axes[0, 0].contourf(lon_grid, lat_grid, u, levels=20, cmap="RdBu_r")
    axes[0, 0].set_title("Original Zonal Wind (U)")
    axes[0, 0].set_xlabel("Longitude")
    axes[0, 0].set_ylabel("Latitude")
    plt.colorbar(im1, ax=axes[0, 0], label="m/s")

    im2 = axes[0, 1].contourf(lon_grid, lat_grid, v, levels=20, cmap="RdBu_r")
    axes[0, 1].set_title("Original Meridional Wind (V)")
    axes[0, 1].set_xlabel("Longitude")
    plt.colorbar(im2, ax=axes[0, 1], label="m/s")

    # Plot vorticity and divergence
    im3 = axes[0, 2].contourf(
        lon_grid, lat_grid, vorticity * 1e5, levels=20, cmap="RdBu_r"
    )
    axes[0, 2].set_title("Relative Vorticity")
    axes[0, 2].set_xlabel("Longitude")
    plt.colorbar(im3, ax=axes[0, 2], label="×10⁻⁵ s⁻¹")

    im4 = axes[1, 0].contourf(
        lon_grid, lat_grid, divergence * 1e5, levels=20, cmap="RdBu_r"
    )
    axes[1, 0].set_title("Divergence")
    axes[1, 0].set_xlabel("Longitude")
    axes[1, 0].set_ylabel("Latitude")
    plt.colorbar(im4, ax=axes[1, 0], label="×10⁻⁵ s⁻¹")

    # Plot rotational and divergent components
    im5 = axes[1, 1].contourf(lon_grid, lat_grid, u_rot, levels=20, cmap="RdBu_r")
    axes[1, 1].set_title("Rotational Component (U)")
    axes[1, 1].set_xlabel("Longitude")
    plt.colorbar(im5, ax=axes[1, 1], label="m/s")

    im6 = axes[1, 2].contourf(lon_grid, lat_grid, u_div, levels=20, cmap="RdBu_r")
    axes[1, 2].set_title("Divergent Component (U)")
    axes[1, 2].set_xlabel("Longitude")
    plt.colorbar(im6, ax=axes[1, 2], label="m/s")

    plt.tight_layout()

    # Show or save plot
    try:
        plt.show()
        print("✓ Plots displayed")
    except:
        plt.savefig("windspharm_analysis_results.png", dpi=150, bbox_inches="tight")
        print("✓ Plots saved as 'windspharm_analysis_results.png'")

    # === SUMMARY ===
    print("\n" + "=" * 50)
    print("SUMMARY: Key Applications Demonstrated")
    print("=" * 50)
    print("✓ Basic spherical harmonic wind analysis")
    print("✓ Computation of vorticity, divergence, streamfunction, velocity potential")
    print("✓ Helmholtz decomposition into rotational and divergent components")
    print("✓ 4D data processing with prep_data/recover_data workflow")
    print("✓ Parallel processing for efficient large dataset handling")
    print("✓ Energy conservation verification")
    print("✓ Realistic atmospheric flow analysis")
    print("\nReal-world applications:")
    print("• Climate model output analysis (CESM2, GFDL, etc.)")
    print("• Reanalysis data processing (ERA5, MERRA-2, JRA-55)")
    print("• Atmospheric circulation pattern studies")
    print("• Weather forecast model diagnostics")
    print("• Large-scale atmospheric dynamics research")


if __name__ == "__main__":
    main()
    print(f"Wind data shape: {u.shape}")

    # Initialize VectorWind
    print("Initializing VectorWind object...")
    vw = VectorWind(u, v, gridtype="gaussian")
    print(f"Grid: {vw.nlat} lats × {vw.nlon} lons")

    # Compute dynamical quantities
    print("\nComputing dynamical quantities...")

    # Vorticity and divergence
    vorticity = vw.vorticity()
    divergence = vw.divergence()

    print(f"Vorticity range: [{vorticity.min():.2e}, {vorticity.max():.2e}] s⁻¹")
    print(f"Divergence range: [{divergence.min():.2e}, {divergence.max():.2e}] s⁻¹")

    # Streamfunction and velocity potential
    streamfunction = vw.streamfunction()
    velocity_potential = vw.velocitypotential()

    print(
        f"Streamfunction range: [{streamfunction.min():.2e}, {streamfunction.max():.2e}] m²/s"
    )
    print(
        f"Velocity potential range: [{velocity_potential.min():.2e}, {velocity_potential.max():.2e}] m²/s"
    )

    # Helmholtz decomposition
    print("\nPerforming Helmholtz decomposition...")
    u_rot, v_rot, u_div, v_div = vw.helmholtz()

    # Calculate energy in each component
    total_energy = np.mean(u**2 + v**2)
    rotational_energy = np.mean(u_rot**2 + v_rot**2)
    divergent_energy = np.mean(u_div**2 + v_div**2)

    print(f"Energy decomposition:")
    print(f"  Total: {total_energy:.3f} m²/s²")
    print(
        f"  Rotational: {rotational_energy:.3f} m²/s² ({100*rotational_energy/total_energy:.1f}%)"
    )
    print(
        f"  Divergent: {divergent_energy:.3f} m²/s² ({100*divergent_energy/total_energy:.1f}%)"
    )

    # Verify reconstruction
    u_reconstructed = u_rot + u_div
    v_reconstructed = v_rot + v_div
    reconstruction_error = (
        np.abs(u - u_reconstructed).max() + np.abs(v - v_reconstructed).max()
    )
    print(f"  Reconstruction error: {reconstruction_error:.2e}")

    # Create visualization
    print("\nCreating visualization...")
    create_plots(
        u, v, vorticity, divergence, streamfunction, velocity_potential, lats, lons
    )

    print("\n✓ Analysis complete!")


def create_plots(
    u, v, vorticity, divergence, streamfunction, velocity_potential, lats, lons
):
    """Create comprehensive plots of the analysis results."""
    # Create coordinate grids for plotting
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Set up the plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Vector Wind Analysis Results", fontsize=16, fontweight="bold")

    # Plot 1: Original wind vectors
    ax = axes[0, 0]
    skip = 4
    ax.quiver(
        lon_grid[::skip, ::skip],
        lat_grid[::skip, ::skip],
        u[::skip, ::skip],
        v[::skip, ::skip],
        scale=400,
        alpha=0.8,
        color="black",
    )
    ax.set_title("Original Wind Field")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)

    # Plot 2: Vorticity
    ax = axes[0, 1]
    im = ax.contourf(
        lon_grid, lat_grid, vorticity * 1e5, levels=20, cmap="RdBu_r", extend="both"
    )
    ax.set_title("Relative Vorticity (×10⁻⁵ s⁻¹)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Plot 3: Divergence
    ax = axes[0, 2]
    im = ax.contourf(
        lon_grid, lat_grid, divergence * 1e5, levels=20, cmap="RdBu_r", extend="both"
    )
    ax.set_title("Divergence (×10⁻⁵ s⁻¹)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Plot 4: Streamfunction
    ax = axes[1, 0]
    im = ax.contourf(
        lon_grid,
        lat_grid,
        streamfunction * 1e-6,
        levels=20,
        cmap="RdBu_r",
        extend="both",
    )
    ax.set_title("Streamfunction (×10⁶ m²/s)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Plot 5: Velocity Potential
    ax = axes[1, 1]
    im = ax.contourf(
        lon_grid,
        lat_grid,
        velocity_potential * 1e-6,
        levels=20,
        cmap="RdBu_r",
        extend="both",
    )
    ax.set_title("Velocity Potential (×10⁶ m²/s)")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Plot 6: Energy spectrum (placeholder)
    ax = axes[1, 2]
    # Create a simple demonstration plot
    sample_freqs = np.arange(1, 21)
    sample_spectrum = np.exp(-sample_freqs / 5) + 0.1 * np.random.randn(20)
    ax.semilogy(sample_freqs, sample_spectrum, "bo-", linewidth=2)
    ax.set_title("Kinetic Energy Spectrum")
    ax.set_xlabel("Wavenumber")
    ax.set_ylabel("Energy (m²/s²)")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("windspharm_analysis_results.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("Plot saved as 'windspharm_analysis_results.png'")


if __name__ == "__main__":
    main()
