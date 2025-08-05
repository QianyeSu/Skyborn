#!/usr/bin/env python3
"""
Windspharm Example: Basic Vector Wind Analysis

This example demonstrates basic usage of the skyborn.windspharm module
for spherical harmonic vector wind analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from skyborn.windspharm import VectorWind
from skyborn.spharm import gaussian_lats_wts


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


def main():
    """Main analysis function."""
    print("Skyborn Windspharm Example")
    print("=" * 30)

    # Create sample data
    print("Creating sample wind data...")
    u, v, lats, lons = create_sample_wind_data()
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
