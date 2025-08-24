"""
ECS (Equilibrium Climate Sensitivity) Emergent Constraint Analysis

This script analyzes the ECS data using emergent constraint methods.
ECS represents the equilibrium global surface temperature change following
a doubling of atmospheric CO2 concentration.

IPCC best estimate: ECS = 3.0°C

Data source: hot_model_ECS.xlsx
Method reference: https://github.com/blackcata/Emergent_Constraints/tree/master
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from skyborn.calc import (
    emergent_constraint_posterior,
    gaussian_pdf,
    pearson_correlation,
)

# Add Skyborn to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))


def load_ecs_data(filepath: str) -> pd.DataFrame:
    """
    Load ECS data from Excel file.

    Parameters
    ----------
    filepath : str
        Path to the Excel file containing ECS data

    Returns
    -------
    pd.DataFrame
        DataFrame with model names and ECS values
    """
    try:
        # Read the Excel file
        df = pd.read_excel(filepath)
        print(f"✅ Successfully loaded ECS data from {filepath}")
        print(f"📊 Data shape: {df.shape}")
        print(f"📋 Columns: {list(df.columns)}")

        return df
    except Exception as e:
        print(f"❌ Error loading ECS data: {e}")
        return None


def analyze_ecs_data():
    """
    Perform emergent constraint analysis on ECS data.
    """
    # Set up plotting
    plt.style.use("default")
    sns.set_palette("husl")
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "figure.figsize": (16, 12),
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 16,
        }
    )

    # Load ECS data
    ecs_file = os.path.join("DATA", "hot_model_ECS.xlsx")
    df = load_ecs_data(ecs_file)

    if df is None:
        return

    # Display first few rows to understand the structure
    print("\n📋 First few rows of the data:")
    print(df.head())

    # Extract ECS values (assuming they're in columns like 'ECS150', 'ECS130', etc.)
    ecs_columns = [col for col in df.columns if "ECS" in col.upper()]
    print(f"\n🔍 Found ECS columns: {ecs_columns}")

    if not ecs_columns:
        print("❌ No ECS columns found in the data")
        return

    # Use the first ECS column for demonstration
    ecs_col = ecs_columns[0]
    print(f"📊 Using column: {ecs_col}")

    # Clean the data - remove NaN values
    ecs_data = df[ecs_col].dropna()
    model_names = (
        df.loc[ecs_data.index, "Model Name"]
        if "Model Name" in df.columns
        else df.loc[ecs_data.index].index
    )

    print(f"\n📈 ECS Statistics:")
    print(f"   Number of models: {len(ecs_data)}")
    print(f"   Mean ECS: {ecs_data.mean():.2f}°C")
    print(f"   Std ECS: {ecs_data.std():.2f}°C")
    print(f"   Range: {ecs_data.min():.2f} - {ecs_data.max():.2f}°C")

    # IPCC best estimate
    ipcc_best_estimate = 3.0
    ipcc_uncertainty = 0.5  # Approximate uncertainty

    print(f"\n🎯 IPCC Best Estimate: {ipcc_best_estimate}°C ± {ipcc_uncertainty}°C")

    # For demonstration, we'll create a synthetic constraint variable
    # In real analysis, this would be an observable quantity correlated with ECS
    np.random.seed(42)
    n_models = len(ecs_data)

    # Create synthetic constraint variable (e.g., cloud feedback parameter)
    # This simulates a present-day observable that correlates with ECS
    constraint_var = np.random.normal(0.5, 0.2, n_models)

    # Add some correlation with ECS for realistic emergent constraint
    correlation_strength = 0.7
    noise_level = np.sqrt(1 - correlation_strength**2)

    # Normalize ECS for correlation calculation
    ecs_normalized = (ecs_data - ecs_data.mean()) / ecs_data.std()
    constraint_var = (
        correlation_strength * ecs_normalized
        + noise_level * np.random.normal(0, 1, n_models)
    )

    print(
        f"\n🔗 Synthetic constraint correlation with ECS: {pearson_correlation(constraint_var, ecs_data.values):.3f}"
    )

    # Create grids for PDF calculations
    constraint_grid = np.linspace(
        constraint_var.min() - 0.5, constraint_var.max() + 0.5, 80
    )
    ecs_grid = np.linspace(1.5, 6.0, 80)

    # Observational constraint PDF (synthetic observation)
    obs_mean = constraint_var.mean()  # Simulated observed value
    obs_std = 0.1  # Observational uncertainty
    obs_pdf = gaussian_pdf(obs_mean, obs_std, constraint_grid)

    # Convert to xarray for compatibility
    import xarray as xr

    constraint_data = xr.DataArray(constraint_var, dims=["model"])
    ecs_data_xr = xr.DataArray(ecs_data.values, dims=["model"])

    # Apply emergent constraint
    posterior_pdf, posterior_std, posterior_mean = emergent_constraint_posterior(
        constraint_data, ecs_data_xr, constraint_grid, ecs_grid, obs_pdf
    )

    # Calculate uncertainty reduction
    unconstrained_std = ecs_data.std()
    uncertainty_reduction = (1 - posterior_std / unconstrained_std) * 100

    # Create comprehensive visualization
    fig = plt.figure(figsize=(20, 15))

    # Subplot 1: ECS distribution and IPCC comparison
    ax1 = plt.subplot(2, 3, 1)

    # Histogram of model ECS values
    ax1.hist(
        ecs_data,
        bins=15,
        alpha=0.7,
        color="skyblue",
        edgecolor="black",
        label=f"Climate Models (n={len(ecs_data)})",
        density=True,
    )

    # IPCC estimate
    ax1.axvline(
        ipcc_best_estimate,
        color="red",
        linewidth=3,
        label=f"IPCC Best Estimate: {ipcc_best_estimate}°C",
    )
    ax1.axvspan(
        ipcc_best_estimate - ipcc_uncertainty,
        ipcc_best_estimate + ipcc_uncertainty,
        alpha=0.3,
        color="red",
        label=f"IPCC Range: ±{ipcc_uncertainty}°C",
    )

    # Model mean
    ax1.axvline(
        ecs_data.mean(),
        color="blue",
        linestyle="--",
        label=f"Model Mean: {ecs_data.mean():.2f}°C",
    )

    ax1.set_xlabel("ECS (°C)", fontsize=12)
    ax1.set_ylabel("Probability Density", fontsize=12)
    ax1.set_title("ECS Distribution: Models vs IPCC", fontsize=14, fontweight="bold")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Subplot 2: Constraint relationship
    ax2 = plt.subplot(2, 3, 2)

    scatter = ax2.scatter(
        constraint_var,
        ecs_data,
        c=range(len(ecs_data)),
        cmap="viridis",
        s=60,
        alpha=0.7,
        edgecolors="black",
    )

    # Regression line
    z = np.polyfit(constraint_var, ecs_data, 1)
    p = np.poly1d(z)
    ax2.plot(
        constraint_grid,
        p(constraint_grid),
        "r--",
        alpha=0.8,
        linewidth=2,
        label=f"R² = {pearson_correlation(constraint_var, ecs_data.values)**2:.3f}",
    )

    # Observational constraint
    ax2.axvline(
        obs_mean,
        color="orange",
        linewidth=3,
        alpha=0.8,
        label=f"Observation: {obs_mean:.2f}±{obs_std:.2f}",
    )
    ax2.fill_betweenx(
        [ecs_data.min(), ecs_data.max()],
        obs_mean - obs_std,
        obs_mean + obs_std,
        alpha=0.3,
        color="orange",
    )

    ax2.set_xlabel("Constraint Variable (Observable)", fontsize=12)
    ax2.set_ylabel("ECS (°C)", fontsize=12)
    ax2.set_title("Emergent Relationship", fontsize=14, fontweight="bold")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Subplot 3: Observational constraint PDF
    ax3 = plt.subplot(2, 3, 3)
    ax3.plot(constraint_grid, obs_pdf, "orange", linewidth=3, label="Observational PDF")
    ax3.fill_between(constraint_grid, obs_pdf, alpha=0.3, color="orange")
    ax3.axvline(obs_mean, color="red", linestyle="--", alpha=0.8)
    ax3.set_xlabel("Constraint Variable", fontsize=12)
    ax3.set_ylabel("Probability Density", fontsize=12)
    ax3.set_title("Observational Constraint", fontsize=14, fontweight="bold")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Subplot 4: Before vs After constraint
    ax4 = plt.subplot(2, 3, 4)

    # Unconstrained distribution
    ecs_dense = np.linspace(ecs_grid.min(), ecs_grid.max(), 200)
    unconstrained_pdf = gaussian_pdf(ecs_data.mean(), unconstrained_std, ecs_dense)

    ax4.plot(
        ecs_dense,
        unconstrained_pdf,
        "blue",
        linewidth=2,
        alpha=0.7,
        label=f"Unconstrained\n(μ={ecs_data.mean():.2f}, σ={unconstrained_std:.2f})",
    )
    ax4.plot(
        ecs_grid,
        posterior_pdf / posterior_pdf.max() * unconstrained_pdf.max(),
        "red",
        linewidth=3,
        label=f"Constrained\n(μ={posterior_mean:.2f}, σ={posterior_std:.2f})",
    )

    ax4.fill_between(ecs_dense, unconstrained_pdf, alpha=0.3, color="blue")
    ax4.fill_between(
        ecs_grid,
        posterior_pdf / posterior_pdf.max() * unconstrained_pdf.max(),
        alpha=0.4,
        color="red",
    )

    # IPCC estimate
    ax4.axvline(
        ipcc_best_estimate,
        color="green",
        linewidth=2,
        linestyle=":",
        label=f"IPCC: {ipcc_best_estimate}°C",
    )

    ax4.set_xlabel("ECS (°C)", fontsize=12)
    ax4.set_ylabel("Normalized Probability", fontsize=12)
    ax4.set_title("Constraint Effect on ECS", fontsize=14, fontweight="bold")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # Subplot 5: Uncertainty reduction
    ax5 = plt.subplot(2, 3, 5)

    categories = ["Unconstrained\nModels", "Constrained\nEstimate", "IPCC\nEstimate"]
    means = [ecs_data.mean(), posterior_mean, ipcc_best_estimate]
    stds = [unconstrained_std, posterior_std, ipcc_uncertainty]
    colors = ["blue", "red", "green"]

    bars = ax5.bar(
        categories,
        means,
        yerr=stds,
        capsize=5,
        color=colors,
        alpha=0.7,
        edgecolor="black",
    )

    for i, (mean, std) in enumerate(zip(means, stds)):
        ax5.annotate(
            f"{mean:.2f}±{std:.2f}",
            xy=(i, mean),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontweight="bold",
        )

    ax5.set_ylabel("ECS (°C)", fontsize=12)
    ax5.set_title(
        f"Uncertainty Reduction: {uncertainty_reduction:.1f}%",
        fontsize=14,
        fontweight="bold",
    )
    ax5.grid(True, alpha=0.3, axis="y")

    # Subplot 6: Summary statistics
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis("off")

    # Model statistics
    highest_ecs_idx = ecs_data.idxmax()
    lowest_ecs_idx = ecs_data.idxmin()

    summary_text = f"""
    ECS ANALYSIS SUMMARY

    Dataset: {len(ecs_data)} Climate Models
    Method: Emergent Constraints
    IPCC Best Estimate: {ipcc_best_estimate}°C

    MODEL RESULTS:
    • Original Mean: {ecs_data.mean():.2f} ± {unconstrained_std:.2f}°C
    • Constrained Mean: {posterior_mean:.2f} ± {posterior_std:.2f}°C
    • Uncertainty Reduction: {uncertainty_reduction:.1f}%

    MODEL RANGE:
    • Lowest ECS: {ecs_data.min():.2f}°C
    • Highest ECS: {ecs_data.max():.2f}°C
    • IPCC Compatibility: {'Compatible' if abs(posterior_mean - ipcc_best_estimate) < 0.5 else 'Check'}

    Reference:
    github.com/blackcata/Emergent_Constraints
    """

    ax6.text(
        0.05,
        0.95,
        summary_text,
        transform=ax6.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.8),
    )

    plt.tight_layout()
    plt.suptitle(
        "ECS Emergent Constraint Analysis Dashboard",
        fontsize=18,
        fontweight="bold",
        y=0.98,
    )
    plt.subplots_adjust(top=0.93)

    # Save the figure
    output_path = os.path.join(
        "docs", "source", "images", "ecs_emergent_constraints.png"
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")

    plt.show()

    # Print detailed results
    print(f"\nECS Emergent Constraint Analysis Complete!")
    print(f"Original ECS range: {ecs_data.mean():.2f} ± {unconstrained_std:.2f}°C")
    print(f"Constrained ECS: {posterior_mean:.2f} ± {posterior_std:.2f}°C")
    print(f"Uncertainty reduced by: {uncertainty_reduction:.1f}%")
    print(
        f"IPCC comparison: Constrained = {posterior_mean:.2f}°C, IPCC = {ipcc_best_estimate}°C"
    )
    print(f"Difference from IPCC: {abs(posterior_mean - ipcc_best_estimate):.2f}°C")

    return {
        "original_mean": ecs_data.mean(),
        "original_std": unconstrained_std,
        "constrained_mean": posterior_mean,
        "constrained_std": posterior_std,
        "uncertainty_reduction": uncertainty_reduction,
        "ipcc_estimate": ipcc_best_estimate,
    }


if __name__ == "__main__":
    print("Starting ECS Emergent Constraint Analysis...")
    results = analyze_ecs_data()
