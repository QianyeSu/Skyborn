"""
Generate documentation images for the Skyborn gallery.

This script runs the emergent constraint analysis and saves high-quality
figures for use in the documentation gallery.

Usage:
    python generate_docs_images.py

Output:
    - docs/source/images/emergent_constraints_dashboard.png
    - docs/source/images/method_comparison.png
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add src to path for development
current_dir = os.path.dirname(os.path.abspath(__file__))
skyborn_src = os.path.join(os.path.dirname(current_dir), "src")
sys.path.insert(0, skyborn_src)

try:
    import skyborn as skb
    from skyborn.calc import gaussian_pdf, pearson_correlation

    print("Using installed Skyborn functions")
except ImportError:
    print("Using fallback functions")

    def gaussian_pdf(mu, sigma, x):
        """Fallback Gaussian PDF calculation."""
        return (
            1 / np.sqrt(2 * np.pi * sigma**2) * np.exp(-1 / 2 * ((x - mu) / sigma) ** 2)
        )

    def pearson_correlation(x, y):
        """Fallback Pearson correlation."""
        return np.corrcoef(x.flatten(), y.flatten())[0, 1]


def setup_plotting():
    """Configure matplotlib for high-quality documentation figures."""
    plt.style.use("default")
    sns.set_style("whitegrid")
    sns.set_palette("husl")

    # High-quality figure settings
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "figure.figsize": (20, 15),
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 16,
            "legend.fontsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
            "font.family": "DejaVu Sans",
        }
    )


def generate_dashboard_image():
    """Generate the main emergent constraints dashboard image."""
    setup_plotting()

    # Generate synthetic ECS data
    np.random.seed(42)
    n_models = 30

    # Generate ECS values with realistic distribution
    ecs_values = np.random.normal(3.2, 0.9, n_models)
    ecs_values = np.clip(ecs_values, 1.5, 5.5)  # Realistic ECS range

    # Add some models with higher/lower ECS for diversity
    ecs_values[0:3] = [2.1, 4.8, 5.2]  # Add extreme values

    model_names = np.array([f"CMIP_Model_{i+1:02d}" for i in range(n_models)])

    # Create synthetic constraint variable
    np.random.seed(123)
    constraint_strength = 0.7  # Correlation strength
    noise_level = np.sqrt(1 - constraint_strength**2)

    # Normalize ECS for constraint generation
    ecs_normalized = (ecs_values - ecs_values.mean()) / ecs_values.std()
    constraint_var = (
        constraint_strength * ecs_normalized + noise_level * np.random.randn(n_models)
    )

    # Add some realistic offset and scaling
    constraint_var = constraint_var * 0.3 + 0.5  # Scale to reasonable physical range

    # Calculate correlation
    correlation = pearson_correlation(constraint_var, ecs_values)

    # Set up grids for PDF calculations
    constraint_grid = np.linspace(
        constraint_var.min() - 0.3, constraint_var.max() + 0.3, 80
    )
    ecs_grid = np.linspace(1.5, 5.5, 80)

    # Observational constraint (synthetic but realistic)
    obs_mean = constraint_var.mean() + 0.05  # Slight offset from model mean
    obs_std = 0.08  # Observational uncertainty
    obs_pdf = gaussian_pdf(obs_mean, obs_std, constraint_grid)

    # Apply emergent constraint using simplified method
    # Linear regression
    slope, intercept = np.polyfit(constraint_var, ecs_values, 1)
    predicted_ecs = slope * constraint_var + intercept
    residuals = ecs_values - predicted_ecs
    prediction_error = np.std(residuals)

    # Calculate constrained distribution
    constrained_mean = slope * obs_mean + intercept
    constrained_std = prediction_error * obs_std / np.std(constraint_var)

    # Calculate uncertainty reduction
    original_std = ecs_values.std()
    uncertainty_reduction = (1 - constrained_std / original_std) * 100

    # IPCC AR6 best estimate and likely range
    ipcc_best_estimate = 3.0
    ipcc_likely_range = (2.5, 4.0)
    ipcc_uncertainty = (ipcc_likely_range[1] - ipcc_likely_range[0]) / 4

    # Create comprehensive visualization dashboard
    fig, axes = plt.subplots(2, 3, figsize=(20, 14))
    axes = axes.flatten()

    # Color scheme (no emojis)
    colors = {
        "models": "skyblue",
        "constrained": "red",
        "ipcc": "green",
        "observation": "orange",
    }

    # 1. ECS Distribution Comparison
    ax = axes[0]
    ax.hist(
        ecs_values,
        bins=15,
        alpha=0.7,
        color=colors["models"],
        density=True,
        label=f"Climate Models (n={n_models})",
        edgecolor="black",
        linewidth=0.5,
    )

    # Add IPCC reference
    ax.axvline(
        ipcc_best_estimate,
        color=colors["ipcc"],
        linewidth=3,
        label=f"IPCC Best Estimate: {ipcc_best_estimate}°C",
    )
    ax.axvspan(
        ipcc_likely_range[0],
        ipcc_likely_range[1],
        alpha=0.3,
        color=colors["ipcc"],
        label=f"IPCC Likely Range: {ipcc_likely_range[0]}-{ipcc_likely_range[1]}°C",
    )

    # Add model mean
    ax.axvline(
        ecs_values.mean(),
        color="blue",
        linestyle="--",
        linewidth=2,
        label=f"Model Mean: {ecs_values.mean():.2f}°C",
    )

    ax.set_xlabel("Equilibrium Climate Sensitivity (°C)", fontsize=14)
    ax.set_ylabel("Probability Density", fontsize=14)
    ax.set_title("ECS Distribution in Climate Models", fontsize=16, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # 2. Emergent Relationship
    ax = axes[1]
    scatter = ax.scatter(
        constraint_var,
        ecs_values,
        c=range(n_models),
        cmap="viridis",
        s=80,
        alpha=0.8,
        edgecolors="black",
        linewidth=0.5,
    )

    # Add regression line
    ax.plot(
        constraint_grid,
        slope * constraint_grid + intercept,
        "r--",
        linewidth=3,
        label=f"Linear Fit: R² = {correlation**2:.3f}",
    )

    # Add observational constraint
    ax.axvline(
        obs_mean,
        color=colors["observation"],
        linewidth=3,
        label=f"Observation: {obs_mean:.3f}±{obs_std:.3f}",
    )
    ax.fill_betweenx(
        [ecs_values.min() - 0.5, ecs_values.max() + 0.5],
        obs_mean - obs_std,
        obs_mean + obs_std,
        alpha=0.3,
        color=colors["observation"],
    )

    ax.set_xlabel("Constraint Variable", fontsize=14)
    ax.set_ylabel("Equilibrium Climate Sensitivity (°C)", fontsize=14)
    ax.set_title("Emergent Relationship", fontsize=16, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label("Model Index", fontsize=12)

    # 3. Observational Constraint PDF
    ax = axes[2]
    ax.plot(
        constraint_grid,
        obs_pdf,
        color=colors["observation"],
        linewidth=4,
        label="Observational PDF",
    )
    ax.fill_between(constraint_grid, obs_pdf, alpha=0.4, color=colors["observation"])
    ax.axvline(obs_mean, color="darkred", linestyle="--", linewidth=2)

    ax.set_xlabel("Constraint Variable", fontsize=14)
    ax.set_ylabel("Probability Density", fontsize=14)
    ax.set_title("Observational Constraint", fontsize=16, fontweight="bold")
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    # 4. Before and After Constraint Comparison
    ax = axes[3]
    ecs_dense = np.linspace(1.0, 6.0, 200)

    # Original distribution
    original_pdf = gaussian_pdf(ecs_values.mean(), original_std, ecs_dense)
    ax.plot(
        ecs_dense,
        original_pdf,
        "b-",
        linewidth=3,
        alpha=0.8,
        label=f"Unconstrained: {ecs_values.mean():.2f}±{original_std:.2f}°C",
    )
    ax.fill_between(ecs_dense, original_pdf, alpha=0.3, color="blue")

    # Constrained distribution
    constrained_pdf = gaussian_pdf(constrained_mean, constrained_std, ecs_dense)
    ax.plot(
        ecs_dense,
        constrained_pdf,
        color=colors["constrained"],
        linewidth=4,
        label=f"Constrained: {constrained_mean:.2f}±{constrained_std:.2f}°C",
    )
    ax.fill_between(ecs_dense, constrained_pdf, alpha=0.4, color=colors["constrained"])

    # IPCC reference
    ax.axvline(
        ipcc_best_estimate,
        color=colors["ipcc"],
        linewidth=3,
        linestyle=":",
        label=f"IPCC: {ipcc_best_estimate}°C",
    )

    ax.set_xlabel("Equilibrium Climate Sensitivity (°C)", fontsize=14)
    ax.set_ylabel("Probability Density", fontsize=14)
    ax.set_title("Constraint Effect", fontsize=16, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # 5. Uncertainty Reduction Bar Chart
    ax = axes[4]
    categories = ["Unconstrained", "Constrained", "IPCC\nBest Estimate"]
    means = [ecs_values.mean(), constrained_mean, ipcc_best_estimate]
    stds = [original_std, constrained_std, ipcc_uncertainty]
    bar_colors = ["blue", colors["constrained"], colors["ipcc"]]

    bars = ax.bar(
        categories,
        means,
        yerr=stds,
        capsize=8,
        color=bar_colors,
        alpha=0.7,
        edgecolor="black",
        linewidth=1,
    )

    # Add value annotations
    for i, (mean, std) in enumerate(zip(means, stds)):
        ax.annotate(
            f"{mean:.2f}±{std:.2f}°C",
            xy=(i, mean),
            xytext=(0, 15),
            textcoords="offset points",
            ha="center",
            fontweight="bold",
            fontsize=11,
        )

    ax.set_ylabel("Equilibrium Climate Sensitivity (°C)", fontsize=14)
    ax.set_title(
        f"Uncertainty Reduction: {uncertainty_reduction:.1f}%",
        fontsize=16,
        fontweight="bold",
    )
    ax.grid(True, alpha=0.3, axis="y")

    # 6. Summary Statistics Panel
    ax = axes[5]
    ax.axis("off")

    summary_text = f"""ECS EMERGENT CONSTRAINT ANALYSIS

DATA:
• Climate Models: {n_models}
• ECS Range: {ecs_values.min():.2f} - {ecs_values.max():.2f}°C

EMERGENT RELATIONSHIP:
• Correlation (r): {correlation:.3f}
• R-squared: {correlation**2:.3f}
• Regression slope: {slope:.2f}

OBSERVATIONAL CONSTRAINT:
• Observed value: {obs_mean:.3f} ± {obs_std:.3f}
• Constraint strength: Strong

RESULTS:
• Original ECS: {ecs_values.mean():.2f} ± {original_std:.2f}°C
• Constrained ECS: {constrained_mean:.2f} ± {constrained_std:.2f}°C
• Uncertainty reduction: {uncertainty_reduction:.1f}%

COMPARISON WITH IPCC:
• IPCC best estimate: {ipcc_best_estimate}°C
• Difference: {abs(constrained_mean - ipcc_best_estimate):.2f}°C
• Agreement: {'Excellent' if abs(constrained_mean - ipcc_best_estimate) < 0.3 else 'Good' if abs(constrained_mean - ipcc_best_estimate) < 0.5 else 'Moderate'}

REFERENCE:
github.com/blackcata/Emergent_Constraints
"""

    ax.text(
        0.05,
        0.95,
        summary_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8),
    )

    plt.tight_layout()
    plt.suptitle(
        "ECS Emergent Constraints Analysis Dashboard",
        fontsize=20,
        fontweight="bold",
        y=0.98,
    )
    plt.subplots_adjust(top=0.93)

    # Save the figure
    output_dir = Path("docs/source/images")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "emergent_constraints_dashboard.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")

    print(f"Generated: {output_file}")
    plt.close()


def generate_method_comparison():
    """Generate a comparison figure showing methodology improvements."""
    setup_plotting()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Left panel: Traditional approach
    ax1.text(
        0.5,
        0.9,
        "Traditional Approach",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="navy",
        transform=ax1.transAxes,
    )

    traditional_steps = [
        "1. Use all model data equally",
        "2. Calculate simple ensemble mean",
        "3. Large uncertainty range",
        "4. No observational constraints",
        "5. High projection uncertainty",
    ]

    for i, step in enumerate(traditional_steps):
        ax1.text(0.1, 0.75 - i * 0.12, step, fontsize=15, transform=ax1.transAxes)

    ax1.text(
        0.5,
        0.1,
        "Result: Wide uncertainty range\nLimited confidence in projections",
        ha="center",
        va="center",
        fontsize=15,
        color="darkred",
        transform=ax1.transAxes,
        style="italic",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="mistyrose"),
    )

    # Right panel: Emergent constraint approach
    ax2.text(
        0.5,
        0.9,
        "Emergent Constraint Approach",
        ha="center",
        va="center",
        fontsize=18,
        fontweight="bold",
        color="darkgreen",
        transform=ax2.transAxes,
    )

    emergent_steps = [
        "1. Find inter-model relationships",
        "2. Apply observational constraints",
        "3. Weight models by performance",
        "4. Reduce projection uncertainty",
        "5. Improved confidence bounds",
    ]

    for i, step in enumerate(emergent_steps):
        ax2.text(0.1, 0.75 - i * 0.12, step, fontsize=15, transform=ax2.transAxes)

    ax2.text(
        0.5,
        0.1,
        "Result: Narrower uncertainty range\nHigher confidence in projections",
        ha="center",
        va="center",
        fontsize=15,
        color="darkgreen",
        transform=ax2.transAxes,
        style="italic",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen"),
    )

    # Remove axes
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Add frames
    for spine in ax1.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2)
    for spine in ax2.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2)

    plt.tight_layout()
    # plt.suptitle("Climate Projection Methods Comparison",
    #              fontsize=20, fontweight="bold", y=1.05)

    # Save the figure
    output_dir = Path("docs/source/images")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "method_comparison.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")

    print(f"Generated: {output_file}")
    plt.close()


def main():
    """Generate all documentation images."""
    print("Generating documentation images...")

    try:
        generate_dashboard_image()
        generate_method_comparison()
        print("\nAll documentation images generated successfully!")

    except Exception as e:
        print(f"Error generating images: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
