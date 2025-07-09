"""
Example usage of the reorganized Skyborn calc module.

This example demonstrates how to use the emergent constraint methods
and statistical calculations in the new calc module structure.
"""

from skyborn.calc import (
    calc_GAUSSIAN_PDF,
    calc_PDF_EC_PRIOR,
    linear_regression,
    pearson_correlation,
)
import numpy as np
import xarray as xr
import sys
import os

# Add the src directory to Python path for development
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# Import from the new calc module


def main():
    """
    Example usage of emergent constraint and statistical functions.
    """
    print("🧮 Skyborn Calc Module Example")
    print("=" * 40)

    # Example 1: Gaussian PDF calculation
    print("\n1. Gaussian PDF Calculation:")
    x = np.linspace(-3, 3, 100)
    pdf = calc_GAUSSIAN_PDF(mu=0, sigma=1, x=x)
    print(f"   PDF max value: {pdf.max():.4f}")
    print(f"   PDF at x=0: {pdf[len(x)//2]:.4f}")

    # Example 2: Linear regression with synthetic data
    print("\n2. Linear Regression Example:")
    np.random.seed(42)
    n_samples = 50

    # Create synthetic 3D data
    predictor = np.random.randn(n_samples)
    data_3d = np.zeros((n_samples, 10, 10))

    for i in range(n_samples):
        # Create data with linear relationship plus noise
        data_3d[i, :, :] = 2 * predictor[i] + np.random.randn(10, 10) * 0.5

    # Perform regression
    coeff, p_values = linear_regression(data_3d, predictor)
    print(f"   Mean regression coefficient: {coeff.mean():.4f}")
    print(f"   Expected coefficient: 2.0")

    # Example 3: Correlation calculation
    print("\n3. Correlation Calculation:")
    x_data = np.random.randn(100)
    y_data = 0.8 * x_data + 0.2 * np.random.randn(100)
    correlation = pearson_correlation(x_data, y_data)
    print(f"   Pearson correlation: {correlation:.4f}")

    print("\n✅ All examples completed successfully!")


if __name__ == "__main__":
    main()
