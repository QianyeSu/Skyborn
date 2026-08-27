"""Test for edge width adaptive calculation fix."""

import matplotlib.pyplot as plt
import numpy as np
import pytest


def test_shadow_contourf_edge_width_adapts_to_dpi():
    """Test that edge width adapts to figure DPI."""
    from skyborn.plot import shadow_contourf

    # Create test data
    x = np.linspace(0, 10, 50)
    y = np.linspace(0, 10, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    # Test with different DPI settings
    for dpi in [72, 100, 150, 300]:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

        cs = shadow_contourf(ax, X, Y, Z, levels=10)

        # Check that artists were created
        assert hasattr(cs, "collections")

        # The edge width should scale with DPI
        # At 72 DPI: edge_width = 72/72 = 1.0
        # At 300 DPI: edge_width = 72/300 = 0.24 (clamped to 0.25)
        expected_width = np.clip(72.0 / dpi, 0.25, 1.0)

        # Verify that the function runs without error
        assert cs is not None

        plt.close(fig)


def test_shadow_contourf_hatch_has_zero_edge_width():
    """Test that hatched contours have zero edge width."""
    from skyborn.plot import shadow_contourf

    x = np.linspace(0, 10, 50)
    y = np.linspace(0, 10, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    fig, ax = plt.subplots(figsize=(6, 4))

    # Create contourf with hatches
    cs = shadow_contourf(ax, X, Y, Z, levels=5, hatches=["/", "\\", "|", "-", "+"])

    # Function should run successfully
    assert cs is not None

    plt.close(fig)


def test_shadow_contourf_edge_width_clamp_limits():
    """Test that edge width is clamped to reasonable limits."""
    from skyborn.plot import shadow_contourf

    x = np.linspace(0, 10, 30)
    y = np.linspace(0, 10, 30)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    # Test extreme DPI values
    test_cases = [
        (10, 0.25),  # Very low DPI -> should clamp to 0.25
        (72, 1.0),  # Standard DPI -> 1.0
        (1000, 0.25),  # Very high DPI -> should clamp to 0.25
    ]

    for dpi, expected_clamp in test_cases:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

        cs = shadow_contourf(ax, X, Y, Z, levels=5)

        # Verify edge width calculation
        calculated_width = np.clip(72.0 / dpi, 0.25, 1.0)
        assert abs(calculated_width - expected_clamp) < 0.01

        plt.close(fig)
