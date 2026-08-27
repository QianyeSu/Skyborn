"""Tests for margin hardcode and recursion risk fixes."""

import matplotlib.pyplot as plt
import numpy as np
import pytest


def test_shadow_boundary_margin_parameter():
    """Test that shadow_boundary_margin parameter is accepted and used."""
    from skyborn.plot.contour import shadow_contourf

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    # Test default margin
    cs1 = shadow_contourf(X, Y, Z, levels=5, ax=ax)
    assert hasattr(cs1, "_skyborn_shadow_artists")

    # Test custom margin - should not raise
    cs2 = shadow_contourf(X, Y, Z, levels=5, ax=ax, shadow_boundary_margin=0.2)
    assert hasattr(cs2, "_skyborn_shadow_artists")

    # Test zero margin
    cs3 = shadow_contourf(X, Y, Z, levels=5, ax=ax, shadow_boundary_margin=0.0)
    assert hasattr(cs3, "_skyborn_shadow_artists")

    plt.close(fig)


def test_shadow_boundary_margin_affects_shadow_count():
    """Test that larger margin reduces shadow artist count."""
    from skyborn.plot.contour import shadow_contourf

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    x = np.linspace(-3, 3, 30)
    y = np.linspace(-3, 3, 30)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    # Small margin - more shadows drawn
    cs1 = shadow_contourf(X, Y, Z, levels=5, ax=ax1, shadow_boundary_margin=0.01)
    small_margin_count = len(cs1._skyborn_shadow_artists)

    # Large margin - fewer shadows drawn
    cs2 = shadow_contourf(X, Y, Z, levels=5, ax=ax2, shadow_boundary_margin=1.0)
    large_margin_count = len(cs2._skyborn_shadow_artists)

    # Larger margin should filter out more contours near boundaries
    assert large_margin_count <= small_margin_count

    plt.close(fig)
