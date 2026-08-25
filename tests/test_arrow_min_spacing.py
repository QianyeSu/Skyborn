"""Comprehensive test suite for arrow_min_spacing parameter."""

import matplotlib.pyplot as plt
import numpy as np
import pytest

from skyborn.plot import arrow_contour


def test_arrow_min_spacing_basic():
    """Test basic arrow_min_spacing parameter."""
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    Z = X**2 + Y**2

    fig, ax = plt.subplots()
    cs = arrow_contour(ax, X, Y, Z, levels=5, arrow_count=3, arrow_min_spacing=30)

    artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
    assert len(artists) > 0, "Should create arrow artists"
    plt.close(fig)


def test_arrow_min_spacing_large_value():
    """Test with large arrow_min_spacing that drastically reduces arrows."""
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X * 2) * np.cos(Y * 2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Default spacing
    cs1 = arrow_contour(ax1, X, Y, Z, levels=8, arrow_count=3)
    artists1 = getattr(cs1, "_skyborn_arrow_contour_artists", [])

    # Large spacing
    cs2 = arrow_contour(ax2, X, Y, Z, levels=8, arrow_count=3, arrow_min_spacing=100)
    artists2 = getattr(cs2, "_skyborn_arrow_contour_artists", [])

    # Both should produce artists
    assert len(artists1) > 0
    assert len(artists2) > 0

    plt.close(fig)


def test_arrow_min_spacing_none():
    """Test that arrow_min_spacing=None uses default behavior."""
    x = np.linspace(-2, 2, 50)
    y = np.linspace(-2, 2, 50)
    X, Y = np.meshgrid(x, y)
    Z = X * Y

    fig, ax = plt.subplots()
    cs = arrow_contour(ax, X, Y, Z, levels=6, arrow_count=2, arrow_min_spacing=None)

    artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
    assert len(artists) > 0
    plt.close(fig)


def test_arrow_min_spacing_with_different_styles():
    """Test arrow_min_spacing works with all arrow styles."""
    x = np.linspace(-2, 2, 80)
    y = np.linspace(-2, 2, 80)
    X, Y = np.meshgrid(x, y)
    Z = X**2 - Y**2

    for style in ["swept", "line", "filled"]:
        fig, ax = plt.subplots()
        cs = arrow_contour(
            ax,
            X,
            Y,
            Z,
            levels=5,
            arrow_count=2,
            arrow_min_spacing=40,
            arrow_style=style,
        )
        artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
        assert len(artists) > 0, f"Should create arrows for style={style}"
        plt.close(fig)


def test_arrow_min_spacing_validation():
    """Test that invalid arrow_min_spacing values are rejected."""
    x = np.linspace(-1, 1, 30)
    y = np.linspace(-1, 1, 30)
    X, Y = np.meshgrid(x, y)
    Z = X + Y

    fig, ax = plt.subplots()

    # Negative spacing should raise ValueError
    with pytest.raises(ValueError, match="arrow_min_spacing"):
        arrow_contour(ax, X, Y, Z, levels=3, arrow_min_spacing=-10)

    # Zero spacing should raise ValueError
    with pytest.raises(ValueError, match="arrow_min_spacing"):
        arrow_contour(ax, X, Y, Z, levels=3, arrow_min_spacing=0)

    plt.close(fig)


def test_arrow_min_spacing_dense_contours():
    """Test arrow_min_spacing with densely packed contours."""
    x = np.linspace(-3, 3, 150)
    y = np.linspace(-3, 3, 150)
    X, Y = np.meshgrid(x, y)
    # Create tight gradient
    Z = np.sin(X * 3) * np.cos(Y * 3) * 5

    fig, ax = plt.subplots()
    cs = arrow_contour(
        ax, X, Y, Z, levels=20, arrow_count=4, arrow_min_spacing=50  # Many levels
    )

    artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
    assert len(artists) > 0, "Should handle dense contours"
    plt.close(fig)


def test_arrow_min_spacing_single_arrow():
    """Test arrow_min_spacing with arrow_count=1."""
    x = np.linspace(-2, 2, 60)
    y = np.linspace(-2, 2, 60)
    X, Y = np.meshgrid(x, y)
    Z = X**2 + Y**2

    fig, ax = plt.subplots()
    cs = arrow_contour(
        ax,
        X,
        Y,
        Z,
        levels=4,
        arrow_count=1,  # Single arrow per contour
        arrow_min_spacing=30,
    )

    artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
    assert len(artists) > 0
    plt.close(fig)


def test_arrow_min_spacing_with_absolute_length():
    """Test arrow_min_spacing combined with arrow_length_points."""
    x = np.linspace(-2, 2, 80)
    y = np.linspace(-2, 2, 80)
    X, Y = np.meshgrid(x, y)
    Z = np.exp(-(X**2) - Y**2) * 10

    fig, ax = plt.subplots()
    cs = arrow_contour(
        ax,
        X,
        Y,
        Z,
        levels=6,
        arrow_count=3,
        arrow_length_points=15,  # Absolute length
        arrow_min_spacing=50,
    )

    artists = getattr(cs, "_skyborn_arrow_contour_artists", [])
    assert len(artists) > 0, "Should work with absolute arrow length"
    plt.close(fig)


if __name__ == "__main__":
    # Run tests
    test_arrow_min_spacing_basic()
    print("PASS test_arrow_min_spacing_basic")

    test_arrow_min_spacing_large_value()
    print("PASS test_arrow_min_spacing_large_value")

    test_arrow_min_spacing_none()
    print("PASS test_arrow_min_spacing_none")

    test_arrow_min_spacing_with_different_styles()
    print("PASS test_arrow_min_spacing_with_different_styles")

    test_arrow_min_spacing_validation()
    print("PASS test_arrow_min_spacing_validation")

    test_arrow_min_spacing_dense_contours()
    print("PASS test_arrow_min_spacing_dense_contours")

    test_arrow_min_spacing_single_arrow()
    print("PASS test_arrow_min_spacing_single_arrow")

    test_arrow_min_spacing_with_absolute_length()
    print("PASS test_arrow_min_spacing_with_absolute_length")

    print("\nAll tests passed!")
