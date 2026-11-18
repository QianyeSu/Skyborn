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

        # Check that layered artists were created without deprecated .collections
        assert hasattr(cs, "_skyborn_shadow_artists")
        assert len(cs._skyborn_shadow_artists) > 0

        # The edge width should scale with DPI
        # At 72 DPI: edge_width = 72/72 = 1.0
        # At 300 DPI: edge_width = 72/300 = 0.24 (clamped to 0.25)
        expected_width = np.clip(72.0 / dpi, 0.25, 1.0)
        fill_collections = [
            artist
            for artist in cs._skyborn_shadow_artists
            if not artist.get_antialiaseds()[0]
        ]
        assert len(fill_collections) > 0
        for col in fill_collections:
            assert col.get_linewidths()[0] == pytest.approx(expected_width)

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

    # Function should run successfully with zero edge width on hatched fills
    assert cs is not None
    fill_collections = [
        artist
        for artist in cs._skyborn_shadow_artists
        if not artist.get_antialiaseds()[0]
    ]
    assert len(fill_collections) > 0
    for col in fill_collections:
        assert col.get_linewidths()[0] == pytest.approx(0.0)

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
        (10, 1.0),  # Very low DPI -> 72/10=7.2 clamps to 1.0
        (72, 1.0),  # Standard DPI -> 1.0
        (1000, 0.25),  # Very high DPI -> 72/1000=0.072 clamps to 0.25
    ]

    for dpi, expected_clamp in test_cases:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)

        cs = shadow_contourf(ax, X, Y, Z, levels=5)

        # Verify edge width calculation and rendered fill collection widths
        calculated_width = np.clip(72.0 / dpi, 0.25, 1.0)
        assert abs(calculated_width - expected_clamp) < 0.01
        fill_collections = [
            artist
            for artist in cs._skyborn_shadow_artists
            if not artist.get_antialiaseds()[0]
        ]
        assert len(fill_collections) > 0
        for col in fill_collections:
            assert col.get_linewidths()[0] == pytest.approx(expected_clamp)

        plt.close(fig)
