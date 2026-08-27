"""Tests for arrow_contour_clabel exception handling."""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest


def test_arrow_contour_clabel_with_valid_contours():
    """Test that arrow_contour_clabel works with normal contours."""
    from skyborn.plot.contour import arrow_contour, arrow_contour_clabel

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    cs = arrow_contour(X, Y, Z, levels=5, ax=ax)

    # Should not raise exception
    labels = arrow_contour_clabel(cs)

    # clabel returns list of text objects, not dict
    assert isinstance(labels, list)
    assert len(labels) > 0
    plt.close(fig)


def test_arrow_contour_clabel_with_degenerate_geometry():
    """Test that arrow_contour_clabel handles degenerate geometry gracefully."""
    from skyborn.plot.contour import arrow_contour, arrow_contour_clabel

    fig, ax = plt.subplots()

    # Create data that might produce degenerate contours
    x = np.linspace(0, 1, 20)
    y = np.linspace(0, 1, 20)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)  # Flat field - no contours
    Z[10, 10] = 1.0  # Single spike

    cs = arrow_contour(X, Y, Z, levels=[0.5], ax=ax)

    # Should handle gracefully, possibly with warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", RuntimeWarning)
        labels = arrow_contour_clabel(cs)

        # clabel returns list of text objects
        assert isinstance(labels, list)

        # May have warnings but should not crash
        if len(w) > 0:
            assert any(
                "Failed to compute label position" in str(warning.message)
                for warning in w
                if issubclass(warning.category, RuntimeWarning)
            )

    plt.close(fig)


def test_arrow_contour_clabel_without_arrows():
    """Test that arrow_contour_clabel handles contours without arrow data."""
    from skyborn.plot.contour import arrow_contour_clabel

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    # Regular contour without arrows
    cs = ax.contour(X, Y, Z, levels=5)

    # Should handle missing arrow data gracefully
    labels = arrow_contour_clabel(cs)

    # clabel returns list of text objects
    assert isinstance(labels, list)
    plt.close(fig)


def test_arrow_contour_positions_function_directly():
    """Test the internal _arrow_contour_label_positions function."""
    from skyborn.plot._core.contour_arrows import _arrow_contour_label_positions
    from skyborn.plot.contour import arrow_contour

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    cs = arrow_contour(X, Y, Z, levels=5, ax=ax)

    # Should return list of tuples
    positions = _arrow_contour_label_positions(cs)
    assert isinstance(positions, list)

    # Each position should be a tuple of two floats
    for pos in positions:
        assert isinstance(pos, tuple)
        assert len(pos) == 2
        assert isinstance(pos[0], float)
        assert isinstance(pos[1], float)

    plt.close(fig)


def test_arrow_contour_positions_with_corrupted_artist():
    """Test that _arrow_contour_label_positions handles corrupted artist data."""
    from unittest.mock import Mock

    from skyborn.plot._core.contour_arrows import _arrow_contour_label_positions
    from skyborn.plot.contour import arrow_contour

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    cs = arrow_contour(X, Y, Z, levels=5, ax=ax)

    # Corrupt one of the arrow artists to trigger exception
    if (
        hasattr(cs, "_skyborn_arrow_contour_artists")
        and cs._skyborn_arrow_contour_artists
    ):
        # Create a mock artist that raises an exception
        bad_artist = Mock()
        bad_artist.get_segments.side_effect = RuntimeError("Simulated geometry error")

        # Insert the bad artist into the list
        cs._skyborn_arrow_contour_artists.insert(0, bad_artist)

        # Should handle gracefully with warning
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", RuntimeWarning)
            positions = _arrow_contour_label_positions(cs)

            # Should still return a list (possibly shorter)
            assert isinstance(positions, list)

            # Should have emitted a warning about the failed line
            assert len(w) > 0
            assert any(
                "Failed to compute label position" in str(warning.message)
                for warning in w
                if issubclass(warning.category, RuntimeWarning)
            )

    plt.close(fig)


def test_arrow_contour_positions_with_transform_failure():
    """Test exception handling when coordinate transform fails."""
    from unittest.mock import Mock, patch

    from skyborn.plot._core.contour_arrows import _arrow_contour_label_positions
    from skyborn.plot.contour import arrow_contour

    fig, ax = plt.subplots()
    x = np.linspace(-3, 3, 50)
    y = np.linspace(-3, 3, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y)

    cs = arrow_contour(X, Y, Z, levels=5, ax=ax)

    # Mock the transform to fail
    bad_transform = Mock()
    bad_transform.transform.side_effect = ValueError("Transform failure")
    bad_transform.inverted.return_value = bad_transform

    original_transform = cs.get_transform
    cs.get_transform = lambda: bad_transform

    # Should handle gracefully with warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always", RuntimeWarning)
        positions = _arrow_contour_label_positions(cs)

        # Should return empty list or partial list
        assert isinstance(positions, list)

        # Should have warnings
        assert len(w) > 0
        assert any(
            "Failed to compute label position" in str(warning.message)
            for warning in w
            if issubclass(warning.category, RuntimeWarning)
        )

    # Restore original transform
    cs.get_transform = original_transform
    plt.close(fig)
