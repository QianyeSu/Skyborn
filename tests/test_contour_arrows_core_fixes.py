"""Test fixes for contour_arrows_core dtype and memory safety issues."""

import numpy as np
import pytest


def test_point_at_distance_requires_float64():
    """Test that point_at_distance rejects non-float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import point_at_distance

    # Create float32 array (should be rejected)
    vertices_f32 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float32)

    with pytest.raises(TypeError, match="float64"):
        point_at_distance(vertices_f32, 1.0)


def test_point_at_distance_accepts_float64():
    """Test that point_at_distance accepts float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import point_at_distance

    # Create float64 array (should work)
    vertices_f64 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float64)

    result = point_at_distance(vertices_f64, 1.0)
    assert result is not None
    assert result.dtype == np.float64


def test_point_at_distance_handles_non_contiguous():
    """Test that point_at_distance handles non-contiguous arrays."""
    from skyborn.plot._core.contour_arrows_core import point_at_distance

    # Create non-contiguous array via slicing
    vertices_full = np.array(
        [[0.0, 0.0], [0.5, 0.5], [1.0, 1.0], [1.5, 0.5], [2.0, 0.0]], dtype=np.float64
    )
    vertices_strided = vertices_full[::2]  # Non-contiguous

    assert not vertices_strided.flags["C_CONTIGUOUS"]
    assert not vertices_strided.flags["F_CONTIGUOUS"]

    # Should still work by internally converting to contiguous
    result = point_at_distance(vertices_strided, 1.0)
    assert result is not None
    assert result.dtype == np.float64


def test_select_arrow_end_distances_requires_float64():
    """Test that select_arrow_end_distances rejects non-float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import select_arrow_end_distances

    vertices_f32 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float32)

    with pytest.raises(TypeError, match="float64"):
        select_arrow_end_distances(vertices_f32, 3.0, 2, 0.5, None)


def test_select_arrow_end_distances_accepts_float64():
    """Test that select_arrow_end_distances accepts float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import select_arrow_end_distances

    vertices_f64 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float64)

    result = select_arrow_end_distances(vertices_f64, 3.0, 2, 0.5, None)
    assert result is not None
    assert result.dtype == np.float64
    assert len(result) <= 2


def test_select_arrow_end_distances_handles_non_contiguous():
    """Test that select_arrow_end_distances handles non-contiguous arrays."""
    from skyborn.plot._core.contour_arrows_core import select_arrow_end_distances

    vertices_full = np.array(
        [[0.0, 0.0], [0.5, 0.5], [1.0, 1.0], [1.5, 0.5], [2.0, 0.0]], dtype=np.float64
    )
    vertices_strided = vertices_full[::2]

    assert not vertices_strided.flags["C_CONTIGUOUS"]

    result = select_arrow_end_distances(vertices_strided, 3.0, 2, 0.5, None)
    assert result is not None
    assert result.dtype == np.float64


def test_local_straightness_score_requires_float64():
    """Test that local_straightness_score rejects non-float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import local_straightness_score

    vertices_f32 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float32)

    # local_straightness_score calls point_at_distance internally,
    # so it should also fail with float32
    with pytest.raises(TypeError, match="float64"):
        local_straightness_score(vertices_f32, 1.5, 0.5)


def test_local_straightness_score_accepts_float64():
    """Test that local_straightness_score accepts float64 arrays."""
    from skyborn.plot._core.contour_arrows_core import local_straightness_score

    vertices_f64 = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]], dtype=np.float64)

    result = local_straightness_score(vertices_f64, 1.5, 0.5)
    assert isinstance(result, float)
    assert np.isfinite(result)
