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


def test_local_tangent_at_distance_matches_display_space_geometry():
    """The tangent API should return the normalized forward local direction."""
    from skyborn.plot._core.contour_arrows_core import local_tangent_at_distance

    vertices = np.array(
        [[0.0, 0.0], [5.0, 0.0], [6.0, 2.0], [7.0, 0.0], [12.0, 0.0]],
        dtype=np.float64,
    )
    total_length = float(np.sum(np.hypot(*np.diff(vertices, axis=0).T)))
    distance = 2.5

    tangent = local_tangent_at_distance(vertices, distance, total_length, 0.5)

    np.testing.assert_allclose(tangent, [1.0, 0.0], atol=1e-12)
    assert tangent.dtype == np.float64


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


def test_select_single_arrow_prefers_straight_section_over_midpoint():
    """A single arrow should use the best local tangent, not a fixed midpoint."""
    from skyborn.plot._core.contour_arrows_core import select_arrow_end_distances

    # The path midpoint lies on a sharp bend. Long straight sections exist on
    # both sides, so a straightness-aware selector should avoid that midpoint.
    vertices = np.array(
        [[0.0, 0.0], [5.0, 0.0], [6.0, 2.0], [7.0, 0.0], [12.0, 0.0]],
        dtype=np.float64,
    )
    total_length = float(np.sum(np.hypot(*np.diff(vertices, axis=0).T)))

    result = select_arrow_end_distances(vertices, total_length, 1, 1.5, None)

    assert result.shape == (1,)
    assert abs(float(result[0]) - total_length / 2.0) > 1.0


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
