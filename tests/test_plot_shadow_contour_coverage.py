"""
Add explicit tests for shadow_contourf helper functions to achieve 100% coverage.
"""

import matplotlib

matplotlib.use("Agg")

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pytest

from skyborn.plot.contour import (
    _normalize_shadow_backend,
    _validate_shadow_offset,
    shadow_contourf,
)


def test_validate_shadow_offset_valid():
    """Test _validate_shadow_offset with valid input."""
    result = _validate_shadow_offset((2.0, -3.0))
    assert result == (2.0, -3.0)

    result = _validate_shadow_offset([1.5, -2.5])
    assert result == (1.5, -2.5)


def test_validate_shadow_offset_invalid():
    """Test _validate_shadow_offset with invalid input."""
    with pytest.raises(ValueError, match="shadow_offset must be a two-item"):
        _validate_shadow_offset((1.0,))

    with pytest.raises(ValueError, match="shadow_offset must be a two-item"):
        _validate_shadow_offset((1.0, 2.0, 3.0))

    with pytest.raises(ValueError, match="shadow_offset must be a two-item"):
        _validate_shadow_offset(1.0)


def test_normalize_shadow_backend_valid():
    """Test _normalize_shadow_backend with valid aliases."""
    assert _normalize_shadow_backend("standard", "test") == "standard"
    assert _normalize_shadow_backend("matplotlib", "test") == "standard"
    assert _normalize_shadow_backend("fast", "test") == "fast"
    assert _normalize_shadow_backend("contourpy", "test") == "fast"
    assert _normalize_shadow_backend("auto", "test") == "auto"

    # Test case insensitive
    assert _normalize_shadow_backend("STANDARD", "test") == "standard"
    assert _normalize_shadow_backend("  Fast  ", "test") == "fast"


def test_normalize_shadow_backend_invalid():
    """Test _normalize_shadow_backend with invalid input."""
    with pytest.raises(ValueError, match="test must be one of"):
        _normalize_shadow_backend("invalid", "test")

    with pytest.raises(ValueError, match="test must be one of"):
        _normalize_shadow_backend("slow", "test")


def test_shadow_contourf_with_blur():
    """Test shadow_contourf with blur to cover _blur_filter_factory."""
    pytest.importorskip("scipy")  # Ensure scipy is available for gaussian_filter

    x = np.linspace(-2.0, 2.0, 16)
    y = np.linspace(-1.5, 1.5, 12)
    xx, yy = np.meshgrid(x, y)
    z = np.sin(xx * 2.0) + np.cos(yy * 3.0)

    fig, ax = plt.subplots()

    # Test with blur > 0 (creates filter and actually uses it when drawing)
    result = shadow_contourf(
        x,
        y,
        z,
        levels=5,
        ax=ax,
        shadow_blur=2.5,
    )
    assert len(result._skyborn_shadow_artists) > 0

    # Trigger rendering to ensure the blur filter is actually called
    fig.canvas.draw()

    plt.close(fig)

    # Test with blur = 0 (no filter)
    fig, ax = plt.subplots()
    result = shadow_contourf(
        x,
        y,
        z,
        levels=5,
        ax=ax,
        shadow_blur=0.0,
    )
    assert len(result._skyborn_shadow_artists) > 0

    plt.close(fig)


def test_shadow_contourf_with_hatches():
    """Test shadow_contourf with hatches to cover _iter_contour_layers hatch logic."""
    x = np.linspace(-2.0, 2.0, 16)
    y = np.linspace(-1.5, 1.5, 12)
    xx, yy = np.meshgrid(x, y)
    z = np.sin(xx * 2.0) + np.cos(yy * 3.0)

    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 5),
        hatches=["/", "\\", "x", "."],
        ax=ax,
    )

    assert len(result._skyborn_shadow_artists) > 0

    plt.close(fig)


def test_shadow_contourf_empty_path():
    """Test shadow_contourf handles empty paths in _path_touches_view_boundary."""
    x = np.linspace(-2.0, 2.0, 16)
    y = np.linspace(-1.5, 1.5, 12)
    xx, yy = np.meshgrid(x, y)
    # Create data that might produce empty paths
    z = np.ones_like(xx) * 5.0

    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=[0, 1, 2, 3, 4],
        ax=ax,
    )

    # Should not crash even if some paths are empty
    assert isinstance(result, plt.matplotlib.contour.QuadContourSet)

    plt.close(fig)


def test_shadow_contourf_z_only_with_extent():
    """Test shadow_contourf with z-only input and extent to cover _initialize_contour_xy."""
    z = np.random.randn(10, 15)

    fig, ax = plt.subplots()

    # Test with extent
    result = shadow_contourf(
        z,
        levels=5,
        extent=[0, 10, 0, 5],
        ax=ax,
    )

    assert isinstance(result, plt.matplotlib.contour.QuadContourSet)

    plt.close(fig)

    # Test without extent
    fig, ax = plt.subplots()

    result = shadow_contourf(
        z,
        levels=5,
        ax=ax,
    )

    assert isinstance(result, plt.matplotlib.contour.QuadContourSet)

    plt.close(fig)


def test_shadow_contourf_z_only_with_origin():
    """Test shadow_contourf with origin parameter to cover _initialize_contour_xy."""
    z = np.random.randn(10, 15)

    fig, ax = plt.subplots()

    result = shadow_contourf(
        z,
        levels=5,
        origin="upper",
        ax=ax,
    )

    assert isinstance(result, plt.matplotlib.contour.QuadContourSet)

    plt.close(fig)


from types import SimpleNamespace
from unittest.mock import Mock

from matplotlib.collections import LineCollection
from matplotlib.path import Path
from matplotlib.transforms import Bbox, IdentityTransform

import skyborn.plot._core.contour_arrows as arrow_module
import skyborn.plot.contour as contour_module
from skyborn.plot._core.contour_arrows import (
    _arrow_contour_label_positions,
    _arrow_end_distances,
    _build_arrow_segments,
    _build_arrow_segments_python,
    _build_arrow_segments_swept_python,
    _build_arrow_triangles_python,
    _copy_line_collection_properties,
    _is_closed_path,
    _label_anchor_away_from_arrows,
    _line_collection_linestyle_from_contour,
    _local_straightness_score,
    _orient_closed_vertices,
    _signed_area,
    _validate_contour_direction,
    _validate_positive_float,
    _validate_positive_int,
)
from skyborn.plot._core.vector_engine import (
    Grid,
    _build_ncl_curve,
    _clip_display_step_to_viewport,
    _corrected_ncl_display_origin,
    _default_ncl_box_center_candidates,
    _resolve_curly_anchor,
    _resolve_ncl_length_scale,
    _resolve_ncl_reference_length_px,
    _select_ncl_centers,
    _trace_ncl_curve,
    _warn_legacy_streamline_controls,
)
from skyborn.plot.contour import (
    _add_layered_shadow_artists,
    _apply_artist_filter,
    _auto_contour_levels,
    _check_contour_xyz,
    _contourpy_contourf,
    _contourpy_generator_kwargs,
    _contourpy_supported,
    _ContourpyCall,
    _hide_contour_artists,
    _initialize_contour_xy,
    _install_layered_remove_hook,
    _iter_contour_layers,
    _path_touches_view_boundary,
    _read_xyz_levels_contourpy_call,
    _read_z_levels_contourpy_call,
    _resolve_contourpy_input,
    _resolve_contourpy_levels,
    _validate_contourpy_levels,
)


def _fake_contour_set(*, segments=None, levels=None, edgecolors=None):
    return SimpleNamespace(
        levels=[1.0] if levels is None else levels,
        allsegs=[[[0.0, 0.0], [1.0, 0.0]]] if segments is None else segments,
        get_edgecolors=lambda: [] if edgecolors is None else edgecolors,
        get_transform=lambda: IdentityTransform(),
        get_zorder=lambda: 2.0,
        get_linewidths=lambda: [],
        get_linestyles=lambda: [],
    )


def test_arrow_validation_and_geometry_fallbacks():
    assert _validate_positive_int("2", "count") == 2
    assert _validate_positive_float("1.5", "size") == pytest.approx(1.5)
    assert _validate_contour_direction("CLOCKWISE", "direction") == "clockwise"
    with pytest.raises(ValueError):
        _validate_positive_int(0, "count")
    with pytest.raises(ValueError):
        _validate_positive_float(0, "size")
    with pytest.raises(ValueError):
        _validate_contour_direction("sideways", "direction")

    assert _line_collection_linestyle_from_contour("solid", 1.0) == "solid"
    style = _line_collection_linestyle_from_contour((2.0, (4.0, 8.0)), 2.0)
    assert style == (1.0, [2.0, 4.0])
    assert _line_collection_linestyle_from_contour((2.0, None), 2.0) == (2.0, None)
    assert _line_collection_linestyle_from_contour((2.0, (4.0,)), 0.0) == (
        2.0,
        (4.0,),
    )

    assert _signed_area(np.ones((2, 2))) == 0.0
    assert _is_closed_path(np.ones((2, 2))) is False
    assert _is_closed_path(np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]]))
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]])
    assert _signed_area(vertices) != 0.0
    assert (
        _orient_closed_vertices(
            vertices, clockwise=False, data_to_display_transform=IdentityTransform()
        ).shape
        == vertices.shape
    )
    flat = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    np.testing.assert_allclose(
        _orient_closed_vertices(
            flat, clockwise=True, data_to_display_transform=IdentityTransform()
        ),
        flat,
    )
    np.testing.assert_allclose(_arrow_end_distances(10.0, 2), [2.5, 7.5])
    assert (
        _local_straightness_score(
            np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]]), 1.0, 0.5
        )
        >= 0.0
    )


@pytest.mark.parametrize(
    "builder",
    [
        _build_arrow_triangles_python,
        _build_arrow_segments_python,
        _build_arrow_segments_swept_python,
    ],
)
def test_arrow_builders_skip_zero_length_vector(monkeypatch, builder):
    monkeypatch.setattr(
        arrow_module, "_select_arrow_end_distances", lambda *args: [1.0]
    )
    monkeypatch.setattr(
        arrow_module,
        "_point_at_distance",
        Mock(side_effect=[np.array([0.0, 0.0]), np.array([1.0, 0.0])]),
    )
    monkeypatch.setattr(arrow_module.np, "hypot", lambda *args: 0.0)
    result, metadata = builder(np.array([[0.0, 0.0], [1.0, 0.0]]), 1.0, 1, 0.5, 0.4)
    assert np.asarray(result).size == 0
    assert np.asarray(metadata).size == 0


def test_copy_properties_supports_missing_and_fallback_attributes():
    target = LineCollection([])
    clip = plt.Circle((0.0, 0.0), 1.0)
    source = SimpleNamespace(
        get_alpha=lambda: None,
        get_antialiaseds=lambda: np.array([], dtype=bool),
        get_clip_on=lambda: False,
        get_clip_box=lambda: None,
        get_clip_path=lambda: clip,
        get_path_effects=lambda: [],
        get_rasterized=lambda: True,
    )
    _copy_line_collection_properties(source, target)
    assert target.get_clip_on() is False
    assert target.get_rasterized() is True

    fallback = SimpleNamespace(
        get_antialiased=lambda: np.array([True]),
        get_clip_on=lambda: True,
        get_clip_box=lambda: None,
        get_clip_path=lambda: None,
        get_path_effects=lambda: [],
        get_rasterized=lambda: False,
    )
    _copy_line_collection_properties(fallback, target)
    assert target.get_antialiaseds()[0]


@pytest.mark.parametrize(
    "builder",
    [
        _build_arrow_triangles_python,
        _build_arrow_segments_python,
        _build_arrow_segments_swept_python,
    ],
)
def test_arrow_builders_skip_invalid_and_use_tangent_fallback(monkeypatch, builder):
    monkeypatch.setattr(
        arrow_module, "_select_arrow_end_distances", lambda *args: np.array([1.0, 2.0])
    )
    monkeypatch.setattr(
        arrow_module,
        "_point_at_distance",
        Mock(
            side_effect=[
                None,
                np.array([1.0, 0.0]),
                np.array([0.0, 0.0]),
                np.array([1.0, 0.0]),
                np.array([1.0, 0.0]),
                np.array([2.0, 0.0]),
            ]
        ),
    )
    monkeypatch.setattr(
        arrow_module, "_local_tangent_at_distance_c", lambda *args: None
    )
    result, metadata = builder(
        np.array([[0.0, 0.0], [2.0, 0.0]]),
        2.0,
        2,
        0.5,
        0.4,
    )
    result = np.asarray(result)
    assert result.size > 0
    assert np.asarray(metadata).ndim in {2, 3}


def test_arrow_native_wrapper_and_label_helpers_cover_empty_paths(monkeypatch):
    monkeypatch.setattr(
        arrow_module,
        "_native_build_arrow_segments",
        lambda *args: (np.empty((0, 2, 2)), np.empty((0, 2, 2))),
    )
    heads, metadata = _build_arrow_segments(
        np.array([[0.0, 0.0], [1.0, 0.0]]), 1.0, 1, 0.2, 0.4
    )
    assert heads.shape == (0, 2, 2)
    assert metadata.shape == (0, 2, 2)

    path = np.array([[0.0, 0.0], [0.0, 0.0]])
    assert _label_anchor_away_from_arrows(path, [], closed=False) is None
    assert (
        _label_anchor_away_from_arrows(
            np.array([[0.0, 0.0], [1.0, 0.0]]), [], closed=False
        )
        is not None
    )
    assert (
        _label_anchor_away_from_arrows(
            np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]]),
            [np.array([1.0, 0.0])],
            closed=True,
        )
        is not None
    )

    contour_set = SimpleNamespace(
        _skyborn_arrow_contour_artists=[
            SimpleNamespace(get_segments=lambda: []),
            SimpleNamespace(get_segments=lambda: [np.array([[0.0, 0.0]])]),
        ],
        get_transform=lambda: IdentityTransform(),
    )
    assert _arrow_contour_label_positions(contour_set) == []
    valid_line = SimpleNamespace(
        get_segments=lambda: [np.array([[0.0, 0.0], [1.0, 0.0]])],
        _skyborn_contour_arrow_segments=[],
    )
    valid_set = SimpleNamespace(
        _skyborn_arrow_contour_artists=[valid_line],
        get_transform=lambda: IdentityTransform(),
    )
    monkeypatch.setattr(
        arrow_module, "_label_anchor_away_from_arrows", lambda *args: None
    )
    assert _arrow_contour_label_positions(valid_set) == []


def test_arrow_add_artist_fallbacks_and_remove_hook(monkeypatch):
    fig, ax = plt.subplots()
    try:
        contour_set = _fake_contour_set(
            segments=[
                [],
                [[[0.0, 0.0], [0.0, 0.0]]],
                [[[0.0, 0.0], [1.0, 0.0]]],
            ],
            levels=[0.0, 1.0, 2.0],
        )
        monkeypatch.setattr(
            arrow_module,
            "_build_arrow_segments_python",
            lambda *args: (np.empty((0, 2, 2)), np.empty((0, 2, 2))),
        )
        monkeypatch.setattr(
            arrow_module,
            "_build_arrow_triangles_python",
            lambda *args: ([], np.empty((0, 2, 2))),
        )
        arrows = arrow_module._add_contour_arrows(
            contour_set,
            ax,
            arrow_count=1,
            arrow_size=0.4,
            arrow_length_fraction=0.1,
            arrow_length_points=None,
            arrow_max_length=10.0,
            positive_direction="clockwise",
            arrow_color=None,
            arrow_linewidth=None,
            zorder=None,
            arrow_style="line",
        )
        assert len(arrows) == 1
        explicit_width = arrow_module._add_contour_arrows(
            contour_set,
            ax,
            arrow_count=1,
            arrow_size=0.4,
            arrow_length_fraction=0.1,
            arrow_length_points=None,
            arrow_max_length=10.0,
            positive_direction="clockwise",
            arrow_color="red",
            arrow_linewidth=2.0,
            zorder=3.0,
            arrow_style="line",
        )
        assert len(explicit_width) == 1

        filled = arrow_module._add_contour_arrows(
            contour_set,
            ax,
            arrow_count=1,
            arrow_size=0.4,
            arrow_length_fraction=0.1,
            arrow_length_points=3.0,
            arrow_max_length=10.0,
            positive_direction="counterclockwise",
            arrow_color=None,
            arrow_linewidth=None,
            zorder=5.0,
            arrow_style="filled",
        )
        assert len(filled) == 1

        removed = {"count": 0}

        class _Removable:
            def remove(self):
                removed["count"] += 1
                raise ValueError("already removed")

        contour = SimpleNamespace(
            remove=lambda: removed.__setitem__("count", removed["count"] + 1)
        )
        arrow_module._install_arrow_remove_hook(contour, [_Removable()])
        contour.remove()
        assert removed["count"] == 2
    finally:
        plt.close(fig)


def test_contour_coordinate_and_layer_helpers_cover_fallbacks():
    z = np.arange(6.0).reshape(2, 3)
    x0, y0 = _initialize_contour_xy(z, None, None)
    assert x0.shape == z.shape and y0.shape == z.shape
    x1, y1 = _initialize_contour_xy(z, None, (0.0, 2.0, 0.0, 1.0))
    assert x1.shape == z.shape and y1.shape == z.shape
    x2, y2 = _initialize_contour_xy(z, "lower", None)
    assert x2.shape == z.shape and y2.shape == z.shape
    x3, y3 = _initialize_contour_xy(z, "upper", (0.0, 2.0, 0.0, 1.0))
    np.testing.assert_allclose(y3[0], y3[0, 0] + (y3[0, 1] - y3[0, 0]) * 0)
    with pytest.raises(TypeError):
        _initialize_contour_xy(np.ones(3), None, None)
    with pytest.raises(TypeError):
        _initialize_contour_xy(np.ones((1, 3)), None, None)

    class _Artist:
        def __init__(self):
            self.visible = True

        def set_visible(self, value):
            self.visible = value

    contour_set = _Artist()
    contour_set.collections = [_Artist()]
    _hide_contour_artists(contour_set)
    assert contour_set.visible is False and contour_set.collections[0].visible is False
    _hide_contour_artists(SimpleNamespace())

    empty_path = Path(np.empty((0, 2)))
    fig, ax = plt.subplots()
    try:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        assert _path_touches_view_boundary(empty_path, ax) is False
        assert bool(_path_touches_view_boundary(Path([[0.0, 0.5]]), ax)) is True
        assert bool(_path_touches_view_boundary(Path([[0.5, 0.5]]), ax)) is False
    finally:
        plt.close(fig)


def test_contour_layer_iteration_filter_and_remove_hook():
    path = Path([[0.2, 0.2], [0.3, 0.3]])
    direct = SimpleNamespace(
        hatches=["/"],
        get_paths=lambda: [path],
        get_facecolors=lambda: [np.array([1.0, 0.0, 0.0, 1.0])],
    )
    assert list(_iter_contour_layers(direct))[0][2] == "/"

    collection = SimpleNamespace(
        get_facecolors=lambda: [],
        get_hatch=lambda: "x",
        get_paths=lambda: [path],
    )
    fallback = SimpleNamespace(collections=[collection], hatches=[])
    assert list(_iter_contour_layers(fallback))[0][2] == "x"

    artist = SimpleNamespace(
        set_agg_filter=lambda value: setattr(artist, "filter", value)
    )
    _apply_artist_filter(artist, object())
    _apply_artist_filter(SimpleNamespace(), None)
    assert hasattr(artist, "filter")

    fig, ax = plt.subplots()
    try:
        cs = SimpleNamespace(
            get_transform=lambda: IdentityTransform(),
            get_zorder=lambda: 1.0,
            get_alpha=lambda: 0.5,
            hatches=[],
            get_paths=lambda: [path],
            get_facecolors=lambda: [np.array([0.2, 0.3, 0.4, 1.0])],
            set_visible=lambda value: None,
            collections=[],
            remove=lambda: None,
        )
        artists = _add_layered_shadow_artists(
            cs, ax, (1.0, -1.0), "black", 0.3, 0.0, boundary_margin=0.0
        )
        assert len(artists) == 2
        removed = {"count": 0}

        class _AlreadyRemoved:
            def remove(self):
                raise ValueError("gone")

        hook_target = SimpleNamespace(
            remove=lambda: removed.__setitem__("count", removed["count"] + 1)
        )
        _install_layered_remove_hook(hook_target, [_AlreadyRemoved()])
        hook_target.remove()
        assert removed["count"] == 1
    finally:
        plt.close(fig)


def test_contour_xyz_levels_and_contourpy_dispatch(monkeypatch):
    fig, ax = plt.subplots()
    try:
        z = np.arange(6.0).reshape(2, 3)
        x, y, checked = _check_contour_xyz(ax, [0, 1, 2], [0, 1], z, {})
        assert x.shape == y.shape == checked.shape
        z_call = _read_z_levels_contourpy_call(ax, [z, 3], {})
        assert z_call.levels_arg == 3
        xyz_call = _read_xyz_levels_contourpy_call(ax, [[0, 1, 2], [0, 1], z, 3], {})
        assert xyz_call.levels_arg == 3
        for bad in (
            np.ones(3),
            np.ones((1, 3)),
        ):
            with pytest.raises(TypeError):
                _check_contour_xyz(ax, [0, 1, 2], [0, 1], bad, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, [0, 1], [0, 1], z, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, [0, 1, 2], [0], z, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, np.ones((2, 2)), [0, 1], z, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, np.ones((2, 2)), np.ones((2, 2)), z, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, np.ones((2, 3)), np.ones((3, 2)), z, {})
        with pytest.raises(TypeError):
            _check_contour_xyz(ax, np.ones((2, 3, 1)), np.ones((2, 3, 1)), z, {})

        assert _resolve_contourpy_levels(
            np.array([[True, False]]), None, False
        ).tolist() == [0.0, 0.5, 1.0]
        assert _resolve_contourpy_levels(np.array([[1.0, 2.0]]), 3, False).ndim == 1
        assert _resolve_contourpy_levels(
            np.array([[1.0, 2.0]]), [1, 2, 3], False
        ).tolist() == [1.0, 2.0, 3.0]
        with pytest.raises(ValueError):
            _validate_contourpy_levels(np.array([1.0]))
        with pytest.raises(ValueError):
            _validate_contourpy_levels(np.array([1.0, 1.0]))

        assert _contourpy_supported({}) is True
        assert _contourpy_supported({"extend": "both"}) is False
        assert _contourpy_supported({"locator": object()}) is False
        generated = _contourpy_generator_kwargs({"algorithm": "mpl2005", "nchunk": 2})
        assert generated["corner_mask"] is False
        assert generated["chunk_size"] == 2
        assert (
            _contourpy_generator_kwargs({"algorithm": "mpl2014", "corner_mask": True})[
                "corner_mask"
            ]
            is True
        )
        assert contour_module._precomputed_contour_kwargs(
            {
                "levels": 3,
                "algorithm": "mpl2014",
                "corner_mask": True,
                "nchunk": 2,
                "cmap": "viridis",
            }
        ) == {"cmap": "viridis"}

        call = _ContourpyCall(
            x=np.arange(3.0),
            y=np.arange(2.0),
            z=np.ma.asarray(z),
            levels_arg=3,
        )
        monkeypatch.setattr(
            contour_module, "_CONTOURPY_CALL_READERS", {2: lambda *args: call}
        )
        assert _resolve_contourpy_input(ax, [z, 3], {}) is not None
        assert _resolve_contourpy_input(ax, [z, 3], {"levels": 4}) is None
        monkeypatch.setattr(contour_module, "_CONTOURPY_CALL_READERS", {})
        assert _resolve_contourpy_input(ax, [z], {}) is None

        assert _contourpy_contourf(ax, [z], {"locator": object()}) is None
        monkeypatch.setattr(
            contour_module, "_resolve_contourpy_input", lambda *args: None
        )
        assert _contourpy_contourf(ax, [z], {}) is None
    finally:
        plt.close(fig)


def test_contourpy_empty_geometry_and_shadow_public_fallbacks(monkeypatch):
    fig, ax = plt.subplots()
    try:
        z = np.arange(12.0).reshape(3, 4)
        monkeypatch.setattr(
            contour_module.contourpy,
            "contour_generator",
            lambda *args, **kwargs: SimpleNamespace(
                filled=lambda lower, upper: ([], [])
            ),
        )
        assert _contourpy_contourf(ax, [z], {}) is None
        assert _auto_contour_levels(1.0, 1.0, 2).size >= 2

        monkeypatch.undo()
        from skyborn.plot.contour import (
            arrow_contour,
            arrow_contour_clabel,
            shadow_contourf,
        )

        cs = arrow_contour(ax, z, levels=3, arrows=False)
        assert cs._skyborn_arrow_contour_artists == []
        manual = arrow_contour_clabel(cs, manual=[])
        assert isinstance(manual, list)
        sf = shadow_contourf(ax, z, levels=3, shadow=False)
        assert sf._skyborn_shadow_artists == []
    finally:
        plt.close(fig)


def test_contour_remaining_public_and_locator_branches(monkeypatch):
    class _Locator:
        def __init__(self, *args, **kwargs):
            del args, kwargs

        def tick_values(self, zmin, zmax):
            del zmin, zmax
            return np.array([0.0, 2.0])

    monkeypatch.setattr(contour_module.mpl.ticker, "MaxNLocator", _Locator)
    assert _auto_contour_levels(1.0, 1.0, 2).size == 2

    fig, ax = plt.subplots()
    try:
        z = np.arange(12.0).reshape(3, 4)
        plt.sca(ax)
        from skyborn.plot.contour import arrow_contour, shadow_contourf

        assert arrow_contour(z, levels=3, arrows=False).axes is ax
        assert shadow_contourf(z, levels=3, shadow=False).axes is ax
    finally:
        plt.close(fig)


def test_vector_engine_remaining_fallback_branches():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _warn_legacy_streamline_controls("bad", False)
    assert any(item.category is FutureWarning for item in caught)

    assert _resolve_ncl_reference_length_px(
        1.0, 5.0, 1.0, 0.0, 1.0, 10.0
    ) == pytest.approx(2.0)
    scale = _resolve_ncl_length_scale(1.0, 5.0, 1.0, 0.0, 1.0, 10.0)
    assert scale["adjust_min"] is True
    assert _resolve_curly_anchor(None, "forward") == "tail"

    assert _corrected_ncl_display_origin([1, 2], None).tolist() == [1.0, 2.0]
    assert _corrected_ncl_display_origin([2, 5], [1, 2]).tolist() == pytest.approx(
        [1.6666667, 4.0]
    )

    viewport = Bbox.from_extents(0.0, 0.0, 10.0, 10.0)
    assert _clip_display_step_to_viewport([5, 5], [6, 6], viewport)[1] is False
    assert _clip_display_step_to_viewport([5, 5], [np.nan, 6], viewport)[1] is True
    assert _clip_display_step_to_viewport([5, 5], [20, 5], viewport)[0][
        0
    ] == pytest.approx(10.0)
    assert _clip_display_step_to_viewport([5, 5], [-20, 5], viewport)[0][
        0
    ] == pytest.approx(0.0)
    assert _clip_display_step_to_viewport([5, 5], [5, 20], viewport)[0][
        1
    ] == pytest.approx(10.0)
    assert _clip_display_step_to_viewport([5, 5], [5, -20], viewport)[0][
        1
    ] == pytest.approx(0.0)

    tiny_grid = SimpleNamespace(nx=1, ny=1)
    assert _default_ncl_box_center_candidates(tiny_grid).shape == (0, 2)


def test_vector_engine_selection_and_trace_retry_branches():
    grid = Grid(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]))
    axes = SimpleNamespace(bbox=Bbox.from_extents(0.0, 0.0, 10.0, 10.0))
    selected = _select_ncl_centers(
        grid=grid,
        magnitude=np.ones(grid.shape),
        transform=IdentityTransform(),
        axes=axes,
        density=1.0,
        start_points=np.array([[0.5, 0.5], [1.5, 1.5]]),
        min_distance=0.1,
        sample_grid_field_array=lambda grid, field, points: np.ones(len(points)),
        thin_ncl_display_candidates=lambda points, bbox, spacing: np.array([1]),
        thin_ncl_mapped_candidates=lambda points, spacing: np.array([0]),
    )
    assert len(selected) == 1
    assert selected[0][0].tolist() == [1.5, 1.5]

    assert (
        _trace_ncl_curve(
            np.array([0.0, 0.0]),
            0.0,
            "tail",
            grid,
            np.ones(grid.shape),
            np.ones(grid.shape),
            IdentityTransform(),
            1.0,
            1.0,
            axes.bbox,
            trace_ncl_direction_fn=lambda *args, **kwargs: np.array(
                [[0.0, 0.0], [1.0, 0.0]]
            ),
        )
        is None
    )
    short = _trace_ncl_curve(
        np.array([0.0, 0.0]),
        4.0,
        "tail",
        grid,
        np.ones(grid.shape),
        np.ones(grid.shape),
        IdentityTransform(),
        1.0,
        1.0,
        axes.bbox,
        trace_ncl_direction_fn=lambda *args, **kwargs: np.array([[0.0, 0.0]]),
    )
    assert short.shape == (1, 2)

    attempts = {"count": 0}

    def trace_retry(**kwargs):
        attempts["count"] += 1
        return (
            np.array([[0.0, 0.0]])
            if attempts["count"] == 1
            else np.array([[0.0, 0.0], [1.0, 0.0]])
        )

    curve = _build_ncl_curve(
        start_point=np.array([0.0, 0.0]),
        total_length_px=4.0,
        anchor="tail",
        grid=grid,
        u=np.ones(grid.shape),
        v=np.ones(grid.shape),
        transform=IdentityTransform(),
        step_px=1.0,
        speed_scale=1.0,
        viewport=axes.bbox,
        trace_ncl_curve_fn=trace_retry,
        evaluate_ncl_display_curve_fn=lambda curve, transform, viewport=None: (
            np.asarray(curve, dtype=float),
            False,
        ),
    )
    assert curve is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
