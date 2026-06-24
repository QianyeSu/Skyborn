"""Tests for shadowed filled-contour plotting helpers."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pytest
from matplotlib.collections import LineCollection, PathCollection
from matplotlib.colors import LogNorm
from matplotlib.contour import QuadContourSet

from skyborn.plot import arrow_contour as public_arrow_contour
from skyborn.plot import scatter as public_scatter
from skyborn.plot import shadow_contourf as public_shadow_contourf
from skyborn.plot._core.contour_arrows import _build_arrow_segments_python
from skyborn.plot.contour import arrow_contour, shadow_contourf
from skyborn.plot.contour_core import (
    build_arrow_segments as native_build_arrow_segments,
)


def _sample_field():
    x = np.linspace(-2.0, 2.0, 16)
    y = np.linspace(-1.5, 1.5, 12)
    xx, yy = np.meshgrid(x, y)
    z = np.sin(xx * 2.0) + np.cos(yy * 3.0)
    return x, y, z


def _sample_signed_circles():
    x = np.linspace(-2.0, 2.0, 101)
    y = np.linspace(-2.0, 2.0, 101)
    xx, yy = np.meshgrid(x, y)
    z = xx**2 + yy**2 - 1.0
    return x, y, z


def _cross2d(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _arrow_cross_in_display(ax, start, end):
    center = ax.transData.transform([0.0, 0.0])
    start = ax.transData.transform(start)
    end = ax.transData.transform(end)
    midpoint = 0.5 * (start + end)
    return _cross2d(midpoint - center, end - start)


def _assert_linestyle_equal(actual, expected):
    actual_offset, actual_dashes = actual
    expected_offset, expected_dashes = expected
    np.testing.assert_allclose(actual_offset, expected_offset)
    if actual_dashes is None or expected_dashes is None:
        assert actual_dashes is expected_dashes
    else:
        np.testing.assert_allclose(actual_dashes, expected_dashes)


def _assert_linestyle_solid(actual):
    actual_offset, actual_dashes = actual
    np.testing.assert_allclose(actual_offset, 0.0)
    assert actual_dashes is None


def test_shadow_contourf_is_exported():
    assert public_shadow_contourf is shadow_contourf


def test_arrow_contour_is_exported():
    assert public_arrow_contour is arrow_contour


def test_shadow_scatter_and_arrow_contour_can_be_combined():
    x, y, z = _sample_field()
    significance_mask = z > 1.0
    fig, ax = plt.subplots(figsize=(5, 4))

    filled = public_shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        ax=ax,
        cmap="viridis",
        shadow_offset=(2.5, -2.5),
        shadow_alpha=0.25,
        shadow_blur=0.6,
        zorder=1,
    )
    stipple = public_scatter(
        x,
        y,
        where=significance_mask,
        ax=ax,
        placement="cells",
        density=1.2,
        s=7.0,
        c="black",
        alpha=0.45,
        linewidths=0.0,
        zorder=4,
    )
    lines = public_arrow_contour(
        x,
        y,
        z,
        levels=[-1.0, 0.0, 1.0],
        ax=ax,
        colors="white",
        linewidths=1.2,
        arrow_count=2,
        arrow_size=0.45,
        arrow_length_fraction=0.04,
        zorder=5,
    )

    assert isinstance(filled, QuadContourSet)
    assert isinstance(stipple, PathCollection)
    assert isinstance(lines, QuadContourSet)
    assert len(filled._skyborn_shadow_artists) > 0
    assert len(stipple.get_offsets()) > 0
    assert len(lines._skyborn_contour_arrows) > 0
    assert not filled.get_visible()
    assert not lines.get_visible()
    assert all(artist.get_zorder() == 5 for artist in lines._skyborn_contour_arrows)
    assert stipple.get_zorder() == 4

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_adds_sign_directed_arrows_to_closed_levels():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    result = arrow_contour(
        x,
        y,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        colors="black",
        arrow_count=2,
        arrow_size=0.45,
    )

    line_artists = result._skyborn_contour_arrows
    assert isinstance(result, QuadContourSet)
    assert len(line_artists) == 2
    assert all(isinstance(artist, LineCollection) for artist in line_artists)
    assert all(
        artist._skyborn_contour_kind == "arrow_contour" for artist in line_artists
    )
    assert not result.get_visible()

    positive_lines = [
        artist
        for artist in line_artists
        if artist._skyborn_contour_level == pytest.approx(0.5)
    ]
    negative_lines = [
        artist
        for artist in line_artists
        if artist._skyborn_contour_level == pytest.approx(-0.5)
    ]
    assert len(positive_lines) == 1
    assert len(negative_lines) == 1
    assert all(
        artist._skyborn_contour_direction == "clockwise" for artist in positive_lines
    )
    assert all(
        artist._skyborn_contour_direction == "counterclockwise"
        for artist in negative_lines
    )
    assert all(len(artist.get_segments()) == 5 for artist in line_artists)

    for start, end in positive_lines[0]._skyborn_contour_arrow_segments:
        start = np.asarray(start)
        end = np.asarray(end)
        midpoint = 0.5 * (start + end)
        assert _cross2d(midpoint, end - start) < 0.0

    for start, end in negative_lines[0]._skyborn_contour_arrow_segments:
        start = np.asarray(start)
        end = np.asarray(end)
        midpoint = 0.5 * (start + end)
        assert _cross2d(midpoint, end - start) > 0.0

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_can_reverse_positive_closed_direction():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    result = arrow_contour(
        x,
        y,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        colors="black",
        arrow_count=1,
        positive_direction="counterclockwise",
    )

    artists_by_level = {
        artist._skyborn_contour_level: artist
        for artist in result._skyborn_contour_arrows
    }
    positive_line = artists_by_level[0.5]
    negative_line = artists_by_level[-0.5]

    assert positive_line._skyborn_contour_direction == "counterclockwise"
    assert negative_line._skyborn_contour_direction == "clockwise"

    start, end = positive_line._skyborn_contour_arrow_segments[0]
    start = np.asarray(start)
    end = np.asarray(end)
    midpoint = 0.5 * (start + end)
    assert _cross2d(midpoint, end - start) > 0.0

    start, end = negative_line._skyborn_contour_arrow_segments[0]
    start = np.asarray(start)
    end = np.asarray(end)
    midpoint = 0.5 * (start + end)
    assert _cross2d(midpoint, end - start) < 0.0

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_uses_visual_direction_on_inverted_axes():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(2.0, -2.0)

    result = arrow_contour(
        x,
        y,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        colors="black",
        arrow_count=1,
    )

    artists_by_level = {
        artist._skyborn_contour_level: artist
        for artist in result._skyborn_contour_arrows
    }
    positive_line = artists_by_level[0.5]
    negative_line = artists_by_level[-0.5]

    assert positive_line._skyborn_contour_direction == "clockwise"
    assert negative_line._skyborn_contour_direction == "counterclockwise"

    start, end = positive_line._skyborn_contour_arrow_segments[0]
    assert _arrow_cross_in_display(ax, start, end) < 0.0

    start, end = negative_line._skyborn_contour_arrow_segments[0]
    assert _arrow_cross_in_display(ax, start, end) > 0.0

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_accepts_positional_axes_and_disable_arrows():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    result = arrow_contour(ax, x, y, z, levels=[0.5], arrows=False)

    assert isinstance(result, QuadContourSet)
    assert result.axes is ax
    assert result._skyborn_contour_arrows == []

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_remove_cleans_generated_arrows():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    result = arrow_contour(x, y, z, levels=[0.5], ax=ax, arrow_count=2)
    arrows = list(result._skyborn_contour_arrows)

    assert len(arrows) == 1
    result.remove()

    assert result._skyborn_contour_arrows == []
    assert all(arrow.axes is None for arrow in arrows)

    plt.close(fig)


def test_arrow_contour_accepts_cartopy_transform_when_available():
    cartopy = pytest.importorskip("cartopy.crs")
    x, y, z = _sample_signed_circles()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=cartopy.PlateCarree())

    result = arrow_contour(
        x,
        y,
        z,
        levels=[0.5],
        ax=ax,
        transform=cartopy.PlateCarree(),
        arrow_count=1,
    )

    assert isinstance(result, QuadContourSet)
    assert len(result._skyborn_contour_arrows) == 1

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_validates_options():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    with pytest.raises(TypeError, match="Axes both positionally"):
        arrow_contour(ax, x, y, z, ax=ax)

    with pytest.raises(ValueError, match="arrow_count"):
        arrow_contour(x, y, z, ax=ax, arrow_count=0)

    with pytest.raises(ValueError, match="arrow_size"):
        arrow_contour(x, y, z, ax=ax, arrow_size=0.0)

    with pytest.raises(ValueError, match="positive_direction"):
        arrow_contour(x, y, z, ax=ax, positive_direction="reverse")

    plt.close(fig)


def test_native_arrow_segment_builder_matches_python():
    theta = np.linspace(0.0, 2.0 * np.pi, 80)
    vertices = np.column_stack([np.cos(theta), np.sin(theta)])
    total_length = float(np.sum(np.hypot(*np.diff(vertices, axis=0).T)))
    arrow_length = 0.2

    py_heads, py_meta = _build_arrow_segments_python(
        vertices,
        total_length,
        arrow_count=2,
        arrow_length=arrow_length,
        arrow_size=0.45,
    )
    native_heads, native_meta = native_build_arrow_segments(
        vertices,
        2,
        arrow_length,
        0.45,
    )

    assert native_heads.shape == py_heads.shape
    assert native_meta.shape == py_meta.shape
    np.testing.assert_allclose(native_heads, py_heads, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(native_meta, py_meta, rtol=1e-12, atol=1e-12)


def test_arrow_contour_accepts_projected_cartopy_axes_when_available():
    cartopy = pytest.importorskip("cartopy.crs")
    lon = np.linspace(-120.0, 120.0, 81)
    lat = np.linspace(-60.0, 60.0, 61)
    lon2d, lat2d = np.meshgrid(lon, lat)
    z = np.sin(np.deg2rad(lon2d)) + 0.6 * np.cos(np.deg2rad(lat2d * 2.0))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=cartopy.Robinson())
    ax.set_global()

    result = arrow_contour(
        lon,
        lat,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        transform=cartopy.PlateCarree(),
        arrow_count=1,
        colors="black",
    )

    assert isinstance(result, QuadContourSet)
    assert len(result._skyborn_contour_arrows) >= 1

    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_preserves_contour_linestyles():
    x, y, z = _sample_signed_circles()

    fig, ax = plt.subplots()
    result = arrow_contour(
        x,
        y,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        colors="black",
        linewidths=[1.5, 2.5],
        negative_linestyles="dashed",
    )

    artists_by_level = {
        artist._skyborn_contour_level: artist
        for artist in result._skyborn_contour_arrows
    }
    _assert_linestyle_equal(
        artists_by_level[-0.5].get_linestyles()[0],
        result.get_linestyles()[0],
    )
    for arrow_linestyle in artists_by_level[-0.5].get_linestyles()[1:]:
        _assert_linestyle_solid(arrow_linestyle)
    _assert_linestyle_equal(
        artists_by_level[0.5].get_linestyles()[0],
        result.get_linestyles()[1],
    )
    for arrow_linestyle in artists_by_level[0.5].get_linestyles()[1:]:
        _assert_linestyle_solid(arrow_linestyle)
    fig.canvas.draw()
    plt.close(fig)

    fig, ax = plt.subplots()
    result = arrow_contour(
        x,
        y,
        z,
        levels=[-0.5, 0.5],
        ax=ax,
        colors="black",
        linewidths=[2.0, 3.0],
        linestyles=[(3.0, (2.0, 4.0)), "dashdot"],
    )

    artists_by_level = {
        artist._skyborn_contour_level: artist
        for artist in result._skyborn_contour_arrows
    }
    _assert_linestyle_equal(
        artists_by_level[-0.5].get_linestyles()[0],
        result.get_linestyles()[0],
    )
    for arrow_linestyle in artists_by_level[-0.5].get_linestyles()[1:]:
        _assert_linestyle_solid(arrow_linestyle)
    _assert_linestyle_equal(
        artists_by_level[0.5].get_linestyles()[0],
        result.get_linestyles()[1],
    )
    for arrow_linestyle in artists_by_level[0.5].get_linestyles()[1:]:
        _assert_linestyle_solid(arrow_linestyle)
    fig.canvas.draw()
    plt.close(fig)


def test_arrow_contour_preserves_alpha_and_antialiasing():
    x, y, z = _sample_signed_circles()
    fig, ax = plt.subplots()

    result = arrow_contour(
        x,
        y,
        z,
        levels=[0.5],
        ax=ax,
        colors="black",
        alpha=0.35,
        antialiased=False,
        zorder=8,
    )

    arrow = result._skyborn_contour_arrows[0]
    assert arrow.get_alpha() == pytest.approx(0.35)
    assert arrow.get_edgecolors()[0, 3] == pytest.approx(0.35)
    assert not bool(arrow.get_antialiaseds()[0])
    assert arrow.get_zorder() == 8

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_accepts_matplotlib_contourf_arguments():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        cmap="viridis",
        extend="both",
        ax=ax,
        shadow_offset=(3.0, -4.0),
        shadow_alpha=0.4,
    )

    assert isinstance(result, QuadContourSet)
    assert result.extend == "both"
    assert result._skyborn_shadow_contour_set is None
    assert result._skyborn_shadow_method == "layered"
    assert len(result._skyborn_shadow_artists) >= 2 * (len(result.levels) - 1)
    assert len(result._skyborn_shadow_artists) % 2 == 0
    assert not result.get_visible()

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_accepts_positional_axes():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(ax, x, y, z, levels=5, colors=["red", "blue"])

    assert isinstance(result, QuadContourSet)
    assert result.axes is ax
    assert not result.get_visible()

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_fast_backend_returns_quad_contour_set():
    x, y, z = _sample_field()
    levels = np.linspace(-2.0, 2.0, 9)
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=levels,
        cmap="viridis",
        ax=ax,
        shadow_blur=0.0,
        shadow_backend="fast",
    )

    assert isinstance(result, QuadContourSet)
    assert result._skyborn_shadow_backend == "fast"
    assert result._skyborn_shadow_engine == "contourpy"
    assert result._skyborn_shadow_method == "layered"
    assert len(result.get_paths()) == len(levels) - 1
    assert len(result._skyborn_shadow_artists) >= 2 * (len(levels) - 1)
    assert not result.get_visible()

    fig.colorbar(result, ax=ax)
    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_fast_backend_falls_back_for_extend():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        cmap="viridis",
        extend="both",
        ax=ax,
        shadow_backend="fast",
    )

    assert isinstance(result, QuadContourSet)
    assert result.extend == "both"
    assert result._skyborn_shadow_backend == "standard"
    assert result._skyborn_shadow_engine == "matplotlib"

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_fast_backend_accepts_algorithm():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        cmap="viridis",
        ax=ax,
        shadow_blur=0.0,
        shadow_backend="fast",
        algorithm="serial",
    )

    assert result._skyborn_shadow_backend == "fast"
    assert result._skyborn_shadow_engine == "contourpy"

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_fast_backend_supports_bool_default_levels():
    x = np.linspace(-2.0, 2.0, 16)
    y = np.linspace(-1.5, 1.5, 12)
    xx, yy = np.meshgrid(x, y)
    z = xx**2 + yy**2 > 1.0
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        colors=["white", "black"],
        ax=ax,
        shadow_blur=0.0,
        shadow_backend="fast",
    )

    assert result._skyborn_shadow_backend == "fast"
    assert result._skyborn_shadow_engine == "contourpy"
    np.testing.assert_allclose(result.levels, [0.0, 0.5, 1.0])

    fig.canvas.draw()
    plt.close(fig)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"locator": mticker.MaxNLocator(5)},
        {"norm": LogNorm()},
        {"data": {"x": np.linspace(-2.0, 2.0, 16)}},
    ],
)
def test_shadow_contourf_fast_backend_falls_back_for_unsupported_kwargs(kwargs):
    x, y, z = _sample_field()
    z = z - z.min() + 0.1
    args = (x, y, z)
    call_kwargs = dict(kwargs)
    if "data" in call_kwargs:
        args = ("x", y, z)
    fig, ax = plt.subplots()

    result = shadow_contourf(
        *args,
        levels=5,
        cmap="viridis",
        ax=ax,
        shadow_backend="fast",
        **call_kwargs,
    )

    assert result._skyborn_shadow_backend == "standard"
    assert result._skyborn_shadow_engine == "matplotlib"

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_accepts_legacy_shadow_engine_alias():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        cmap="viridis",
        ax=ax,
        shadow_blur=0.0,
        shadow_engine="contourpy",
    )

    assert result._skyborn_shadow_backend == "fast"
    assert result._skyborn_shadow_engine == "contourpy"

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_path_effect_method_draws_with_path_effects():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 9),
        ax=ax,
        shadow_method="path_effect",
        shadow_offset=(3.0, -4.0),
        shadow_alpha=0.4,
    )

    assert result._skyborn_shadow_method == "path_effect"
    assert result._skyborn_shadow_artists == []

    effects = result.get_path_effects()
    assert isinstance(effects[0], mpatheffects.SimplePatchShadow)
    assert isinstance(effects[1], mpatheffects.Normal)

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_layered_preserves_alpha_and_hatches():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 5),
        hatches=["/", "x", ".", "+"],
        alpha=0.7,
        ax=ax,
    )

    fill_artists = result._skyborn_shadow_artists[1::2]
    assert any(artist.get_hatch() == "/" for artist in fill_artists)
    assert all(artist.get_alpha() == 0.7 for artist in fill_artists)

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_layered_remove_cleans_generated_artists():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(x, y, z, levels=5, ax=ax)
    artists = list(result._skyborn_shadow_artists)

    assert artists
    result.remove()

    assert result._skyborn_shadow_artists == []
    assert all(artist.axes is None for artist in artists)

    plt.close(fig)


def test_shadow_contourf_rejects_duplicate_axes():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    with pytest.raises(TypeError, match="Axes both positionally"):
        shadow_contourf(ax, x, y, z, ax=ax)

    plt.close(fig)


def test_shadow_contourf_overlay_uses_single_shadow_contour_set():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(
        x,
        y,
        z,
        levels=np.linspace(-2.0, 2.0, 7),
        cmap="plasma",
        ax=ax,
        shadow_method="overlay",
        shadow_offset=(4.0, -3.0),
        shadow_alpha=0.25,
        shadow_blur=1.5,
    )

    shadow_set = result._skyborn_shadow_contour_set
    assert isinstance(result, QuadContourSet)
    assert isinstance(shadow_set, QuadContourSet)
    assert shadow_set is not result
    assert shadow_set.get_agg_filter() is not None
    assert result._skyborn_shadow_artists == []

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_accepts_cartopy_transform_when_available():
    cartopy = pytest.importorskip("cartopy.crs")
    x, y, z = _sample_field()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=cartopy.PlateCarree())

    result = shadow_contourf(
        x,
        y,
        z,
        levels=6,
        ax=ax,
        transform=cartopy.PlateCarree(),
        shadow_offset=(2.0, -2.0),
    )

    assert isinstance(result, QuadContourSet)
    assert result._skyborn_shadow_method == "layered"

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_can_disable_shadow():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    result = shadow_contourf(x, y, z, levels=4, ax=ax, shadow=False)

    assert result._skyborn_shadow_method is None
    assert result._skyborn_shadow_contour_set is None
    assert result._skyborn_shadow_artists == []
    assert result.get_path_effects() in (None, [])

    fig.canvas.draw()
    plt.close(fig)


def test_shadow_contourf_validates_shadow_options():
    x, y, z = _sample_field()
    fig, ax = plt.subplots()

    with pytest.raises(ValueError, match="shadow_method"):
        shadow_contourf(x, y, z, ax=ax, shadow_method="bad")

    with pytest.raises(ValueError, match="shadow_offset"):
        shadow_contourf(x, y, z, ax=ax, shadow_offset=(1.0, 2.0, 3.0))

    with pytest.raises(ValueError, match="shadow_backend"):
        shadow_contourf(x, y, z, ax=ax, shadow_backend="contourpy-fast")

    with pytest.raises(TypeError, match="shadow_backend"):
        shadow_contourf(
            x,
            y,
            z,
            ax=ax,
            shadow_backend="fast",
            shadow_engine="contourpy",
        )

    plt.close(fig)
