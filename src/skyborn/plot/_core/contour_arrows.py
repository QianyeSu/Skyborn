"""Arrowed-contour rendering helpers."""

from __future__ import annotations

import types
from typing import Any, Iterable, List, Optional

import matplotlib as mpl
import numpy as np
from matplotlib.collections import LineCollection

from ..contour_core import build_arrow_segments as _native_build_arrow_segments


def _validate_positive_int(value: Any, name: str) -> int:
    integer = int(value)
    if integer < 1:
        raise ValueError(f"{name} must be at least 1")
    return integer


def _validate_positive_float(value: Any, name: str) -> float:
    scalar = float(value)
    if scalar <= 0.0:
        raise ValueError(f"{name} must be positive")
    return scalar


def _validate_contour_direction(value: Any, name: str) -> str:
    direction = str(value).lower()
    if direction not in {"clockwise", "counterclockwise"}:
        raise ValueError(f"{name} must be 'clockwise' or 'counterclockwise'")
    return direction


def _copy_line_collection_properties(source: Any, target: LineCollection) -> None:
    if hasattr(source, "get_alpha"):
        alpha = source.get_alpha()
        if alpha is not None:
            target.set_alpha(alpha)
    if hasattr(source, "get_antialiaseds"):
        antialiaseds = source.get_antialiaseds()
        if len(antialiaseds):
            target.set_antialiaseds([antialiaseds[0]])
    elif hasattr(source, "get_antialiased"):
        antialiased = source.get_antialiased()
        if len(antialiased):
            target.set_antialiaseds([antialiased[0]])
    if hasattr(source, "get_clip_on"):
        target.set_clip_on(source.get_clip_on())
    if hasattr(source, "get_clip_box"):
        target.set_clip_box(source.get_clip_box())
    if hasattr(source, "get_clip_path"):
        clip_path = source.get_clip_path()
        if clip_path is not None:
            target.set_clip_path(clip_path)
    if hasattr(source, "get_path_effects"):
        target.set_path_effects(source.get_path_effects())
    if hasattr(source, "get_rasterized"):
        target.set_rasterized(source.get_rasterized())


def _line_collection_linestyle_from_contour(linestyle: Any, linewidth: Any) -> Any:
    if not (
        isinstance(linestyle, tuple)
        and len(linestyle) == 2
        and linestyle[1] is not None
    ):
        return linestyle

    linewidth_value = float(np.asarray(linewidth).reshape(-1)[0])
    if linewidth_value <= 0.0:
        return linestyle

    offset, dashes = linestyle
    return (
        float(offset) / linewidth_value,
        [float(dash) / linewidth_value for dash in dashes],
    )


def _signed_area(vertices: np.ndarray) -> float:
    if len(vertices) < 3:
        return 0.0
    x = vertices[:, 0]
    y = vertices[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _is_closed_path(vertices: np.ndarray) -> bool:
    if len(vertices) < 3:
        return False
    return bool(np.allclose(vertices[0], vertices[-1]))


def _orient_closed_vertices(vertices: np.ndarray, clockwise: bool) -> np.ndarray:
    area = _signed_area(vertices)
    if area == 0.0:
        return vertices
    is_clockwise = area < 0.0
    if is_clockwise == clockwise:
        return vertices
    return vertices[::-1]


def _iter_level_segments(
    contour_set: Any,
) -> Iterable[tuple[int, float, np.ndarray]]:
    levels = list(getattr(contour_set, "levels", ()))
    allsegs = list(getattr(contour_set, "allsegs", ()))
    for level_index, (level, segments) in enumerate(zip(levels, allsegs)):
        for segment in segments:
            vertices = np.asarray(segment, dtype=float)
            if vertices.ndim == 2 and vertices.shape[0] >= 2 and vertices.shape[1] >= 2:
                yield level_index, float(level), vertices[:, :2]


def _point_at_distance(vertices: np.ndarray, distance: float) -> Optional[np.ndarray]:
    deltas = np.diff(vertices, axis=0)
    lengths = np.hypot(deltas[:, 0], deltas[:, 1])
    total = float(np.sum(lengths))
    if total <= 0.0:
        return None

    target = float(np.clip(distance, 0.0, total))
    cumulative = np.cumsum(lengths)
    index = int(np.searchsorted(cumulative, target, side="right"))
    index = min(index, len(lengths) - 1)
    previous = 0.0 if index == 0 else float(cumulative[index - 1])
    segment_length = float(lengths[index])
    if segment_length <= 0.0:
        return vertices[index].copy()
    fraction = (target - previous) / segment_length
    return vertices[index] + fraction * deltas[index]


def _arrow_end_distances(total_length: float, arrow_count: int) -> np.ndarray:
    spacing = total_length / float(arrow_count)
    return (
        np.linspace(spacing, total_length, arrow_count, endpoint=True) - spacing * 0.5
    )


def _local_straightness_score(
    vertices: np.ndarray,
    distance: float,
    arrow_length: float,
) -> float:
    start = _point_at_distance(vertices, max(0.0, distance - arrow_length))
    middle = _point_at_distance(vertices, max(0.0, distance - arrow_length * 0.5))
    end = _point_at_distance(vertices, distance)
    if start is None or middle is None or end is None:
        return -np.inf

    first = middle - start
    second = end - middle
    first_length = float(np.hypot(first[0], first[1]))
    second_length = float(np.hypot(second[0], second[1]))
    chord_length = float(np.hypot(*(end - start)))
    if first_length <= 0.0 or second_length <= 0.0 or chord_length <= 0.0:
        return -np.inf

    cosine = float(
        np.clip(np.dot(first, second) / (first_length * second_length), -1.0, 1.0)
    )
    turn_angle = float(np.arccos(cosine))
    straight_ratio = min(chord_length / max(arrow_length, 1e-12), 1.0)
    return straight_ratio - 0.65 * (turn_angle / np.pi)


def _select_arrow_end_distances(
    vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
) -> np.ndarray:
    if arrow_count <= 1:
        return _arrow_end_distances(total_length, arrow_count)

    lower = min(total_length, arrow_length * 1.1)
    upper = max(lower, total_length - arrow_length * 0.5)
    if upper <= lower:
        return _arrow_end_distances(total_length, arrow_count)

    sample_count = max(arrow_count * 28, 64)
    candidates = np.linspace(lower, upper, sample_count)
    scored = [
        (
            _local_straightness_score(vertices, float(distance), arrow_length),
            float(distance),
        )
        for distance in candidates
    ]
    scored = [item for item in scored if np.isfinite(item[0])]
    if not scored:
        return _arrow_end_distances(total_length, arrow_count)

    min_spacing = max(arrow_length * 2.5, total_length / max(arrow_count * 2.2, 1.0))
    selected: List[float] = []
    for _score, distance in sorted(scored, key=lambda item: item[0], reverse=True):
        if all(abs(distance - previous) >= min_spacing for previous in selected):
            selected.append(distance)
            if len(selected) == arrow_count:
                break

    if len(selected) < arrow_count:
        for distance in _arrow_end_distances(total_length, arrow_count):
            if all(
                abs(float(distance) - previous) >= arrow_length for previous in selected
            ):
                selected.append(float(distance))
                if len(selected) == arrow_count:
                    break

    return np.asarray(sorted(selected[:arrow_count]), dtype=float)


def _segment_total_length(vertices: np.ndarray) -> float:
    deltas = np.diff(vertices, axis=0)
    return float(np.sum(np.hypot(deltas[:, 0], deltas[:, 1])))


def _build_arrow_segments_python(
    display_vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    arrow_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    head_segments = []
    head_metadata = []
    for end_distance in _select_arrow_end_distances(
        display_vertices,
        total_length,
        arrow_count,
        arrow_length,
    ):
        start = _point_at_distance(
            display_vertices,
            max(0.0, float(end_distance) - arrow_length),
        )
        end = _point_at_distance(display_vertices, float(end_distance))
        if start is None or end is None or np.allclose(start, end):
            continue

        vector = end - start
        vector_length = float(np.hypot(vector[0], vector[1]))
        if vector_length <= 0.0:
            continue
        tangent = vector / vector_length
        normal = np.array([-tangent[1], tangent[0]])
        width = vector_length * arrow_size
        base = start
        head_segments.extend(
            [
                [end, base + normal * width * 0.5],
                [end, base - normal * width * 0.5],
            ]
        )
        head_metadata.append([start, end])

    return np.asarray(head_segments, dtype=float), np.asarray(
        head_metadata, dtype=float
    )


def _build_arrow_segments(
    display_vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    arrow_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    del total_length
    return _native_build_arrow_segments(
        display_vertices,
        arrow_count,
        arrow_length,
        arrow_size,
    )


def _add_contour_arrows(
    contour_set: Any,
    ax: Any,
    arrow_count: int,
    arrow_size: float,
    arrow_length_fraction: float,
    arrow_max_length: float,
    positive_direction: str,
    arrow_color: Any,
    arrow_linewidth: Any,
    zorder: Any,
) -> List[Any]:
    arrows: List[Any] = []
    edgecolors = getattr(contour_set, "get_edgecolors", lambda: [])()
    transform = contour_set.get_transform()
    inverted_transform = transform.inverted()
    base_zorder = contour_set.get_zorder() if zorder is None else zorder

    for level_index, level, segment in _iter_level_segments(contour_set):
        closed = _is_closed_path(segment)
        direction = "forward"
        vertices = segment
        display_vertices = np.asarray(transform.transform(vertices), dtype=float)
        if closed and level != 0.0:
            positive_clockwise = positive_direction == "clockwise"
            clockwise = positive_clockwise if level > 0.0 else not positive_clockwise
            original_display_vertices = display_vertices
            display_vertices = _orient_closed_vertices(
                original_display_vertices,
                clockwise=clockwise,
            )
            if display_vertices is not original_display_vertices:
                vertices = vertices[::-1]
            direction = "clockwise" if clockwise else "counterclockwise"

        total_length = _segment_total_length(display_vertices)
        if total_length <= 0.0:
            continue

        color = arrow_color
        if color is None and len(edgecolors):
            color = edgecolors[level_index % len(edgecolors)]
        if color is None:
            color = mpl.rcParams["lines.color"]

        linewidth = arrow_linewidth
        if linewidth is None:
            linewidths = getattr(contour_set, "get_linewidths", lambda: [])()
            linewidth = (
                linewidths[level_index % len(linewidths)] if len(linewidths) else 1.0
            )
        linestyles = getattr(contour_set, "get_linestyles", lambda: [])()
        linestyle = (
            linestyles[level_index % len(linestyles)] if len(linestyles) else "solid"
        )
        linestyle = _line_collection_linestyle_from_contour(linestyle, linewidth)

        max_length_pixels = arrow_max_length * ax.figure.dpi / 72.0
        min_length_pixels = min(4.0 * ax.figure.dpi / 72.0, total_length * 0.25)
        local_arrow_length = min(
            total_length * arrow_length_fraction, max_length_pixels
        )
        local_arrow_length = max(local_arrow_length, min_length_pixels)
        line_segments = [vertices]
        display_head_segments, display_head_metadata = _build_arrow_segments(
            display_vertices,
            total_length,
            arrow_count,
            local_arrow_length,
            arrow_size,
        )
        if len(display_head_segments):
            data_head_segments = inverted_transform.transform(
                display_head_segments.reshape(-1, 2)
            ).reshape(display_head_segments.shape)
            line_segments.extend(data_head_segments)
        if len(display_head_metadata):
            data_head_metadata = inverted_transform.transform(
                display_head_metadata.reshape(-1, 2)
            ).reshape(display_head_metadata.shape)
            head_metadata = [
                (tuple(start), tuple(end)) for start, end in data_head_metadata
            ]
        else:
            head_metadata = []

        segment_linestyles = [linestyle] + ["solid"] * (len(line_segments) - 1)
        line = LineCollection(
            line_segments,
            colors=[color],
            linewidths=[linewidth],
            linestyles=segment_linestyles,
            transform=transform,
            zorder=base_zorder,
        )
        _copy_line_collection_properties(contour_set, line)
        line._skyborn_contour_level = level
        line._skyborn_contour_direction = direction
        line._skyborn_contour_kind = "arrow_contour"
        line._skyborn_contour_arrow_segments = head_metadata
        ax.add_collection(line)
        arrows.append(line)

    return arrows


def _install_arrow_remove_hook(contour_set: Any, arrows: List[Any]) -> None:
    original_remove = contour_set.remove

    def _remove(self):
        for arrow in list(arrows):
            try:
                arrow.remove()
            except ValueError:
                pass
        arrows.clear()
        return original_remove()

    contour_set.remove = types.MethodType(_remove, contour_set)
