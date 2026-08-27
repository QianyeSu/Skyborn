"""Arrowed-contour rendering helpers."""

from __future__ import annotations

import types
from typing import Any, Iterable, List, Optional

import matplotlib as mpl
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection

from ..contour_core import build_arrow_segments as _native_build_arrow_segments

# Import C-accelerated functions
from .contour_arrows_core import local_straightness_score as _local_straightness_score_c
from .contour_arrows_core import point_at_distance as _point_at_distance_c
from .contour_arrows_core import (
    select_arrow_end_distances as _select_arrow_end_distances_c,
)


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


def _orient_closed_vertices(
    vertices: np.ndarray, clockwise: bool, data_to_display_transform: Any
) -> np.ndarray:
    """Orient closed path vertices to the requested direction in display space.

    Parameters
    ----------
    vertices : ndarray
        Path vertices in data coordinates
    clockwise : bool
        Target orientation (True for clockwise in display coordinates)
    data_to_display_transform : matplotlib.transforms.Transform
        Transform from data to display coordinates

    Returns
    -------
    ndarray
        Vertices in data coordinates, potentially reversed to match target orientation
    """
    # Transform to display coordinates for orientation check
    display_vertices = data_to_display_transform.transform(vertices)
    area = _signed_area(display_vertices)
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
    """Find point at given arc-length distance along path (C-accelerated)."""
    return _point_at_distance_c(vertices, distance)


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
    """Score straightness at given distance for arrow placement (C-accelerated)."""
    return _local_straightness_score_c(vertices, distance, arrow_length)


def _select_arrow_end_distances(
    vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    min_spacing_override: float | None = None,
) -> np.ndarray:
    """Select optimal arrow placement distances (C-accelerated)."""
    return _select_arrow_end_distances_c(
        vertices, total_length, arrow_count, arrow_length, min_spacing_override
    )


def _segment_total_length(vertices: np.ndarray) -> float:
    deltas = np.diff(vertices, axis=0)
    return float(np.sum(np.hypot(deltas[:, 0], deltas[:, 1])))


def _build_arrow_triangles_python(
    display_vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    arrow_size: float,
    min_spacing_override: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build filled triangle arrowheads instead of V-shaped line segments."""
    triangles = []
    head_metadata = []
    for end_distance in _select_arrow_end_distances(
        display_vertices,
        total_length,
        arrow_count,
        arrow_length,
        min_spacing_override,
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

        # Create triangle vertices: tip, left barb, right barb
        triangle = np.array(
            [
                end,
                base + normal * width * 0.5,
                base - normal * width * 0.5,
            ]
        )
        triangles.append(triangle)
        head_metadata.append([start, end])

    return triangles, np.asarray(head_metadata, dtype=float)


def _build_arrow_segments_python(
    display_vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    arrow_size: float,
    min_spacing_override: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    head_segments = []
    head_metadata = []
    for end_distance in _select_arrow_end_distances(
        display_vertices,
        total_length,
        arrow_count,
        arrow_length,
        min_spacing_override,
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


def _build_arrow_segments_swept_python(
    display_vertices: np.ndarray,
    total_length: float,
    arrow_count: int,
    arrow_length: float,
    arrow_size: float,
    arrow_angle: float = 22.5,
    min_spacing_override: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build NCL-style swept-back arrowhead line segments.

    Mirrors the geometry in NCAR Graphics' ``drwvec.f``: each barb starts at
    the arrow tip and sweeps back at ``arrow_angle`` degrees from the
    reversed shaft direction, rather than fanning out perpendicular to the
    shaft from its base. This produces a narrower, sharper head than the
    flat-backed V used by the default "line" style.
    """
    theta = np.radians(arrow_angle)
    cos_t = float(np.cos(theta))
    sin_t = float(np.sin(theta))

    head_segments = []
    head_metadata = []
    for end_distance in _select_arrow_end_distances(
        display_vertices,
        total_length,
        arrow_count,
        arrow_length,
        min_spacing_override,
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
        barb_length = vector_length * arrow_size

        barb_left = end - barb_length * (cos_t * tangent + sin_t * normal)
        barb_right = end - barb_length * (cos_t * tangent - sin_t * normal)
        head_segments.extend(
            [
                [end, barb_left],
                [end, barb_right],
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
    arrow_length_points: float | None,
    arrow_max_length: float,
    positive_direction: str,
    arrow_color: Any,
    arrow_linewidth: Any,
    zorder: Any,
    arrow_style: str = "swept",
    arrow_angle: float = 22.5,
    arrow_min_spacing: float | None = None,
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
        if closed and level != 0.0:
            positive_clockwise = positive_direction == "clockwise"
            clockwise = positive_clockwise if level > 0.0 else not positive_clockwise
            vertices = _orient_closed_vertices(
                vertices,
                clockwise=clockwise,
                data_to_display_transform=transform,
            )
            direction = "clockwise" if clockwise else "counterclockwise"

        display_vertices = np.asarray(transform.transform(vertices), dtype=float)

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

        # Priority: arrow_length_points (absolute) > arrow_length_fraction (relative)
        if arrow_length_points is not None:
            # Absolute arrow length in display pixels
            local_arrow_length = arrow_length_points * ax.figure.dpi / 72.0
        else:
            # Relative arrow length with max cap
            max_length_pixels = arrow_max_length * ax.figure.dpi / 72.0
            local_arrow_length = min(
                total_length * arrow_length_fraction, max_length_pixels
            )

        # Ensure minimum length for visibility
        min_length_pixels = min(4.0 * ax.figure.dpi / 72.0, total_length * 0.25)
        local_arrow_length = max(local_arrow_length, min_length_pixels)

        # Build arrows based on style
        if arrow_style == "filled":
            # Use filled triangles for better visual appearance
            display_triangles, display_head_metadata = _build_arrow_triangles_python(
                display_vertices,
                total_length,
                arrow_count,
                local_arrow_length,
                arrow_size,
                arrow_min_spacing,
            )

            # Draw contour line
            line = LineCollection(
                [vertices],
                colors=[color],
                linewidths=[linewidth],
                linestyles=[linestyle],
                transform=transform,
                zorder=base_zorder,
            )
            _copy_line_collection_properties(contour_set, line)
            line._skyborn_contour_level = level
            line._skyborn_contour_direction = direction
            line._skyborn_contour_kind = "arrow_contour"
            ax.add_collection(line)
            arrows.append(line)

            # Draw filled triangle arrowheads
            if len(display_triangles):
                data_triangles = [
                    inverted_transform.transform(tri) for tri in display_triangles
                ]
                poly = PolyCollection(
                    data_triangles,
                    facecolors=[color],
                    edgecolors=[color],
                    linewidths=[linewidth * 0.5],
                    transform=transform,
                    zorder=base_zorder + 0.1,
                )
                _copy_line_collection_properties(contour_set, poly)
                poly._skyborn_contour_level = level
                poly._skyborn_contour_direction = direction
                poly._skyborn_contour_kind = "arrow_contour_heads"
                ax.add_collection(poly)
                arrows.append(poly)

            # Store metadata
            if len(display_head_metadata):
                data_head_metadata = inverted_transform.transform(
                    display_head_metadata.reshape(-1, 2)
                ).reshape(display_head_metadata.shape)
                head_metadata = [
                    (tuple(start), tuple(end)) for start, end in data_head_metadata
                ]
            else:
                head_metadata = []
            line._skyborn_contour_arrow_segments = head_metadata
        else:
            # Line-based styles: flat-backed V ("line") or NCL-style
            # swept-back barbs ("swept")
            line_segments = [vertices]
            if arrow_style == "swept":
                display_head_segments, display_head_metadata = (
                    _build_arrow_segments_swept_python(
                        display_vertices,
                        total_length,
                        arrow_count,
                        local_arrow_length,
                        arrow_size,
                        arrow_angle,
                        arrow_min_spacing,
                    )
                )
            else:
                display_head_segments, display_head_metadata = (
                    _build_arrow_segments_python(
                        display_vertices,
                        total_length,
                        arrow_count,
                        local_arrow_length,
                        arrow_size,
                        arrow_min_spacing,
                    )
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


def _path_cumulative_lengths(vertices: np.ndarray) -> np.ndarray:
    deltas = np.diff(vertices, axis=0)
    lengths = np.hypot(deltas[:, 0], deltas[:, 1])
    return np.concatenate([[0.0], np.cumsum(lengths)])


def _nearest_arclength(
    vertices: np.ndarray, cumulative: np.ndarray, point: np.ndarray
) -> float:
    distances = np.hypot(vertices[:, 0] - point[0], vertices[:, 1] - point[1])
    index = int(np.argmin(distances))
    return float(cumulative[index])


def _label_anchor_away_from_arrows(
    display_path: np.ndarray,
    arrow_midpoints_display: List[np.ndarray],
    closed: bool,
) -> Optional[np.ndarray]:
    total_length = _segment_total_length(display_path)
    if total_length <= 0.0:
        return None
    cumulative = _path_cumulative_lengths(display_path)

    # For open paths, keep candidates away from the endpoints: the point
    # farthest from every arrow midpoint is otherwise almost always right at
    # one of the ends (e.g. where the line exits the axes), which pushes
    # labels off to the plot border instead of routing them around arrows.
    if closed:
        lower, upper = 0.0, total_length
    else:
        margin = total_length * 0.35
        lower, upper = margin, max(margin, total_length - margin)

    if not arrow_midpoints_display:
        target = (lower + upper) / 2.0
    else:
        occupied = [
            _nearest_arclength(display_path, cumulative, mid)
            for mid in arrow_midpoints_display
        ]
        candidates = np.linspace(lower, upper, 200)

        def min_gap(distance: float) -> float:
            gaps = [abs(distance - occupied_distance) for occupied_distance in occupied]
            if closed:
                gaps = [min(gap, total_length - gap) for gap in gaps]
            return min(gaps)

        target = max(candidates, key=min_gap)

    return _point_at_distance(display_path, float(target))


def _arrow_contour_label_positions(contour_set: Any) -> List[tuple[float, float]]:
    """Compute ``clabel`` anchor points that avoid ``arrow_contour`` arrowheads.

    For each contour line drawn by ``arrow_contour``, picks the point along
    that line (by arc length in display coordinates) farthest from the
    midpoints of its arrowheads, so a label placed there won't land on top
    of an arrow.
    """
    arrow_artists = getattr(contour_set, "_skyborn_arrow_contour_artists", None) or []
    if not arrow_artists:
        return []

    transform = contour_set.get_transform()
    inverted_transform = transform.inverted()
    positions: List[tuple[float, float]] = []

    for line in arrow_artists:
        try:
            get_segments = getattr(line, "get_segments", None)
            segments = get_segments() if get_segments is not None else []
            if not segments:
                continue
            main_path = np.asarray(segments[0], dtype=float)
            if main_path.shape[0] < 2:
                continue
            display_path = transform.transform(main_path)
            closed = _is_closed_path(main_path)

            arrow_spans = getattr(line, "_skyborn_contour_arrow_segments", None) or []
            midpoints_display = []
            for start, end in arrow_spans:
                start_arr = np.asarray(start, dtype=float)
                end_arr = np.asarray(end, dtype=float)
                midpoint_data = (start_arr + end_arr) / 2.0
                midpoints_display.append(transform.transform(midpoint_data))

            anchor_display = _label_anchor_away_from_arrows(
                display_path, midpoints_display, closed
            )
            if anchor_display is None:
                continue
            anchor_data = inverted_transform.transform(anchor_display)
            positions.append((float(anchor_data[0]), float(anchor_data[1])))
        except Exception as e:
            # Silently skip problematic contour lines rather than failing
            # the entire label placement. This can happen if the line geometry
            # is degenerate or if coordinate transforms fail.
            import warnings

            warnings.warn(
                f"Failed to compute label position for contour line: {e}",
                RuntimeWarning,
                stacklevel=3,
            )
            continue

    return positions


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
