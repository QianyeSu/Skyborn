"""Low-level vector artist sizing helpers for Skyborn plots."""

from __future__ import annotations

import numpy as np


def _ncl_arrow_edge_size_px(magnitude_value, max_mag, min_edge_px, max_edge_px):
    max_mag = max(float(max_mag), 1e-12)
    vmf = np.clip(float(magnitude_value) / max_mag, 0.0, 1.0)
    return float(min_edge_px) + (float(max_edge_px) - float(min_edge_px)) * vmf


def _resolve_open_arrow_size(edge_size_px):
    edge_size_px = max(float(edge_size_px), 1.0)
    half_angle = 0.5
    shaft_length_px = edge_size_px * np.cos(half_angle)
    head_width_px = 2.0 * edge_size_px * np.sin(half_angle)
    return shaft_length_px, head_width_px
