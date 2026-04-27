"""Native C helpers for Skyborn NCL-like curly-vector tracing.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from ._ncl_curly_core import (
    compute_display_cell_valid,
    generate_cell_candidates,
    sample_display_grid,
    sample_grid_field,
    sample_grid_field_array,
    thin_display_candidates,
    thin_mapped_candidates,
    trace_ncl_direction,
    validate_display_curve,
)

__all__ = [
    "compute_display_cell_valid",
    "generate_cell_candidates",
    "sample_display_grid",
    "trace_ncl_direction",
    "sample_grid_field",
    "sample_grid_field_array",
    "thin_display_candidates",
    "thin_mapped_candidates",
    "validate_display_curve",
]
