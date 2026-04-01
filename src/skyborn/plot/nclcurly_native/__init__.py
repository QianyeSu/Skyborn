"""Native C helpers for NCL-like curved-vector tracing."""

from ._ncl_curly_core import (
    sample_grid_field,
    sample_grid_field_array,
    thin_mapped_candidates,
    trace_ncl_direction,
)

__all__ = [
    "trace_ncl_direction",
    "sample_grid_field",
    "sample_grid_field_array",
    "thin_mapped_candidates",
]
