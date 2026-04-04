"""Compatibility shims for dataset-vector helpers during plot migration."""

from __future__ import annotations

from ._adapters import dataset_vector as _dataset_vector

_ISSUED_PLOT_WARNINGS = _dataset_vector._ISSUED_PLOT_WARNINGS
_warn_plot_once = _dataset_vector._warn_plot_once
_apply_dataset_isel = _dataset_vector._apply_dataset_isel
_get_plot_dataarray = _dataset_vector._get_plot_dataarray
_transpose_2d_dataarray_to_dims = _dataset_vector._transpose_2d_dataarray_to_dims
_filled_scalar_field_array = _dataset_vector._filled_scalar_field_array
_extract_curly_vector_dataset_source = (
    _dataset_vector._extract_curly_vector_dataset_source
)
_prepare_dataset_style_field = _dataset_vector._prepare_dataset_style_field
