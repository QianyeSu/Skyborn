"""Compiled-kernel dispatch helpers for internal MK submodules."""

import sys
from importlib import import_module
from importlib import util as importlib_util
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
from typing import Tuple

import numpy as np

from ._compiled import (
    _as_core_input_2d,
    _load_core_module_from,
    _require_kernel,
    _resolve_kernel_handles,
)


def _load_core_module():
    """Best-effort loader for the compiled Mann-Kendall core extension."""
    return _load_core_module_from(
        file_path=__file__,
        package_name=__package__,
        extension_suffixes=EXTENSION_SUFFIXES,
        path_type=Path,
        importlib_util_module=importlib_util,
        import_module_func=import_module,
        sys_modules=sys.modules,
    )


_core_module = _load_core_module()
_kernel_handles = _resolve_kernel_handles(_core_module)
_score_variance_kernel = _kernel_handles["_score_variance_kernel"]
_score_variance_slope_kernel = _kernel_handles["_score_variance_slope_kernel"]
_sen_slope_kernel = _kernel_handles["_sen_slope_kernel"]
_grouped_sen_slope_kernel = _kernel_handles["_grouped_sen_slope_kernel"]
_grouped_correlated_stats_kernel = _kernel_handles["_grouped_correlated_stats_kernel"]
_partial_stats_kernel = _kernel_handles["_partial_stats_kernel"]
_partial_stats_sen_kernel = _kernel_handles["_partial_stats_sen_kernel"]


def _score_variance_batch(
    data_2d: np.ndarray, modified: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Run the compiled 2D score/variance kernel."""
    kernel = _require_kernel(_score_variance_kernel, "mk_score_var_batch")
    s_values, var_values = kernel(_as_core_input_2d(data_2d), int(modified))
    return np.asarray(s_values, dtype=np.float64), np.asarray(
        var_values, dtype=np.float64
    )


def _sen_slope_batch(data_2d: np.ndarray) -> np.ndarray:
    """Run the compiled 2D Theil-Sen slope kernel."""
    kernel = _require_kernel(_sen_slope_kernel, "sen_slope_batch")
    slopes = kernel(_as_core_input_2d(data_2d))
    return np.asarray(slopes, dtype=np.float64)


def _grouped_sen_slope_batch(data_2d: np.ndarray, period: int) -> np.ndarray:
    """Run the compiled grouped 2D Theil-Sen slope kernel."""
    kernel = _require_kernel(_grouped_sen_slope_kernel, "grouped_sen_slope_batch")
    slopes = kernel(_as_core_input_2d(data_2d), int(period))
    return np.asarray(slopes, dtype=np.float64)


def _grouped_correlated_stats_batch(
    data_2d: np.ndarray, period: int
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Run the compiled grouped correlated-statistics kernel."""
    kernel = _require_kernel(
        _grouped_correlated_stats_kernel, "grouped_correlated_stats_batch"
    )
    s_values, var_values, denom = kernel(_as_core_input_2d(data_2d), int(period))
    return (
        np.asarray(s_values, dtype=np.float64),
        np.asarray(var_values, dtype=np.float64),
        float(denom),
    )


def _partial_stats_batch(
    response_2d: np.ndarray, covariate_2d: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run the compiled partial-MK statistics kernel."""
    kernel = _require_kernel(_partial_stats_kernel, "partial_stats_batch")
    s_values, var_values, tau_values = kernel(
        _as_core_input_2d(response_2d),
        _as_core_input_2d(covariate_2d),
    )
    return (
        np.asarray(s_values, dtype=np.float64),
        np.asarray(var_values, dtype=np.float64),
        np.asarray(tau_values, dtype=np.float64),
    )


def _partial_stats_sen_batch(
    response_2d: np.ndarray, covariate_2d: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run the compiled partial-MK statistics + Sen-slope kernel."""
    kernel = _require_kernel(_partial_stats_sen_kernel, "partial_stats_sen_batch")
    s_values, var_values, tau_values, slopes = kernel(
        _as_core_input_2d(response_2d),
        _as_core_input_2d(covariate_2d),
    )
    return (
        np.asarray(s_values, dtype=np.float64),
        np.asarray(var_values, dtype=np.float64),
        np.asarray(tau_values, dtype=np.float64),
        np.asarray(slopes, dtype=np.float64),
    )


def _score_variance_slope_batch(
    data_2d: np.ndarray, modified: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run the compiled 2D combined score/variance/slope kernel."""
    kernel = _require_kernel(_score_variance_slope_kernel, "mk_score_var_sen_batch")
    s_values, var_values, slopes = kernel(_as_core_input_2d(data_2d), int(modified))
    return (
        np.asarray(s_values, dtype=np.float64),
        np.asarray(var_values, dtype=np.float64),
        np.asarray(slopes, dtype=np.float64),
    )
