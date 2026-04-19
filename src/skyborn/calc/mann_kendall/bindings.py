"""Shared compiled-kernel dispatch helpers for internal MK submodules."""

import sys
from importlib import import_module
from importlib import util as importlib_util
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
from typing import Tuple

import numpy as np


def _load_core_module():
    """Best-effort loader for the compiled Mann-Kendall core extension."""
    backend_dir = Path(__file__).resolve().parent
    for probe_dir in (backend_dir / "build", backend_dir):
        for suffix in EXTENSION_SUFFIXES:
            candidate = probe_dir / f"mann_kendall_core{suffix}"
            if not candidate.exists():
                continue

            try:
                previous = sys.modules.pop("mann_kendall_core", None)
                spec = importlib_util.spec_from_file_location(
                    "mann_kendall_core", candidate
                )
                if spec is None or spec.loader is None:
                    if previous is not None:
                        sys.modules["mann_kendall_core"] = previous
                    continue
                core = importlib_util.module_from_spec(spec)
                spec.loader.exec_module(core)
                return core
            except Exception:
                if previous is not None:
                    sys.modules["mann_kendall_core"] = previous
                continue

    candidate_names = []
    if __package__:
        candidate_names.append(f"{__package__}.mann_kendall_core")
    candidate_names.append("skyborn.calc.mann_kendall.mann_kendall_core")

    for module_name in candidate_names:
        try:
            return import_module(module_name)
        except Exception:
            continue

    return None


_core_module = _load_core_module()
_score_variance_kernel = getattr(_core_module, "mk_score_var_batch", None)
_score_variance_slope_kernel = getattr(_core_module, "mk_score_var_sen_batch", None)
_sen_slope_kernel = getattr(_core_module, "sen_slope_batch", None)
_grouped_sen_slope_kernel = getattr(_core_module, "grouped_sen_slope_batch", None)
_grouped_correlated_stats_kernel = getattr(
    _core_module, "grouped_correlated_stats_batch", None
)
_partial_stats_kernel = getattr(_core_module, "partial_stats_batch", None)
_partial_stats_sen_kernel = getattr(_core_module, "partial_stats_sen_batch", None)


def _as_core_input_2d(data_2d: np.ndarray) -> np.ndarray:
    """Convert 2D data to the float64 Fortran layout expected by the core."""
    array = np.asarray(data_2d, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("Expected a 2D array for Mann-Kendall core kernels.")
    if array.flags.f_contiguous:
        return array
    return np.asfortranarray(array)


def _require_kernel(function, function_name: str):
    """Return a compiled kernel entry point or raise a clear import error."""
    if function is None:
        raise ImportError(
            "skyborn.calc.mann_kendall requires the compiled "
            f"mann_kendall_core entry '{function_name}'. "
            "Install a prebuilt wheel or build the Skyborn extensions first."
        )
    return function


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
