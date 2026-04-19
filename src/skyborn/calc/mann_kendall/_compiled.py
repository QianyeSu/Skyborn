"""Shared compiled-core helpers for internal Mann-Kendall modules."""

from pathlib import Path
from typing import Callable, Dict, Iterable, MutableMapping, Optional

import numpy as np

_KERNEL_ENTRY_POINTS = {
    "_score_variance_kernel": "mk_score_var_batch",
    "_score_variance_slope_kernel": "mk_score_var_sen_batch",
    "_yue_wang_slope_kernel": "mk_yue_wang_score_var_sen_batch",
    "_hamed_rao_slope_kernel": "mk_hamed_rao_score_var_sen_batch",
    "_sen_slope_kernel": "sen_slope_batch",
    "_grouped_sen_slope_kernel": "grouped_sen_slope_batch",
    "_grouped_correlated_stats_kernel": "grouped_correlated_stats_batch",
    "_partial_stats_kernel": "partial_stats_batch",
    "_partial_stats_sen_kernel": "partial_stats_sen_batch",
}


def _load_core_module_from(
    *,
    file_path: str,
    package_name: Optional[str],
    extension_suffixes: Iterable[str],
    path_type: type[Path],
    importlib_util_module,
    import_module_func: Callable[[str], object],
    sys_modules: MutableMapping[str, object],
):
    """Best-effort loader for the compiled ``mann_kendall_core`` extension."""
    backend_dir = path_type(file_path).resolve().parent
    for probe_dir in (backend_dir / "build", backend_dir):
        for suffix in extension_suffixes:
            candidate = probe_dir / f"mann_kendall_core{suffix}"
            if not candidate.exists():
                continue

            previous = sys_modules.pop("mann_kendall_core", None)
            try:
                spec = importlib_util_module.spec_from_file_location(
                    "mann_kendall_core", candidate
                )
                if spec is None or spec.loader is None:
                    if previous is not None:
                        sys_modules["mann_kendall_core"] = previous
                    continue
                core = importlib_util_module.module_from_spec(spec)
                spec.loader.exec_module(core)
                return core
            except Exception:
                if previous is not None:
                    sys_modules["mann_kendall_core"] = previous
                continue

    candidate_names = []
    if package_name:
        candidate_names.append(f"{package_name}.mann_kendall_core")
    candidate_names.append("skyborn.calc.mann_kendall.mann_kendall_core")

    for module_name in candidate_names:
        try:
            return import_module_func(module_name)
        except Exception:
            continue

    return None


def _resolve_kernel_handles(core_module) -> Dict[str, object]:
    """Return module-level kernel handles from the compiled extension."""
    return {
        variable_name: getattr(core_module, entry_name, None)
        for variable_name, entry_name in _KERNEL_ENTRY_POINTS.items()
    }


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
