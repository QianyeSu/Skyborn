"""Top-level package for Skyborn.

The package root intentionally stays lightweight. Public convenience exports are
resolved lazily so that ``import skyborn`` does not import optional plotting,
xarray, or compiled-extension stacks before they are needed.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

__version__ = "0.4.1"

_SUBMODULE_EXPORTS = {
    "ROF",
    "calc",
    "gridfill",
    "interp",
    "plot",
    "spharm",
    "windspharm",
}

_OBJECT_EXPORTS = {
    # Calculation convenience exports kept for backward compatibility.
    "calc_GAUSSIAN_PDF": (
        "skyborn.calc.emergent_constraints",
        "calc_GAUSSIAN_PDF",
    ),
    "calc_PDF_EC": ("skyborn.calc.emergent_constraints", "calc_PDF_EC"),
    "calc_PDF_EC_PRIOR": (
        "skyborn.calc.emergent_constraints",
        "calc_PDF_EC_PRIOR",
    ),
    "calculate_potential_temperature": (
        "skyborn.calc.calculations",
        "calculate_potential_temperature",
    ),
    "convert_longitude_range": (
        "skyborn.calc.calculations",
        "convert_longitude_range",
    ),
    "emergent_constraint_posterior": (
        "skyborn.calc.emergent_constraints",
        "emergent_constraint_posterior",
    ),
    "emergent_constraint_prior": (
        "skyborn.calc.emergent_constraints",
        "emergent_constraint_prior",
    ),
    "find_std_from_PDF": ("skyborn.calc.emergent_constraints", "find_std_from_PDF"),
    "gaussian_pdf": ("skyborn.calc.emergent_constraints", "gaussian_pdf"),
    "geostrophic": ("skyborn.calc", "geostrophic"),
    "kendall_correlation": ("skyborn.calc.calculations", "kendall_correlation"),
    "linear_regression": ("skyborn.calc.calculations", "linear_regression"),
    "mann_kendall_multidim": (
        "skyborn.calc.mann_kendall",
        "mann_kendall_multidim",
    ),
    "mann_kendall_test": ("skyborn.calc.mann_kendall", "mann_kendall_test"),
    "mann_kendall_xarray": ("skyborn.calc.mann_kendall", "mann_kendall_xarray"),
    "mk_multidim": ("skyborn.calc.mann_kendall", "mk_multidim"),
    "mk_test": ("skyborn.calc.mann_kendall", "mk_test"),
    "pearson_correlation": ("skyborn.calc.calculations", "pearson_correlation"),
    "potential_intensity": ("skyborn.calc.GPI.core", "potential_intensity"),
    "spearman_correlation": ("skyborn.calc.calculations", "spearman_correlation"),
    "trend_analysis": ("skyborn.calc.mann_kendall", "trend_analysis"),
    "troposphere": ("skyborn.calc", "troposphere"),
    # Top-level science utilities.
    "granger_causality": ("skyborn.causality", "granger_causality"),
    "liang_causality": ("skyborn.causality", "liang_causality"),
    "calculate_gradient": ("skyborn.gradients", "calculate_gradient"),
    "calculate_meridional_gradient": (
        "skyborn.gradients",
        "calculate_meridional_gradient",
    ),
    "calculate_vertical_gradient": (
        "skyborn.gradients",
        "calculate_vertical_gradient",
    ),
    "calculate_zonal_gradient": ("skyborn.gradients", "calculate_zonal_gradient"),
    # Historical gridfill aliases.
    "fill": ("skyborn.gridfill.gridfill", "fill"),
    "gridfill_fill": ("skyborn.gridfill.gridfill", "fill"),
}

__all__ = sorted((*_SUBMODULE_EXPORTS, *_OBJECT_EXPORTS, "__version__"))


def __getattr__(name: str) -> Any:
    """Resolve legacy top-level exports on first access."""

    if name in _SUBMODULE_EXPORTS:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    if name in _OBJECT_EXPORTS:
        module_name, object_name = _OBJECT_EXPORTS[name]
        value = getattr(import_module(module_name), object_name)
        globals()[name] = value
        return value

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Include lazy exports in interactive introspection."""

    return sorted(set(globals()) | set(__all__))
