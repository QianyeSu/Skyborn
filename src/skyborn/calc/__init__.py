"""Calculation routines for Skyborn.

The calculation package exposes historical convenience imports lazily. This
keeps ``import skyborn.calc`` and imports of leaf modules such as
``skyborn.calc.calculations`` from pulling in compiled extension stacks that are
only needed by specific diagnostics.

This module contains various calculation functions including:
- Statistical calculations and linear regression
- Emergent constraint methods
- PDF calculations and analysis
- Mann-Kendall trend analysis
- WMO tropopause calculation (with Fortran extensions)
- Geostrophic wind calculation (with SIMD-optimized Fortran extensions)
- Tropical cyclone potential intensity calculation (with Fortran extensions)
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_SUBMODULE_EXPORTS = {
    "GPI",
    "dcape",
    "geostrophic",
    "growth_rate",
    "mann_kendall",
    "troposphere",
    "ventilation",
}

_OBJECT_EXPORTS = {
    # Statistical calculations.
    "calculate_potential_temperature": (
        "skyborn.calc.calculations",
        "calculate_potential_temperature",
    ),
    "convert_longitude_range": (
        "skyborn.calc.calculations",
        "convert_longitude_range",
    ),
    "kendall_correlation": ("skyborn.calc.calculations", "kendall_correlation"),
    "linear_regression": ("skyborn.calc.calculations", "linear_regression"),
    "pearson_correlation": ("skyborn.calc.calculations", "pearson_correlation"),
    "spatial_correlation": ("skyborn.calc.calculations", "spatial_correlation"),
    "spearman_correlation": ("skyborn.calc.calculations", "spearman_correlation"),
    "calculate_dcape": ("skyborn.calc.dcape", "calculate_dcape"),
    # Emergent constraints.
    "calc_GAUSSIAN_PDF": (
        "skyborn.calc.emergent_constraints",
        "calc_GAUSSIAN_PDF",
    ),
    "calc_PDF_EC": ("skyborn.calc.emergent_constraints", "calc_PDF_EC"),
    "calc_PDF_EC_PRIOR": (
        "skyborn.calc.emergent_constraints",
        "calc_PDF_EC_PRIOR",
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
    # Compiled or extension-backed diagnostics.
    "GeostrophicWind": ("skyborn.calc.geostrophic", "GeostrophicWind"),
    "geostrophic_speed": ("skyborn.calc.geostrophic", "geostrophic_speed"),
    "geostrophic_uv": ("skyborn.calc.geostrophic", "geostrophic_uv"),
    "geostrophic_wind": ("skyborn.calc.geostrophic", "geostrophic_wind"),
    "potential_intensity": ("skyborn.calc.GPI.core", "potential_intensity"),
    # Tropical cyclone ventilation index (Chavas et al. 2025).
    "absolute_vorticity_850": (
        "skyborn.calc.ventilation",
        "absolute_vorticity_850",
    ),
    "air_sea_disequilibrium": (
        "skyborn.calc.ventilation",
        "air_sea_disequilibrium",
    ),
    "entropy_deficit": ("skyborn.calc.ventilation", "entropy_deficit"),
    "genesis_potential_index": (
        "skyborn.calc.ventilation",
        "genesis_potential_index",
    ),
    "ventilation_components": (
        "skyborn.calc.ventilation",
        "ventilation_components",
    ),
    "ventilated_pi": ("skyborn.calc.ventilation", "ventilated_pi"),
    "ventilation_index": ("skyborn.calc.ventilation", "ventilation_index"),
    "vertical_wind_shear": ("skyborn.calc.ventilation", "vertical_wind_shear"),
    "baroc_growth_rate": ("skyborn.calc.growth_rate", "baroc_growth_rate"),
    "barot_growth_rate": ("skyborn.calc.growth_rate", "barot_growth_rate"),
    "trop_wmo": ("skyborn.calc.troposphere", "trop_wmo"),
    "trop_wmo_profile": ("skyborn.calc.troposphere", "trop_wmo_profile"),
    # Mann-Kendall trend analysis.
    "mann_kendall_multidim": (
        "skyborn.calc.mann_kendall",
        "mann_kendall_multidim",
    ),
    "mann_kendall_test": ("skyborn.calc.mann_kendall", "mann_kendall_test"),
    "mann_kendall_xarray": ("skyborn.calc.mann_kendall", "mann_kendall_xarray"),
    "mk_multidim": ("skyborn.calc.mann_kendall", "mk_multidim"),
    "mk_test": ("skyborn.calc.mann_kendall", "mk_test"),
    "partial_mann_kendall_multidim": (
        "skyborn.calc.mann_kendall",
        "partial_mann_kendall_multidim",
    ),
    "partial_mann_kendall_test": (
        "skyborn.calc.mann_kendall",
        "partial_mann_kendall_test",
    ),
    "partial_mann_kendall_xarray": (
        "skyborn.calc.mann_kendall",
        "partial_mann_kendall_xarray",
    ),
    "partial_multidim": ("skyborn.calc.mann_kendall", "partial_multidim"),
    "partial_test": ("skyborn.calc.mann_kendall", "partial_test"),
    "trend_analysis": ("skyborn.calc.mann_kendall", "trend_analysis"),
}

__all__ = sorted((*_SUBMODULE_EXPORTS, *_OBJECT_EXPORTS))


def __getattr__(name: str) -> Any:
    """Resolve calculation package compatibility exports on first access."""

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
