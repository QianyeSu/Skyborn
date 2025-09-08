"""
Skyborn Tropical Cyclone Potential Intensity Module

This module provides optimized Fortran implementations for calculating
tropical cyclone potential intensity (PI), with high-level Python
interfaces for multi-dimensional meteorological data handling.

Features:
- High-performance Fortran implementation with OpenMP support
- Support for 1D profiles, 3D gridded, and 4D time series data
- Automatic dimension detection and unit conversion
- Direct processing of atmospheric profile data
- XArray integration with metadata preservation

Physical Basis:
The potential intensity calculation is based on Emanuel's thermodynamic theory:
- Thermodynamic disequilibrium between ocean and atmosphere
- Environmental temperature and humidity profiles
- Sea surface temperature constraints

Equations:
- PI² = (Ck/Cd) * (SST/T₀) * (CAPE* - CAPE)
- Where Ck/Cd is exchange coefficient ratio, CAPE* is saturation CAPE

Input Requirements:
- sst: Sea surface temperature [K]
- psl: Sea level pressure [Pa]
- pressure_levels: Atmospheric pressure levels [mb or hPa]
- temperature: Temperature profiles [K]
- mixing_ratio: Water vapor mixing ratio [kg/kg]

Output:
- min_pressure: Minimum central pressure [mb]
- pi: Potential intensity (maximum wind speed) [m/s]
- error_flag: Error status (1 = success, other values = error)

Examples
--------
# NumPy interface
>>> from skyborn.calc.GPI import potential_intensity
>>> min_p, pi, err = potential_intensity(sst_val, psl_val, p_levels, temp_profile, mixr_profile)

# XArray interface
>>> from skyborn.calc.GPI.xarray import potential_intensity
>>> result = potential_intensity(sst, psl, pressure_levels, temperature, mixing_ratio)
>>> print(f"PI: {result.pi.values} m/s, Min P: {result.min_pressure.values} mb")

References
----------
Emanuel, K. (1995). Sensitivity of tropical cyclones to surface exchange
coefficients and a revised steady-state model incorporating eye dynamics.
Journal of the Atmospheric Sciences, 52(22), 3969-3976.

Bister, M., & Emanuel, K. A. (2002). Low frequency variability of tropical
cyclone potential intensity 1. Interannual to interdecadal variability.
Journal of Geophysical Research, 107(D24), 4801.
"""

# Import xarray submodule for user access
from . import xarray

# Import high-level Python interface (user-facing API)
from .interface import potential_intensity

__version__ = "1.0.0"

# Define public API - only high-level numpy interface
__all__ = [
    "potential_intensity",
]
