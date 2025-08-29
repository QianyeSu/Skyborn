"""
Skyborn Geostrophic Wind Calculation Module

This module provides optimized Fortran implementations for calculating
geostrophic wind components from geopotential height fields, with high-level
Python interfaces for multi-dimensional data handling.

Features:
- SIMD-optimized Fortran implementation with OpenMP support
- Support for 2D, 3D, and 4D data arrays with shape (nlat, mlon, ...)
- Automatic handling of both south-to-north and north-to-south latitude ordering
- Direct processing without unnecessary array transposes
- Compatible with various longitude boundary conditions (cyclic/non-cyclic)
- Integration with windspharm prep_data/recover_data for flexible dimension handling

Low-level Fortran Functions:
- z2geouv: Main 2D geostrophic wind calculation
- z2geouv_3d: 3D version for multiple levels/times
- z2geouv_4d: 4D version for level×time datasets
- zuvnew: Internal function for north-to-south data reordering
- z2guv: Core SIMD-optimized calculation routine

High-level Python Interface:
- geostrophic_wind: Main interface function with automatic dimension handling
- GeostrophicWind: Class-based interface with derived quantities
- geostrophic_uv, geostrophic_speed, geostrophic_direction: Convenience functions

Physical Equations:
- Geostrophic balance: f * Vg = -∇(Z) x k
- ug = -(g/f) * dZ/dy   (zonal wind component)
- vg =  (g/f) * dZ/dx   (meridional wind component)

where:
- g = 9.80616 m/s² (gravity)
- f = 2Ω sin(lat) (Coriolis parameter)
- Ω = 7.292×10⁻⁵ s⁻¹ (Earth rotation rate)
- Z = geopotential height [gpm]

Examples
--------
# Simple 2D calculation
>>> ug, vg = geostrophic_wind(z, glon, glat, 'yx')

# Multi-dimensional data
>>> z4d = np.random.randn(12, 17, 73, 144)  # (time, level, lat, lon)
>>> ug, vg = geostrophic_wind(z4d, glon, glat, 'tzyx')

# Class-based interface with derived quantities
>>> gw = GeostrophicWind(z, glon, glat, 'tzyx')
>>> speed = gw.speed()
"""

# Import high-level Python interface (user-facing API)
from .interface import (
    GeostrophicWind,
    geostrophic_speed,
    geostrophic_uv,
    geostrophic_wind,
)

__version__ = "1.0.0"

# Define public API - only high-level functions for users
__all__ = [
    "geostrophic_wind",
    "GeostrophicWind",
    "geostrophic_uv",
    "geostrophic_speed",
]
