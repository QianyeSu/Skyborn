"""
Troposphere calculation module.

This module provides WMO tropopause height calculation functionality
for atmospheric science research.

Main Functions:
    trop_wmo : WMO tropopause calculation for numpy arrays
    trop_wmo_profile : WMO tropopause calculation for 1D profiles

For xarray interface:
    >>> from skyborn.calc.troposphere.xarray import trop_wmo
"""

from . import xarray
from .tropopause import trop_wmo, trop_wmo_profile

__all__ = ["trop_wmo", "trop_wmo_profile"]
