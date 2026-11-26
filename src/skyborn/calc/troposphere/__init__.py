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

from ..._package_alias import register_module_alias
from . import xarray
from .core import trop_wmo, trop_wmo_profile, tropopause_height

register_module_alias(
    "skyborn.calc.troposphere.tropopause_height",
    "skybornbackend.calc.troposphere.tropopause_height",
)

__all__ = ["trop_wmo", "trop_wmo_profile"]
