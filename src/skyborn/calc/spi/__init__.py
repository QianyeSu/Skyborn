"""
Standardized Precipitation Index (SPI) calculation module.

This module provides efficient, vectorized calculation of the Standardized 
Precipitation Index for multi-dimensional climate datasets.
"""

from .core import standardized_precipitation_index, spi
from .xarray import spi_xarray

__all__ = ["standardized_precipitation_index", "spi", "spi_xarray"]