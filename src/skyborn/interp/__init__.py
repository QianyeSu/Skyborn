"""
Interpolation and regridding utilities for skyborn.

This module provides various interpolation methods including:
- Nearest neighbor interpolation
- Bilinear interpolation
- Conservative interpolation
- Hybrid-sigma to pressure level interpolation
- Curvilinear grid interpolation (RCM/WRF/NARR)
- Grid to/from triple format conversion
"""

from .interpolation import (
    interp_hybrid_to_pressure,
    interp_multidim,
    interp_sigma_to_hybrid,
)

# Import curvilinear grid interpolation functions (require compiled Fortran modules)
from .rcm2points import rcm2points
from .rcm2rgrid import rcm2rgrid, rgrid2rcm
from .regridding import (
    BilinearRegridder,
    ConservativeRegridder,
    Grid,
    NearestRegridder,
    Regridder,
    nearest_neighbor_indices,
    regrid_dataset,
)
from .triple_to_grid import grid_to_triple, triple_to_grid

__all__ = [
    "interp_hybrid_to_pressure",
    "interp_multidim",
    "interp_sigma_to_hybrid",
    "BilinearRegridder",
    "ConservativeRegridder",
    "Grid",
    "NearestRegridder",
    "Regridder",
    "nearest_neighbor_indices",
    "regrid_dataset",
    "rcm2points",
    "rcm2rgrid",
    "rgrid2rcm",
    "grid_to_triple",
    "triple_to_grid",
]
