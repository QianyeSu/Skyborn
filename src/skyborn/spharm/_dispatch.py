"""Top-level workflow dispatch for spharm backends."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from . import spherical_harmonics as _rectangular_backend
from .reduced_gaussian import ReducedGaussianSpharmt
from .reduced_gaussian import regrid as _reduced_regrid
from .reduced_gaussian import regriduv as _reduced_regriduv


def regrid(grdin, grdout, datagrid, ntrunc=None, smooth=None):
    """Dispatch scalar spectral regridding by backend type."""
    if isinstance(grdin, ReducedGaussianSpharmt) or isinstance(
        grdout, ReducedGaussianSpharmt
    ):
        if not isinstance(grdin, ReducedGaussianSpharmt) or not isinstance(
            grdout, ReducedGaussianSpharmt
        ):
            raise TypeError(
                "regrid needs both grdin and grdout to be ReducedGaussianSpharmt "
                "when using the reduced-grid backend"
            )
        return _reduced_regrid(
            grdin,
            grdout,
            datagrid,
            ntrunc=ntrunc,
            smooth=smooth,
        )
    return _rectangular_backend.regrid(
        grdin,
        grdout,
        datagrid,
        ntrunc=ntrunc,
        smooth=smooth,
    )


def regriduv(
    grdin, grdout, ugrid, vgrid, ntrunc: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """Dispatch vector-wind spectral regridding by backend type."""
    if isinstance(grdin, ReducedGaussianSpharmt) or isinstance(
        grdout, ReducedGaussianSpharmt
    ):
        if not isinstance(grdin, ReducedGaussianSpharmt) or not isinstance(
            grdout, ReducedGaussianSpharmt
        ):
            raise TypeError(
                "regriduv needs both grdin and grdout to be ReducedGaussianSpharmt "
                "when using the reduced-grid backend"
            )
        return _reduced_regriduv(
            grdin,
            grdout,
            ugrid,
            vgrid,
            ntrunc=ntrunc,
        )
    return _rectangular_backend.regriduv(
        grdin,
        grdout,
        ugrid,
        vgrid,
        ntrunc=ntrunc,
    )
