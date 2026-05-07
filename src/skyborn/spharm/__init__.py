"""
spharm - Spherical harmonic transforms for atmospheric/oceanic modeling.

This package exposes the main rectangular-grid `Spharmt` interface together
with the experimental packed reduced-Gaussian backend.
"""

from .reduced_gaussian import (
    ReducedGaussianBackendError,
    ReducedGaussianGrid,
    ReducedGaussianSpharmt,
    ReducedGaussianValidationError,
    regriduv,
)
from .spherical_harmonics import *

__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"
