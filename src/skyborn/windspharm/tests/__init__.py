"""Tests for the `windspharm` package."""

from __future__ import absolute_import

import pytest

from skyborn import windspharm
from .utils import error


# Create a mapping from interface name to VectorWind class.
solvers = {"standard": windspharm.standard.VectorWind}
try:
    solvers["iris"] = windspharm.iris.VectorWind
except AttributeError:
    pass
try:
    solvers["xarray"] = windspharm.xarray.VectorWind
except AttributeError:
    pass


class VectorWindTest(object):
    """Base class for vector wind tests."""

    def assert_error_is_zero(self, f1, f2):
        assert error(f1, f2) == pytest.approx(0.0, abs=1e-5)
