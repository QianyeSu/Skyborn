"""Utilities for constructing tests."""

from __future__ import absolute_import

import numpy as np

try:
    from iris.cube import Cube
except ImportError:
    pass
try:
    import xarray as xr
except ImportError:
    try:
        import xray as xr
    except ImportError:
        pass


def __tomasked(*args):
    """Convert supported data types to masked arrays.

    The conversion is safe, so anything not recognised is just returned.

    """

    def __asma(a):
        try:
            if isinstance(a, Cube):
                # Retrieve the data from the cube.
                a = a.data
        except NameError:
            pass
        try:
            if isinstance(a, xr.DataArray):
                a = a.values
        except NameError:
            pass
        return a

    return [__asma(a) for a in args]


def error(a, b):
    """Compute the error between two arrays.

    Computes RMSD normalized by the range of the second input.

    """
    a, b = __tomasked(a, b)

    return np.sqrt(((a - b) ** 2).mean()) / (np.max(b) - np.min(b))


if __name__ == "__main__":
    pass
