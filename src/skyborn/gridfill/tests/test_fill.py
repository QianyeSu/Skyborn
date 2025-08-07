"""Test the `gridfill` package."""

# Copyright (c) 2013 Andrew Dawson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
import os

import pytest
import numpy as np
import numpy.ma as ma

from gridfill import fill


def _test_data_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))


def _test_data_file(filename):
    return os.path.join(_test_data_dir(), filename)


def reference_solution(cyclic, initzonal):
    grid_data = ma.masked_invalid(np.load(_test_data_file("grid.npy")))
    cyclic_str = "_cyclic" if cyclic else ""
    zonal_str = "_zonal" if initzonal else ""
    soln_data = np.load(
        _test_data_file("soln{!s}{!s}.npy".format(cyclic_str, zonal_str))
    )
    return grid_data, soln_data


class GridFillTest(object):
    """Base class for `gridfill` tests."""

    # parameters to be set by each test subclass
    cyclic = None
    initzonal = None

    # parameters for relaxation scheme, same for all tests
    eps = 1e-4
    relax = 0.6
    itermax = 2000

    @classmethod
    def setup_class(cls):
        cls.grid, cls.soln = reference_solution(cls.cyclic, cls.initzonal)

    def test_single_grid(self):
        filled, c = fill(
            self.grid[0],
            1,
            0,
            self.eps,
            relax=self.relax,
            itermax=self.itermax,
            initzonal=self.initzonal,
            cyclic=self.cyclic,
            verbose=False,
        )
        self.assert_array_almost_equal(filled, self.soln[0])

    def test_multi_grid(self):
        filled, c = fill(
            self.grid,
            2,
            1,
            self.eps,
            relax=self.relax,
            itermax=self.itermax,
            initzonal=self.initzonal,
            cyclic=self.cyclic,
            verbose=False,
        )
        self.assert_array_almost_equal(filled, self.soln)

    def test_not_masked(self):
        with pytest.raises(TypeError):
            fill(
                self.grid.filled(fill_value=np.nan),
                self.eps,
                relax=self.relax,
                itermax=self.itermax,
                initzonal=self.initzonal,
                cyclic=self.cyclic,
                verbose=False,
            )

    def assert_array_almost_equal(self, a, b):
        np.testing.assert_array_almost_equal(a, b)


class TestFloat32(GridFillTest):
    """Test with 32-bit float input."""

    cyclic = False
    initzonal = False

    @classmethod
    def setup_class(cls):
        cls.grid, cls.soln = reference_solution(cls.cyclic, cls.initzonal)
        cls.grid = cls.grid.astype(np.float32)


class TestFillNonCyclicInitZero(GridFillTest):
    """Non-cyclic, initialized with zeros."""

    cyclic = False
    initzonal = False


class TestFillNonCyclicInitZonal(GridFillTest):
    """Non-cyclic, initialized with zonal mean."""

    cyclic = False
    initzonal = True


class TestFillCyclicInitZero(GridFillTest):
    """Cyclic, initialized with zeros."""

    cyclic = True
    initzonal = False


class TestFillCyclicInitZonal(GridFillTest):
    """Cyclic, initialized with zonal mean."""

    cyclic = True
    initzonal = True
