"""Test the `fill_cube` function from the `gridfill` package."""

import warnings

import pytest

try:
    import iris.coords
    import iris.cube
except ImportError:
    pytest.skip(
        "Cannot import iris, fill_cube() will not be tested", allow_module_level=True
    )
import numpy as np

from skyborn.gridfill import fill_cube
from .test_fill import reference_solution


def make_cube(grid_data, cyclic):
    shape = grid_data.shape
    z_coord = iris.coords.DimCoord(np.arange(shape[0]), "altitude")
    x_coord = iris.coords.DimCoord(np.arange(shape[2]), "longitude", circular=cyclic)
    y_coord = iris.coords.DimCoord(np.arange(shape[1]), "latitude")
    dim_coords_and_dims = [(z_coord, 0), (y_coord, 1), (x_coord, 2)]
    cube = iris.cube.Cube(grid_data, dim_coords_and_dims=dim_coords_and_dims)
    return cube


class CubeFillTest(object):
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
        grid, cls.soln = reference_solution(cls.cyclic, cls.initzonal)
        cls.cube = make_cube(grid, cls.cyclic)

    def test_single_grid(self):
        filled = fill_cube(
            self.cube[0],
            self.eps,
            relax=self.relax,
            itermax=self.itermax,
            initzonal=self.initzonal,
            verbose=False,
        )
        self.assert_array_almost_equal(filled.data, self.soln[0])

    def test_multi_grid(self):
        filled = fill_cube(
            self.cube,
            self.eps,
            relax=self.relax,
            itermax=self.itermax,
            initzonal=self.initzonal,
            verbose=False,
        )
        self.assert_array_almost_equal(filled.data, self.soln)

    def test_multi_grid_inplace(self):
        cube = self.cube.copy()
        filled = fill_cube(
            cube,
            self.eps,
            relax=self.relax,
            itermax=self.itermax,
            initzonal=self.initzonal,
            verbose=False,
            inplace=True,
        )
        self.assert_array_almost_equal(filled.data, self.soln)
        self.assert_array_almost_equal(filled.data, cube.data)
        assert not hasattr(filled.data, "mask")
        assert not hasattr(cube.data, "mask")

    def test_not_masked(self):
        cube = self.cube.copy(data=self.cube.data.filled(fill_value=np.nan))
        with pytest.raises(TypeError):
            fill_cube(
                cube,
                self.eps,
                relax=self.relax,
                itermax=self.itermax,
                initzonal=self.initzonal,
                verbose=False,
            )

    def test_not_converged_warning(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # always trigger all warnings
            fill_cube(
                self.cube,
                self.eps / 1e7,
                relax=self.relax,
                itermax=self.itermax,
                initzonal=self.initzonal,
                verbose=False,
            )
            assert str(w[-1].message) == (
                "gridfill did not converge on 3 out " "of 3 slices"
            )
            assert issubclass(w[-1].category, UserWarning)

    def test_not_converged_result(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # always trigger all warnings
            filled, notconverged = fill_cube(
                self.cube,
                self.eps / 1e7,
                relax=self.relax,
                itermax=self.itermax,
                initzonal=self.initzonal,
                full_output=True,
                verbose=False,
            )
        self.assert_array_almost_equal(notconverged, np.ones(3, dtype=bool))

    def assert_array_almost_equal(self, a, b):
        np.testing.assert_array_almost_equal(a, b)


class TestFillCubeNonCyclicInitZonal(CubeFillTest):
    """Non-cyclic, initialized with zonal mean."""

    cyclic = False
    initzonal = True


class TestFillCubeCyclicInitZonal(CubeFillTest):
    """Cyclic, initialized with zonal mean."""

    cyclic = True
    initzonal = True
