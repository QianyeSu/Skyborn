"""
Introduction
============

This module provides a python interface to the NCAR
U{SPHEREPACK<https://www2.cisl.ucar.edu/resources/legacy/spherepack>} library.
It is not a one-to-one wrapper for the SPHEREPACK routines, rather
it provides a simple interface oriented toward working with
atmospheric general circulation model (GCM) data.

Requirements
============
 - U{numpy<http://numeric.scipy.org>}, and a fortran compiler
 supported by numpy.f2py.

Installation
============
 - U{Download<http://code.google.com/p/pyspharm/downloads/list>} module source,
 untar.
 - run C{python setup.py install} (as root if necessary).
 The SPHEREPACK fortran source files will be downloaded automatically by the
 setup.py script, since the SPHEREPACK license prohibits redistribution.
 To specify the fortran compiler to use (e.g. g95) run
 C{python setup.py config_fc --fcompiler=g95 install}. C{f2py -c --help-fcompiler}
 will show you what fortran compilers are available.

Usage
=====

>>> import spharm
>>> x=spharm.Spharmt(144,72,rsphere=8e6,gridtype='gaussian',legfunc='computed')

creates a class instance for spherical harmonic calculations on a 144x72
gaussian grid on a sphere with radius 8000 km. The associated legendre
functions are recomputed on the fly (instead of pre-computed and stored).
Default values of rsphere, gridtype and legfunc are 6.3712e6, 'regular'
and 'stored'. Real-world examples are included in the source distribution.

Class methods
=============
 - grdtospec: grid to spectral transform (spherical harmonic analysis).
 - spectogrd: spectral to grid transform (spherical harmonic synthesis).
 - getuv:  compute u and v winds from spectral coefficients of vorticity
 and divergence.
 - getvrtdivspec: get spectral coefficients of vorticity and divergence
 from u and v winds.
 - getgrad: compute the vector gradient given spectral coefficients.
 - getpsichi: compute streamfunction and velocity potential from winds.
 - specsmooth:  isotropic spectral smoothing.

Functions
=========
 - regrid:  spectral re-gridding, with optional spectral smoothing and/or
 truncation.
 - gaussian_lats_wts: compute gaussian latitudes and weights.
 - getspecindx: compute indices of zonal wavenumber and degree
 for complex spherical harmonic coefficients.
 - legendre: compute associated legendre functions.
 - getgeodesicpts: computes the points on the surface of the sphere
 corresponding to a twenty-sided (icosahedral) geodesic.
 - specintrp: spectral interpolation to an arbitrary point on the sphere.

Conventions
===========

The gridded data is assumed to be oriented such that i=1 is the
Greenwich meridian and j=1 is the northernmost point. Grid indices
increase eastward and southward. If nlat is odd the equator is included.
If nlat is even the equator will lie half way between points nlat/2
and (nlat/2)+1. nlat must be at least 3. For regular grids
(gridtype='regular') the poles will be included when nlat is odd.
The grid increment in longitude is 2*pi/nlon radians. For example,
nlon = 72 for a five degree grid. nlon must be greater than or
equal to 4. The efficiency of the computation is improved when nlon
is a product of small prime numbers.

The spectral data is assumed to be in a complex array of dimension
(ntrunc+1)*(ntrunc+2)/2. ntrunc is the triangular truncation limit
(ntrunc = 42 for T42). ntrunc must be <= nlat-1. Coefficients are
ordered so that first (nm=0) is m=0,n=0, second is m=0,n=1,
nm=ntrunc is m=0,n=ntrunc, nm=ntrunc+1 is m=1,n=1, etc.
The values of m (degree) and n (order) as a function of the index
nm are given by the arrays indxm, indxn returned by getspecindx.

The associated legendre polynomials are normalized so that the
integral (pbar(n,m,theta)**2)*sin(theta) on the interval theta=0 to pi
is 1, where pbar(m,n,theta)=sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m)))*
sin(theta)**m/(2**n*factorial(n)) times the (n+m)th derivative of
(x**2-1)**n with respect to x=cos(theta).
theta = pi/2 - phi, where phi is latitude and theta is colatitude.
Therefore, cos(theta) = sin(phi) and sin(theta) = cos(phi).
Note that pbar(0,0,theta)=sqrt(2)/2, and pbar(1,0,theta)=.5*sqrt(6)*sin(lat).

The default grid type is regular (equally spaced latitude points).
Set gridtype='gaussian' when creating a class instance
for gaussian latitude points.

Quantities needed to compute spherical harmonics are precomputed and stored
when the class instance is created with legfunc='stored' (the default).
If legfunc='computed', they are recomputed on the fly on each method call.
The storage requirements for legfunc="stored" increase like nlat**2, while
those for legfunc='stored' increase like nlat**3.  However, for
repeated method invocations on a single class instance, legfunc="stored"
will always be faster.

@contact: U{Jeff Whitaker<mailto:jeffrey.s.whitaker@noaa.gov>}

@version: 1.0.7

@license: Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that copyright notice and this permission notice appear in
supporting documentation.
THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.
"""

import math
from typing import Any, Callable, Dict, Literal, Optional, Tuple, Union

__all__ = [
    "Spharmt",
    "SpheremackError",
    "ValidationError",
    "regrid",
    "gaussian_lats_wts",
    "getspecindx",
    "getgeodesicpts",
    "legendre",
    "specintrp",
]
import sys

import numpy as np
from numpy.typing import NDArray

# Try relative import first (when imported as part of skyborn.spharm)
from . import _spherepack

# Type aliases for better readability
GridType = Literal["regular", "gaussian"]
LegendreFunc = Literal["stored", "computed"]
FloatArray = NDArray[np.floating]
ComplexArray = NDArray[np.complexfloating]
IntArray = NDArray[np.integer]

# define a list of instance variables that cannot be rebound or unbound.
_private_vars = ["nlon", "nlat", "gridtype", "legfunc", "rsphere", "precision"]

# Constants for better code clarity
DEFAULT_EARTH_RADIUS = 6.3712e6
MIN_NLON = 4
MIN_NLAT = 3
DEGREES_TO_RADIANS = math.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / math.pi


class SpheremackError(Exception):
    """Custom exception for SPHEREPACK related errors."""

    pass


class ValidationError(ValueError):
    """Custom exception for input validation errors."""

    pass


class Spharmt:
    """
    Spherical harmonic transform class with optimized performance and type safety.

    Attributes:
        nlat: Number of latitudes (read-only after initialization)
        nlon: Number of longitudes (read-only after initialization)
        rsphere: Sphere radius in meters (read-only after initialization)
        gridtype: Grid type - 'regular' or 'gaussian' (read-only after initialization)
        legfunc: Legendre function handling - 'stored' or 'computed' (read-only after initialization)
    """

    def __init__(
        self,
        nlon: int,
        nlat: int,
        rsphere: float = DEFAULT_EARTH_RADIUS,
        gridtype: GridType = "regular",
        legfunc: LegendreFunc = "stored",
        precision: Literal["auto", "single", "double"] = "auto",
    ) -> None:
        """
        Create a Spharmt class instance with input validation and optimization.

        Args:
            nlon: Number of longitudes (must be >= 4)
            nlat: Number of latitudes (must be >= 3)
            rsphere: Sphere radius in meters (default: Earth radius)
            gridtype: Grid type - 'regular' or 'gaussian'
            legfunc: Legendre function handling - 'stored' or 'computed'

        Raises:
            ValidationError: If input parameters are invalid
            SpheremackError: If SPHEREPACK initialization fails
        """
        # Validate inputs with better error messages
        self._validate_inputs(nlon, nlat, rsphere, gridtype, legfunc)

        # Set read-only attributes
        self.nlon = nlon
        self.nlat = nlat
        self.rsphere = rsphere
        self.gridtype = gridtype
        self.legfunc = legfunc
        self.precision = self._validate_precision(precision)

        # Pre-compute common values to avoid repeated calculations
        self._n1, self._n2 = self._compute_transform_sizes()

        # Initialize SPHEREPACK work arrays using strategy pattern
        self._initialize_spherepack_arrays()

    def _validate_inputs(
        self, nlon: int, nlat: int, rsphere: float, gridtype: str, legfunc: str
    ) -> None:
        """Validate all input parameters with descriptive error messages."""
        validators = [
            (rsphere > 0.0, f"rsphere must be positive, got {rsphere}"),
            (nlon >= MIN_NLON, f"nlon must be >= {MIN_NLON}, got {nlon}"),
            (nlat >= MIN_NLAT, f"nlat must be >= {MIN_NLAT}, got {nlat}"),
            (
                gridtype in ["regular", "gaussian"],
                f'gridtype must be "regular" or "gaussian", got "{gridtype}"',
            ),
            (
                legfunc in ["stored", "computed"],
                f'legfunc must be "stored" or "computed", got "{legfunc}"',
            ),
        ]

        for condition, error_msg in validators:
            if not condition:
                raise ValidationError(error_msg)

    @staticmethod
    def _validate_precision(precision: str) -> str:
        """Validate the requested public-output precision mode."""
        if precision not in ("auto", "single", "double"):
            raise ValidationError(
                f'precision must be "auto", "single", or "double", got "{precision}"'
            )
        return precision

    @staticmethod
    def _input_public_real_dtype(array: np.ndarray) -> np.dtype:
        """Map an input array dtype to the matching public real dtype."""
        dtype = np.dtype(np.asarray(array).dtype)
        if dtype == np.float32 or dtype == np.complex64:
            return np.dtype(np.float32)
        return np.dtype(np.float64)

    def _public_real_dtype(self, *arrays: np.ndarray) -> np.dtype:
        """Resolve the public real output dtype for this instance."""
        if self.precision == "single":
            return np.dtype(np.float32)
        if self.precision == "double":
            return np.dtype(np.float64)

        dtype = np.dtype(np.float32)
        for array in arrays:
            dtype = np.result_type(dtype, self._input_public_real_dtype(array))
        return np.dtype(np.float64 if dtype == np.float64 else np.float32)

    def _public_complex_dtype(self, *arrays: np.ndarray) -> np.dtype:
        """Resolve the public complex output dtype for this instance."""
        return (
            np.dtype(np.complex128)
            if self._public_real_dtype(*arrays) == np.dtype(np.float64)
            else np.dtype(np.complex64)
        )

    def _compute_transform_sizes(self) -> Tuple[int, int]:
        """Pre-compute transform sizes to avoid repeated calculations."""
        nlat, nlon = self.nlat, self.nlon

        n1 = min(nlat, (nlon + 1) // 2 if nlon % 2 else (nlon + 2) // 2)
        n2 = (nlat + 1) // 2 if nlat % 2 else nlat // 2

        return n1, n2

    def _initialize_spherepack_arrays(self) -> None:
        """Initialize SPHEREPACK work arrays using strategy pattern to reduce code duplication."""
        initializers = {
            ("regular", "stored"): self._init_regular_stored,
            ("regular", "computed"): self._init_regular_computed,
            ("gaussian", "stored"): self._init_gaussian_stored,
            ("gaussian", "computed"): self._init_gaussian_computed,
        }

        key = (self.gridtype, self.legfunc)
        initializer = initializers.get(key)

        if initializer is None:
            raise ValidationError(
                f"Unsupported combination: gridtype='{self.gridtype}', legfunc='{self.legfunc}'"
            )

        initializer()

    def _call_spherepack_safely(
        self, func: Callable, *args, operation_name: str, **kwargs
    ) -> Any:
        """Safely call SPHEREPACK functions with consistent error handling."""
        try:
            result = func(*args, **kwargs)
            if isinstance(result, tuple) and len(result) > 1:
                *data, ierror = result
                if ierror != 0:
                    raise SpheremackError(
                        f"{operation_name} failed with error code {ierror}"
                    )
                return data[0] if len(data) == 1 else data
            return result
        except Exception as e:
            if isinstance(e, SpheremackError):
                raise
            raise SpheremackError(f"{operation_name} failed: {str(e)}") from e

    def _init_regular_stored(self) -> None:
        """Initialize arrays for regular grid with stored Legendre functions."""
        nlat, nlon = self.nlat, self.nlon
        n1, n2 = self._n1, self._n2

        # Scalar harmonic analysis initialization
        lshaes = (n1 * n2 * (nlat + nlat - n1 + 1)) // 2 + nlon + 15
        lwork = 5 * nlat * n2 + 3 * ((n1 - 2) * (nlat + nlat - n1 - 1)) // 2

        self.wshaes = self._call_spherepack_safely(
            _spherepack.shaesi,
            nlat,
            nlon,
            lshaes,
            lwork,
            nlat + 1,
            operation_name="shaesi initialization",
        )

        # Scalar harmonic synthesis initialization
        self.wshses = self._call_spherepack_safely(
            _spherepack.shsesi,
            nlat,
            nlon,
            lshaes,
            lwork,
            nlat + 1,
            operation_name="shsesi initialization",
        )

        # Vector harmonic analysis initialization
        lvhaes = n1 * n2 * (nlat + nlat - n1 + 1) + nlon + 15
        lwork_vec = 3 * (max(n1 - 2, 0) * (nlat + nlat - n1 - 1)) // 2 + 5 * n2 * nlat

        self.wvhaes = self._call_spherepack_safely(
            _spherepack.vhaesi,
            nlat,
            nlon,
            lvhaes,
            lwork_vec,
            2 * (nlat + 1),
            operation_name="vhaesi initialization",
        )

        # Vector harmonic synthesis initialization
        self.wvhses = self._call_spherepack_safely(
            _spherepack.vhsesi,
            nlat,
            nlon,
            lvhaes,
            lwork_vec,
            2 * (nlat + 1),
            operation_name="vhsesi initialization",
        )

    def _init_regular_computed(self) -> None:
        """Initialize arrays for regular grid with computed Legendre functions."""
        nlat, nlon = self.nlat, self.nlon
        n1, n2 = self._n1, self._n2

        # Scalar harmonic analysis initialization
        lshaec = (
            2 * nlat * n2 + 3 * ((n1 - 2) * (nlat + nlat - n1 - 1)) // 2 + nlon + 15
        )

        self.wshaec = self._call_spherepack_safely(
            _spherepack.shaeci,
            nlat,
            nlon,
            lshaec,
            2 * (nlat + 1),
            operation_name="shaeci initialization",
        )

        # Scalar harmonic synthesis initialization
        self.wshsec = self._call_spherepack_safely(
            _spherepack.shseci,
            nlat,
            nlon,
            lshaec,
            2 * (nlat + 1),
            operation_name="shseci initialization",
        )

        # Vector harmonic analysis initialization
        lvhaec = 4 * nlat * n2 + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1) + nlon + 15

        self.wvhaec = self._call_spherepack_safely(
            _spherepack.vhaeci,
            nlat,
            nlon,
            lvhaec,
            2 * (nlat + 1),
            operation_name="vhaeci initialization",
        )

        # Vector harmonic synthesis initialization
        self.wvhsec = self._call_spherepack_safely(
            _spherepack.vhseci,
            nlat,
            nlon,
            lvhaec,
            2 * (nlat + 1),
            operation_name="vhseci initialization",
        )

    def _init_gaussian_stored(self) -> None:
        """Initialize arrays for Gaussian grid with stored Legendre functions."""
        nlat, nlon = self.nlat, self.nlon
        n1, n2 = self._n1, self._n2

        # Scalar harmonic analysis initialization
        lshags = (
            nlat * (3 * (n1 + n2) - 2)
            + (n1 - 1) * (n2 * (2 * nlat - n1) - 3 * n1) // 2
            + nlon
            + 15
        )
        lwork = 4 * nlat * (nlat + 2) + 2
        ldwork = nlat * (nlat + 4)

        self.wshags = self._call_spherepack_safely(
            _spherepack.shagsi,
            nlat,
            nlon,
            lshags,
            lwork,
            ldwork,
            operation_name="shagsi initialization",
        )

        # Scalar harmonic synthesis initialization
        self.wshsgs = self._call_spherepack_safely(
            _spherepack.shsgsi,
            nlat,
            nlon,
            lshags,
            lwork,
            ldwork,
            operation_name="shsgsi initialization",
        )

        # Vector harmonic analysis initialization
        lvhags = (nlat + 1) * (nlat + 1) * nlat // 2 + nlon + 15
        ldwork_vec = (3 * nlat * (nlat + 3) + 2) // 2

        self.wvhags = self._call_spherepack_safely(
            _spherepack.vhagsi,
            nlat,
            nlon,
            lvhags,
            ldwork_vec,
            operation_name="vhagsi initialization",
        )

        # Vector harmonic synthesis initialization
        lvhsgs = n1 * n2 * (nlat + nlat - n1 + 1) + nlon + 15 + 2 * nlat

        self.wvhsgs = self._call_spherepack_safely(
            _spherepack.vhsgsi,
            nlat,
            nlon,
            lvhsgs,
            ldwork_vec,
            operation_name="vhsgsi initialization",
        )

    def _init_gaussian_computed(self) -> None:
        """Initialize arrays for Gaussian grid with computed Legendre functions."""
        nlat, nlon = self.nlat, self.nlon
        n1, n2 = self._n1, self._n2

        # Scalar harmonic analysis initialization
        lshagc = nlat * (2 * n2 + 3 * n1 - 2) + 3 * n1 * (1 - n1) // 2 + nlon + 15
        ldwork = nlat * (nlat + 4)

        self.wshagc = self._call_spherepack_safely(
            _spherepack.shagci,
            nlat,
            nlon,
            lshagc,
            ldwork,
            operation_name="shagci initialization",
        )

        # Scalar harmonic synthesis initialization
        self.wshsgc = self._call_spherepack_safely(
            _spherepack.shsgci,
            nlat,
            nlon,
            lshagc,
            ldwork,
            operation_name="shsgci initialization",
        )

        # Vector harmonic analysis initialization
        lvhagc = (
            4 * nlat * n2 + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1) + nlon + n2 + 15
        )
        ldwork_vec = 2 * nlat * (nlat + 1) + 1

        self.wvhagc = self._call_spherepack_safely(
            _spherepack.vhagci,
            nlat,
            nlon,
            lvhagc,
            ldwork_vec,
            operation_name="vhagci initialization",
        )

        # Vector harmonic synthesis initialization
        lvhsgc = 4 * nlat * n2 + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1) + nlon + 15

        self.wvhsgc = self._call_spherepack_safely(
            _spherepack.vhsgci,
            nlat,
            nlon,
            lvhsgc,
            ldwork_vec,
            operation_name="vhsgci initialization",
        )

    def __setattr__(self, key: str, val: Any) -> None:
        """Prevent modification of read-only instance variables."""
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError(f"Attempt to rebind read-only instance variable {key}")
        else:
            self.__dict__[key] = val

    def __delattr__(self, key: str) -> None:
        """Prevent deletion of read-only instance variables."""
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError(f"Attempt to unbind read-only instance variable {key}")
        else:
            del self.__dict__[key]

    def __repr__(self) -> str:
        """Return string representation of Spharmt instance."""
        return f"Spharmt({self.nlon:d}, {self.nlat:d}, {self.rsphere:e}, '{self.gridtype}', '{self.legfunc}')"

    def _validate_grid_data(
        self, data: FloatArray, operation_name: str
    ) -> Tuple[int, FloatArray, Tuple[int, ...]]:
        """Validate grid data and flatten non-horizontal dimensions."""
        data = np.asarray(data)
        if data.ndim < 2:
            raise ValidationError(
                f"{operation_name} needs at least a rank 2 array, got rank {data.ndim}"
            )

        if data.shape[0] != self.nlat or data.shape[1] != self.nlon:
            raise ValidationError(
                f"{operation_name} needs an array of size {self.nlat} by {self.nlon}, "
                f"got {data.shape[0]} by {data.shape[1]}"
            )

        if data.ndim == 2:
            nt = 1
            normalized_data = np.expand_dims(data, 2)
            extra_shape: Tuple[int, ...] = ()
        else:
            extra_shape = tuple(data.shape[2:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized_data = data.reshape(self.nlat, self.nlon, nt)

        return nt, normalized_data, extra_shape

    def _validate_spectral_data(
        self, data: ComplexArray, operation_name: str
    ) -> Tuple[int, int, ComplexArray, Tuple[int, ...]]:
        """Validate spectral data and flatten non-coefficient dimensions."""
        data = np.asarray(data)
        if data.ndim < 1:
            raise ValidationError(
                f"{operation_name} needs at least a rank 1 array, got rank {data.ndim}"
            )

        # Infer ntrunc from spectral data size
        ntrunc = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - data.shape[0])))

        if ntrunc > self.nlat - 1:
            raise ValidationError(
                f"ntrunc too large - can be max of {self.nlat - 1}, got {ntrunc}"
            )

        if data.ndim == 1:
            nt = 1
            normalized_data = np.expand_dims(data, 1)
            extra_shape: Tuple[int, ...] = ()
        else:
            extra_shape = tuple(data.shape[1:])
            nt = int(np.prod(extra_shape, dtype=int))
            normalized_data = data.reshape(data.shape[0], nt)

        return nt, ntrunc, normalized_data, extra_shape

    def _restore_grid_shape(
        self,
        datagrid: FloatArray,
        extra_shape: Tuple[int, ...],
        output_dtype: Optional[np.dtype] = None,
    ) -> FloatArray:
        """Restore flattened grid fields to their original extra dimensions."""
        array = np.asarray(datagrid)
        if output_dtype is None:
            if self.precision == "double":
                output_dtype = np.float64
            elif self.precision == "single":
                output_dtype = np.float32
        if output_dtype is not None:
            array = array.astype(np.dtype(output_dtype), copy=False)
        if not extra_shape:
            return array[:, :, 0]
        return array.reshape((self.nlat, self.nlon) + extra_shape)

    def _restore_spectral_shape(
        self,
        dataspec: ComplexArray,
        extra_shape: Tuple[int, ...],
        output_dtype: Optional[np.dtype] = None,
    ) -> ComplexArray:
        """Restore flattened spectral fields to their original extra dimensions."""
        array = np.asarray(dataspec)
        if output_dtype is None:
            if self.precision == "double":
                output_dtype = np.complex128
            elif self.precision == "single":
                output_dtype = np.complex64
        if output_dtype is not None:
            array = array.astype(np.dtype(output_dtype), copy=False)
        if not extra_shape:
            return array[:, 0]
        return array.reshape((array.shape[0],) + extra_shape)

    def _spectogrd_pair(
        self, spec_a: ComplexArray, spec_b: ComplexArray
    ) -> Tuple[FloatArray, FloatArray]:
        """Synthesize two scalar spectra together through one batched call."""
        if spec_a.shape != spec_b.shape:
            raise ValidationError("paired spectra must have the same shape")

        paired_specs = np.stack((spec_a, spec_b), axis=-1)
        paired_grids = self.spectogrd(paired_specs)
        return paired_grids[..., 0], paired_grids[..., 1]

    def _invlapspec(
        self, dataspec: ComplexArray, operation_name: str = "_invlapspec"
    ) -> ComplexArray:
        """Apply the inverse Laplacian to a scalar spectrum and preserve shape."""
        _, _, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, operation_name
        )
        invlap_spec = _spherepack.invlap(normalized_spec, self.rsphere)
        return self._restore_spectral_shape(invlap_spec, extra_shape)

    def _analyze_vector_harmonics(
        self,
        ugrid: FloatArray,
        vgrid: FloatArray,
        ntrunc: Optional[int],
        operation_name: str,
        ityp: int = 0,
    ) -> Tuple[
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
        Tuple[int, ...],
        int,
    ]:
        """Run one vector harmonic analysis and return raw coefficient arrays."""
        if ugrid.shape != vgrid.shape:
            raise ValidationError("ugrid and vgrid must have the same shape")

        nt, normalized_u, extra_shape = self._validate_grid_data(ugrid, operation_name)
        _, normalized_v, v_extra_shape = self._validate_grid_data(vgrid, operation_name)
        if extra_shape != v_extra_shape:
            raise ValidationError("ugrid and vgrid must have consistent dimensions")
        return self._analyze_vector_harmonics_flattened(
            normalized_u,
            normalized_v,
            extra_shape,
            ntrunc,
            operation_name,
            ityp=ityp,
        )

    def _analyze_vector_harmonics_flattened(
        self,
        normalized_u: FloatArray,
        normalized_v: FloatArray,
        extra_shape: Tuple[int, ...],
        ntrunc: Optional[int],
        operation_name: str,
        ityp: int = 0,
    ) -> Tuple[
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
        Tuple[int, ...],
        int,
    ]:
        """Run one vector harmonic analysis from wrapper-owned flattened 3D arrays."""
        if normalized_u.shape != normalized_v.shape:
            raise ValidationError("ugrid and vgrid must have consistent dimensions")
        if normalized_u.ndim != 3:
            raise ValidationError(
                f"{operation_name} flattened inputs must be rank 3, got rank {normalized_u.ndim}"
            )
        if normalized_u.shape[:2] != (self.nlat, self.nlon):
            raise ValidationError(
                f"{operation_name} flattened inputs must have shape ({self.nlat}, {self.nlon}, nt)"
            )

        nt = normalized_u.shape[2]
        ntrunc = self._validate_ntrunc(ntrunc, self.nlat - 1)
        if ityp not in (0, 1, 2):
            raise ValidationError(
                f"{operation_name} invalid ityp {ityp}; expected 0, 1, or 2"
            )

        # Convert from geographical to mathematical coordinates.
        w = normalized_u
        v = -normalized_v

        vha_functions = {
            ("regular", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaes,
                    v,
                    w,
                    self.wvhaes,
                    (2 * nt + 1) * self.nlat * self.nlon,
                    ityp=ityp,
                    operation_name="vhaes transform",
                )
            ),
            ("regular", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaec,
                    v,
                    w,
                    self.wvhaec,
                    self.nlat * (2 * nt * self.nlon + max(6 * self._n2, self.nlon)),
                    ityp=ityp,
                    operation_name="vhaec transform",
                )
            ),
            ("gaussian", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhags,
                    v,
                    w,
                    self.wvhags,
                    (2 * nt + 1) * self.nlat * self.nlon,
                    ityp=ityp,
                    operation_name="vhags transform",
                )
            ),
            ("gaussian", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhagc,
                    v,
                    w,
                    self.wvhagc,
                    2 * self.nlat * (2 * self.nlon * nt + 3 * self._n2),
                    ityp=ityp,
                    operation_name="vhagc transform",
                )
            ),
        }

        br, bi, cr, ci = vha_functions[(self.gridtype, self.legfunc)]()
        return br, bi, cr, ci, extra_shape, ntrunc

    def _analyze_vector_harmonic_component(
        self,
        ugrid: FloatArray,
        vgrid: FloatArray,
        ntrunc: Optional[int],
        operation_name: str,
        component: str,
    ) -> Tuple[FloatArray, FloatArray, Tuple[int, ...], int]:
        """Run one vector harmonic analysis and return only one coefficient family."""
        if ugrid.shape != vgrid.shape:
            raise ValidationError("ugrid and vgrid must have the same shape")

        nt, normalized_u, extra_shape = self._validate_grid_data(ugrid, operation_name)
        _, normalized_v, v_extra_shape = self._validate_grid_data(vgrid, operation_name)
        if extra_shape != v_extra_shape:
            raise ValidationError("ugrid and vgrid must have consistent dimensions")
        return self._analyze_vector_harmonic_component_flattened(
            normalized_u,
            normalized_v,
            extra_shape,
            ntrunc,
            operation_name,
            component,
        )

    def _analyze_vector_harmonic_component_flattened(
        self,
        normalized_u: FloatArray,
        normalized_v: FloatArray,
        extra_shape: Tuple[int, ...],
        ntrunc: Optional[int],
        operation_name: str,
        component: str,
    ) -> Tuple[FloatArray, FloatArray, Tuple[int, ...], int]:
        """Run one single-family vector harmonic analysis from flattened 3D arrays."""
        if normalized_u.shape != normalized_v.shape:
            raise ValidationError("ugrid and vgrid must have consistent dimensions")
        if normalized_u.ndim != 3:
            raise ValidationError(
                f"{operation_name} flattened inputs must be rank 3, got rank {normalized_u.ndim}"
            )
        if normalized_u.shape[:2] != (self.nlat, self.nlon):
            raise ValidationError(
                f"{operation_name} flattened inputs must have shape ({self.nlat}, {self.nlon}, nt)"
            )

        nt = normalized_u.shape[2]
        ntrunc = self._validate_ntrunc(ntrunc, self.nlat - 1)
        if component not in {"div", "vrt"}:
            raise ValidationError(
                f"{operation_name} invalid component {component!r}; expected 'div' or 'vrt'"
            )

        # Convert from geographical to mathematical coordinates.
        w = normalized_u
        v = -normalized_v

        work_regular = (2 * nt + 1) * self.nlat * self.nlon
        work_regular_computed = self.nlat * (
            2 * nt * self.nlon + max(6 * self._n2, self.nlon)
        )
        work_gaussian_computed = 2 * self.nlat * (2 * self.nlon * nt + 3 * self._n2)

        vha_functions = {
            ("regular", "stored", "div"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaesdiv,
                    v,
                    w,
                    self.wvhaes,
                    operation_name="vhaes divergence transform",
                )
            ),
            ("regular", "stored", "vrt"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaesvrt,
                    v,
                    w,
                    self.wvhaes,
                    work_regular,
                    operation_name="vhaes vorticity transform",
                )
            ),
            ("regular", "computed", "div"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaecdiv,
                    v,
                    w,
                    self.wvhaec,
                    work_regular_computed,
                    operation_name="vhaec divergence transform",
                )
            ),
            ("regular", "computed", "vrt"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhaecvrt,
                    v,
                    w,
                    self.wvhaec,
                    work_regular_computed,
                    operation_name="vhaec vorticity transform",
                )
            ),
            ("gaussian", "stored", "div"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhagsdiv,
                    v,
                    w,
                    self.wvhags,
                    work_regular,
                    operation_name="vhags divergence transform",
                )
            ),
            ("gaussian", "stored", "vrt"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhagsvrt,
                    v,
                    w,
                    self.wvhags,
                    work_regular,
                    operation_name="vhags vorticity transform",
                )
            ),
            ("gaussian", "computed", "div"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhagcdiv,
                    v,
                    w,
                    self.wvhagc,
                    work_gaussian_computed,
                    operation_name="vhagc divergence transform",
                )
            ),
            ("gaussian", "computed", "vrt"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhagcvrt,
                    v,
                    w,
                    self.wvhagc,
                    work_gaussian_computed,
                    operation_name="vhagc vorticity transform",
                )
            ),
        }

        coeff_r, coeff_i = vha_functions[(self.gridtype, self.legfunc, component)]()
        return coeff_r, coeff_i, extra_shape, ntrunc

    def _validate_ntrunc(self, ntrunc: Optional[int], max_allowed: int) -> int:
        """Validate and return appropriate ntrunc value."""
        if ntrunc is None:
            return self.nlat - 1

        if ntrunc < 0 or ntrunc > max_allowed:
            raise ValidationError(
                f"ntrunc must be between 0 and {max_allowed}, got {ntrunc}"
            )

        return ntrunc

    def grdtospec(
        self, datagrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """
        Grid to spectral transform (spherical harmonic analysis) with type safety.

        Args:
            datagrid: Grid data with shape (nlat, nlon) or (nlat, nlon, *extra_dims)
            ntrunc: Optional spectral truncation limit (default: nlat-1)

        Returns:
            Complex spherical harmonic coefficients

        Raises:
            ValidationError: If input data has invalid dimensions
            SpheremackError: If SPHEREPACK transform fails
        """
        # Validate inputs
        nt, normalized_data, extra_shape = self._validate_grid_data(
            datagrid, "grdtospec"
        )
        ntrunc = self._validate_ntrunc(ntrunc, self.nlat - 1)

        # Select transform function based on configuration
        transform_functions = {
            ("regular", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shaes,
                    normalized_data,
                    self.wshaes,
                    (nt + 1) * self.nlat * self.nlon,
                    operation_name="shaes transform",
                )
            ),
            ("regular", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shaec,
                    normalized_data,
                    self.wshaec,
                    self.nlat * (nt * self.nlon + max(3 * self._n2, self.nlon)),
                    operation_name="shaec transform",
                )
            ),
            ("gaussian", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shags,
                    normalized_data,
                    self.wshags,
                    self.nlat * self.nlon * (nt + 1),
                    operation_name="shags transform",
                )
            ),
            ("gaussian", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shagc,
                    normalized_data,
                    self.wshagc,
                    self.nlat * (self.nlon * nt + max(3 * self._n2, self.nlon)),
                    operation_name="shagc transform",
                )
            ),
        }

        transform_func = transform_functions[(self.gridtype, self.legfunc)]
        a, b = transform_func()

        # Convert 2D real and imaginary arrays to 1D complex array
        dataspec = _spherepack.twodtooned(a, b, ntrunc)

        return self._restore_spectral_shape(
            dataspec, extra_shape, self._public_complex_dtype(datagrid)
        )

    def spectogrd(self, dataspec: ComplexArray) -> FloatArray:
        """
        Spectral to grid transform (spherical harmonic synthesis) with type safety.

        Args:
            dataspec: Complex spectral coefficients

        Returns:
            Grid data with shape (nlat, nlon) or (nlat, nlon, *extra_dims)

        Raises:
            ValidationError: If input data has invalid dimensions
            SpheremackError: If SPHEREPACK transform fails
        """
        # Validate inputs and get normalized data
        nt, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, "spectogrd"
        )

        # Convert 1D complex array to 2D real and imaginary arrays
        a, b = _spherepack.onedtotwod(normalized_spec, self.nlat)

        # Select transform function based on configuration
        transform_functions = {
            ("regular", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shses,
                    self.nlon,
                    a,
                    b,
                    self.wshses,
                    (nt + 1) * self.nlat * self.nlon,
                    operation_name="shses transform",
                )
            ),
            ("regular", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shsec,
                    self.nlon,
                    a,
                    b,
                    self.wshsec,
                    self.nlat * (nt * self.nlon + max(3 * self._n2, self.nlon)),
                    operation_name="shsec transform",
                )
            ),
            ("gaussian", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shsgs,
                    self.nlon,
                    a,
                    b,
                    self.wshsgs,
                    self.nlat * self.nlon * (nt + 1),
                    operation_name="shsgs transform",
                )
            ),
            ("gaussian", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.shsgc,
                    self.nlon,
                    a,
                    b,
                    self.wshsgc,
                    self.nlat * (self.nlon * nt + max(3 * self._n2, self.nlon)),
                    operation_name="shsgc transform",
                )
            ),
        }

        transform_func = transform_functions[(self.gridtype, self.legfunc)]
        datagrid = transform_func()

        return self._restore_grid_shape(
            datagrid, extra_shape, self._public_real_dtype(dataspec)
        )

    def getvrtdivspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> Tuple[ComplexArray, ComplexArray]:
        """
        Compute spectral coefficients of vorticity and divergence from vector wind.

        Args:
            ugrid: Zonal wind grid data
            vgrid: Meridional wind grid data
            ntrunc: Optional spectral truncation limit

        Returns:
            Tuple of (vorticity_spec, divergence_spec) coefficients

        Raises:
            ValidationError: If input arrays have mismatched shapes or invalid dimensions
            SpheremackError: If SPHEREPACK transform fails
        """
        br, bi, cr, ci, extra_shape, ntrunc = self._analyze_vector_harmonics(
            ugrid, vgrid, ntrunc, "getvrtdivspec"
        )
        vrtspec, divspec = _spherepack.twodtooned_vrtdiv(
            br, bi, cr, ci, ntrunc, self.rsphere
        )

        return (
            self._restore_spectral_shape(
                vrtspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
            ),
            self._restore_spectral_shape(
                divspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
            ),
        )

    def getvrtspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute only vorticity spectral coefficients from vector wind."""
        cr, ci, extra_shape, ntrunc = self._analyze_vector_harmonic_component(
            ugrid, vgrid, ntrunc, "getvrtspec", component="vrt"
        )
        vrtspec = _spherepack.twodtooned_vrt(cr, ci, ntrunc, self.rsphere)
        return self._restore_spectral_shape(
            vrtspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
        )

    def getdivspec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute only divergence spectral coefficients from vector wind."""
        br, bi, extra_shape, ntrunc = self._analyze_vector_harmonic_component(
            ugrid, vgrid, ntrunc, "getdivspec", component="div"
        )
        divspec = _spherepack.twodtooned_div(br, bi, ntrunc, self.rsphere)
        return self._restore_spectral_shape(
            divspec, extra_shape, self._public_complex_dtype(ugrid, vgrid)
        )

    def _synthesize_vector_harmonics(
        self,
        br: FloatArray,
        bi: FloatArray,
        cr: FloatArray,
        ci: FloatArray,
        nt: int,
        extra_shape: Tuple[int, ...],
        operation_name: str,
        ityp: int = 0,
    ) -> Tuple[FloatArray, FloatArray]:
        """Run one vector harmonic synthesis and restore grid dimensions."""
        if ityp not in (0, 1, 2):
            raise ValidationError(
                f"{operation_name} invalid ityp {ityp}; expected 0, 1, or 2"
            )

        vhs_functions = {
            ("regular", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhses,
                    self.nlon,
                    br,
                    bi,
                    cr,
                    ci,
                    self.wvhses,
                    (2 * nt + 1) * self.nlat * self.nlon,
                    ityp=ityp,
                    operation_name="vhses transform",
                )
            ),
            ("regular", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhsec,
                    self.nlon,
                    br,
                    bi,
                    cr,
                    ci,
                    self.wvhsec,
                    self.nlat * (2 * nt * self.nlon + max(6 * self._n2, self.nlon)),
                    ityp=ityp,
                    operation_name="vhsec transform",
                )
            ),
            ("gaussian", "stored"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhsgs,
                    self.nlon,
                    br,
                    bi,
                    cr,
                    ci,
                    self.wvhsgs,
                    (2 * nt + 1) * self.nlat * self.nlon,
                    ityp=ityp,
                    operation_name="vhsgs transform",
                )
            ),
            ("gaussian", "computed"): (
                lambda: self._call_spherepack_safely(
                    _spherepack.vhsgc,
                    self.nlon,
                    br,
                    bi,
                    cr,
                    ci,
                    self.wvhsgc,
                    self.nlat * (2 * nt * self.nlon + max(6 * self._n2, self.nlon)),
                    ityp=ityp,
                    operation_name="vhsgc transform",
                )
            ),
        }

        vhs_func = vhs_functions[(self.gridtype, self.legfunc)]
        v, w = vhs_func()
        public_dtype = self._public_real_dtype(vrtspec, divspec)
        return self._restore_grid_shape(
            w, extra_shape, public_dtype
        ), -self._restore_grid_shape(v, extra_shape, public_dtype)

    def _getuv_spec_component(
        self, dataspec: ComplexArray, component: str, operation_name: str
    ) -> Tuple[FloatArray, FloatArray]:
        """Synthesize one vector family from a single spectral field."""
        nt, _, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, operation_name
        )

        work_regular = (2 * nt + 1) * self.nlat * self.nlon
        work_regular_computed = self.nlat * (
            2 * nt * self.nlon + max(6 * self._n2, self.nlon)
        )

        if component == "div":
            br, bi = _spherepack.onedtotwod_div(
                normalized_spec, self.nlat, self.rsphere
            )
            synthesis_functions = {
                ("regular", "stored"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsesdiv,
                        self.nlon,
                        br,
                        bi,
                        self.wvhses,
                        work_regular,
                        operation_name="vhses divergence transform",
                    )
                ),
                ("regular", "computed"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsecdiv,
                        self.nlon,
                        br,
                        bi,
                        self.wvhsec,
                        work_regular_computed,
                        operation_name="vhsec divergence transform",
                    )
                ),
                ("gaussian", "stored"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsgsdiv,
                        self.nlon,
                        br,
                        bi,
                        self.wvhsgs,
                        work_regular,
                        operation_name="vhsgs divergence transform",
                    )
                ),
                ("gaussian", "computed"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsgcdiv,
                        self.nlon,
                        br,
                        bi,
                        self.wvhsgc,
                        work_regular_computed,
                        operation_name="vhsgc divergence transform",
                    )
                ),
            }
            v, w = synthesis_functions[(self.gridtype, self.legfunc)]()
            public_dtype = self._public_real_dtype(normalized_spec)
            return self._restore_grid_shape(
                w, extra_shape, public_dtype
            ), -self._restore_grid_shape(v, extra_shape, public_dtype)

        if component == "vrt":
            cr, ci = _spherepack.onedtotwod_vrt(
                normalized_spec, self.nlat, self.rsphere
            )
            synthesis_functions = {
                ("regular", "stored"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsesvrt,
                        self.nlon,
                        cr,
                        ci,
                        self.wvhses,
                        work_regular,
                        operation_name="vhses vorticity transform",
                    )
                ),
                ("regular", "computed"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsecvrt,
                        self.nlon,
                        cr,
                        ci,
                        self.wvhsec,
                        work_regular_computed,
                        operation_name="vhsec vorticity transform",
                    )
                ),
                ("gaussian", "stored"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsgsvrt,
                        self.nlon,
                        cr,
                        ci,
                        self.wvhsgs,
                        work_regular,
                        operation_name="vhsgs vorticity transform",
                    )
                ),
                ("gaussian", "computed"): (
                    lambda: self._call_spherepack_safely(
                        _spherepack.vhsgcvrt,
                        self.nlon,
                        cr,
                        ci,
                        self.wvhsgc,
                        work_regular_computed,
                        operation_name="vhsgc vorticity transform",
                    )
                ),
            }
            v, w = synthesis_functions[(self.gridtype, self.legfunc)]()
            public_dtype = self._public_real_dtype(normalized_spec)
            return self._restore_grid_shape(
                w, extra_shape, public_dtype
            ), -self._restore_grid_shape(v, extra_shape, public_dtype)

        raise ValidationError(
            f'{operation_name} invalid component "{component}"; expected "div" or "vrt"'
        )

    def getuv(
        self, vrtspec: ComplexArray, divspec: ComplexArray
    ) -> Tuple[FloatArray, FloatArray]:
        """
        Compute vector wind from spectral coefficients of vorticity and divergence.

        Args:
            vrtspec: Vorticity spectral coefficients
            divspec: Divergence spectral coefficients

        Returns:
            Tuple of (ugrid, vgrid) wind components

        Raises:
            ValidationError: If input arrays have mismatched shapes
            SpheremackError: If SPHEREPACK transform fails
        """
        # Validate inputs
        if vrtspec.shape != divspec.shape:
            raise ValidationError("vrtspec and divspec must have the same shape")

        nt_vrt, ntrunc_vrt, normalized_vrt, extra_shape = self._validate_spectral_data(
            vrtspec, "getuv"
        )
        (
            nt_div,
            ntrunc_div,
            normalized_div,
            div_extra_shape,
        ) = self._validate_spectral_data(divspec, "getuv")

        if (
            nt_vrt != nt_div
            or ntrunc_vrt != ntrunc_div
            or extra_shape != div_extra_shape
        ):
            raise ValidationError("vrtspec and divspec must have consistent dimensions")

        # Convert 1D complex arrays to 2D vector harmonic arrays
        br, bi, cr, ci = _spherepack.onedtotwod_vrtdiv(
            normalized_vrt, normalized_div, self.nlat, self.rsphere
        )
        return self._synthesize_vector_harmonics(
            br, bi, cr, ci, nt_vrt, extra_shape, "getuv", ityp=0
        )

    def getpsichi_component(
        self,
        ugrid: FloatArray,
        vgrid: FloatArray,
        ntrunc: Optional[int],
        component: str,
        operation_name: str,
    ) -> FloatArray:
        """Compute one streamfunction/velocity-potential component safely."""
        fieldspec = self.getpsichi_spec_component(
            ugrid, vgrid, ntrunc, component, operation_name
        )
        return self.spectogrd(fieldspec)

    def getpsichi_spec_component(
        self,
        ugrid: FloatArray,
        vgrid: FloatArray,
        ntrunc: Optional[int],
        component: str,
        operation_name: str,
    ) -> ComplexArray:
        """Compute one streamfunction/velocity-potential spectrum safely."""
        source_spec = (
            self.getvrtspec(ugrid, vgrid, ntrunc)
            if component == "psi"
            else self.getdivspec(ugrid, vgrid, ntrunc)
        )
        return self._invlapspec(source_spec, operation_name)

    def getpsispec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute streamfunction spectral coefficients without gridding."""
        return self.getpsichi_spec_component(
            ugrid, vgrid, ntrunc, component="psi", operation_name="getpsispec"
        )

    def getchispec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> ComplexArray:
        """Compute velocity-potential spectral coefficients without gridding."""
        return self.getpsichi_spec_component(
            ugrid, vgrid, ntrunc, component="chi", operation_name="getchispec"
        )

    def getpsichispec(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> Tuple[ComplexArray, ComplexArray]:
        """Compute streamfunction and velocity-potential spectra safely."""
        if ugrid.shape != vgrid.shape:
            raise ValidationError("ugrid and vgrid must have the same shape")

        _, _, extra_shape = self._validate_grid_data(ugrid, "getpsichispec")
        ntrunc = self._validate_ntrunc(ntrunc, self.nlat - 1)

        vrtspec, divspec = self.getvrtdivspec(ugrid, vgrid, ntrunc)
        _, _, normalized_vrt, spec_extra_shape = self._validate_spectral_data(
            vrtspec, "getpsichispec"
        )
        _, _, normalized_div, div_spec_extra_shape = self._validate_spectral_data(
            divspec, "getpsichispec"
        )
        if spec_extra_shape != extra_shape or div_spec_extra_shape != extra_shape:
            raise ValidationError(
                "vorticity and divergence spectra have inconsistent dimensions"
            )

        psispec = _spherepack.invlap(normalized_vrt, self.rsphere)
        chispec = _spherepack.invlap(normalized_div, self.rsphere)
        return (
            self._restore_spectral_shape(psispec, extra_shape),
            self._restore_spectral_shape(chispec, extra_shape),
        )

    def getpsi(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> FloatArray:
        """
        Compute streamfunction from vector wind without materializing chi.

        Args:
            ugrid: Zonal wind grid data
            vgrid: Meridional wind grid data
            ntrunc: Optional spectral truncation limit

        Returns:
            Streamfunction grid

        Raises:
            ValidationError: If input arrays have invalid dimensions
        """
        return self.getpsichi_component(
            ugrid, vgrid, ntrunc, component="psi", operation_name="getpsi"
        )

    def getchi(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> FloatArray:
        """
        Compute velocity potential from vector wind without materializing psi.

        Args:
            ugrid: Zonal wind grid data
            vgrid: Meridional wind grid data
            ntrunc: Optional spectral truncation limit

        Returns:
            Velocity-potential grid

        Raises:
            ValidationError: If input arrays have invalid dimensions
        """
        return self.getpsichi_component(
            ugrid, vgrid, ntrunc, component="chi", operation_name="getchi"
        )

    def getpsichi(
        self, ugrid: FloatArray, vgrid: FloatArray, ntrunc: Optional[int] = None
    ) -> Tuple[FloatArray, FloatArray]:
        """
        Compute streamfunction and velocity potential from vector wind.

        Args:
            ugrid: Zonal wind grid data
            vgrid: Meridional wind grid data
            ntrunc: Optional spectral truncation limit

        Returns:
            Tuple of (streamfunction_grid, velocity_potential_grid)

        Raises:
            ValidationError: If input arrays have invalid dimensions
        """
        # Validate inputs
        if ugrid.shape != vgrid.shape:
            raise ValidationError("ugrid and vgrid must have the same shape")

        _, _, extra_shape = self._validate_grid_data(ugrid, "getpsichi")
        ntrunc = self._validate_ntrunc(ntrunc, self.nlat - 1)

        # Compute spectral coefficients of vorticity and divergence
        vrtspec, divspec = self.getvrtdivspec(ugrid, vgrid, ntrunc)

        _, _, normalized_vrt, spec_extra_shape = self._validate_spectral_data(
            vrtspec, "getpsichi"
        )
        _, _, normalized_div, div_spec_extra_shape = self._validate_spectral_data(
            divspec, "getpsichi"
        )
        if spec_extra_shape != extra_shape or div_spec_extra_shape != extra_shape:
            raise ValidationError(
                "vorticity and divergence spectra have inconsistent dimensions"
            )

        # Convert to spectral coefficients of streamfunction and velocity potential
        psispec = _spherepack.invlap(normalized_vrt, self.rsphere)
        chispec = _spherepack.invlap(normalized_div, self.rsphere)

        # Transform both scalar spectra back to grid in one batched synthesis.
        return self._spectogrd_pair(
            self._restore_spectral_shape(psispec, extra_shape),
            self._restore_spectral_shape(chispec, extra_shape),
        )

    def getgrad(self, chispec: ComplexArray) -> Tuple[FloatArray, FloatArray]:
        """
        Compute vector gradient from spectral coefficients.

        Args:
            chispec: Spectral coefficients of scalar field

        Returns:
            Tuple of (u_gradient, v_gradient) components

        Raises:
            ValidationError: If input data has invalid dimensions
        """
        nt, ntrunc, normalized_spec, extra_shape = self._validate_spectral_data(
            chispec, "getgrad"
        )

        # Convert chispec to divergence spec using Laplacian
        divspec = _spherepack.lap(normalized_spec, self.rsphere)
        return self._getuv_spec_component(
            self._restore_spectral_shape(divspec, extra_shape),
            component="div",
            operation_name="getgrad",
        )

    def specsmooth(self, datagrid: FloatArray, smooth: FloatArray) -> FloatArray:
        """
        Apply isotropic spectral smoothing to grid data.

        Args:
            datagrid: Grid data to smooth
            smooth: Smoothing factors as function of total wavenumber

        Returns:
            Smoothed grid data

        Raises:
            ValidationError: If input dimensions are invalid
        """
        # Validate inputs
        self._validate_grid_data(datagrid, "specsmooth")

        if smooth.ndim != 1 or smooth.shape[0] != self.nlat:
            raise ValidationError(
                f"smooth must be rank 1 with size {self.nlat}, got shape {smooth.shape}"
            )

        # Grid to spectral transform
        dataspec = self.grdtospec(datagrid, self.nlat - 1)

        # Apply smoothing in spectral space
        _, _, normalized_spec, extra_shape = self._validate_spectral_data(
            dataspec, "specsmooth"
        )
        smoothed_spec = _spherepack.multsmoothfact(normalized_spec, smooth)

        # Spectral to grid transform
        return self.spectogrd(self._restore_spectral_shape(smoothed_spec, extra_shape))


# Standalone functions with type annotations


def regrid(
    grdin: Spharmt,
    grdout: Spharmt,
    datagrid: FloatArray,
    ntrunc: Optional[int] = None,
    smooth: Optional[FloatArray] = None,
) -> FloatArray:
    """
    Regrid data using spectral interpolation with optional smoothing and truncation.

    Args:
        grdin: Input grid Spharmt instance
        grdout: Output grid Spharmt instance
        datagrid: Data on input grid
        ntrunc: Optional spectral truncation limit
        smooth: Optional smoothing factors

    Returns:
        Interpolated data on output grid

    Raises:
        ValidationError: If input parameters are invalid
    """
    datagrid = np.asarray(datagrid)

    # Validate input data dimensions
    if datagrid.ndim < 2:
        raise ValidationError(
            f"regrid needs at least a rank 2 array, got rank {datagrid.ndim}"
        )

    if datagrid.shape[0] != grdin.nlat or datagrid.shape[1] != grdin.nlon:
        raise ValidationError(
            f"regrid needs input array of size {grdin.nlat} by {grdin.nlon}, "
            f"got {datagrid.shape[0]} by {datagrid.shape[1]}"
        )

    if smooth is not None and (smooth.ndim != 1 or smooth.shape[0] != grdout.nlat):
        raise ValidationError("smooth must be rank 1 with size grdout.nlat")

    # Set default truncation
    if ntrunc is None:
        ntrunc = min(grdout.nlat - 1, grdin.nlat - 1)

    # Perform spectral interpolation
    dataspec = grdin.grdtospec(datagrid, ntrunc)

    if smooth is not None:
        _, _, normalized_spec, extra_shape = grdin._validate_spectral_data(
            dataspec, "regrid"
        )
        dataspec = grdin._restore_spectral_shape(
            _spherepack.multsmoothfact(normalized_spec, smooth), extra_shape
        )

    return grdout.spectogrd(dataspec)


def gaussian_lats_wts(nlat: int) -> Tuple[FloatArray, FloatArray]:
    """
    Compute Gaussian latitudes and quadrature weights.

    Args:
        nlat: Number of Gaussian latitudes desired

    Returns:
        Tuple of (latitudes_in_degrees, quadrature_weights)

    Raises:
        SpheremackError: If SPHEREPACK computation fails
    """
    try:
        colats, wts, ierror = _spherepack.gaqd(nlat)
        if ierror != 0:
            raise SpheremackError(f"gaqd failed with error code {ierror}")
    except Exception as e:
        raise SpheremackError(f"Failed to compute Gaussian latitudes: {str(e)}") from e

    # Convert colatitudes to degrees north latitude
    lats = 90.0 - colats * RADIANS_TO_DEGREES

    return lats, wts


def getspecindx(ntrunc: int) -> Tuple[IntArray, IntArray]:
    """
    Compute indices of zonal wavenumber and degree for spherical harmonic coefficients.

    Args:
        ntrunc: Spherical harmonic triangular truncation limit

    Returns:
        Tuple of (zonal_wavenumber_indices, degree_indices)
    """
    if ntrunc < 0:
        raise ValidationError(f"ntrunc must be non-negative, got {ntrunc}")

    indexn = np.indices((ntrunc + 1, ntrunc + 1))[1, :, :]
    indexm = np.indices((ntrunc + 1, ntrunc + 1))[0, :, :]
    indices = np.nonzero(np.greater(indexn, indexm - 1).flatten())
    indxn = np.take(indexn.flatten(), indices)
    indxm = np.take(indexm.flatten(), indices)

    return np.squeeze(indxm), np.squeeze(indxn)


def getgeodesicpts(m: int) -> Tuple[FloatArray, FloatArray]:
    """
    Compute lat/lon values for icosahedral geodesic points.

    Args:
        m: Number of points on edge of single geodesic triangle

    Returns:
        Tuple of (latitudes, longitudes) in degrees

    Raises:
        ValidationError: If m is invalid
    """
    if m < 2:
        raise ValidationError(f"m must be >= 2, got {m}")

    try:
        x, y, z = _spherepack.ihgeod(m)
    except Exception as e:
        raise SpheremackError(f"Failed to compute geodesic points: {str(e)}") from e

    # Convert Cartesian coordinates to lat/lon
    r1 = x * x + y * y
    r = np.sqrt(r1 + z * z)
    r1 = np.sqrt(r1)

    xtmp = np.where(np.logical_or(x, y), x, np.ones(x.shape, np.float32))
    ztmp = np.where(np.logical_or(r1, z), z, np.ones(z.shape, np.float32))

    lons = RADIANS_TO_DEGREES * np.arctan2(y, xtmp) + 180.0
    lats = RADIANS_TO_DEGREES * np.arctan2(r1, ztmp) - 90.0

    # Initialize output arrays
    total_points = 10 * (m - 1) ** 2 + 2
    lat = np.zeros(total_points, np.float32)
    lon = np.zeros(total_points, np.float32)

    # First two points are poles
    lat[0] = 90.0
    lat[1] = -90.0
    lon[0] = 0.0
    lon[1] = 0.0

    # Fill remaining points
    lat[2:] = lats[0 : 2 * (m - 1), 0 : m - 1, :].flatten()
    lon[2:] = lons[0 : 2 * (m - 1), 0 : m - 1, :].flatten()

    return lat, lon


def legendre(lat: float, ntrunc: int) -> FloatArray:
    """
    Calculate associated Legendre functions for triangular truncation.

    Args:
        lat: Latitude in degrees
        ntrunc: Triangular truncation limit

    Returns:
        Associated Legendre functions array

    Raises:
        ValidationError: If parameters are invalid
    """
    if not -90.0 <= lat <= 90.0:
        raise ValidationError(f"lat must be between -90 and 90 degrees, got {lat}")

    if ntrunc < 0:
        raise ValidationError(f"ntrunc must be non-negative, got {ntrunc}")

    try:
        return _spherepack.getlegfunc(lat, ntrunc)
    except Exception as e:
        raise SpheremackError(f"Failed to compute Legendre functions: {str(e)}") from e


def specintrp(
    lon: float, dataspec: ComplexArray, legfuncs: FloatArray
) -> Union[float, complex]:
    """
    Spectral interpolation to arbitrary point on sphere.

    Args:
        lon: Longitude in degrees
        dataspec: Spectral coefficients
        legfuncs: Associated Legendre functions

    Returns:
        Interpolated value

    Raises:
        ValidationError: If inputs have inconsistent truncation limits
        SpheremackError: If interpolation fails
    """
    # Validate truncation consistency
    ntrunc1 = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - dataspec.shape[0])))
    ntrunc2 = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - legfuncs.shape[0])))

    if ntrunc1 != ntrunc2:
        raise ValidationError(
            f"dataspec and legfuncs have inconsistent spectral truncations: {ntrunc1} vs {ntrunc2}"
        )

    try:
        return _spherepack.specintrp(
            lon * DEGREES_TO_RADIANS, ntrunc1, dataspec, legfuncs
        )
    except Exception as e:
        raise SpheremackError(f"Spectral interpolation failed: {str(e)}") from e
