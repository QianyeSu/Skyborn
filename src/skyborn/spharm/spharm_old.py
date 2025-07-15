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
 The SPHERPACK fortran source files will be downloaded automatically by the
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

import sys
import math
import numpy

from . import _spherepack

# define a list of instance variables that cannot be rebound
# or unbound.
_private_vars = ["nlon", "nlat", "gridtype", "legfunc", "rsphere"]
__version__ = "1.0.9"


class Spharmt:
    """
    spherical harmonic transform class.

    @ivar nlat: number of latitudes (set when class instance is created,
    cannot be changed).

    @ivar nlon: number of longitudes (set when class instance is created,
    cannot be changed).

    @ivar rsphere: The radius of the sphere in meters (set when class
    instance is created, cannot be changed).

    @ivar legfunc: 'stored' or 'computed'.  If 'stored',
    associated legendre functions are precomputed and stored when the
    class instance is created.  If 'computed', associated
    legendre functions are computed on the fly when transforms are
    requested. Set when class instance is created, cannot be changed.

    @ivar gridtype: 'regular' (equally spaced in longitude and latitude)
    or 'gaussian' (equally spaced in longitude, latitudes located at
    roots of ordinary Legendre polynomial of degree nlat). Set when class
    instance is created, cannot be changed.
    """

    def __setattr__(self, key, val):
        """
        prevent modification of read-only instance variables.
        """
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError("Attempt to rebind read-only instance variable " + key)
        else:
            self.__dict__[key] = val

    def __delattr__(self, key):
        """
        prevent deletion of read-only instance variables.
        """
        if key in self.__dict__ and key in _private_vars:
            raise AttributeError("Attempt to unbind read-only instance variable " + key)
        else:
            del self.__dict__[key]

    def __init__(
        self, nlon, nlat, rsphere=6.3712e6, gridtype="regular", legfunc="stored"
    ):
        """
        create a Spharmt class instance.

        @param nlon: Number of longitudes. The grid must be oriented from
        east to west, with the first point at the Greenwich meridian
        and the last point at 360-delta degrees east
        (where delta = 360/nlon degrees). Must be >= 4. Transforms will
        be faster when nlon is the product of small primes.

        @param nlat: Number of latitudes.  The grid must be oriented from north
        to south. If nlat is odd the equator is included.
        If nlat is even the equator will lie half way between points
        points nlat/2 and (nlat/2)+1. Must be >=3.

        @keyword rsphere: The radius of the sphere in meters.
        Default 6371200 (the value for Earth).

        @keyword legfunc: 'stored' (default) or 'computed'.  If 'stored',
        associated legendre functions are precomputed and stored when the
        class instance is created.  This uses O(nlat**3) memory, but
        speeds up the spectral transforms.  If 'computed', associated
        legendre functions are computed on the fly when transforms are
        requested.  This uses O(nlat**2) memory, but slows down the spectral
        transforms a bit.

        @keyword gridtype: 'regular' (default) or 'gaussian'. Regular grids
        will include the poles and equator if nlat is odd.  Gaussian
        grids never include the poles, but will include the equator if
        nlat is odd.

        """
        # sanity checks.
        if rsphere > 0.0:
            self.rsphere = rsphere
        else:
            msg = (
                "Spharmt.__init__ illegal value of rsphere (%s) - must be postitive"
                % (rsphere)
            )
            raise ValueError(msg)
        if nlon > 3:
            self.nlon = nlon
        else:
            msg = "Spharmt.__init__ illegal value of nlon (%s) - must be at least 4" % (
                nlon,
            )
            raise ValueError(msg)
        if nlat > 2:
            self.nlat = nlat
        else:
            msg = "Spharmt.__init__ illegal value of nlat (%s) - must be at least 3" % (
                nlat,
            )
            raise ValueError(msg)
        if gridtype != "regular" and gridtype != "gaussian":
            msg = (
                'Spharmt.__init__ illegal value of gridtype (%s) - must be either "gaussian" or "regular"'
                % gridtype
            )
            raise ValueError(msg)
        else:
            self.gridtype = gridtype

        if legfunc != "computed" and legfunc != "stored":
            msg = (
                'Spharmt.__init__ illegal value of legfunc (%s) - must be either "computed" or "stored"'
                % legfunc
            )
            raise ValueError(msg)
        else:
            self.legfunc = legfunc

        if nlon % 2:  # nlon is odd
            n1 = min(nlat, (nlon + 1) / 2)
        else:
            n1 = min(nlat, (nlon + 2) / 2)
        if nlat % 2:  # nlat is odd
            n2 = (nlat + 1) / 2
        else:
            n2 = nlat / 2

        if gridtype == "regular":
            if legfunc == "stored":
                lshaes = (n1 * n2 * (nlat + nlat - n1 + 1)) / 2 + nlon + 15
                lwork = 5 * nlat * n2 + 3 * ((n1 - 2) * (nlat + nlat - n1 - 1)) / 2
                wshaes, ierror = _spherepack.shaesi(nlat, nlon, lshaes, lwork, nlat + 1)
                if ierror != 0:
                    msg = (
                        "In return from call to shaesi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshaes = wshaes
                lshses = lshaes
                wshses, ierror = _spherepack.shsesi(nlat, nlon, lshses, lwork, nlat + 1)
                if ierror != 0:
                    msg = (
                        "In return from call to shsesi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshses = wshses
                lvhaes = n1 * n2 * (nlat + nlat - n1 + 1) + nlon + 15
                lwork = (
                    3 * (max(n1 - 2, 0) * (nlat + nlat - n1 - 1)) / 2 + 5 * n2 * nlat
                )
                wvhaes, ierror = _spherepack.vhaesi(
                    nlat, nlon, lvhaes, lwork, 2 * (nlat + 1)
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhaesi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhaes = wvhaes
                lwork = (
                    3 * (max(n1 - 2, 0) * (nlat + nlat - n1 - 1)) / 2 + 5 * n2 * nlat
                )
                lvhses = n1 * n2 * (nlat + nlat - n1 + 1) + nlon + 15
                wvhses, ierror = _spherepack.vhsesi(
                    nlat, nlon, lvhses, lwork, 2 * (nlat + 1)
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhsesi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhses = wvhses
            else:
                lshaec = (
                    2 * nlat * n2
                    + 3 * ((n1 - 2) * (nlat + nlat - n1 - 1)) / 2
                    + nlon
                    + 15
                )
                wshaec, ierror = _spherepack.shaeci(nlat, nlon, lshaec, 2 * (nlat + 1))
                if ierror != 0:
                    msg = (
                        "In return from call to shaeci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshaec = wshaec
                lshsec = lshaec
                wshsec, ierror = _spherepack.shseci(nlat, nlon, lshsec, 2 * (nlat + 1))
                if ierror != 0:
                    msg = (
                        "In return from call to shseci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshsec = wshsec
                lvhaec = (
                    4 * nlat * n2 + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1) + nlon + 15
                )
                wvhaec, ierror = _spherepack.vhaeci(nlat, nlon, lvhaec, 2 * (nlat + 1))
                if ierror != 0:
                    msg = (
                        "In return from call to vhaeci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhaec = wvhaec
                lvhsec = lvhaec
                wvhsec, ierror = _spherepack.vhseci(nlat, nlon, lvhsec, 2 * (nlat + 1))
                if ierror != 0:
                    msg = (
                        "In return from call to vhseci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhsec = wvhsec

        elif gridtype == "gaussian":
            if legfunc == "stored":
                lshags = (
                    nlat * (3 * (n1 + n2) - 2)
                    + (n1 - 1) * (n2 * (2 * nlat - n1) - 3 * n1) / 2
                    + nlon
                    + 15
                )
                lwork = 4 * nlat * (nlat + 2) + 2
                ldwork = nlat * (nlat + 4)
                wshags, ierror = _spherepack.shagsi(nlat, nlon, lshags, lwork, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shagsi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshags = wshags
                lshsgs = lshags
                wshsgs, ierror = _spherepack.shsgsi(nlat, nlon, lshsgs, lwork, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shsgsi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshsgs = wshsgs
                lvhags = (nlat + 1) * (nlat + 1) * nlat / 2 + nlon + 15
                ldwork = (3 * nlat * (nlat + 3) + 2) / 2
                wvhags, ierror = _spherepack.vhagsi(nlat, nlon, lvhags, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhagsi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhags = wvhags
                lvhsgs = n1 * n2 * (nlat + nlat - n1 + 1) + nlon + 15 + 2 * nlat
                ldwork = (3 * nlat * (nlat + 3) + 2) / 2
                wvhsgs, ierror = _spherepack.vhsgsi(nlat, nlon, lvhsgs, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhsgsi in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhsgs = wvhsgs
            else:
                lshagc = (
                    nlat * (2 * n2 + 3 * n1 - 2) + 3 * n1 * (1 - n1) / 2 + nlon + 15
                )
                wshagc, ierror = _spherepack.shagci(
                    nlat, nlon, lshagc, nlat * (nlat + 4)
                )
                if ierror != 0:
                    msg = (
                        "In return from call to shagci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshagc = wshagc
                lshsgc = lshagc
                wshsgc, ierror = _spherepack.shsgci(
                    nlat, nlon, lshsgc, nlat * (nlat + 4)
                )
                if ierror != 0:
                    msg = (
                        "In return from call to shsgci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wshsgc = wshsgc
                lvhagc = (
                    4 * nlat * n2
                    + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1)
                    + nlon
                    + n2
                    + 15
                )
                ldwork = 2 * nlat * (nlat + 1) + 1
                wvhagc, ierror = _spherepack.vhagci(nlat, nlon, lvhagc, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhagci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhagc = wvhagc
                lvhsgc = (
                    4 * nlat * n2 + 3 * max(n1 - 2, 0) * (2 * nlat - n1 - 1) + nlon + 15
                )
                wvhsgc, ierror = _spherepack.vhsgci(nlat, nlon, lvhsgc, ldwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhsgci in Spharmt.__init__ ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
                self.wvhsgc = wvhsgc

    def grdtospec(self, datagrid, ntrunc=None):
        """
        grid to spectral transform (spherical harmonic analysis).

        @param datagrid: rank 2 or 3 numpy float32 array with shape (nlat,nlon) or
        (nlat,nlon,nt), where nt is the number of grids to be transformed.  If
        datagrid is rank 2, nt is assumed to be 1.

        @keyword ntrunc:  optional spectral truncation limit.
        (default self.nlat-1)

        @return: C{B{dataspec}} - rank 1 or 2 numpy complex array with shape
        (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing
        complex spherical harmonic coefficients resulting from the spherical
        harmonic analysis of datagrid.
        """

        # check that datagrid is rank 2 or 3 with size (self.nlat, self.nlon) or
        # (self.nlat, self.nlon, nt) where nt is number of grids to transform.

        idim = datagrid.ndim

        if len(datagrid.shape) > 3:
            msg = "grdtospec needs a rank two or three array, got %d" % (
                len(datagrid.shape),
            )
            raise ValueError(msg)

        if datagrid.shape[0] != self.nlat or datagrid.shape[1] != self.nlon:
            msg = "grdtospec needs an array of size %d by %d, got %d by %d" % (
                self.nlat,
                self.nlon,
                datagrid.shape[0],
                datagrid.shape[1],
            )
            raise ValueError(msg)

        # check ntrunc.

        if ntrunc is None:
            ntrunc = self.nlat - 1

        if ntrunc < 0 or ntrunc + 1 > datagrid.shape[0]:
            msg = "ntrunc must be between 0 and %d" % (datagrid.shape[0] - 1,)
            raise ValueError(msg)

        nlat = self.nlat
        nlon = self.nlon
        if nlat % 2:  # nlat is odd
            n2 = (nlat + 1) / 2
        else:
            n2 = nlat / 2

        if len(datagrid.shape) == 2:
            nt = 1
            datagrid = numpy.reshape(datagrid, (nlat, nlon, 1))
        else:
            nt = datagrid.shape[2]

        if self.gridtype == "regular":

            # do grid to spectral transform.

            if self.legfunc == "stored":
                lwork = (nt + 1) * nlat * nlon
                a, b, ierror = _spherepack.shaes(datagrid, self.wshaes, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shaes in Spharmt.grdtospec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (nt * nlon + max(3 * n2, nlon))
                a, b, ierror = _spherepack.shaec(datagrid, self.wshaec, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shaec in Spharmt.grdtospec ierror =  %d"
                        % ierror
                    )

                    raise ValueError(msg)

        # gaussian grid.

        elif self.gridtype == "gaussian":

            # do grid to spectral transform.

            if self.legfunc == "stored":
                lwork = nlat * nlon * (nt + 1)
                a, b, ierror = _spherepack.shags(datagrid, self.wshags, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shags in Spharmt.grdtospec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (nlon * nt + max(3 * n2, nlon))
                a, b, ierror = _spherepack.shagc(datagrid, self.wshagc, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shagc in Spharmt.grdtospec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # convert 2d real and imag spectral arrays into 1d complex array.

        dataspec = _spherepack.twodtooned(a, b, ntrunc)

        if idim == 2:
            return numpy.squeeze(dataspec)
        else:
            return dataspec

    def spectogrd(self, dataspec):
        """
        spectral to grid transform (spherical harmonic synthesis).

        @param dataspec: rank 1 or 2 numpy complex array with shape
        (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing
        complex spherical harmonic coefficients (where ntrunc is the
        triangular truncation limit and nt is the number of spectral arrays
        to be transformed). If dataspec is rank 1, nt is assumed to be 1.

        @return: C{B{datagrid}} - rank 2 or 3 numpy float32 array with shape
        (nlat,nlon) or (nlat,nlon,nt) containing the gridded data resulting from
        the spherical harmonic synthesis of dataspec.
        """

        # make sure dataspec is rank 1 or 2.

        idim = dataspec.ndim

        if len(dataspec.shape) > 2:
            msg = "spectogrd needs a rank one or two array, got %d" % (
                len(dataspec.shape),
            )
            raise ValueError(msg)

        nlat = self.nlat
        nlon = self.nlon
        if nlat % 2:  # nlat is odd
            n2 = (nlat + 1) / 2
        else:
            n2 = nlat / 2

        if len(dataspec.shape) == 1:
            nt = 1
            dataspec = numpy.reshape(dataspec, (dataspec.shape[0], 1))
        else:
            nt = dataspec.shape[1]

        ntrunc = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - dataspec.shape[0])))
        if ntrunc > nlat - 1:
            msg = "ntrunc too large - can be max of %d, got %d" % (nlat - 1, ntrunc)
            raise ValueError(msg)

        a, b = _spherepack.onedtotwod(dataspec, nlat)

        # regular grid.

        if self.gridtype == "regular":

            # do spectral to grid transform.

            if self.legfunc == "stored":
                lwork = (nt + 1) * nlat * nlon
                datagrid, ierror = _spherepack.shses(nlon, a, b, self.wshses, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shses in Spharmt.spectogrd ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (nt * nlon + max(3 * n2, nlon))
                datagrid, ierror = _spherepack.shsec(nlon, a, b, self.wshsec, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shsec in Spharmt.spectogrd ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # gaussian grid.

        elif self.gridtype == "gaussian":

            # do spectral to grid transform.

            if self.legfunc == "stored":
                lwork = nlat * nlon * (nt + 1)
                datagrid, ierror = _spherepack.shsgs(nlon, a, b, self.wshsgs, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shsgs in Spharmt.spectogrd ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (nlon * nt + max(3 * n2, nlon))
                datagrid, ierror = _spherepack.shsgc(nlon, a, b, self.wshsgc, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to shsgc in Spharmt.spectogrd ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        if idim == 1:
            return numpy.squeeze(datagrid)
        else:
            return datagrid

    def getvrtdivspec(self, ugrid, vgrid, ntrunc=None):
        """
        compute spectral coefficients of vorticity and divergence given vector wind.

        @param ugrid: rank 2 or 3 numpy float32 array containing grid of zonal
        winds.  Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number
        of grids to be transformed.  If ugrid is rank 2, nt is assumed to be 1.

        @param vgrid: rank 2 or 3 numpy float32 array containing grid of meridional
        winds.  Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number
        of grids to be transformed.  Both ugrid and vgrid must have the same shape.

        @keyword ntrunc:  optional spectral truncation limit.
        (default self.nlat-1)

        @return: C{B{vrtspec, divspec}} - rank 1 or 2 numpy complex arrays
        of vorticity and divergence spherical harmonic coefficients with shape
        shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt).
        """

        # make sure ugrid,vgrid are rank 2 or 3 and same shape.

        idim = ugrid.ndim

        shapeu = ugrid.shape
        shapev = vgrid.shape

        if ntrunc is None:
            ntrunc = self.nlat - 1

        if shapeu != shapev:
            msg = "getvrtdivspec input arrays must be same shape!"
            raise ValueError(msg)

        if len(shapeu) != 2 and len(shapeu) != 3:
            msg = "getvrtdivspec needs rank two or three arrays!"
            raise ValueError(msg)

        if shapeu[0] != self.nlat or shapeu[1] != self.nlon:
            msg = (
                "getvrtdivspec needs input arrays whose first two dimensions are si%d and %d, got %d and %d"
                % (
                    self.nlat,
                    self.nlon,
                    ugrid.shape[0],
                    ugrid.shape[1],
                )
            )
            raise ValueError(msg)

        # check ntrunc.

        if ntrunc < 0 or ntrunc + 1 > shapeu[0]:
            msg = "ntrunc must be between 0 and %d" % (ugrid.shape[0] - 1,)
            raise ValueError(msg)

        nlat = self.nlat
        nlon = self.nlon
        if nlat % 2:  # nlat is odd
            n2 = (nlat + 1) / 2
        else:
            n2 = nlat / 2
        rsphere = self.rsphere

        # convert from geographical to math coordinates, add extra dimension
        # if necessary.

        if len(shapeu) == 2:
            nt = 1
            w = numpy.reshape(ugrid, (nlat, nlon, 1))
            v = -numpy.reshape(vgrid, (nlat, nlon, 1))
        else:
            nt = shapeu[2]
            w = ugrid
            v = -vgrid

        # regular grid.

        if self.gridtype == "regular":

            # vector harmonic analysis.

            if self.legfunc == "stored":
                lwork = (2 * nt + 1) * nlat * nlon
                br, bi, cr, ci, ierror = _spherepack.vhaes(v, w, self.wvhaes, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhaes in Spharmt.getvrtdivspec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (2 * nt * nlon + max(6 * n2, nlon))
                br, bi, cr, ci, ierror = _spherepack.vhaec(v, w, self.wvhaec, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhaec in Spharmt.getvrtdivspec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # gaussian grid.

        elif self.gridtype == "gaussian":

            # vector harmonic analysis.

            if self.legfunc == "stored":
                lwork = (2 * nt + 1) * nlat * nlon
                br, bi, cr, ci, ierror = _spherepack.vhags(v, w, self.wvhags, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhags in Spharmt.getvrtdivspec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = 2 * nlat * (2 * nlon * nt + 3 * n2)
                br, bi, cr, ci, ierror = _spherepack.vhagc(v, w, self.wvhagc, lwork)
                if ierror != 0:
                    msg = (
                        "In return from call to vhagc in Spharmt.getvrtdivspec ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # convert vector harmonic coeffs to 1d complex coefficients
        # of vorticity and divergence.

        vrtspec, divspec = _spherepack.twodtooned_vrtdiv(
            br, bi, cr, ci, ntrunc, rsphere
        )

        if idim == 2:
            return numpy.squeeze(vrtspec), numpy.squeeze(divspec)
        else:
            return vrtspec, divspec

    def getuv(self, vrtspec, divspec):
        """
        compute vector wind on grid given complex spectral coefficients
        of vorticity and divergence.

        @param vrtspec: rank 1 or 2 numpy complex array of vorticity spectral
        coefficients, with shape (ntrunc+1)*(ntrunc+2)/2 or
        ((ntrunc+1)*(ntrunc+2)/2,nt) (where ntrunc is the triangular truncation
        and nt is the number of spectral arrays to be transformed).
        If vrtspec is rank 1, nt is assumed to be 1.

        @param divspec: rank 1 or 2 numpy complex array of divergence spectral
        coefficients, with shape (ntrunc+1)*(ntrunc+2)/2 or
        ((ntrunc+1)*(ntrunc+2)/2,nt) (where ntrunc is the triangular truncation
        and nt is the number of spectral arrays to be transformed).
        Both vrtspec and divspec must have the same shape.

        @return: C{B{ugrid, vgrid}} - rank 2 or 3 numpy float32 arrays containing
        gridded zonal and meridional winds. Shapes are either (nlat,nlon) or
        (nlat,nlon,nt).
        """

        idim = vrtspec.ndim

        shapevrt = vrtspec.shape
        shapediv = divspec.shape

        # make sure vrtspec, divspec are rank 1 or 2, and have the same shape.

        if shapevrt != shapediv:
            msg = "vrtspec, divspec must be same size in getuv!"
            raise ValueError(msg)

        if len(shapevrt) != 1 and len(shapevrt) != 2:
            msg = "getuv needs rank one or two input arrays!"
            raise ValueError(msg)

        # infer ntrunc from size of dataspec (dataspec must be rank 1!)
        # dataspec is assumed to have size (ntrunc+1)*(ntrunc+2)/2

        ntrunc = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - vrtspec.shape[0])))

        nlat = self.nlat
        nlon = self.nlon
        if nlat % 2:  # nlat is odd
            n2 = (nlat + 1) / 2
        else:
            n2 = nlat / 2
        rsphere = self.rsphere

        if len(vrtspec.shape) == 1:
            nt = 1
            vrtspec = numpy.reshape(vrtspec, (vrtspec.shape[0], 1))
            divspec = numpy.reshape(divspec, (divspec.shape[0], 1))
        else:
            nt = vrtspec.shape[1]

        # convert 1d complex arrays of vort, div to 2d vector harmonic arrays.

        br, bi, cr, ci = _spherepack.onedtotwod_vrtdiv(vrtspec, divspec, nlat, rsphere)

        # regular grid.

        if self.gridtype == "regular":

            # vector harmonic synthesis.

            if self.legfunc == "stored":
                lwork = (2 * nt + 1) * nlat * nlon
                v, w, ierror = _spherepack.vhses(
                    nlon, br, bi, cr, ci, self.wvhses, lwork
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhses in Spharmt.getuv ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (2 * nt * nlon + max(6 * n2, nlon))
                v, w, ierror = _spherepack.vhsec(
                    nlon, br, bi, cr, ci, self.wvhsec, lwork
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhsec in Spharmt.getuv ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # gaussian grid.

        elif self.gridtype == "gaussian":

            # vector harmonic synthesis.

            if self.legfunc == "stored":
                lwork = (2 * nt + 1) * nlat * nlon
                v, w, ierror = _spherepack.vhsgs(
                    nlon, br, bi, cr, ci, self.wvhsgs, lwork
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhsgs in Spharmt.getuv ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)
            else:
                lwork = nlat * (2 * nt * nlon + max(6 * n2, nlon))
                v, w, ierror = _spherepack.vhsgc(
                    nlon, br, bi, cr, ci, self.wvhsgc, lwork
                )
                if ierror != 0:
                    msg = (
                        "In return from call to vhsgc in Spharmt.getuv ierror =  %d"
                        % ierror
                    )
                    raise ValueError(msg)

        # convert to u and v in geographical coordinates.

        if idim == 1:
            return numpy.reshape(w, (nlat, nlon)), -numpy.reshape(v, (nlat, nlon))
        else:
            return w, -v

    def getpsichi(self, ugrid, vgrid, ntrunc=None):
        """
        compute streamfunction and velocity potential on grid given vector wind.

        @param ugrid: rank 2 or 3 numpy float32 array containing grid of zonal
        winds.  Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number
        of grids to be transformed.  If ugrid is rank 2, nt is assumed to be 1.

        @param vgrid: rank 2 or 3 numpy float32 array containing grid of meridional
        winds.  Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number
        of grids to be transformed.  Both ugrid and vgrid must have the same shape.

        @keyword ntrunc:  optional spectral truncation limit.
        (default self.nlat-1)

        @return: C{B{psigrid, chigrid}} - rank 2 or 3 numpy float32 arrays
        of gridded streamfunction and velocity potential. Shapes are either
        (nlat,nlon) or (nlat,nlon,nt).
        """

        # make sure ugrid,vgrid are rank 2 or 3 and same shape.

        idim = ugrid.ndim

        shapeu = ugrid.shape
        shapev = vgrid.shape

        if ntrunc is None:
            ntrunc = self.nlat - 1

        if shapeu != shapev:
            msg = "getvrtdivspec input arrays must be same shape!"
            raise ValueError(msg)

        if len(shapeu) != 2 and len(shapeu) != 3:
            msg = "getvrtdivspec needs rank two or three arrays!"
            raise ValueError(msg)

        if shapeu[0] != self.nlat or shapeu[1] != self.nlon:
            msg = (
                "getpsichi needs input arrays whose first two dimensions are si%d and %d, got %d and %d"
                % (
                    self.nlat,
                    self.nlon,
                    ugrid.shape[0],
                    ugrid.shape[1],
                )
            )
            raise ValueError(msg)

        # check ntrunc.

        if ntrunc < 0 or ntrunc + 1 > ugrid.shape[0]:
            msg = "ntrunc must be between 0 and %d" % (ugrid.shape[0] - 1,)
            raise ValueError(msg)

        # compute spectral coeffs of vort, div.

        vrtspec, divspec = self.getvrtdivspec(ugrid, vgrid, ntrunc)

        # number of grids to compute.

        if len(vrtspec.shape) == 1:
            nt = 1
            vrtspec = numpy.reshape(vrtspec, ((ntrunc + 1) * (ntrunc + 2) // 2, 1))
            divspec = numpy.reshape(divspec, ((ntrunc + 1) * (ntrunc + 2) // 2, 1))
        else:
            nt = vrtspec.shape[1]

        # convert to spectral coeffs of psi, chi.

        psispec = _spherepack.invlap(vrtspec, self.rsphere)
        chispec = _spherepack.invlap(divspec, self.rsphere)

        # inverse transform to grid.

        psigrid = self.spectogrd(psispec)
        chigrid = self.spectogrd(chispec)

        if idim == 2:
            return numpy.squeeze(psigrid), numpy.squeeze(chigrid)
        else:
            return psigrid, chigrid

    def getgrad(self, chispec):
        """
        compute vector gradient on grid given complex spectral coefficients.

        @param chispec: rank 1 or 2 numpy complex array with shape
        (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing
        complex spherical harmonic coefficients (where ntrunc is the
        triangular truncation limit and nt is the number of spectral arrays
        to be transformed). If chispec is rank 1, nt is assumed to be 1.

        @return: C{B{uchi, vchi}} - rank 2 or 3 numpy float32 arrays containing
        gridded zonal and meridional components of the vector gradient.
        Shapes are either (nlat,nlon) or (nlat,nlon,nt).
        """

        # make sure chispec is rank 1 or 2.

        idim = chispec.ndim

        if len(chispec.shape) != 1 and len(chispec.shape) != 2:
            msg = "getgrad needs rank one or two arrays!"
            raise ValueError(msg)

        # infer ntrunc from size of chispec (chispec must be rank 1!)
        # chispec is assumed to have size (ntrunc+1)*(ntrunc+2)/2

        ntrunc = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - chispec.shape[0])))

        # number of grids to compute.

        if len(chispec.shape) == 1:
            nt = 1
            chispec = numpy.reshape(chispec, ((ntrunc + 1) * (ntrunc + 2) // 2, 1))
        else:
            nt = chispec.shape[1]

        # convert chispec to divspec.

        divspec = _spherepack.lap(chispec, self.rsphere)

        # call getuv, with vrtspec=0, to get uchi,vchi.

        uchi, vchi = self.getuv(numpy.zeros(chispec.shape, chispec.dtype), divspec)

        if idim == 1:
            return numpy.squeeze(uchi), numpy.squeeze(vchi)
        else:
            return uchi, vchi

    def specsmooth(self, datagrid, smooth):
        """
        isotropic spectral smoothing on a sphere.

        @param datagrid: rank 2 or 3 numpy float32 array with shape (nlat,nlon) or
        (nlat,nlon,nt), where nt is the number of grids to be smoothed.  If
        datagrid is rank 2, nt is assumed to be 1.

        @param smooth: rank 1 array of length nlat containing smoothing factors
        as a function of total wavenumber.

        @return: C{B{datagrid}} - rank 2 or 3 numpy float32 array with shape
        (nlat,nlon) or (nlat,nlon,nt) containing the smoothed grids.
        """

        # check that datagrid is rank 2 or 3 with size (self.nlat, self.nlon) or
        # (self.nlat, self.nlon, nt) where nt is number of grids to transform.

        if len(datagrid.shape) > 3:
            msg = "specsmooth needs a rank two or three array, got %d" % (
                len(datagrid.shape),
            )
            raise ValueError(msg)

        if datagrid.shape[0] != self.nlat or datagrid.shape[1] != self.nlon:
            msg = "specsmooth needs an array of size %d by %d, got %d by %d" % (
                self.nlat,
                self.nlon,
                datagrid.shape[0],
                datagrid.shape[1],
            )
            raise ValueError(msg)

        # make sure smooth is rank 1, same size as datagrid.shape[0]

        if len(smooth.shape) != 1 or smooth.shape[0] != datagrid.shape[0]:
            msg = "smooth must be rank 1 and same size as datagrid.shape[0] in specsmooth!"
            raise ValueError(msg)

        # grid to spectral transform.

        nlat = self.nlat
        dataspec = self.grdtospec(datagrid, nlat - 1)

        # multiply spectral coeffs. by smoothing factor.

        dataspec = _spherepack.multsmoothfact(dataspec, smooth)

        # spectral to grid transform.

        datagrid = self.spectogrd(dataspec)

        return datagrid


def regrid(grdin, grdout, datagrid, ntrunc=None, smooth=None):
    """
    regrid data using spectral interpolation, while performing
    optional spectral smoothing and/or truncation.

    @param grdin: Spharmt class instance describing input grid.

    @param grdout: Spharmt class instance describing output grid.

    @param datagrid: data on input grid (grdin.nlat x grdin.nlon). If
    datagrid is rank 3, last dimension is the number of grids to interpolate.

    @keyword ntrunc:  optional spectral truncation limit for datagrid
    (default min(grdin.nlat-1,grdout.nlat-1)).

    @keyword smooth: rank 1 array of length grdout.nlat containing smoothing
    factors as a function of total wavenumber (default is no smoothing).

    @return: C{B{datagrid}} - interpolated (and optionally smoothed) array(s)
    on grdout.nlon x grdout.nlat grid.
    """

    # check that datagrid is rank 2 or 3 with size (grdin.nlat, grdin.nlon) or
    # (grdin.nlat, grdin.nlon, nt) where nt is number of grids to transform.

    if len(datagrid.shape) > 3:
        msg = "regrid needs a rank two or three array, got %d" % (len(datagrid.shape),)
        raise ValueError(msg)

    if datagrid.shape[0] != grdin.nlat or datagrid.shape[1] != grdin.nlon:
        msg = "grdtospec needs an array of size %d by %d, got %d by %d" % (
            grdin.nlat,
            grdin.nlon,
            datagrid.shape[0],
            datagrid.shape[1],
        )
        raise ValueError(msg)

    if smooth is not None and (
        len(smooth.shape) != 1 or smooth.shape[0] != grdout.nlat
    ):
        msg = "smooth must be rank 1 size grdout.nlat in regrid!"
        raise ValueError(msg)

    if ntrunc is None:
        ntrunc = min(grdout.nlat - 1, grdin.nlat - 1)

    dataspec = grdin.grdtospec(datagrid, ntrunc)

    if smooth is not None:
        dataspec = _spherepack.multsmoothfact(dataspec, smooth)

    return grdout.spectogrd(dataspec)


def gaussian_lats_wts(nlat):
    """
    compute the gaussian latitudes (in degrees) and quadrature weights.

    @param nlat: number of gaussian latitudes desired.

    @return: C{B{lats, wts}} - rank 1 numpy float64 arrays containing
    gaussian latitudes (in degrees north) and gaussian quadrature weights.
    """

    # get the gaussian colatitudes and weights using gaqd.

    colats, wts, ierror = _spherepack.gaqd(nlat)

    if ierror != 0:
        msg = "In return from call to gaqd ierror =  %d" % ierror
        raise ValueError(msg)

    # convert to degrees north latitude.

    lats = 90.0 - colats * 180.0 / math.pi

    return lats, wts


def getspecindx(ntrunc):
    """
    compute indices of zonal wavenumber (indxm) and degree (indxn)
    for complex spherical harmonic coefficients.

    @param ntrunc: spherical harmonic triangular truncation limit.

    @return: C{B{indxm, indxn}} - rank 1 numpy Int32 arrays
    containing zonal wavenumber (indxm) and degree (indxn) of
    spherical harmonic coefficients.
    """

    indexn = numpy.indices((ntrunc + 1, ntrunc + 1))[1, :, :]
    indexm = numpy.indices((ntrunc + 1, ntrunc + 1))[0, :, :]
    indices = numpy.nonzero(numpy.greater(indexn, indexm - 1).flatten())
    indxn = numpy.take(indexn.flatten(), indices)
    indxm = numpy.take(indexm.flatten(), indices)

    return numpy.squeeze(indxm), numpy.squeeze(indxn)


def getgeodesicpts(m):
    """
    computes the lat/lon values of the points on the surface of the sphere
    corresponding to a twenty-sided (icosahedral) geodesic.

    @param m: the number of points on the edge of a single geodesic triangle.
    There are 10*(m-1)**2+2 total geodesic points, including the poles.

    @return: C{B{lats, lons}} - rank 1 numpy float32 arrays containing
    the latitudes and longitudes of the geodesic points (in degrees). These
    points are nearly evenly distributed on the surface of the sphere.
    """
    x, y, z = _spherepack.ihgeod(m)
    # convert cartesian coords to lat/lon.
    rad2dg = 180.0 / math.pi
    r1 = x * x + y * y
    r = numpy.sqrt(r1 + z * z)
    r1 = numpy.sqrt(r1)
    xtmp = numpy.where(numpy.logical_or(x, y), x, numpy.ones(x.shape, numpy.float32))
    ztmp = numpy.where(numpy.logical_or(r1, z), z, numpy.ones(z.shape, numpy.float32))
    lons = rad2dg * numpy.arctan2(y, xtmp) + 180.0
    lats = rad2dg * numpy.arctan2(r1, ztmp) - 90.0
    lat = numpy.zeros(10 * (m - 1) ** 2 + 2, numpy.float32)
    lon = numpy.zeros(10 * (m - 1) ** 2 + 2, numpy.float32)
    # first two points are poles.
    lat[0] = 90
    lat[1] = -90.0
    lon[0] = 0.0
    lon[1] = 0.0
    lat[2:] = lats[0 : 2 * (m - 1), 0 : m - 1, :].flatten()
    lon[2:] = lons[0 : 2 * (m - 1), 0 : m - 1, :].flatten()
    return lat, lon


def legendre(lat, ntrunc):
    """
    calculate associated legendre functions for triangular truncation T(ntrunc),
    at a given latitude.

    @param lat:  the latitude (in degrees) to compute the associate legendre
    functions.

    @param ntrunc:  the triangular truncation limit.

    @return: C{B{pnm}} - rank 1 numpy float32 array containing the
    (C{B{ntrunc}}+1)*(C{B{ntrunc}}+2)/2 associated legendre functions at
    latitude C{B{lat}}.
    """
    return _spherepack.getlegfunc(lat, ntrunc)


def specintrp(lon, dataspec, legfuncs):
    """
    spectral interpolation given spherical harmonic coefficients.

    @param lon: longitude (in degrees) of point on a sphere to interpolate to.

    @param dataspec:  spectral coefficients of function to interpolate.

    @param legfuncs: associated legendre functions with same triangular
    truncation as C{B{dataspec}} (computed using L{legendre}), computed
    at latitude of interpolation point.

    @return: C{B{ob}} - interpolated value.
    """
    ntrunc1 = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - dataspec.shape[0])))
    ntrunc2 = int(-1.5 + 0.5 * math.sqrt(9.0 - 8.0 * (1.0 - legfuncs.shape[0])))
    if ntrunc1 != ntrunc2:
        raise ValueError(
            "first dimensions of dataspec and legfuncs in Spharmt.specintrp imply inconsistent spectral truncations - they must be the same!"
        )
    return _spherepack.specintrp((math.pi / 180.0) * lon, ntrunc1, dataspec, legfuncs)
