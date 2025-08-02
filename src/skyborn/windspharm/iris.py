"""
Spherical harmonic vector wind computations with Iris interface.

This module provides a VectorWind class that works with Iris Cubes,
preserving coordinate information and metadata throughout the computation process.
It serves as a high-level interface to the standard VectorWind implementation.

Main Class:
    VectorWind: Iris-aware interface for wind field analysis

Example:
    >>> import iris
    >>> from skyborn.windspharm.iris import VectorWind
    >>>
    >>> # Load wind data as Iris cubes
    >>> u = iris.load_cube('u_wind.nc')
    >>> v = iris.load_cube('v_wind.nc')
    >>>
    >>> # Create VectorWind instance
    >>> vw = VectorWind(u, v)
    >>>
    >>> # Compute with preserved metadata
    >>> vorticity = vw.vorticity()
    >>> streamfunction = vw.streamfunction()
"""

from __future__ import annotations
from typing import Optional, Tuple, Union, Any

try:
    from iris.cube import Cube
    from iris.util import reverse
except ImportError:
    raise ImportError(
        "Iris is required for the iris interface. "
        "Install with: conda install -c conda-forge iris"
    )

from . import standard
from ._common import get_apiorder, inspect_gridtype, to3d

# Type aliases for better readability
IrisCube = Cube
LegFunc = str  # 'stored' or 'computed'
GridType = str  # 'regular' or 'gaussian'


class VectorWind:
    """
    Vector wind analysis using Iris cubes.

    This class provides a high-level interface for spherical harmonic wind analysis
    that preserves Iris coordinate information and metadata. It wraps the standard
    VectorWind implementation while maintaining CF-compliant attributes.

    Parameters
    ----------
    u, v : iris.cube.Cube
        Zonal and meridional wind components. Must have the same dimension
        coordinates and contain no missing values. Should include latitude
        and longitude dimensions with appropriate coordinate information.
    rsphere : float, default 6.3712e6
        Earth radius in meters for spherical harmonic computations.
    legfunc : {'stored', 'computed'}, default 'stored'
        Legendre function computation method:
        - 'stored': precompute and store (faster, more memory)
        - 'computed': compute on-the-fly (slower, less memory)

    Attributes
    ----------
    _api : standard.VectorWind
        Underlying standard VectorWind instance
    _reorder : list
        Dimension reordering for output reconstruction
    _ishape : tuple
        Original data shape
    _coords : list
        Original coordinate information

    Examples
    --------
    >>> import iris
    >>> from skyborn.windspharm.iris import VectorWind
    >>>
    >>> # Load wind components
    >>> u = iris.load_cube('u850.nc')
    >>> v = iris.load_cube('v850.nc')
    >>>
    >>> # Create VectorWind instance
    >>> vw = VectorWind(u, v)
    >>>
    >>> # Compute vorticity with preserved metadata
    >>> vorticity = vw.vorticity()
    >>> print(vorticity.attributes)  # CF-compliant attributes
    >>>
    >>> # Helmholtz decomposition
    >>> u_chi, v_chi, u_psi, v_psi = vw.helmholtz()
    """

    def __init__(
        self,
        u: IrisCube,
        v: IrisCube,
        rsphere: float = 6.3712e6,
        legfunc: LegFunc = "stored",
    ) -> None:
        """
        Initialize VectorWind instance with comprehensive validation.

        This method performs thorough validation of input wind components including
        checks for cube compatibility, coordinate consistency, and proper formatting.

        Parameters
        ----------
        u, v : iris.cube.Cube
            Zonal and meridional wind components with matching coordinates
        rsphere : float, default 6.3712e6
            Earth radius in meters
        legfunc : {'stored', 'computed'}, default 'stored'
            Legendre function computation method

        Raises
        ------
        TypeError
            If u or v are not Iris cubes
        ValueError
            If u and v don't have matching dimensions or coordinates
        """
        # Validate input types
        if not isinstance(u, Cube):
            raise TypeError(f"u must be iris.cube.Cube, got {type(u).__name__}")
        if not isinstance(v, Cube):
            raise TypeError(f"v must be iris.cube.Cube, got {type(v).__name__}")

        # Validate coordinate compatibility
        self._validate_cube_compatibility(u, v)

        # Extract and validate latitude/longitude coordinates
        lat, lat_dim = _dim_coord_and_dim(u, "latitude")
        lon, lon_dim = _dim_coord_and_dim(v, "longitude")

        # Ensure north-to-south latitude ordering
        if lat.points[0] < lat.points[1]:
            u = reverse(u, lat_dim)
            v = reverse(v, lat_dim)
            lat, lat_dim = _dim_coord_and_dim(u, "latitude")

        # Determine grid type
        gridtype = inspect_gridtype(lat.points)

        # Calculate dimension ordering for API compatibility
        apiorder, self._reorder = get_apiorder(u.ndim, lat_dim, lon_dim)

        # Prepare cubes for processing (make copies to avoid modifying originals)
        u = u.copy()
        v = v.copy()
        u.transpose(apiorder)
        v.transpose(apiorder)

        # Store original structure for output reconstruction
        self._ishape = u.shape
        self._coords = u.dim_coords

        # Convert to 3D format and initialize standard API
        u_data = to3d(u.data)
        v_data = to3d(v.data)

        self._api = standard.VectorWind(
            u_data, v_data, gridtype=gridtype, rsphere=rsphere, legfunc=legfunc
        )

    def _validate_cube_compatibility(self, u: IrisCube, v: IrisCube) -> None:
        """
        Validate that u and v cubes have compatible coordinates.

        Parameters
        ----------
        u, v : iris.cube.Cube
            Wind components to validate

        Raises
        ------
        ValueError
            If cubes don't have matching dimension coordinates
        """
        if u.dim_coords != v.dim_coords:
            raise ValueError(
                f"u and v must have identical dimension coordinates. "
                f"u coords: {[c.name() for c in u.dim_coords]}, "
                f"v coords: {[c.name() for c in v.dim_coords]}"
            )

    def _metadata(self, var: Any, **attributes: Any) -> IrisCube:
        """
        Create Iris cube with proper metadata and coordinate information.

        Parameters
        ----------
        var : array_like
            Data to wrap in Iris cube
        **attributes
            Additional attributes to set on the cube

        Returns
        -------
        iris.cube.Cube
            Properly formatted cube with coordinates and metadata
        """
        # Reshape to original structure
        var = var.reshape(self._ishape)

        # Create cube with coordinates
        cube = Cube(var, dim_coords_and_dims=list(zip(self._coords, range(var.ndim))))

        # Restore original dimension order
        cube.transpose(self._reorder)

        # Set attributes
        for attribute, value in attributes.items():
            setattr(cube, attribute, value)

        return cube

    def u(self) -> IrisCube:
        """
        Get zonal component of vector wind.

        Returns
        -------
        iris.cube.Cube
            Zonal (eastward) wind component with CF-compliant attributes

        Examples
        --------
        >>> u_wind = vw.u()
        >>> print(u_wind.standard_name)  # 'eastward_wind'
        """
        u = self._api.u
        return self._metadata(
            u,
            standard_name="eastward_wind",
            units="m s**-1",
            long_name="eastward component of wind",
        )

    def v(self) -> IrisCube:
        """
        Get meridional component of vector wind.

        Returns
        -------
        iris.cube.Cube
            Meridional (northward) wind component with CF-compliant attributes

        Examples
        --------
        >>> v_wind = vw.v()
        >>> print(v_wind.standard_name)  # 'northward_wind'
        """
        v = self._api.v
        return self._metadata(
            v,
            standard_name="northward_wind",
            units="m s**-1",
            long_name="northward component of wind",
        )

    def magnitude(self) -> IrisCube:
        """
        Calculate wind speed (magnitude of vector wind).

        Returns
        -------
        iris.cube.Cube
            Wind speed with CF-compliant attributes

        Examples
        --------
        >>> wind_speed = vw.magnitude()
        >>> print(wind_speed.standard_name)  # 'wind_speed'
        """
        m = self._api.magnitude()
        return self._metadata(
            m, standard_name="wind_speed", units="m s**-1", long_name="wind speed"
        )

    def vrtdiv(self, truncation: Optional[int] = None) -> Tuple[IrisCube, IrisCube]:
        """
        Calculate relative vorticity and horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        vorticity : iris.cube.Cube
            Relative vorticity with CF-compliant attributes
        divergence : iris.cube.Cube
            Horizontal divergence with CF-compliant attributes

        See Also
        --------
        vorticity : Calculate only vorticity
        divergence : Calculate only divergence

        Examples
        --------
        >>> vrt, div = vw.vrtdiv()
        >>> vrt_t13, div_t13 = vw.vrtdiv(truncation=13)
        """
        vrt, div = self._api.vrtdiv(truncation=truncation)

        vrt_cube = self._metadata(
            vrt,
            units="s**-1",
            standard_name="atmosphere_relative_vorticity",
            long_name="relative vorticity",
        )

        div_cube = self._metadata(
            div,
            units="s**-1",
            standard_name="divergence_of_wind",
            long_name="horizontal divergence",
        )

        return vrt_cube, div_cube

    def vorticity(self, truncation=None):
        """Relative vorticity.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *vrt*
            The relative vorticity.

        **See also:**

        `~VectorWind.vrtdiv`, `~VectorWind.absolutevorticity`.

        **Examples:**

        Compute the relative vorticity::

            vrt = w.vorticity()

        Compute the relative vorticity and apply spectral truncation at
        triangular T13::

            vrtT13 = w.vorticity(truncation=13)

        """
        vrt = self._api.vorticity(truncation=truncation)
        vrt = self._metadata(
            vrt,
            units="s**-1",
            standard_name="atmosphere_relative_vorticity",
            long_name="relative vorticity",
        )
        return vrt

    def divergence(self, truncation=None):
        """Horizontal divergence.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *div*
            The divergence.

        **See also:**

        `~VectorWind.vrtdiv`.

        **Examples:**

        Compute the divergence::

            div = w.divergence()

        Compute the divergence and apply spectral truncation at
        triangular T13::

            divT13 = w.divergence(truncation=13)

        """
        div = self._api.divergence(truncation=truncation)
        div = self._metadata(
            div,
            units="s**-1",
            standard_name="divergence_of_wind",
            long_name="horizontal divergence",
        )
        return div

    def planetaryvorticity(self, omega=None):
        """Planetary vorticity (Coriolis parameter).

        **Optional argument:**

        *omega*
            Earth's angular velocity. The default value if not specified
            is 7.292x10**-5 s**-1.

        **Returns:**

        *pvorticity*
            The planetary vorticity.

        **See also:**

        `~VectorWind.absolutevorticity`.

        **Example:**

        Compute planetary vorticity using default values::

            pvrt = w.planetaryvorticity()

        Override the default value for Earth's angular velocity::

            pvrt = w.planetaryvorticity(omega=7.2921150)

        """
        f = self._api.planetaryvorticity(omega=omega)
        f = self._metadata(
            f,
            units="s**-1",
            standard_name="coriolis_parameter",
            long_name="planetary vorticity (coriolis parameter)",
        )
        return f

    def absolutevorticity(self, omega=None, truncation=None):
        """Absolute vorticity (sum of relative and planetary vorticity).

        **Optional arguments:**

        *omega*
            Earth's angular velocity. The default value if not specified
            is 7.292x10**-5 s**-1.

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *avorticity*
            The absolute (relative + planetary) vorticity.

        **See also:**

        `~VectorWind.vorticity`, `~VectorWind.planetaryvorticity`.

        **Examples:**

        Compute absolute vorticity::

            avrt = w.absolutevorticity()

        Compute absolute vorticity and apply spectral truncation at
        triangular T13, also override the default value for Earth's
        angular velocity::

            avrt = w.absolutevorticity(omega=7.2921150, truncation=13)

        """
        avrt = self._api.absolutevorticity(omega=omega, truncation=truncation)
        avrt = self._metadata(
            avrt,
            units="s**-1",
            standard_name="atmosphere_absolute_vorticity",
            long_name="absolute vorticity",
        )
        return avrt

    def sfvp(self, truncation=None):
        """Streamfunction and velocity potential.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *sf*, *vp*
            The streamfunction and velocity potential respectively.

        **See also:**

        `~VectorWind.streamfunction`, `~VectorWind.velocitypotential`.

        **Examples:**

        Compute streamfunction and velocity potential::

            sf, vp = w.sfvp()

        Compute streamfunction and velocity potential and apply spectral
        truncation at triangular T13::

            sfT13, vpT13 = w.sfvp(truncation=13)

        """
        sf, vp = self._api.sfvp(truncation=truncation)
        sf = self._metadata(
            sf,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_streamfunction",
            long_name="streamfunction",
        )
        vp = self._metadata(
            vp,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_velocity_potential",
            long_name="velocity potential",
        )
        return sf, vp

    def streamfunction(self, truncation=None):
        """Streamfunction.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *sf*
            The streamfunction.

        **See also:**

        `~VectorWind.sfvp`.

        **Examples:**

        Compute streamfunction::

            sf = w.streamfunction()

        Compute streamfunction and apply spectral truncation at
        triangular T13::

            sfT13 = w.streamfunction(truncation=13)

        """
        sf = self._api.streamfunction(truncation=truncation)
        sf = self._metadata(
            sf,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_streamfunction",
            long_name="streamfunction",
        )
        return sf

    def velocitypotential(self, truncation=None):
        """Velocity potential.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *vp*
            The velocity potential.

        **See also:**

        `~VectorWind.sfvp`.

        **Examples:**

        Compute velocity potential::

            vp = w.velocity potential()

        Compute velocity potential and apply spectral truncation at
        triangular T13::

            vpT13 = w.velocity potential(truncation=13)

        """
        vp = self._api.velocitypotential(truncation=truncation)
        vp = self._metadata(
            vp,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_velocity_potential",
            long_name="velocity potential",
        )
        return vp

    def helmholtz(self, truncation=None):
        """Irrotational and non-divergent components of the vector wind.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *uchi*, *vchi*, *upsi*, *vpsi*
            Zonal and meridional components of irrotational and
            non-divergent wind components respectively.

        **See also:**

        `~VectorWind.irrotationalcomponent`,
        `~VectorWind.nondivergentcomponent`.

        **Examples:**

        Compute the irrotational and non-divergent components of the
        vector wind::

            uchi, vchi, upsi, vpsi = w.helmholtz()

        Compute the irrotational and non-divergent components of the
        vector wind and apply spectral truncation at triangular T13::

            uchiT13, vchiT13, upsiT13, vpsiT13 = w.helmholtz(truncation=13)

        """
        uchi, vchi, upsi, vpsi = self._api.helmholtz(truncation=truncation)
        uchi = self._metadata(
            uchi, units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi = self._metadata(
            vchi, units="m s**-1", long_name="irrotational_northward_wind"
        )
        upsi = self._metadata(
            upsi, units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi = self._metadata(
            vpsi, units="m s**-1", long_name="non_divergent_northward_wind"
        )
        return uchi, vchi, upsi, vpsi

    def irrotationalcomponent(self, truncation=None):
        """Irrotational (divergent) component of the vector wind.

        .. note::

           If both the irrotational and non-divergent components are
           required then `~VectorWind.helmholtz` should be used instead.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *uchi*, *vchi*
            The zonal and meridional components of the irrotational wind
            respectively.

        **See also:**

        `~VectorWind.helmholtz`.

        **Examples:**

        Compute the irrotational component of the vector wind::

            uchi, vchi = w.irrotationalcomponent()

        Compute the irrotational component of the vector wind and apply
        spectral truncation at triangular T13::

            uchiT13, vchiT13 = w.irrotationalcomponent(truncation=13)

        """
        uchi, vchi = self._api.irrotationalcomponent(truncation=truncation)
        uchi = self._metadata(
            uchi, units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi = self._metadata(
            vchi, units="m s**-1", long_name="irrotational_northward_wind"
        )
        return uchi, vchi

    def nondivergentcomponent(self, truncation=None):
        """Non-divergent (rotational) component of the vector wind.

        .. note::

           If both the non-divergent and irrotational components are
           required then `~VectorWind.helmholtz` should be used instead.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *upsi*, *vpsi*
            The zonal and meridional components of the non-divergent
            wind respectively.

        **See also:**

        `~VectorWind.helmholtz`.

        **Examples:**

        Compute the non-divergent component of the vector wind::

            upsi, vpsi = w.nondivergentcomponent()

        Compute the non-divergent component of the vector wind and apply
        spectral truncation at triangular T13::

            upsiT13, vpsiT13 = w.nondivergentcomponent(truncation=13)

        """
        upsi, vpsi = self._api.nondivergentcomponent(truncation=truncation)
        upsi = self._metadata(
            upsi, units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi = self._metadata(
            vpsi, units="m s**-1", long_name="non_divergent_northward_wind"
        )
        return upsi, vpsi

    def gradient(self, chi, truncation=None):
        """Computes the vector gradient of a scalar field on the sphere.

        **Argument:**

        *chi*
            A scalar field. It must be a `~iris.cube.Cube`
            with the same latitude and longitude dimensions as the
            vector wind components that initialized the `VectorWind`
            instance.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation.

        **Returns:**

        *uchi*, *vchi*
            The zonal and meridional components of the vector gradient
            respectively.

        **Examples:**

        Compute the vector gradient of absolute vorticity::

            avrt = w.absolutevorticity()
            avrt_zonal, avrt_meridional = w.gradient(avrt)

        Compute the vector gradient of absolute vorticity and apply
        spectral truncation at triangular T13::

            avrt = w.absolutevorticity()
            avrt_zonalT13, avrt_meridionalT13 = w.gradient(avrt, truncation=13)

        """
        if type(chi) is not Cube:
            raise TypeError("scalar field must be an iris cube")
        name = chi.name()
        lat, lat_dim = _dim_coord_and_dim(chi, "latitude")
        lon, lon_dim = _dim_coord_and_dim(chi, "longitude")
        if lat.points[0] < lat.points[1]:
            # need to reverse latitude dimension
            chi = reverse(chi, lat_dim)
            lat, lat_dim = _dim_coord_and_dim(chi, "latitude")
        apiorder, reorder = get_apiorder(chi.ndim, lat_dim, lon_dim)
        chi = chi.copy()
        chi.transpose(apiorder)
        ishape = chi.shape
        coords = chi.dim_coords
        chi = to3d(chi.data)
        uchi, vchi = self._api.gradient(chi, truncation=truncation)
        uchi = uchi.reshape(ishape)
        vchi = vchi.reshape(ishape)
        uchi = Cube(uchi, dim_coords_and_dims=list(zip(coords, range(uchi.ndim))))
        vchi = Cube(vchi, dim_coords_and_dims=list(zip(coords, range(vchi.ndim))))
        uchi.transpose(reorder)
        vchi.transpose(reorder)
        uchi.long_name = "zonal_gradient_of_{!s}".format(name)
        vchi.long_name = "meridional_gradient_of_{!s}".format(name)
        return uchi, vchi

    def truncate(self, field, truncation=None):
        """Apply spectral truncation to a scalar field.

        This is useful to represent other fields in a way consistent
        with the output of other `VectorWind` methods.

        **Argument:**

        *field*
            A scalar field. It must be a `~iris.cube.Cube`
            with the same latitude and longitude dimensions as the
            vector wind components that initialized the `VectorWind`
            instance.

        **Optional argument:**

        *truncation*
            Truncation limit (triangular truncation) for the spherical
            harmonic computation. If not specified it will default to
            *nlats - 1* where *nlats* is the number of latitudes.

        **Returns:**

        *truncated_field*
            The field with spectral truncation applied.

        **Examples:**

        Truncate a scalar field to the computational resolution of the
        `VectorWind` instance::

            scalar_field_truncated = w.truncate(scalar_field)

        Truncate a scalar field to T21::

            scalar_field_T21 = w.truncate(scalar_field, truncation=21)

        """
        if type(field) is not Cube:
            raise TypeError("scalar field must be an iris cube")
        lat, lat_dim = _dim_coord_and_dim(field, "latitude")
        lon, lon_dim = _dim_coord_and_dim(field, "longitude")
        if lat.points[0] < lat.points[1]:
            # need to reverse latitude dimension
            field = reverse(field, lat_dim)
            lat, lat_dim = _dim_coord_and_dim(field, "latitude")
        apiorder, reorder = get_apiorder(field.ndim, lat_dim, lon_dim)
        field = field.copy()
        field.transpose(apiorder)
        ishape = field.shape
        fielddata = to3d(field.data)
        fieldtrunc = self._api.truncate(fielddata, truncation=truncation)
        field.data = fieldtrunc.reshape(ishape)
        field.transpose(reorder)
        return field


def _dim_coord_and_dim(cube, coord):
    """
    Retrieve a given dimension coordinate from a `~iris.cube.Cube` and
    the dimension number it corresponds to.

    """
    coords = [c for c in cube.dim_coords if coord in c.name()]
    if len(coords) > 1:
        raise ValueError(
            "multiple {!s} coordinates not " "allowed: {!r}".format(coord, cube)
        )
    try:
        c = coords[0]
    except IndexError:
        raise ValueError(
            "cannot get {!s} coordinate " "from cube {!r}".format(coord, cube)
        )
    c_dim = cube.coord_dims(c)
    if len(c_dim) != 1:
        raise ValueError(
            "multiple dimensions with {!s} coordinate "
            "not allowed: {!r}".format(coord, cube)
        )
    c_dim = c_dim[0]
    return c, c_dim
