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

    def vorticity(self, truncation: Optional[int] = None) -> IrisCube:
        """
        Calculate relative vorticity.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        iris.cube.Cube
            Relative vorticity field with CF-compliant attributes

        See Also
        --------
        vrtdiv : Calculate both vorticity and divergence
        absolutevorticity : Calculate absolute vorticity

        Examples
        --------
        >>> vrt = vw.vorticity()
        >>> vrt_t13 = vw.vorticity(truncation=13)
        """
        vrt = self._api.vorticity(truncation=truncation)
        return self._metadata(
            vrt,
            units="s**-1",
            standard_name="atmosphere_relative_vorticity",
            long_name="relative vorticity",
        )

    def divergence(self, truncation: Optional[int] = None) -> IrisCube:
        """
        Calculate horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        iris.cube.Cube
            Horizontal divergence field with CF-compliant attributes

        See Also
        --------
        vrtdiv : Calculate both vorticity and divergence

        Examples
        --------
        >>> div = vw.divergence()
        >>> div_t13 = vw.divergence(truncation=13)
        """
        div = self._api.divergence(truncation=truncation)
        return self._metadata(
            div,
            units="s**-1",
            standard_name="divergence_of_wind",
            long_name="horizontal divergence",
        )

    def planetaryvorticity(self, omega: Optional[float] = None) -> IrisCube:
        """
        Calculate planetary vorticity (Coriolis parameter).

        Parameters
        ----------
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹

        Returns
        -------
        iris.cube.Cube
            Planetary vorticity (Coriolis parameter) with CF-compliant attributes

        See Also
        --------
        absolutevorticity : Calculate absolute vorticity

        Examples
        --------
        >>> f = vw.planetaryvorticity()
        >>> f_custom = vw.planetaryvorticity(omega=7.2921150e-5)
        """
        f = self._api.planetaryvorticity(omega=omega)
        return self._metadata(
            f,
            units="s**-1",
            standard_name="coriolis_parameter",
            long_name="planetary vorticity (coriolis parameter)",
        )

    def absolutevorticity(
        self, omega: Optional[float] = None, truncation: Optional[int] = None
    ) -> IrisCube:
        """
        Calculate absolute vorticity (relative + planetary vorticity).

        Parameters
        ----------
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        iris.cube.Cube
            Absolute vorticity field with CF-compliant attributes

        See Also
        --------
        vorticity : Calculate relative vorticity
        planetaryvorticity : Calculate planetary vorticity

        Examples
        --------
        >>> abs_vrt = vw.absolutevorticity()
        >>> abs_vrt_t13 = vw.absolutevorticity(omega=7.2921150e-5, truncation=13)
        """
        avrt = self._api.absolutevorticity(omega=omega, truncation=truncation)
        return self._metadata(
            avrt,
            units="s**-1",
            standard_name="atmosphere_absolute_vorticity",
            long_name="absolute vorticity",
        )

    def sfvp(self, truncation: Optional[int] = None) -> Tuple[IrisCube, IrisCube]:
        """
        Calculate streamfunction and velocity potential.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        streamfunction : iris.cube.Cube
            Streamfunction field with CF-compliant attributes
        velocity_potential : iris.cube.Cube
            Velocity potential field with CF-compliant attributes

        See Also
        --------
        streamfunction : Calculate only streamfunction
        velocitypotential : Calculate only velocity potential

        Examples
        --------
        >>> psi, chi = vw.sfvp()
        >>> psi_t13, chi_t13 = vw.sfvp(truncation=13)
        """
        sf, vp = self._api.sfvp(truncation=truncation)

        sf_cube = self._metadata(
            sf,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_streamfunction",
            long_name="streamfunction",
        )

        vp_cube = self._metadata(
            vp,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_velocity_potential",
            long_name="velocity potential",
        )

        return sf_cube, vp_cube

    def streamfunction(self, truncation: Optional[int] = None) -> IrisCube:
        """
        Calculate streamfunction.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        iris.cube.Cube
            Streamfunction field with CF-compliant attributes

        See Also
        --------
        sfvp : Calculate both streamfunction and velocity potential

        Examples
        --------
        >>> psi = vw.streamfunction()
        >>> psi_t13 = vw.streamfunction(truncation=13)
        """
        sf = self._api.streamfunction(truncation=truncation)
        return self._metadata(
            sf,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_streamfunction",
            long_name="streamfunction",
        )

    def velocitypotential(self, truncation: Optional[int] = None) -> IrisCube:
        """
        Calculate velocity potential.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        iris.cube.Cube
            Velocity potential field with CF-compliant attributes

        See Also
        --------
        sfvp : Calculate both streamfunction and velocity potential

        Examples
        --------
        >>> chi = vw.velocitypotential()
        >>> chi_t13 = vw.velocitypotential(truncation=13)
        """
        vp = self._api.velocitypotential(truncation=truncation)
        return self._metadata(
            vp,
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_velocity_potential",
            long_name="velocity potential",
        )

    def helmholtz(
        self, truncation: Optional[int] = None
    ) -> Tuple[IrisCube, IrisCube, IrisCube, IrisCube]:
        """
        Perform Helmholtz decomposition of vector wind.

        Decomposes the wind field into irrotational (divergent) and
        non-divergent (rotational) components.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        u_chi : iris.cube.Cube
            Zonal component of irrotational wind
        v_chi : iris.cube.Cube
            Meridional component of irrotational wind
        u_psi : iris.cube.Cube
            Zonal component of non-divergent wind
        v_psi : iris.cube.Cube
            Meridional component of non-divergent wind

        See Also
        --------
        irrotationalcomponent : Get only irrotational component
        nondivergentcomponent : Get only non-divergent component

        Examples
        --------
        >>> u_chi, v_chi, u_psi, v_psi = vw.helmholtz()
        >>> u_chi_t13, v_chi_t13, u_psi_t13, v_psi_t13 = vw.helmholtz(truncation=13)
        """
        uchi, vchi, upsi, vpsi = self._api.helmholtz(truncation=truncation)

        uchi_cube = self._metadata(
            uchi, units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi_cube = self._metadata(
            vchi, units="m s**-1", long_name="irrotational_northward_wind"
        )
        upsi_cube = self._metadata(
            upsi, units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi_cube = self._metadata(
            vpsi, units="m s**-1", long_name="non_divergent_northward_wind"
        )

        return uchi_cube, vchi_cube, upsi_cube, vpsi_cube

    def irrotationalcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[IrisCube, IrisCube]:
        """
        Calculate irrotational (divergent) component of vector wind.

        Note
        ----
        If both irrotational and non-divergent components are needed,
        use `helmholtz()` method for efficiency.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        u_chi : iris.cube.Cube
            Zonal component of irrotational wind
        v_chi : iris.cube.Cube
            Meridional component of irrotational wind

        See Also
        --------
        helmholtz : Complete Helmholtz decomposition
        nondivergentcomponent : Non-divergent component

        Examples
        --------
        >>> u_chi, v_chi = vw.irrotationalcomponent()
        >>> u_chi_t13, v_chi_t13 = vw.irrotationalcomponent(truncation=13)
        """
        uchi, vchi = self._api.irrotationalcomponent(truncation=truncation)

        uchi_cube = self._metadata(
            uchi, units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi_cube = self._metadata(
            vchi, units="m s**-1", long_name="irrotational_northward_wind"
        )

        return uchi_cube, vchi_cube

    def nondivergentcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[IrisCube, IrisCube]:
        """
        Calculate non-divergent (rotational) component of vector wind.

        Note
        ----
        If both non-divergent and irrotational components are needed,
        use `helmholtz()` method for efficiency.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        u_psi : iris.cube.Cube
            Zonal component of non-divergent wind
        v_psi : iris.cube.Cube
            Meridional component of non-divergent wind

        See Also
        --------
        helmholtz : Complete Helmholtz decomposition
        irrotationalcomponent : Irrotational component

        Examples
        --------
        >>> u_psi, v_psi = vw.nondivergentcomponent()
        >>> u_psi_t13, v_psi_t13 = vw.nondivergentcomponent(truncation=13)
        """
        upsi, vpsi = self._api.nondivergentcomponent(truncation=truncation)

        upsi_cube = self._metadata(
            upsi, units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi_cube = self._metadata(
            vpsi, units="m s**-1", long_name="non_divergent_northward_wind"
        )

        return upsi_cube, vpsi_cube

    def gradient(
        self, chi: IrisCube, truncation: Optional[int] = None
    ) -> Tuple[IrisCube, IrisCube]:
        """
        Calculate vector gradient of a scalar field on the sphere.

        Parameters
        ----------
        chi : iris.cube.Cube
            Scalar field with same latitude/longitude dimensions as wind components
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        u_gradient : iris.cube.Cube
            Zonal component of vector gradient
        v_gradient : iris.cube.Cube
            Meridional component of vector gradient

        Examples
        --------
        >>> abs_vrt = vw.absolutevorticity()
        >>> avrt_u, avrt_v = vw.gradient(abs_vrt)
        >>> avrt_u_t13, avrt_v_t13 = vw.gradient(abs_vrt, truncation=13)
        """
        if not isinstance(chi, Cube):
            raise TypeError(
                f"Scalar field must be iris.cube.Cube, got {type(chi).__name__}"
            )

        name = chi.name() or "field"

        # Process coordinate ordering similar to initialization
        lat, lat_dim = _dim_coord_and_dim(chi, "latitude")
        lon, lon_dim = _dim_coord_and_dim(chi, "longitude")

        # Ensure north-to-south latitude ordering
        if lat.points[0] < lat.points[1]:
            chi = reverse(chi, lat_dim)
            lat, lat_dim = _dim_coord_and_dim(chi, "latitude")

        # Reorder for API compatibility
        apiorder, reorder = get_apiorder(chi.ndim, lat_dim, lon_dim)
        chi = chi.copy()
        chi.transpose(apiorder)

        # Store shape and coordinates
        ishape = chi.shape
        coords = chi.dim_coords

        # Compute gradient using standard API
        chi_data = to3d(chi.data)
        uchi, vchi = self._api.gradient(chi_data, truncation=truncation)

        # Reshape and create cubes
        uchi = uchi.reshape(ishape)
        vchi = vchi.reshape(ishape)

        uchi_cube = Cube(uchi, dim_coords_and_dims=list(zip(coords, range(uchi.ndim))))
        vchi_cube = Cube(vchi, dim_coords_and_dims=list(zip(coords, range(vchi.ndim))))

        # Restore original dimension order
        uchi_cube.transpose(reorder)
        vchi_cube.transpose(reorder)

        # Set descriptive names
        uchi_cube.long_name = f"zonal_gradient_of_{name}"
        vchi_cube.long_name = f"meridional_gradient_of_{name}"

        return uchi_cube, vchi_cube

    def truncate(self, field: IrisCube, truncation: Optional[int] = None) -> IrisCube:
        """
        Apply spectral truncation to a scalar field.

        This is useful to represent other fields consistently with the output
        of other VectorWind methods.

        Parameters
        ----------
        field : iris.cube.Cube
            Scalar field with same latitude/longitude dimensions as wind components
        truncation : int, optional
            Triangular truncation limit. If None, defaults to nlat-1

        Returns
        -------
        iris.cube.Cube
            Field with spectral truncation applied

        Examples
        --------
        >>> field_trunc = vw.truncate(scalar_field)
        >>> field_t21 = vw.truncate(scalar_field, truncation=21)
        """
        if not isinstance(field, Cube):
            raise TypeError(f"Field must be iris.cube.Cube, got {type(field).__name__}")

        # Process coordinate ordering
        lat, lat_dim = _dim_coord_and_dim(field, "latitude")
        lon, lon_dim = _dim_coord_and_dim(field, "longitude")

        # Ensure north-to-south latitude ordering
        if lat.points[0] < lat.points[1]:
            field = reverse(field, lat_dim)
            lat, lat_dim = _dim_coord_and_dim(field, "latitude")

        # Reorder for API compatibility
        apiorder, reorder = get_apiorder(field.ndim, lat_dim, lon_dim)
        field = field.copy()
        field.transpose(apiorder)

        # Store shape and apply truncation
        ishape = field.shape
        fielddata = to3d(field.data)
        fieldtrunc = self._api.truncate(fielddata, truncation=truncation)

        # Update field data and restore dimension order
        field.data = fieldtrunc.reshape(ishape)
        field.transpose(reorder)

        return field


def _dim_coord_and_dim(cube: IrisCube, coord: str) -> Tuple[Any, int]:
    """
    Retrieve a dimension coordinate from an Iris cube and its dimension index.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube to search
    coord : str
        Name of coordinate to find (e.g., 'latitude', 'longitude')

    Returns
    -------
    coordinate : iris.coords.DimCoord
        Found dimension coordinate
    dim_index : int
        Dimension index of the coordinate

    Raises
    ------
    ValueError
        If coordinate not found or multiple coordinates found
    """
    # Find coordinates matching the name
    coords = [c for c in cube.dim_coords if coord in c.name()]

    if len(coords) > 1:
        raise ValueError(f"Multiple {coord} coordinates not allowed: {coords}")

    if not coords:
        raise ValueError(f"Cannot find {coord} coordinate in cube {cube}")

    coordinate = coords[0]
    coord_dims = cube.coord_dims(coordinate)

    if len(coord_dims) != 1:
        raise ValueError(
            f"Multiple dimensions with {coord} coordinate not allowed: {coord_dims}"
        )

    dim_index = coord_dims[0]
    return coordinate, dim_index
