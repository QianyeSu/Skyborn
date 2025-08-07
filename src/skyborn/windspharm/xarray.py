"""
Spherical harmonic vector wind computations with xarray interface.

This module provides a VectorWind class that works with xarray DataArrays,
preserving coordinate information and metadata throughout the computation process.
It serves as a high-level interface to the standard VectorWind implementation.

Main Class:
    VectorWind: xarray-aware interface for wind field analysis

Example:
    >>> import xarray as xr
    >>> from skyborn.windspharm.xarray import VectorWind
    >>>
    >>> # Load wind data as xarray DataArrays
    >>> u = xr.open_dataarray('u_wind.nc')
    >>> v = xr.open_dataarray('v_wind.nc')
    >>>
    >>> # Create VectorWind instance
    >>> vw = VectorWind(u, v)
    >>>
    >>> # Compute with preserved metadata
    >>> vorticity = vw.vorticity()
    >>> streamfunction = vw.streamfunction()
"""

from __future__ import annotations
from typing import Optional, Tuple, Union, Callable, Any, List
import warnings

try:
    import xarray as xr
except ImportError:
    try:
        import xray as xr  # Fallback for older versions
    except ImportError:
        raise ImportError(
            "xarray is required for the xarray interface. "
            "Install with: pip install xarray"
        )

from . import standard
from ._common import get_apiorder, inspect_gridtype, to3d

# Type aliases for better readability
DataArray = xr.DataArray
LegFunc = str  # 'stored' or 'computed'


class VectorWind:
    """
    Vector wind analysis using xarray DataArrays.

    This class provides a high-level interface for spherical harmonic wind analysis
    that preserves xarray coordinate information and metadata. It wraps the standard
    VectorWind implementation while maintaining CF-compliant attributes.

    Parameters
    ----------
    u, v : xarray.DataArray
        Zonal and meridional wind components. Must have the same dimensions,
        coordinates, and contain no missing values. Should include latitude
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
    _reorder : tuple
        Original dimension ordering for output reconstruction
    _ishape : tuple
        Original data shape
    _coords : list
        Original coordinate information

    Examples
    --------
    >>> import xarray as xr
    >>> from skyborn.windspharm.xarray import VectorWind
    >>>
    >>> # Load wind components
    >>> u = xr.open_dataarray('u850.nc')
    >>> v = xr.open_dataarray('v850.nc')
    >>>
    >>> # Create VectorWind instance
    >>> vw = VectorWind(u, v)
    >>>
    >>> # Compute vorticity with preserved metadata
    >>> vorticity = vw.vorticity()
    >>> print(vorticity.attrs)  # CF-compliant attributes
    >>>
    >>> # Helmholtz decomposition
    >>> u_chi, v_chi, u_psi, v_psi = vw.helmholtz()
    """

    def __init__(
        self,
        u: DataArray,
        v: DataArray,
        rsphere: float = 6.3712e6,
        legfunc: LegFunc = "stored",
    ) -> None:
        """Initialize VectorWind instance with comprehensive validation."""
        # Validate input types
        if not isinstance(u, xr.DataArray):
            raise TypeError(f"u must be xarray.DataArray, got {type(u).__name__}")
        if not isinstance(v, xr.DataArray):
            raise TypeError(f"v must be xarray.DataArray, got {type(v).__name__}")

        # Validate coordinate compatibility
        self._validate_coordinates(u, v)

        # Find and validate latitude/longitude coordinates
        lat, lat_dim = _find_latitude_coordinate(u)
        lon, lon_dim = _find_longitude_coordinate(u)

        # Ensure north-to-south latitude ordering
        if lat.values[0] < lat.values[1]:
            u = _reverse(u, lat_dim)
            v = _reverse(v, lat_dim)
            lat, lat_dim = _find_latitude_coordinate(u)

        # Determine grid type
        gridtype = inspect_gridtype(lat.values)

        # Prepare data for standard API
        apiorder, _ = get_apiorder(u.ndim, lat_dim, lon_dim)
        apiorder = [u.dims[i] for i in apiorder]

        # Store original structure for output reconstruction
        self._reorder = u.dims

        # Reorder dimensions and prepare data
        u = u.copy().transpose(*apiorder)
        v = v.copy().transpose(*apiorder)

        # Store shape and coordinates for reconstruction
        self._ishape = u.shape
        self._coords = [u.coords[name] for name in u.dims]

        # Convert to 3D and initialize standard API
        u_data = to3d(u.values)
        v_data = to3d(v.values)

        self._api = standard.VectorWind(
            u_data, v_data, gridtype=gridtype, rsphere=rsphere, legfunc=legfunc
        )

    def _validate_coordinates(self, u: DataArray, v: DataArray) -> None:
        """
        Validate that u and v have compatible coordinates.

        Parameters
        ----------
        u, v : DataArray
            Wind components to validate

        Raises
        ------
        ValueError
            If dimensions or coordinate values don't match
        """
        # Check dimension names
        if u.dims != v.dims:
            raise ValueError(
                f"u and v must have identical dimensions. "
                f"Got u: {u.dims}, v: {v.dims}"
            )

        # Check coordinate values
        u_coords = [u.coords[name].values for name in u.dims]
        v_coords = [v.coords[name].values for name in v.dims]

        mismatched_coords = []
        for i, (uc, vc) in enumerate(zip(u_coords, v_coords)):
            try:
                if not (uc == vc).all():
                    mismatched_coords.append(u.dims[i])
            except (ValueError, TypeError):
                # Handle different shapes or types
                mismatched_coords.append(u.dims[i])

        if mismatched_coords:
            raise ValueError(
                f"u and v must have identical coordinate values. "
                f"Mismatched coordinates: {mismatched_coords}"
            )

    def _metadata(self, data: Any, name: str, **attributes: Any) -> DataArray:
        """
        Create DataArray with proper metadata and coordinate information.

        Parameters
        ----------
        data : array_like
            Data to wrap in DataArray
        name : str
            Variable name
        **attributes
            Additional attributes to set

        Returns
        -------
        DataArray
            Properly formatted DataArray with coordinates and metadata
        """
        # Reshape to original structure
        data = data.reshape(self._ishape)

        # Create DataArray with coordinates
        result = xr.DataArray(data, coords=self._coords, name=name)

        # Restore original dimension order
        result = result.transpose(*self._reorder)

        # Set attributes
        for attr, value in attributes.items():
            result.attrs[attr] = value

        return result

    def u(self) -> DataArray:
        """
        Get zonal component of vector wind.

        Returns
        -------
        DataArray
            Zonal (eastward) wind component with CF-compliant attributes

        Examples
        --------
        >>> u_wind = vw.u()
        >>> print(u_wind.attrs['standard_name'])  # 'eastward_wind'
        """
        return self._metadata(
            self._api.u,
            "u",
            units="m s**-1",
            standard_name="eastward_wind",
            long_name="eastward_component_of_wind",
        )

    def v(self) -> DataArray:
        """
        Get meridional component of vector wind.

        Returns
        -------
        DataArray
            Meridional (northward) wind component with CF-compliant attributes

        Examples
        --------
        >>> v_wind = vw.v()
        >>> print(v_wind.attrs['standard_name'])  # 'northward_wind'
        """
        return self._metadata(
            self._api.v,
            "v",
            units="m s**-1",
            standard_name="northward_wind",
            long_name="northward_component_of_wind",
        )

    def magnitude(self) -> DataArray:
        """
        Calculate wind speed (magnitude of vector wind).

        Returns
        -------
        DataArray
            Wind speed with CF-compliant attributes

        Examples
        --------
        >>> wind_speed = vw.magnitude()
        >>> print(wind_speed.attrs['standard_name'])  # 'wind_speed'
        """
        magnitude = self._api.magnitude()
        return self._metadata(
            magnitude,
            "speed",
            units="m s**-1",
            standard_name="wind_speed",
            long_name="wind_speed",
        )

    def vrtdiv(self, truncation: Optional[int] = None) -> Tuple[DataArray, DataArray]:
        """
        Calculate relative vorticity and horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        vorticity : DataArray
            Relative vorticity with CF-compliant attributes
        divergence : DataArray
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

        vrt_da = self._metadata(
            vrt,
            "vorticity",
            units="s**-1",
            standard_name="atmosphere_relative_vorticity",
            long_name="relative_vorticity",
        )

        div_da = self._metadata(
            div,
            "divergence",
            units="s**-1",
            standard_name="divergence_of_wind",
            long_name="horizontal_divergence",
        )

        return vrt_da, div_da

    def vorticity(self, truncation: Optional[int] = None) -> DataArray:
        """
        Calculate relative vorticity.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        DataArray
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
            "vorticity",
            units="s**-1",
            standard_name="atmosphere_relative_vorticity",
            long_name="relative_vorticity",
        )

    def divergence(self, truncation: Optional[int] = None) -> DataArray:
        """
        Calculate horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        DataArray
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
            "divergence",
            units="s**-1",
            standard_name="divergence_of_wind",
            long_name="horizontal_divergence",
        )

    def planetaryvorticity(self, omega: Optional[float] = None) -> DataArray:
        """
        Calculate planetary vorticity (Coriolis parameter).

        Parameters
        ----------
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹

        Returns
        -------
        DataArray
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
            "coriolis",
            units="s**-1",
            standard_name="coriolis_parameter",
            long_name="planetary_vorticity",
        )

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
            "absolute_vorticity",
            units="s**-1",
            standard_name="atmosphere_absolute_vorticity",
            long_name="absolute_vorticity",
        )
        return avrt

    def sfvp(self, truncation: Optional[int] = None) -> Tuple[DataArray, DataArray]:
        """
        Calculate streamfunction and velocity potential.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        streamfunction : DataArray
            Streamfunction field with CF-compliant attributes
        velocity_potential : DataArray
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

        sf_da = self._metadata(
            sf,
            "streamfunction",
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_streamfunction",
            long_name="streamfunction",
        )

        vp_da = self._metadata(
            vp,
            "velocity_potential",
            units="m**2 s**-1",
            standard_name="atmosphere_horizontal_velocity_potential",
            long_name="velocity potential",
        )

        return sf_da, vp_da

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
            "streamfunction",
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
            "velocity_potential",
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
            uchi, "u_chi", units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi = self._metadata(
            vchi, "v_chi", units="m s**-1", long_name="irrotational_northward_wind"
        )
        upsi = self._metadata(
            upsi, "u_psi", units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi = self._metadata(
            vpsi, "v_psi", units="m s**-1", long_name="non_divergent_northward_wind"
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
            uchi, "u_chi", units="m s**-1", long_name="irrotational_eastward_wind"
        )
        vchi = self._metadata(
            vchi, "v_chi", units="m s**-1", long_name="irrotational_northward_wind"
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
            upsi, "u_psi", units="m s**-1", long_name="non_divergent_eastward_wind"
        )
        vpsi = self._metadata(
            vpsi, "v_psi", units="m s**-1", long_name="non_divergent_northward_wind"
        )
        return upsi, vpsi

    def gradient(
        self, chi: DataArray, truncation: Optional[int] = None
    ) -> Tuple[DataArray, DataArray]:
        """
        Calculate vector gradient of a scalar field on the sphere.

        Parameters
        ----------
        chi : DataArray
            Scalar field with same latitude/longitude dimensions as wind components
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation

        Returns
        -------
        u_gradient : DataArray
            Zonal component of vector gradient
        v_gradient : DataArray
            Meridional component of vector gradient

        Examples
        --------
        >>> abs_vrt = vw.absolutevorticity()
        >>> avrt_u, avrt_v = vw.gradient(abs_vrt)
        >>> avrt_u_t13, avrt_v_t13 = vw.gradient(abs_vrt, truncation=13)
        """
        if not isinstance(chi, xr.DataArray):
            raise TypeError(
                f"Scalar field must be xarray.DataArray, got {type(chi).__name__}"
            )

        name = chi.name or "field"

        # Process coordinate ordering similar to initialization
        lat, lat_dim = _find_latitude_coordinate(chi)
        lon, lon_dim = _find_longitude_coordinate(chi)

        # Ensure north-to-south latitude ordering
        if lat.values[0] < lat.values[1]:
            chi = _reverse(chi, lat_dim)
            lat, lat_dim = _find_latitude_coordinate(chi)

        # Reorder for API compatibility
        apiorder, _ = get_apiorder(chi.ndim, lat_dim, lon_dim)
        apiorder = [chi.dims[i] for i in apiorder]
        reorder = chi.dims

        chi = chi.copy().transpose(*apiorder)
        ishape = chi.shape
        coords = [chi.coords[n] for n in chi.dims]

        # Compute gradient using standard API
        chi_data = to3d(chi.values)
        u_grad, v_grad = self._api.gradient(chi_data, truncation=truncation)

        # Reshape and create DataArrays
        u_grad = u_grad.reshape(ishape)
        v_grad = v_grad.reshape(ishape)

        u_name = f"zonal_gradient_of_{name}"
        v_name = f"meridional_gradient_of_{name}"

        u_da = xr.DataArray(
            u_grad, coords=coords, name=u_name, attrs={"long_name": u_name}
        )
        v_da = xr.DataArray(
            v_grad, coords=coords, name=v_name, attrs={"long_name": v_name}
        )

        # Restore original dimension order
        u_da = u_da.transpose(*reorder)
        v_da = v_da.transpose(*reorder)

        return u_da, v_da

    def rossbywavesource(
        self, truncation: Optional[int] = None, omega: Optional[float] = None
    ) -> DataArray:
        """
        Calculate Rossby wave source.

        The Rossby wave source quantifies the generation of Rossby wave activity
        in the atmosphere through the interaction of divergence with absolute
        vorticity and the advection of absolute vorticity by the irrotational wind.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.
            If None, uses the default truncation based on grid resolution.
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹.

        Returns
        -------
        DataArray
            Rossby wave source field with CF-compliant attributes

        See Also
        --------
        absolutevorticity : Calculate absolute vorticity
        divergence : Calculate horizontal divergence
        irrotationalcomponent : Calculate irrotational wind component
        gradient : Calculate vector gradient

        Notes
        -----
        The Rossby wave source is defined as:
        S = -ζₐ∇·v - v_χ·∇ζₐ

        where:
        - ζₐ is absolute vorticity (relative + planetary)
        - ∇·v is horizontal divergence
        - v_χ is the irrotational (divergent) wind component
        - ∇ζₐ is the gradient of absolute vorticity

        Positive values indicate Rossby wave generation, while negative values
        indicate wave absorption or dissipation.

        Examples
        --------
        >>> rws = vw.rossbywavesource()
        >>> rws_t21 = vw.rossbywavesource(truncation=21)
        >>> rws_custom_omega = vw.rossbywavesource(omega=7.2921150e-5)

        # Create a plot of Rossby wave source
        >>> import matplotlib.pyplot as plt
        >>> import cartopy.crs as ccrs
        >>>
        >>> ax = plt.axes(projection=ccrs.PlateCarree())
        >>> rws.plot.contourf(ax=ax, transform=ccrs.PlateCarree(),
        ...                   levels=20, cmap='RdBu_r')
        >>> ax.coastlines()
        >>> ax.gridlines()
        >>> plt.title('Rossby Wave Source')
        >>> plt.show()

        References
        ----------
        Sardeshmukh, P. D., & Hoskins, B. J. (1988). The generation of global
        rotational flow by steady idealized tropical heating. Journal of the
        Atmospheric Sciences, 45(7), 1228-1251.
        """
        rws = self._api.rossbywavesource(truncation=truncation, omega=omega)
        return self._metadata(
            rws,
            "rossby_wave_source",
            units="s**-2",
            standard_name="rossby_wave_source",
            long_name="rossby_wave_source_term",
            description="Generation term for Rossby wave activity",
        )

    def truncate(self, field: DataArray, truncation: Optional[int] = None) -> DataArray:
        """
        Apply spectral truncation to a scalar field.

        Parameters
        ----------
        field : DataArray
            Scalar field with same latitude/longitude dimensions as wind components
        truncation : int, optional
            Triangular truncation limit. If None, defaults to nlat-1

        Returns
        -------
        DataArray
            Field with spectral truncation applied

        Examples
        --------
        >>> field_trunc = vw.truncate(scalar_field)
        >>> field_t21 = vw.truncate(scalar_field, truncation=21)
        """
        if not isinstance(field, xr.DataArray):
            raise TypeError(
                f"Field must be xarray.DataArray, got {type(field).__name__}"
            )

        # Process coordinate ordering
        lat, lat_dim = _find_latitude_coordinate(field)
        lon, lon_dim = _find_longitude_coordinate(field)

        # Ensure north-to-south latitude ordering
        if lat.values[0] < lat.values[1]:
            field = _reverse(field, lat_dim)
            lat, lat_dim = _find_latitude_coordinate(field)

        # Reorder for API compatibility
        apiorder, _ = get_apiorder(field.ndim, lat_dim, lon_dim)
        apiorder = [field.dims[i] for i in apiorder]
        reorder = field.dims

        field = field.copy().transpose(*apiorder)
        ishape = field.shape

        # Apply truncation using standard API
        field_data = to3d(field.values)
        field_trunc = self._api.truncate(field_data, truncation=truncation)

        # Update field values and restore dimension order
        field.values = field_trunc.reshape(ishape)
        field = field.transpose(*reorder)

        return field


def _reverse(array: DataArray, dim: int) -> DataArray:
    """
    Reverse an xarray DataArray along a given dimension.

    Parameters
    ----------
    array : DataArray
        Input DataArray to reverse
    dim : int
        Dimension index to reverse

    Returns
    -------
    DataArray
        Array with specified dimension reversed
    """
    slicers = [slice(None)] * array.ndim
    slicers[dim] = slice(None, None, -1)
    return array[tuple(slicers)]


def _find_coord_and_dim(
    array: DataArray, predicate: Callable[[Any], bool], name: str
) -> Tuple[Any, int]:
    """
    Find a dimension coordinate in DataArray that satisfies a predicate.

    Parameters
    ----------
    array : DataArray
        Input DataArray to search
    predicate : callable
        Function that returns True for the desired coordinate
    name : str
        Name of coordinate type for error messages

    Returns
    -------
    coord : coordinate
        Found coordinate that satisfies predicate
    dim : int
        Dimension index of the coordinate

    Raises
    ------
    ValueError
        If no coordinate or multiple coordinates found
    """
    candidates = [
        coord for coord in [array.coords[n] for n in array.dims] if predicate(coord)
    ]

    if not candidates:
        raise ValueError(f"Cannot find a {name} coordinate")

    if len(candidates) > 1:
        raise ValueError(f"Multiple {name} coordinates are not allowed")

    coord = candidates[0]
    dim = array.dims.index(coord.name)
    return coord, dim


def _find_latitude_coordinate(array: DataArray) -> Tuple[Any, int]:
    """
    Find latitude dimension coordinate in an xarray DataArray.

    Parameters
    ----------
    array : DataArray
        Input DataArray to search for latitude coordinate

    Returns
    -------
    lat_coord : coordinate
        Latitude coordinate
    lat_dim : int
        Latitude dimension index

    Raises
    ------
    ValueError
        If latitude coordinate cannot be found or multiple found
    """

    def is_latitude(coord: Any) -> bool:
        """Check if coordinate represents latitude."""
        return (
            coord.name in ("latitude", "lat")
            or coord.attrs.get("units") == "degrees_north"
            or coord.attrs.get("axis") == "Y"
        )

    return _find_coord_and_dim(array, is_latitude, "latitude")


def _find_longitude_coordinate(array: DataArray) -> Tuple[Any, int]:
    """
    Find longitude dimension coordinate in an xarray DataArray.

    Parameters
    ----------
    array : DataArray
        Input DataArray to search for longitude coordinate

    Returns
    -------
    lon_coord : coordinate
        Longitude coordinate
    lon_dim : int
        Longitude dimension index

    Raises
    ------
    ValueError
        If longitude coordinate cannot be found or multiple found
    """

    def is_longitude(coord: Any) -> bool:
        """Check if coordinate represents longitude."""
        return (
            coord.name in ("longitude", "lon")
            or coord.attrs.get("units") == "degrees_east"
            or coord.attrs.get("axis") == "X"
        )

    return _find_coord_and_dim(array, is_longitude, "longitude")
