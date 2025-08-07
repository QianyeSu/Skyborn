"""
Spherical harmonic vector wind computations.

This module provides the main VectorWind class for analyzing atmospheric wind fields
using spherical harmonic transforms. It enables calculation of vorticity, divergence,
streamfunction, velocity potential, and other dynamical quantities.
"""

from __future__ import annotations
from typing import Optional, Union, Tuple, Literal
import numpy as np
from skyborn.spharm import Spharmt, gaussian_lats_wts

# Type aliases for better readability
ArrayLike = Union[np.ndarray, np.ma.MaskedArray]
GridType = Literal["regular", "gaussian"]
LegFunc = Literal["stored", "computed"]


class VectorWind:
    """
    Vector wind analysis using spherical harmonic transforms.

    This class provides methods for analyzing atmospheric wind fields through
    spherical harmonic decomposition. It can compute vorticity, divergence,
    streamfunction, velocity potential, and perform Helmholtz decomposition.

    Parameters
    ----------
    u : array_like
        Zonal wind component. Shape should be (nlat, nlon) or (nlat, nlon, nt)
        where nlat is latitude points, nlon is longitude points, and nt is
        number of time steps. Latitude dimension should be north-to-south.
    v : array_like
        Meridional wind component. Must have same shape as u.
    gridtype : {'regular', 'gaussian'}, default 'regular'
        Type of input grid:
        - 'regular': evenly-spaced lat/lon grid
        - 'gaussian': Gaussian latitude grid
    rsphere : float, default 6.3712e6
        Earth radius in meters for spherical harmonic computations.
    legfunc : {'stored', 'computed'}, default 'stored'
        Legendre function computation method:
        - 'stored': precompute and store (faster, more memory)
        - 'computed': compute on-the-fly (slower, less memory)

    Attributes
    ----------
    u : ndarray
        Zonal wind component
    v : ndarray
        Meridional wind component
    gridtype : str
        Grid type used
    s : Spharmt
        Spherical harmonic transform object

    Examples
    --------
    >>> import numpy as np
    >>> from skyborn.windspharm import VectorWind
    >>>
    >>> # Create sample wind field
    >>> nlat, nlon = 73, 144
    >>> u = np.random.randn(nlat, nlon)
    >>> v = np.random.randn(nlat, nlon)
    >>>
    >>> # Initialize VectorWind
    >>> vw = VectorWind(u, v, gridtype='gaussian')
    >>>
    >>> # Calculate dynamical quantities
    >>> vorticity = vw.vorticity()
    >>> divergence = vw.divergence()
    >>> streamfunction = vw.streamfunction()
    >>> velocity_potential = vw.velocitypotential()
    """

    def __init__(
        self,
        u: ArrayLike,
        v: ArrayLike,
        gridtype: GridType = "regular",
        rsphere: float = 6.3712e6,
        legfunc: LegFunc = "stored",
    ) -> None:
        """
        Initialize VectorWind instance with comprehensive data validation.

        This method performs thorough validation of input wind components including
        checks for missing values, infinite values, shape compatibility, and
        dimensional requirements.
        """
        # Store original input for reference
        self._input_u = u
        self._input_v = v

        # Step 1: Handle masked arrays and create copies
        self.u = self._process_input_array(u, "u")
        self.v = self._process_input_array(v, "v")

        # Step 2: Comprehensive data validation
        self._validate_input_data()

        # Step 3: Validate dimensions and shapes
        self._validate_dimensions()

        # Step 4: Initialize spherical harmonic transform
        nlat, nlon = self.u.shape[0], self.u.shape[1]
        self._initialize_spharmt(nlon, nlat, gridtype, rsphere, legfunc)

        # Step 5: Create method aliases for backward compatibility
        self._setup_method_aliases()

    def _process_input_array(self, arr: ArrayLike, name: str) -> np.ndarray:
        """
        Process input array handling masked arrays and creating copies.

        Parameters
        ----------
        arr : array_like
            Input array (u or v wind component)
        name : str
            Name of the array ('u' or 'v') for error messages

        Returns
        -------
        ndarray
            Processed numpy array
        """
        try:
            # Handle masked arrays
            if hasattr(arr, "filled"):
                processed = arr.filled(fill_value=np.nan)
            else:
                processed = np.asarray(arr, dtype=np.float64).copy()
        except (ValueError, TypeError) as e:
            raise ValueError(f"Cannot convert {name} to numpy array: {e}")

        return processed

    def _validate_input_data(self) -> None:
        """
        Comprehensive validation of input wind data.

        Checks for:
        - Missing values (NaN)
        - Infinite values (±inf)
        - Shape compatibility
        - Data type compatibility
        """
        # Check for NaN values
        u_has_nan = np.isnan(self.u).any()
        v_has_nan = np.isnan(self.v).any()

        if u_has_nan or v_has_nan:
            nan_locations = []
            if u_has_nan:
                nan_count_u = np.isnan(self.u).sum()
                nan_locations.append(f"u: {nan_count_u} NaN values")
            if v_has_nan:
                nan_count_v = np.isnan(self.v).sum()
                nan_locations.append(f"v: {nan_count_v} NaN values")

            raise ValueError(
                f"Input wind components contain missing values (NaN). "
                f"Found in {', '.join(nan_locations)}. "
                f"Please ensure all wind data is finite and valid."
            )

        # Check for infinite values
        u_has_inf = np.isinf(self.u).any()
        v_has_inf = np.isinf(self.v).any()

        if u_has_inf or v_has_inf:
            inf_locations = []
            if u_has_inf:
                inf_count_u = np.isinf(self.u).sum()
                inf_locations.append(f"u: {inf_count_u} infinite values")
            if v_has_inf:
                inf_count_v = np.isinf(self.v).sum()
                inf_locations.append(f"v: {inf_count_v} infinite values")

            raise ValueError(
                f"Input wind components contain infinite values. "
                f"Found in {', '.join(inf_locations)}. "
                f"Please ensure all wind data is finite."
            )

        # Check for extremely large values that might cause numerical issues
        u_max = np.abs(self.u).max()
        v_max = np.abs(self.v).max()
        max_reasonable = 1e6  # 1000 km/s should be more than enough for any wind

        if u_max > max_reasonable or v_max > max_reasonable:
            raise ValueError(
                f"Input wind components contain extremely large values "
                f"(max |u|={u_max:.2e}, max |v|={v_max:.2e}). "
                f"Please check units and data validity."
            )

    def _validate_dimensions(self) -> None:
        """
        Validate array dimensions and shapes.
        """
        # Check shape compatibility
        if self.u.shape != self.v.shape:
            raise ValueError(
                f"Wind components must have identical shapes. "
                f"Got u: {self.u.shape}, v: {self.v.shape}"
            )

        # Check dimensionality
        if self.u.ndim not in (2, 3):
            raise ValueError(
                f"Wind components must be 2D or 3D arrays. "
                f"Got {self.u.ndim}D arrays with shape {self.u.shape}. "
                f"Expected shapes: (nlat, nlon) or (nlat, nlon, ntime)"
            )

        # Check minimum grid size
        nlat, nlon = self.u.shape[:2]
        if nlat < 3 or nlon < 4:
            raise ValueError(
                f"Grid too small for spherical harmonic analysis. "
                f"Got (nlat={nlat}, nlon={nlon}). "
                f"Minimum requirements: nlat ≥ 3, nlon ≥ 4"
            )

        # Check for reasonable grid sizes
        if nlat > 10000 or nlon > 10000:
            import warnings

            warnings.warn(
                f"Very large grid detected (nlat={nlat}, nlon={nlon}). "
                f"This may require significant memory and computation time.",
                UserWarning,
            )

    def _initialize_spharmt(
        self, nlon: int, nlat: int, gridtype: str, rsphere: float, legfunc: str
    ) -> None:
        """
        Initialize spherical harmonic transform object with validation.

        Parameters
        ----------
        nlon, nlat : int
            Grid dimensions
        gridtype : str
            Grid type
        rsphere : float
            Sphere radius
        legfunc : str
            Legendre function computation method
        """
        # Validate gridtype
        gridtype_lower = gridtype.lower()
        if gridtype_lower not in ("regular", "gaussian"):
            raise ValueError(
                f"Invalid grid type: '{gridtype}'. " f"Must be 'regular' or 'gaussian'"
            )

        # Validate rsphere
        if not isinstance(rsphere, (int, float)) or rsphere <= 0:
            raise ValueError(
                f"Invalid sphere radius: {rsphere}. "
                f"Must be a positive number (in meters)"
            )

        # Validate legfunc
        if legfunc not in ("stored", "computed"):
            raise ValueError(
                f"Invalid legfunc option: '{legfunc}'. "
                f"Must be 'stored' or 'computed'"
            )

        # Store configuration
        self.gridtype = gridtype_lower
        self.rsphere = rsphere
        self.legfunc = legfunc

        try:
            # Create spherical harmonic transform object
            self.s = Spharmt(
                nlon=nlon,
                nlat=nlat,
                gridtype=self.gridtype,
                rsphere=rsphere,
                legfunc=legfunc,
            )
        except Exception as e:
            # Provide more informative error messages
            if "invalid input dimensions" in str(e).lower():
                raise ValueError(
                    f"Invalid grid dimensions for {gridtype_lower} grid: "
                    f"(nlat={nlat}, nlon={nlon}). "
                    f"Please check grid size requirements for spherical harmonics."
                ) from e
            elif "spharm" in str(e).lower():
                raise ValueError(
                    f"Spherical harmonic initialization failed: {e}. "
                    f"Please ensure spharm module is properly compiled."
                ) from e
            else:
                raise ValueError(
                    f"Failed to initialize spherical harmonic transform: {e}"
                ) from e

    def _setup_method_aliases(self) -> None:
        """Set up method aliases for backward compatibility."""
        self.rotationalcomponent = self.nondivergentcomponent
        self.divergentcomponent = self.irrotationalcomponent

    @property
    def grid_info(self) -> dict:
        """
        Get information about the grid configuration.

        Returns
        -------
        dict
            Dictionary containing grid information
        """
        nlat, nlon = self.u.shape[:2]
        return {
            "gridtype": self.gridtype,
            "nlat": nlat,
            "nlon": nlon,
            "shape": self.u.shape,
            "rsphere": self.rsphere,
            "legfunc": self.legfunc,
            "total_points": nlat * nlon,
        }

    def magnitude(self) -> np.ndarray:
        """
        Calculate wind speed (magnitude of vector wind).

        Returns
        -------
        ndarray
            Wind speed field computed as sqrt(u² + v²)

        Examples
        --------
        >>> wind_speed = vw.magnitude()
        """
        return np.hypot(self.u, self.v)

    def vrtdiv(self, truncation: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate relative vorticity and horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.
            If None, uses the default truncation based on grid resolution.

        Returns
        -------
        vorticity : ndarray
            Relative vorticity field
        divergence : ndarray
            Horizontal divergence field

        See Also
        --------
        vorticity : Calculate only vorticity
        divergence : Calculate only divergence

        Examples
        --------
        >>> vrt, div = vw.vrtdiv()
        >>> vrt_t13, div_t13 = vw.vrtdiv(truncation=13)
        """
        vrtspec, divspec = self.s.getvrtdivspec(self.u, self.v, ntrunc=truncation)
        vrtgrid = self.s.spectogrd(vrtspec)
        divgrid = self.s.spectogrd(divspec)
        return vrtgrid, divgrid

    def vorticity(self, truncation: Optional[int] = None) -> np.ndarray:
        """
        Calculate relative vorticity.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        ndarray
            Relative vorticity field (∇ × v)

        See Also
        --------
        vrtdiv : Calculate both vorticity and divergence
        absolutevorticity : Calculate absolute vorticity

        Examples
        --------
        >>> vrt = vw.vorticity()
        >>> vrt_t13 = vw.vorticity(truncation=13)
        """
        vrtspec, _ = self.s.getvrtdivspec(self.u, self.v, ntrunc=truncation)
        return self.s.spectogrd(vrtspec)

    def divergence(self, truncation: Optional[int] = None) -> np.ndarray:
        """
        Calculate horizontal divergence.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        ndarray
            Horizontal divergence field (∇ · v)

        See Also
        --------
        vrtdiv : Calculate both vorticity and divergence

        Examples
        --------
        >>> div = vw.divergence()
        >>> div_t13 = vw.divergence(truncation=13)
        """
        _, divspec = self.s.getvrtdivspec(self.u, self.v, ntrunc=truncation)
        return self.s.spectogrd(divspec)

    def planetaryvorticity(self, omega: Optional[float] = None) -> np.ndarray:
        """
        Calculate planetary vorticity (Coriolis parameter).

        Parameters
        ----------
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹.

        Returns
        -------
        ndarray
            Planetary vorticity field (2Ω sin φ)

        See Also
        --------
        absolutevorticity : Calculate absolute vorticity

        Examples
        --------
        >>> f = vw.planetaryvorticity()
        >>> f_custom = vw.planetaryvorticity(omega=7.2921150e-5)
        """
        if omega is None:
            omega = 7.292e-05  # Earth's angular velocity

        nlat = self.s.nlat
        if self.gridtype == "gaussian":
            lat, _ = gaussian_lats_wts(nlat)
        else:
            if nlat % 2:
                lat = np.linspace(90, -90, nlat)
            else:
                dlat = 180.0 / nlat
                lat = np.arange(90 - dlat / 2.0, -90, -dlat)

        try:
            coriolis = 2.0 * omega * np.sin(np.deg2rad(lat))
        except (TypeError, ValueError):
            raise ValueError(f"invalid value for omega: {omega!r}")

        # Broadcast to match input shape
        indices = [slice(None)] + [np.newaxis] * (len(self.u.shape) - 1)
        return coriolis[tuple(indices)] * np.ones(self.u.shape, dtype=np.float32)

    def absolutevorticity(
        self, omega: Optional[float] = None, truncation: Optional[int] = None
    ) -> np.ndarray:
        """
        Calculate absolute vorticity (relative + planetary vorticity).

        Parameters
        ----------
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹.
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        ndarray
            Absolute vorticity field (ζ + f)

        See Also
        --------
        vorticity : Calculate relative vorticity
        planetaryvorticity : Calculate planetary vorticity

        Examples
        --------
        >>> abs_vrt = vw.absolutevorticity()
        >>> abs_vrt_t13 = vw.absolutevorticity(omega=7.2921150e-5, truncation=13)
        """
        planetary_vrt = self.planetaryvorticity(omega=omega)
        relative_vrt = self.vorticity(truncation=truncation)
        return planetary_vrt + relative_vrt

    def sfvp(self, truncation: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate streamfunction and velocity potential.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        streamfunction : ndarray
            Streamfunction field (ψ)
        velocity_potential : ndarray
            Velocity potential field (χ)

        See Also
        --------
        streamfunction : Calculate only streamfunction
        velocitypotential : Calculate only velocity potential

        Examples
        --------
        >>> psi, chi = vw.sfvp()
        >>> psi_t13, chi_t13 = vw.sfvp(truncation=13)
        """
        return self.s.getpsichi(self.u, self.v, ntrunc=truncation)

    def streamfunction(self, truncation: Optional[int] = None) -> np.ndarray:
        """
        Calculate streamfunction.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        ndarray
            Streamfunction field (ψ)

        See Also
        --------
        sfvp : Calculate both streamfunction and velocity potential

        Examples
        --------
        >>> psi = vw.streamfunction()
        >>> psi_t13 = vw.streamfunction(truncation=13)
        """
        psi, _ = self.sfvp(truncation=truncation)
        return psi

    def velocitypotential(self, truncation: Optional[int] = None) -> np.ndarray:
        """
        Calculate velocity potential.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        ndarray
            Velocity potential field (χ)

        See Also
        --------
        sfvp : Calculate both streamfunction and velocity potential

        Examples
        --------
        >>> chi = vw.velocitypotential()
        >>> chi_t13 = vw.velocitypotential(truncation=13)
        """
        _, chi = self.sfvp(truncation=truncation)
        return chi

    def helmholtz(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Perform Helmholtz decomposition of vector wind.

        Decomposes the wind field into irrotational (divergent) and
        non-divergent (rotational) components.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        u_chi : ndarray
            Zonal component of irrotational wind
        v_chi : ndarray
            Meridional component of irrotational wind
        u_psi : ndarray
            Zonal component of non-divergent wind
        v_psi : ndarray
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
        psi, chi = self.s.getpsichi(self.u, self.v, ntrunc=truncation)
        psi_spec = self.s.grdtospec(psi)
        chi_spec = self.s.grdtospec(chi)
        v_psi, u_psi = self.s.getgrad(psi_spec)
        u_chi, v_chi = self.s.getgrad(chi_spec)
        return u_chi, v_chi, -u_psi, v_psi

    def irrotationalcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate irrotational (divergent) component of vector wind.

        Note
        ----
        If both irrotational and non-divergent components are needed,
        use `helmholtz()` method for efficiency.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        u_chi : ndarray
            Zonal component of irrotational wind
        v_chi : ndarray
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
        _, chi = self.s.getpsichi(self.u, self.v, ntrunc=truncation)
        chi_spec = self.s.grdtospec(chi)
        return self.s.getgrad(chi_spec)

    def nondivergentcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate non-divergent (rotational) component of vector wind.

        Note
        ----
        If both non-divergent and irrotational components are needed,
        use `helmholtz()` method for efficiency.

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        u_psi : ndarray
            Zonal component of non-divergent wind
        v_psi : ndarray
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
        psi, _ = self.s.getpsichi(self.u, self.v, ntrunc=truncation)
        psi_spec = self.s.grdtospec(psi)
        v_psi, u_psi = self.s.getgrad(psi_spec)
        return -u_psi, v_psi

    def gradient(
        self, chi: ArrayLike, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate vector gradient of a scalar field on the sphere.

        Parameters
        ----------
        chi : array_like
            Scalar field with shape (nlat, nlon) or (nlat, nlon, nfields)
            matching the dimensions of the VectorWind instance.
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.

        Returns
        -------
        u_chi : ndarray
            Zonal component of vector gradient (∂χ/∂λ)
        v_chi : ndarray
            Meridional component of vector gradient (∂χ/∂φ)

        Examples
        --------
        >>> abs_vrt = vw.absolutevorticity()
        >>> avrt_u, avrt_v = vw.gradient(abs_vrt)
        >>> avrt_u_t13, avrt_v_t13 = vw.gradient(abs_vrt, truncation=13)
        """
        try:
            chi = chi.filled(fill_value=np.nan) if hasattr(chi, "filled") else chi
        except AttributeError:
            pass

        if np.isnan(chi).any():
            raise ValueError("chi cannot contain missing values")

        try:
            chi_spec = self.s.grdtospec(chi, ntrunc=truncation)
        except ValueError:
            raise ValueError("input field is not compatible")

        return self.s.getgrad(chi_spec)

    def truncate(
        self, field: ArrayLike, truncation: Optional[int] = None
    ) -> np.ndarray:
        """
        Apply spectral truncation to a scalar field.

        This is useful to represent fields consistently with the output
        of other VectorWind methods.

        Parameters
        ----------
        field : array_like
            Scalar field with shape (nlat, nlon) or (nlat, nlon, nfields)
            matching the dimensions of the VectorWind instance.
        truncation : int, optional
            Triangular truncation limit. If None, defaults to nlat-1.

        Returns
        -------
        ndarray
            Field with spectral truncation applied

        Examples
        --------
        >>> field_trunc = vw.truncate(scalar_field)
        >>> field_t21 = vw.truncate(scalar_field, truncation=21)
        """
        try:
            field = (
                field.filled(fill_value=np.nan) if hasattr(field, "filled") else field
            )
        except AttributeError:
            pass

        if np.isnan(field).any():
            raise ValueError("field cannot contain missing values")

        try:
            field_spec = self.s.grdtospec(field, ntrunc=truncation)
        except ValueError:
            raise ValueError("field is not compatible")

        return self.s.spectogrd(field_spec)

    def rossbywavesource(
        self, truncation: Optional[int] = None, omega: Optional[float] = None
    ) -> np.ndarray:
        """
        Calculate Rossby wave source.

        The Rossby wave source is defined as:
        S = -ζₐ∇·v - v_χ·∇ζₐ

        where:
        - ζₐ is absolute vorticity (relative + planetary)
        - ∇·v is horizontal divergence
        - v_χ is the irrotational (divergent) wind component
        - ∇ζₐ is the gradient of absolute vorticity

        Parameters
        ----------
        truncation : int, optional
            Triangular truncation limit for spherical harmonic computation.
            If None, uses the default truncation based on grid resolution.
        omega : float, optional
            Earth's angular velocity in rad/s. Default is 7.292e-5 s⁻¹.

        Returns
        -------
        ndarray
            Rossby wave source field (s⁻²)

        See Also
        --------
        absolutevorticity : Calculate absolute vorticity
        divergence : Calculate horizontal divergence
        irrotationalcomponent : Calculate irrotational wind component
        gradient : Calculate vector gradient

        Notes
        -----
        The Rossby wave source quantifies the generation of Rossby wave activity
        in the atmosphere. Positive values indicate wave generation, while
        negative values indicate wave absorption or dissipation.

        The calculation involves several steps:
        1. Compute absolute vorticity (relative + planetary)
        2. Compute horizontal divergence
        3. Compute irrotational wind components
        4. Compute gradients of absolute vorticity
        5. Combine terms according to the RWS formula

        Examples
        --------
        >>> rws = vw.rossbywavesource()
        >>> rws_t21 = vw.rossbywavesource(truncation=21)
        >>> rws_custom_omega = vw.rossbywavesource(omega=7.2921150e-5)

        References
        ----------
        Sardeshmukh, P. D., & Hoskins, B. J. (1988). The generation of global
        rotational flow by steady idealized tropical heating. Journal of the
        Atmospheric Sciences, 45(7), 1228-1251.
        """
        # Calculate absolute vorticity
        eta = self.absolutevorticity(omega=omega, truncation=truncation)

        # Calculate divergence
        div = self.divergence(truncation=truncation)

        # Calculate irrotational (divergent) wind components
        uchi, vchi = self.irrotationalcomponent(truncation=truncation)

        # Calculate gradients of absolute vorticity
        etax, etay = self.gradient(eta, truncation=truncation)

        # Combine components to form Rossby wave source
        # S = -eta * div - (uchi * etax + vchi * etay)
        rws = -eta * div - (uchi * etax + vchi * etay)

        return rws
