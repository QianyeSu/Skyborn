"""
Spherical harmonic vector wind computations.

This module provides the main VectorWind class for analyzing atmospheric wind fields
using spherical harmonic transforms. It enables calculation of vorticity, divergence,
streamfunction, velocity potential, and other dynamical quantities.
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple, Union

import numpy as np

from skyborn.spharm import Spharmt, gaussian_lats_wts

__all__ = ["VectorWind"]

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
        Zonal wind component. Shape should be (nlat, nlon) or
        (nlat, nlon, *extra_dims), where nlat is latitude points and nlon is
        longitude points. Latitude dimension should be north-to-south.
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
        # Step 1: Handle masked arrays and create copies
        self._output_dtype = self._infer_output_dtype(u, v)
        self.u = self._process_input_array(u, "u")
        self.v = self._process_input_array(v, "v")

        # Step 2: Comprehensive data validation
        self._validate_input_data()

        # Step 3: Validate dimensions and shapes
        self._validate_dimensions()

        # Step 4: Initialize spherical harmonic transform
        nlat, nlon = self.u.shape[0], self.u.shape[1]
        self._initialize_spharmt(nlon, nlat, gridtype, rsphere, legfunc)
        self._setup_cached_arrays()

        # Step 5: Create method aliases for backward compatibility
        self._setup_method_aliases()

    @classmethod
    def _from_prepared(
        cls,
        u: ArrayLike,
        v: ArrayLike,
        gridtype: GridType = "regular",
        rsphere: float = 6.3712e6,
        legfunc: LegFunc = "stored",
    ) -> "VectorWind":
        """Create an instance from validated wrapper-owned arrays."""
        self = cls.__new__(cls)
        self._output_dtype = self._infer_output_dtype(u, v)
        self.u = self._process_input_array(u, "u", copy=False)
        self.v = self._process_input_array(v, "v", copy=False)
        self._validate_input_data()
        self._validate_dimensions()
        nlat, nlon = self.u.shape[0], self.u.shape[1]
        self._initialize_spharmt(nlon, nlat, gridtype, rsphere, legfunc)
        self._setup_cached_arrays()
        self._setup_method_aliases()
        return self

    # ------------------------------------------------------------------
    # Input processing, dtype handling, and validation
    # ------------------------------------------------------------------

    @staticmethod
    def _infer_output_dtype(*arrays: ArrayLike) -> np.dtype:
        """Return the public output dtype for floating inputs."""
        output_dtype = None
        for arr in arrays:
            dtype = np.dtype(np.asarray(arr).dtype)
            if np.issubdtype(dtype, np.floating):
                public_dtype = dtype
            else:
                public_dtype = np.dtype(np.float64)
            output_dtype = (
                public_dtype
                if output_dtype is None
                else np.result_type(output_dtype, public_dtype)
            )
        return np.dtype(np.float64 if output_dtype is None else output_dtype)

    def _restore_output_dtype(
        self, data: np.ndarray, output_dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        """Cast grid-space public outputs back to the input precision."""
        array = np.asarray(data)
        target_dtype = (
            self._output_dtype if output_dtype is None else np.dtype(output_dtype)
        )
        if array.dtype == target_dtype:
            return array
        return array.astype(target_dtype, copy=False)

    def _process_input_array(
        self, arr: ArrayLike, name: str, copy: bool = True
    ) -> np.ndarray:
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
                masked = np.ma.asarray(arr)
                dtype = masked.dtype
                if not np.issubdtype(dtype, np.floating):
                    dtype = np.dtype(np.float64)
                processed = np.ma.asarray(masked, dtype=dtype).filled(fill_value=np.nan)
                if copy:
                    processed = processed.copy()
            else:
                raw = np.asarray(arr)
                dtype = np.dtype(raw.dtype)
                if not np.issubdtype(dtype, np.floating):
                    dtype = np.dtype(np.float64)
                processed = np.asarray(raw, dtype=dtype)
                if copy:
                    processed = processed.copy()
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
        if not np.isfinite(self.u).all() or not np.isfinite(self.v).all():
            nan_locations, inf_locations = self._invalid_value_locations()

            if nan_locations:
                raise ValueError(
                    f"Input wind components contain missing values (NaN). "
                    f"Found in {', '.join(nan_locations)}. "
                    f"Please ensure all wind data is finite and valid."
                )

            raise ValueError(
                f"Input wind components contain infinite values. "
                f"Found in {', '.join(inf_locations)}. "
                f"Please ensure all wind data is finite."
            )

        # Check for extremely large values that might cause numerical issues
        # u_max = np.abs(self.u).max()
        # v_max = np.abs(self.v).max()
        # max_reasonable = 1e6  # 1000 km/s should be more than enough for any wind

        # if u_max > max_reasonable or v_max > max_reasonable:
        #     raise ValueError(
        #         f"Input wind components contain extremely large values "
        #         f"(max |u|={u_max:.2e}, max |v|={v_max:.2e}). "
        #         f"Please check units and data validity."
        #     )

    def _invalid_value_locations(self) -> Tuple[list[str], list[str]]:
        """Count NaN and infinite values after the fast finite check fails."""
        nan_locations: list[str] = []
        inf_locations: list[str] = []

        for name, values in (("u", self.u), ("v", self.v)):
            invalid = ~np.isfinite(values)
            invalid_count = int(invalid.sum())
            if not invalid_count:
                continue

            invalid_values = values[invalid]
            nan_count = int(np.isnan(invalid_values).sum())
            inf_count = invalid_count - nan_count
            if nan_count:
                nan_locations.append(f"{name}: {nan_count} NaN values")
            if inf_count:
                inf_locations.append(f"{name}: {inf_count} infinite values")

        return nan_locations, inf_locations

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
        if self.u.ndim < 2:
            raise ValueError(
                f"Wind components must be at least 2D arrays. "
                f"Got {self.u.ndim}D arrays with shape {self.u.shape}. "
                f"Expected shape: (nlat, nlon) or (nlat, nlon, *extra_dims)"
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
        # if nlat > 10000 or nlon > 10000:
        #     import warnings

        #     warnings.warn(
        #         f"Very large grid detected (nlat={nlat}, nlon={nlon}). "
        #         f"This may require significant memory and computation time.",
        #         UserWarning,
        #     )

    # ------------------------------------------------------------------
    # Transform setup and small per-instance caches
    # ------------------------------------------------------------------

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

    def _setup_cached_arrays(self) -> None:
        """Initialize small per-instance caches for repeated wrapper work."""
        self._latitude_cache: dict[tuple[str, int], np.ndarray] = {}
        self._coriolis_cache: dict[tuple[str, int, float, str], np.ndarray] = {}

    def _latitude_values(self) -> np.ndarray:
        """Return latitude coordinates for this wind grid."""
        key = (self.gridtype, self.s.nlat)
        cache = getattr(self, "_latitude_cache", None)
        if cache is None:
            cache = {}
            self._latitude_cache = cache

        lat = cache.get(key)
        if lat is not None:
            return lat

        nlat = self.s.nlat
        if self.gridtype == "gaussian":
            lat, _ = gaussian_lats_wts(nlat)
        elif nlat % 2:
            lat = np.linspace(90, -90, nlat)
        else:
            dlat = 180.0 / nlat
            lat = np.arange(90 - dlat / 2.0, -90, -dlat)

        cache[key] = np.asarray(lat)
        return cache[key]

    def _coriolis_values(
        self, omega: Optional[float], dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        """Return cached latitude-only Coriolis values."""
        if omega is None:
            omega = 7.292e-05

        try:
            omega_value = float(omega)
        except (TypeError, ValueError):
            raise ValueError(f"invalid value for omega: {omega!r}")

        output_dtype = self._output_dtype if dtype is None else np.dtype(dtype)
        key = (self.gridtype, self.s.nlat, omega_value, output_dtype.str)
        cache = getattr(self, "_coriolis_cache", None)
        if cache is None:
            cache = {}
            self._coriolis_cache = cache

        coriolis = cache.get(key)
        if coriolis is not None:
            return coriolis

        coriolis = 2.0 * omega_value * np.sin(np.deg2rad(self._latitude_values()))
        coriolis = coriolis.astype(output_dtype, copy=False)
        cache[key] = coriolis
        return coriolis

    # ------------------------------------------------------------------
    # Raw numerical helpers. These avoid public dtype restoration between
    # intermediate steps in chained diagnostics.
    # ------------------------------------------------------------------

    def _vrtdiv_raw(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return relative vorticity and divergence before public dtype restore."""
        vrtspec, divspec = self.s.getvrtdivspec(self.u, self.v, ntrunc=truncation)
        return self.s.spectogrd(vrtspec), self.s.spectogrd(divspec)

    def _planetaryvorticity_raw(
        self,
        omega: Optional[float] = None,
        materialize: bool = True,
        dtype: Optional[np.dtype] = None,
    ) -> np.ndarray:
        """Return planetary vorticity before public dtype restore."""
        coriolis = self._coriolis_values(omega, dtype=dtype)
        return self._broadcast_planetaryvorticity(coriolis, materialize=materialize)

    def _broadcast_planetaryvorticity(
        self, coriolis: np.ndarray, materialize: bool
    ) -> np.ndarray:
        """Broadcast latitude-only Coriolis values to the wind field shape."""
        indices = [slice(None)] + [np.newaxis] * (len(self.u.shape) - 1)
        broadcast = np.broadcast_to(coriolis[tuple(indices)], self.u.shape)
        if materialize:
            return np.array(broadcast, copy=True)
        return broadcast

    def _irrotationalcomponent_raw(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return irrotational wind components before public dtype restore."""
        chispec = self.s._getchi_spec(self.u, self.v, ntrunc=truncation)
        return self.s.getgrad(chispec)

    def _gradient_raw(
        self, chi: ArrayLike, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return vector gradient before public dtype restore."""
        try:
            chi = chi.filled(fill_value=np.nan) if hasattr(chi, "filled") else chi
        except AttributeError:
            pass

        if np.any(np.isnan(chi)):
            raise ValueError("chi cannot contain missing values")

        try:
            chi_spec = self.s.grdtospec(chi, ntrunc=truncation)
        except ValueError:
            raise ValueError("input field is not compatible")

        return self.s.getgrad(chi_spec)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

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
        return self._restore_output_dtype(np.hypot(self.u, self.v))

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
        vrt, div = self._vrtdiv_raw(truncation=truncation)
        return (
            self._restore_output_dtype(vrt),
            self._restore_output_dtype(div),
        )

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
        return self._restore_output_dtype(self.s.spectogrd(vrtspec))

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
        return self._restore_output_dtype(self.s.spectogrd(divspec))

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
        return self._planetaryvorticity_raw(
            omega=omega, materialize=True, dtype=self._output_dtype
        )

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
        vrt, _ = self._vrtdiv_raw(truncation=truncation)
        return self._restore_output_dtype(
            vrt + self._planetaryvorticity_raw(omega=omega, materialize=False)
        )

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
        psi, chi = self.s.getpsichi(self.u, self.v, ntrunc=truncation)
        return self._restore_output_dtype(psi), self._restore_output_dtype(chi)

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
        psi = self.s.getpsi(self.u, self.v, ntrunc=truncation)
        return self._restore_output_dtype(psi)

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
        chi = self.s.getchi(self.u, self.v, ntrunc=truncation)
        return self._restore_output_dtype(chi)

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
        psispec, chispec = self.s._getpsichi_spec(self.u, self.v, ntrunc=truncation)
        v_psi, u_psi = self.s.getgrad(psispec)
        u_chi, v_chi = self.s.getgrad(chispec)
        return (
            self._restore_output_dtype(u_chi),
            self._restore_output_dtype(v_chi),
            self._restore_output_dtype(-u_psi),
            self._restore_output_dtype(v_psi),
        )

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
        u_chi, v_chi = self._irrotationalcomponent_raw(truncation=truncation)
        return self._restore_output_dtype(u_chi), self._restore_output_dtype(v_chi)

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
        psispec = self.s._getpsi_spec(self.u, self.v, ntrunc=truncation)
        v_psi, u_psi = self.s.getgrad(psispec)
        return self._restore_output_dtype(-u_psi), self._restore_output_dtype(v_psi)

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
        output_dtype = self._infer_output_dtype(chi)
        u_chi, v_chi = self._gradient_raw(chi, truncation=truncation)
        return (
            self._restore_output_dtype(u_chi, output_dtype),
            self._restore_output_dtype(v_chi, output_dtype),
        )

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

        if np.any(np.isnan(field)):
            raise ValueError("field cannot contain missing values")

        output_dtype = self._infer_output_dtype(field)

        try:
            field_spec = self.s.grdtospec(field, ntrunc=truncation)
        except ValueError:
            raise ValueError("field is not compatible")

        return self._restore_output_dtype(self.s.spectogrd(field_spec), output_dtype)

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
        # Compute relative vorticity and divergence in a single vector analysis.
        vrt, div = self._vrtdiv_raw(truncation=truncation)
        eta = vrt + self._planetaryvorticity_raw(omega=omega, materialize=False)

        # Calculate irrotational (divergent) wind components
        uchi, vchi = self._irrotationalcomponent_raw(truncation=truncation)

        # Calculate gradients of absolute vorticity
        etax, etay = self._gradient_raw(eta, truncation=truncation)

        # Combine components to form Rossby wave source
        # S = -eta * div - (uchi * etax + vchi * etay)
        rws = -eta * div - (uchi * etax + vchi * etay)

        return self._restore_output_dtype(rws)
