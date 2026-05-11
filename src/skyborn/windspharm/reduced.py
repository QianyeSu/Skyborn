"""Reduced-Gaussian vector wind analysis for packed fields."""

from __future__ import annotations

from typing import Dict, Optional, Tuple, Union

import numpy as np

from skyborn.spharm import ReducedGaussianSpharmt, gaussian_lats_wts

ArrayLike = Union[np.ndarray, np.ma.MaskedArray]


class ReducedVectorWind:
    """
    Vector wind analysis on packed reduced Gaussian grids.

    Parameters
    ----------
    u, v:
        Packed zonal and meridional wind components with shape
        ``(sum(pl), ...)``.
    pl:
        Reduced Gaussian row lengths, ordered north-to-south.
    rsphere:
        Sphere radius in meters.
    legfunc:
        Legendre function mode passed to grouped Gaussian ``Spharmt`` backends.
    precision:
        Public output precision mode. ``"auto"`` preserves the promoted input
        floating precision, ``"single"`` returns float32 outputs, and
        ``"double"`` returns float64 outputs.

    Attributes
    ----------
    u, v:
        Packed zonal and meridional wind components.
    pl:
        Reduced Gaussian row lengths ordered north-to-south.
    gridtype:
        Fixed grid type label ``"reduced_gaussian"``.
    rsphere:
        Sphere radius in meters.
    legfunc:
        Legendre function handling mode used by the underlying transform.
    precision:
        Public output precision mode.
    s:
        Underlying :class:`skyborn.spharm.ReducedGaussianSpharmt` instance.
    """

    def __init__(
        self,
        u: ArrayLike,
        v: ArrayLike,
        pl: ArrayLike,
        rsphere: float = 6.3712e6,
        legfunc: str = "stored",
        precision: str = "auto",
    ) -> None:
        """
        Create a reduced-Gaussian vector wind analysis wrapper.

        Args:
            u: Packed zonal wind component with shape ``(sum(pl), ...)``.
            v: Packed meridional wind component with the same shape as ``u``.
            pl: Reduced Gaussian row lengths ordered north-to-south.
            rsphere: Sphere radius in meters.
            legfunc: Legendre function handling mode - ``"stored"`` or
                ``"computed"``.
            precision: Public output precision mode. ``"auto"`` preserves the
                promoted input floating precision, ``"single"`` returns
                float32 outputs, and ``"double"`` returns float64 outputs.

        Raises:
            ValueError: If ``u``, ``v``, or ``pl`` fails validation.
        """
        self.pl = self._validate_pl(pl)
        self._output_dtype = self._infer_output_dtype(u, v)
        self._precision = ReducedGaussianSpharmt._validate_precision(precision)
        self.u = self._process_input_array(u, "u")
        self.v = self._process_input_array(v, "v")
        self._validate_input_data()
        self._extra_shape = tuple(self.u.shape[1:])
        self._nt = (
            int(np.prod(self._extra_shape, dtype=int)) if self._extra_shape else 1
        )
        self._u_backend = self._normalize_grid_backend(self.u)
        self._v_backend = self._normalize_grid_backend(self.v)

        self.s = ReducedGaussianSpharmt(
            self.pl,
            rsphere=rsphere,
            legfunc=legfunc,
            precision=precision,
        )
        self.gridtype = self.s.gridtype
        self.rsphere = self.s.rsphere
        self.legfunc = self.s.legfunc
        self._latitude_cache: Optional[np.ndarray] = None
        self._coriolis_cache: Dict[Tuple[float, str], np.ndarray] = {}
        self._planetary_vorticity_backend_cache: Dict[float, np.ndarray] = {}
        self._planetary_vorticity_spec_cache: Dict[Tuple[float, int], np.ndarray] = {}
        self._spectral_cache: Dict[
            Tuple[str, int], Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]
        ] = {}

    @staticmethod
    def _validate_pl(pl: ArrayLike) -> np.ndarray:
        arr = np.asarray(pl, dtype=np.int32)
        if arr.ndim != 1:
            raise ValueError(f"pl must be rank 1, got shape {arr.shape}")
        if arr.size < 3:
            raise ValueError(f"pl must contain at least 3 rows, got {arr.size}")
        if np.any(arr < 4):
            raise ValueError("each reduced Gaussian row must contain >= 4 points")
        return np.ascontiguousarray(arr)

    @staticmethod
    def _infer_output_dtype(*arrays: ArrayLike) -> np.dtype:
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

    def _resolved_output_dtype(
        self, output_dtype: Optional[np.dtype] = None
    ) -> np.dtype:
        if self._precision == "double":
            return np.dtype(np.float64)
        if self._precision == "single":
            return np.dtype(np.float32)
        return self._output_dtype if output_dtype is None else np.dtype(output_dtype)

    def _restore_output_dtype(
        self, data: np.ndarray, output_dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        array = np.asarray(data)
        target_dtype = self._resolved_output_dtype(output_dtype)
        if array.dtype == target_dtype:
            return array
        return array.astype(target_dtype, copy=False)

    def _normalize_grid_backend(self, data: np.ndarray) -> np.ndarray:
        return np.asfortranarray(
            np.asarray(data).reshape(self.u.shape[0], self._nt), dtype=np.float32
        )

    def _restore_grid(
        self, datagrid: np.ndarray, output_dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        return self.s._restore_grid_shape(
            datagrid, self._extra_shape, self._resolved_output_dtype(output_dtype)
        )

    @staticmethod
    def _process_input_array(arr: ArrayLike, name: str) -> np.ndarray:
        try:
            if hasattr(arr, "filled"):
                masked = np.ma.asarray(arr)
                dtype = masked.dtype
                if not np.issubdtype(dtype, np.floating):
                    dtype = np.dtype(np.float64)
                return np.ma.asarray(masked, dtype=dtype).filled(np.nan).copy()

            raw = np.asarray(arr)
            dtype = np.dtype(raw.dtype)
            if not np.issubdtype(dtype, np.floating):
                dtype = np.dtype(np.float64)
            return np.asarray(raw, dtype=dtype).copy()
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Cannot convert {name} to numpy array: {exc}") from exc

    def _validate_input_data(self) -> None:
        if self.u.shape != self.v.shape:
            raise ValueError(
                f"Wind components must have identical shapes. "
                f"Got u: {self.u.shape}, v: {self.v.shape}"
            )
        if self.u.ndim < 1:
            raise ValueError(
                f"Packed wind components must be at least 1D, got {self.u.ndim}D"
            )
        expected = int(self.pl.sum())
        if self.u.shape[0] != expected:
            raise ValueError(
                f"Packed first dimension must be sum(pl)={expected}, "
                f"got {self.u.shape[0]}"
            )
        if not np.isfinite(self.u).all() or not np.isfinite(self.v).all():
            raise ValueError("Input wind components must be finite")

    def _validate_scalar_field(self, field: ArrayLike, name: str) -> np.ndarray:
        try:
            field = field.filled(np.nan) if hasattr(field, "filled") else field
        except AttributeError:
            pass
        array = np.asarray(field)
        if array.shape != self.u.shape:
            raise ValueError(f"{name} is not compatible")
        if not np.isfinite(array).all():
            raise ValueError(f"{name} cannot contain missing or infinite values")
        return array

    @property
    def grid_info(self) -> dict:
        return {
            "gridtype": self.gridtype,
            "nlat": int(self.pl.size),
            "npoints": int(self.pl.sum()),
            "shape": self.u.shape,
            "rsphere": self.rsphere,
            "legfunc": self.legfunc,
            "pl": self.pl.copy(),
        }

    def _latitude_values(self) -> np.ndarray:
        if self._latitude_cache is None:
            lat, _ = gaussian_lats_wts(int(self.pl.size))
            self._latitude_cache = np.repeat(lat, self.pl)
        return self._latitude_cache

    def _coriolis_values(
        self, omega: Optional[float], dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        omega_value = 7.292e-05 if omega is None else float(omega)
        target_dtype = self._u_backend.dtype if dtype is None else np.dtype(dtype)
        key = (omega_value, target_dtype.str)
        cached = self._coriolis_cache.get(key)
        if cached is not None:
            return cached

        coriolis = np.asarray(
            2.0 * omega_value * np.sin(np.deg2rad(self._latitude_values())),
            dtype=target_dtype,
        )
        self._coriolis_cache[key] = coriolis
        return coriolis

    def _planetary_vorticity(
        self, omega: Optional[float] = None, materialize: bool = True
    ) -> np.ndarray:
        coriolis = self._coriolis_values(
            omega, dtype=self._output_dtype if materialize else self._u_backend.dtype
        )
        indices = [slice(None)] + [np.newaxis] * (self.u.ndim - 1)
        broadcast = np.broadcast_to(coriolis[tuple(indices)], self.u.shape)
        if materialize:
            return np.array(broadcast, copy=True)
        return broadcast

    def _planetary_vorticity_backend(self, omega: Optional[float] = None) -> np.ndarray:
        omega_value = 7.292e-05 if omega is None else float(omega)
        cached = self._planetary_vorticity_backend_cache.get(omega_value)
        if cached is not None:
            return cached

        coriolis = self._coriolis_values(omega_value, dtype=self._u_backend.dtype)
        broadcast = np.broadcast_to(
            coriolis.reshape(self.u.shape[0], 1), self._u_backend.shape
        )
        self._planetary_vorticity_backend_cache[omega_value] = broadcast
        return broadcast

    def _planetary_vorticity_spectral(
        self, truncation: Optional[int] = None, omega: Optional[float] = None
    ) -> np.ndarray:
        ntrunc = self._ntrunc(truncation)
        omega_value = 7.292e-05 if omega is None else float(omega)
        key = (omega_value, ntrunc)
        cached = self._planetary_vorticity_spec_cache.get(key)
        if cached is not None:
            return cached

        coriolis = self._coriolis_values(omega_value, dtype=self._u_backend.dtype)
        grid = np.broadcast_to(
            coriolis.reshape(self.u.shape[0], 1), self._u_backend.shape
        )
        f_spec = self.s._analyze_scalar(
            np.asfortranarray(grid, dtype=np.float32), ntrunc
        )
        self._planetary_vorticity_spec_cache[key] = f_spec
        return f_spec

    def _ntrunc(self, truncation: Optional[int]) -> int:
        return self.s._validate_ntrunc(truncation)

    def _spectogrd_pair(
        self,
        spec_a: np.ndarray,
        spec_b: np.ndarray,
        truncation: Optional[int] = None,
        output_dtype: Optional[np.dtype] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        ntrunc = self._ntrunc(truncation)
        grid_a, grid_b = self.s._synthesize_scalar_pair(spec_a, spec_b, ntrunc)
        return self._restore_grid(grid_a, output_dtype), self._restore_grid(
            grid_b, output_dtype
        )

    def _spectogrd(
        self,
        spec: np.ndarray,
        truncation: Optional[int] = None,
        output_dtype: Optional[np.dtype] = None,
    ) -> np.ndarray:
        ntrunc = self._ntrunc(truncation)
        return self._restore_grid(self.s._synthesize_scalar(spec, ntrunc), output_dtype)

    def _vector_analysis_spectral(
        self, truncation: Optional[int]
    ) -> Tuple[np.ndarray, np.ndarray]:
        ntrunc = self._ntrunc(truncation)
        key = ("vrtdiv", ntrunc)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        vrtspec, divspec = self.s._analyze_wind(
            self._u_backend, self._v_backend, ntrunc
        )
        self._spectral_cache[key] = (vrtspec, divspec)
        self._spectral_cache[("vrt", ntrunc)] = vrtspec
        self._spectral_cache[("div", ntrunc)] = divspec
        return vrtspec, divspec

    def _vorticity_spectral(self, truncation: Optional[int]) -> np.ndarray:
        ntrunc = self._ntrunc(truncation)
        key = ("vrt", ntrunc)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        full = self._spectral_cache.get(("vrtdiv", ntrunc))
        if full is not None:
            vrtspec = full[0]  # type: ignore[index]
        else:
            vrtspec = self.s._analyze_vorticity(
                self._u_backend, self._v_backend, ntrunc
            )
        self._spectral_cache[key] = vrtspec
        return vrtspec

    def _divergence_spectral(self, truncation: Optional[int]) -> np.ndarray:
        ntrunc = self._ntrunc(truncation)
        key = ("div", ntrunc)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        full = self._spectral_cache.get(("vrtdiv", ntrunc))
        if full is not None:
            divspec = full[1]  # type: ignore[index]
        else:
            divspec = self.s._analyze_divergence(
                self._u_backend, self._v_backend, ntrunc
            )
        self._spectral_cache[key] = divspec
        return divspec

    def magnitude(self) -> np.ndarray:
        """Return wind speed."""
        return self._restore_output_dtype(np.hypot(self.u, self.v))

    def vrtdiv(self, truncation: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return relative vorticity and horizontal divergence."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation)
        return self._spectogrd_pair(vrtspec, divspec, truncation)

    def vorticity(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return relative vorticity."""
        vrtspec = self._vorticity_spectral(truncation)
        return self._spectogrd(vrtspec, truncation)

    def divergence(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return horizontal divergence."""
        divspec = self._divergence_spectral(truncation)
        return self._spectogrd(divspec, truncation)

    def planetaryvorticity(self, omega: Optional[float] = None) -> np.ndarray:
        """Return planetary vorticity on the packed reduced grid."""
        return self._restore_output_dtype(self._planetary_vorticity(omega=omega))

    def absolutevorticity(
        self, omega: Optional[float] = None, truncation: Optional[int] = None
    ) -> np.ndarray:
        """Return absolute vorticity."""
        vrt = self.vorticity(truncation=truncation)
        return self._restore_output_dtype(
            vrt + self._planetary_vorticity(omega=omega, materialize=False)
        )

    def sfvp(self, truncation: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return streamfunction and velocity potential."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation)
        return self._spectogrd_pair(
            self.s._invlap(vrtspec),
            self.s._invlap(divspec),
            truncation,
        )

    def streamfunction(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return streamfunction."""
        vrtspec = self._vorticity_spectral(truncation)
        return self._spectogrd(self.s._invlap(vrtspec), truncation)

    def velocitypotential(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return velocity potential."""
        divspec = self._divergence_spectral(truncation)
        return self._spectogrd(self.s._invlap(divspec), truncation)

    def helmholtz(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return irrotational and non-divergent wind components."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation)
        psispec = self.s._invlap(vrtspec)
        chispec = self.s._invlap(divspec)
        u_chi, v_chi, v_psi, u_psi = self.s._synthesize_gradient_pair(
            chispec, psispec, truncation
        )
        return (
            self._restore_grid(u_chi),
            self._restore_grid(v_chi),
            self._restore_grid(-u_psi),
            self._restore_grid(v_psi),
        )

    def irrotationalcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return irrotational wind component."""
        divspec = self._divergence_spectral(truncation)
        u_chi, v_chi = self.s._synthesize_gradient(self.s._invlap(divspec), truncation)
        return self._restore_grid(u_chi), self._restore_grid(v_chi)

    def nondivergentcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return non-divergent wind component."""
        vrtspec = self._vorticity_spectral(truncation)
        psispec = self.s._invlap(vrtspec)
        v_psi, u_psi = self.s._synthesize_gradient(psispec, truncation)
        return self._restore_grid(-u_psi), self._restore_grid(v_psi)

    def gradient(
        self, chi: ArrayLike, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return vector gradient of a packed scalar field."""
        output_dtype = self._infer_output_dtype(chi)
        field = self._validate_scalar_field(chi, "chi")
        chi_spec = self.s._analyze_scalar(
            self._normalize_grid_backend(field), truncation
        )
        u_chi, v_chi = self.s._synthesize_gradient(chi_spec, truncation)
        return (
            self._restore_grid(u_chi, output_dtype),
            self._restore_grid(v_chi, output_dtype),
        )

    def truncate(
        self, field: ArrayLike, truncation: Optional[int] = None
    ) -> np.ndarray:
        """Apply spectral truncation to a packed scalar field."""
        output_dtype = self._infer_output_dtype(field)
        scalar = self._validate_scalar_field(field, "field")
        field_spec = self.s._analyze_scalar(
            self._normalize_grid_backend(scalar), truncation
        )
        return self._spectogrd(field_spec, truncation, output_dtype)

    def rossbywavesource(
        self, truncation: Optional[int] = None, omega: Optional[float] = None
    ) -> np.ndarray:
        """Return Rossby wave source on the packed reduced grid."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation)
        vrt, div = self.s._synthesize_scalar_pair(vrtspec, divspec, truncation)
        eta = vrt + self._planetary_vorticity_backend(omega=omega)
        uchi, vchi = self.s._synthesize_gradient(self.s._invlap(divspec), truncation)
        etaspec = vrtspec + self._planetary_vorticity_spectral(truncation, omega=omega)
        etax, etay = self.s._synthesize_gradient(etaspec, truncation)
        rws = -eta * div - (uchi * etax + vchi * etay)
        return self._restore_grid(rws)
