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
        self.pl = self._validate_pl(pl)
        self._output_dtype = self._infer_output_dtype(u, v)
        self._precision = precision
        self.u = self._process_input_array(u, "u")
        self.v = self._process_input_array(v, "v")
        self._validate_input_data()

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

    def _restore_output_dtype(
        self, data: np.ndarray, output_dtype: Optional[np.dtype] = None
    ) -> np.ndarray:
        array = np.asarray(data)
        if self._precision == "double":
            target_dtype = np.float64
        elif self._precision == "single":
            target_dtype = np.float32
        else:
            target_dtype = (
                self._output_dtype if output_dtype is None else np.dtype(output_dtype)
            )
        if array.dtype == target_dtype:
            return array
        return array.astype(target_dtype, copy=False)

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

    def _planetary_vorticity(
        self, omega: Optional[float] = None, materialize: bool = True
    ) -> np.ndarray:
        omega_value = 7.292e-05 if omega is None else float(omega)
        coriolis = 2.0 * omega_value * np.sin(np.deg2rad(self._latitude_values()))
        indices = [slice(None)] + [np.newaxis] * (self.u.ndim - 1)
        broadcast = np.broadcast_to(coriolis[tuple(indices)], self.u.shape)
        if materialize:
            return np.array(broadcast, copy=True)
        return broadcast

    def _spectogrd_pair(
        self, spec_a: np.ndarray, spec_b: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        if hasattr(self.s, "_spectogrd_pair"):
            return self.s._spectogrd_pair(spec_a, spec_b)
        return self.s.spectogrd(spec_a), self.s.spectogrd(spec_b)

    def _spectral_cache_key(
        self, component: str, truncation: Optional[int]
    ) -> Tuple[str, int]:
        return component, self.s._validate_ntrunc(truncation)

    def _vector_analysis_spectral(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        key = self._spectral_cache_key("vrtdiv", truncation)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        ntrunc = key[1]
        vrtspec, divspec = self.s.getvrtdivspec(self.u, self.v, ntrunc=ntrunc)
        self._spectral_cache[key] = (vrtspec, divspec)
        self._spectral_cache[("vrt", ntrunc)] = vrtspec
        self._spectral_cache[("div", ntrunc)] = divspec
        return vrtspec, divspec

    def _vorticity_spectral(self, truncation: Optional[int] = None) -> np.ndarray:
        key = self._spectral_cache_key("vrt", truncation)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        full = self._spectral_cache.get(("vrtdiv", key[1]))
        if full is not None:
            vrtspec = full[0]  # type: ignore[index]
        else:
            vrtspec = self.s.getvrtspec(self.u, self.v, ntrunc=key[1])
        self._spectral_cache[key] = vrtspec
        return vrtspec

    def _divergence_spectral(self, truncation: Optional[int] = None) -> np.ndarray:
        key = self._spectral_cache_key("div", truncation)
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached  # type: ignore[return-value]

        full = self._spectral_cache.get(("vrtdiv", key[1]))
        if full is not None:
            divspec = full[1]  # type: ignore[index]
        else:
            divspec = self.s.getdivspec(self.u, self.v, ntrunc=key[1])
        self._spectral_cache[key] = divspec
        return divspec

    def magnitude(self) -> np.ndarray:
        """Return wind speed."""
        return self._restore_output_dtype(np.hypot(self.u, self.v))

    def vrtdiv(self, truncation: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Return relative vorticity and horizontal divergence."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation=truncation)
        vrt, div = self._spectogrd_pair(vrtspec, divspec)
        return self._restore_output_dtype(vrt), self._restore_output_dtype(div)

    def vorticity(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return relative vorticity."""
        vrtspec = self._vorticity_spectral(truncation=truncation)
        return self._restore_output_dtype(self.s.spectogrd(vrtspec))

    def divergence(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return horizontal divergence."""
        divspec = self._divergence_spectral(truncation=truncation)
        return self._restore_output_dtype(self.s.spectogrd(divspec))

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
        vrtspec, divspec = self._vector_analysis_spectral(truncation=truncation)
        psi, chi = self._spectogrd_pair(
            self.s._invlapspec(vrtspec, "sfvp"),
            self.s._invlapspec(divspec, "sfvp"),
        )
        return self._restore_output_dtype(psi), self._restore_output_dtype(chi)

    def streamfunction(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return streamfunction."""
        vrtspec = self._vorticity_spectral(truncation=truncation)
        return self._restore_output_dtype(
            self.s.spectogrd(self.s._invlapspec(vrtspec, "streamfunction"))
        )

    def velocitypotential(self, truncation: Optional[int] = None) -> np.ndarray:
        """Return velocity potential."""
        divspec = self._divergence_spectral(truncation=truncation)
        return self._restore_output_dtype(
            self.s.spectogrd(self.s._invlapspec(divspec, "velocitypotential"))
        )

    def helmholtz(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return irrotational and non-divergent wind components."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation=truncation)
        psispec = self.s._invlapspec(vrtspec, "helmholtz")
        chispec = self.s._invlapspec(divspec, "helmholtz")
        u_chi, v_chi, v_psi, u_psi = self.s.getgrad_pair(chispec, psispec)
        return (
            self._restore_output_dtype(u_chi),
            self._restore_output_dtype(v_chi),
            self._restore_output_dtype(-u_psi),
            self._restore_output_dtype(v_psi),
        )

    def irrotationalcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return irrotational wind component."""
        divspec = self._divergence_spectral(truncation=truncation)
        u_chi, v_chi = self.s.getgrad(
            self.s._invlapspec(divspec, "irrotationalcomponent")
        )
        return self._restore_output_dtype(u_chi), self._restore_output_dtype(v_chi)

    def nondivergentcomponent(
        self, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return non-divergent wind component."""
        vrtspec = self._vorticity_spectral(truncation=truncation)
        psispec = self.s._invlapspec(vrtspec, "nondivergentcomponent")
        v_psi, u_psi = self.s.getgrad(psispec)
        return self._restore_output_dtype(-u_psi), self._restore_output_dtype(v_psi)

    def gradient(
        self, chi: ArrayLike, truncation: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return vector gradient of a packed scalar field."""
        output_dtype = self._infer_output_dtype(chi)
        field = self._validate_scalar_field(chi, "chi")
        chi_spec = self.s.grdtospec(field, ntrunc=truncation)
        u_chi, v_chi = self.s.getgrad(chi_spec)
        return (
            self._restore_output_dtype(u_chi, output_dtype),
            self._restore_output_dtype(v_chi, output_dtype),
        )

    def truncate(
        self, field: ArrayLike, truncation: Optional[int] = None
    ) -> np.ndarray:
        """Apply spectral truncation to a packed scalar field."""
        output_dtype = self._infer_output_dtype(field)
        scalar = self._validate_scalar_field(field, "field")
        field_spec = self.s.grdtospec(scalar, ntrunc=truncation)
        return self._restore_output_dtype(self.s.spectogrd(field_spec), output_dtype)

    def rossbywavesource(
        self, truncation: Optional[int] = None, omega: Optional[float] = None
    ) -> np.ndarray:
        """Return Rossby wave source on the packed reduced grid."""
        vrtspec, divspec = self._vector_analysis_spectral(truncation=truncation)
        vrt, div = self._spectogrd_pair(vrtspec, divspec)
        eta = vrt + self._planetary_vorticity(omega=omega, materialize=False)
        uchi, vchi = self.s.getgrad(self.s._invlapspec(divspec, "rossbywavesource"))
        etaspec = self.s.grdtospec(eta, ntrunc=truncation)
        etax, etay = self.s.getgrad(etaspec)
        rws = -eta * div - (uchi * etax + vchi * etay)
        return self._restore_output_dtype(rws)
