"""Focused tests for Chemke-style growth-rate diagnostics."""

from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_allclose

import skyborn.calc.growth_rate.core as growth_rate_core
from skyborn.calc.growth_rate import baroc_growth_rate, barot_growth_rate

OMEGA = 7.292e-5
RADIUS = 6_371_000.0
GAS_CONSTANT_DRY = 287.04
HEAT_CAPACITY = 1004.7
KAPPA = GAS_CONSTANT_DRY / HEAT_CAPACITY


def _theta_from_temperature(
    temperature: np.ndarray, pressure_pa: np.ndarray
) -> np.ndarray:
    """Return potential temperature from temperature and pressure."""

    return (
        np.asarray(temperature, dtype=np.float64)
        * (100000.0 / np.asarray(pressure_pa, dtype=np.float64)) ** KAPPA
    )


def _barot_growth_reference(lat: np.ndarray, u: np.ndarray) -> float:
    """Pure NumPy port of the Chemke barotropic reference algorithm."""

    lat = np.asarray(lat, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    y = np.deg2rad(lat) * RADIUS
    dy1 = (y[1:] + y[:-1]) / 2.0
    dy2 = dy1[1:] - dy1[:-1]
    dy = (y[1:] - y[:-1]) * np.cos(np.deg2rad((lat[1:] + lat[:-1]) / 2.0))
    y1 = (y[1] - y[0]) * np.cos(np.deg2rad((lat[1] + lat[0]) / 2.0))
    y2 = (y[-1] - y[-2]) * np.cos(np.deg2rad((lat[-1] + lat[-2]) / 2.0))
    beta = 2.0 * OMEGA * np.cos(np.deg2rad(lat)) / RADIUS
    waven = np.arange(0.0, 2.0e-5, 1.0e-7, dtype=np.float64)
    growth = np.full(waven.shape, np.nan, dtype=np.float64)

    for k, kval in enumerate(waven):
        a = np.zeros((lat.size, lat.size), dtype=np.float64)
        b = np.zeros((lat.size, lat.size), dtype=np.float64)
        for n in range(lat.size):
            if n == lat.size - 1:
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / y2 * (-1.0 / dy[n - 1] - 1.0 / y2)
                    - kval / y2 * ((u[n - 1] - u[n]) / dy[n - 1] - u[n] / y2)
                )
                a[n, n - 1] = (kval * u[n]) / (y2 * dy[n - 1])
                b[n, n] = -(kval**2) + (1.0 / y2) * (-1.0 / dy[n - 1] - 1.0 / y2)
                b[n, n - 1] = 1.0 / (y2 * dy[n - 1])
            elif n == 0:
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / y1 * (-1.0 / y1 - 1.0 / dy[n])
                    - kval / y1 * ((-u[n]) / y1 - (u[n] - u[n + 1]) / dy[n])
                )
                a[n, n + 1] = (kval * u[n]) / (y1 * dy[n])
                b[n, n] = -(kval**2) + (1.0 / y1) * (-1.0 / y1 - 1.0 / dy[n])
                b[n, n + 1] = 1.0 / (y1 * dy[n])
            else:
                a[n, n - 1] = (kval * u[n]) / (dy2[n - 1] * dy[n - 1])
                a[n, n + 1] = (kval * u[n]) / (dy2[n - 1] * dy[n])
                a[n, n] = (
                    kval * beta[n]
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) / dy2[n - 1] * (-1.0 / dy[n - 1] - 1.0 / dy[n])
                    - kval
                    / dy2[n - 1]
                    * ((u[n - 1] - u[n]) / dy[n - 1] - (u[n] - u[n + 1]) / dy[n])
                )
                b[n, n - 1] = 1.0 / (dy2[n - 1] * dy[n - 1])
                b[n, n + 1] = 1.0 / (dy2[n - 1] * dy[n])
                b[n, n] = -(kval**2) + (1.0 / dy2[n - 1]) * (
                    -1.0 / dy[n - 1] - 1.0 / dy[n]
                )

        eigvals = np.linalg.eigvals(np.linalg.solve(b, a))
        growth[k] = np.max(np.imag(eigvals))

    return float(np.nanmax(growth))


def _baroc_growth_reference(
    u: np.ndarray,
    temperature: np.ndarray,
    pressure_pa: np.ndarray,
    lat: float,
    smooth_window: int = 1,
) -> float:
    """Pure NumPy port of the Chemke baroclinic reference algorithm."""

    u = np.asarray(u, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    pressure_pa = np.asarray(pressure_pa, dtype=np.float64)
    theta = _theta_from_temperature(temperature, pressure_pa)

    dz1 = (pressure_pa[1:] + pressure_pa[:-1]) / 2.0
    dz2 = np.empty_like(pressure_pa)
    dz2[0] = pressure_pa[1] - pressure_pa[0]
    dz2[1:-1] = dz1[1:] - dz1[:-1]
    dz2[-1] = pressure_pa[-1] - pressure_pa[-2]

    f = 2.0 * OMEGA * np.sin(np.deg2rad(lat))
    beta = 2.0 * OMEGA * np.cos(np.deg2rad(lat)) / RADIUS
    rho = pressure_pa / (GAS_CONSTANT_DRY * temperature)
    nz = -((theta[1:] - theta[:-1]) * (1.0 / (rho[:-1] * theta[:-1])))
    waven = np.arange(0.0, 2.0e-5, 1.0e-7, dtype=np.float64)
    growth = np.full(waven.shape, np.nan, dtype=np.float64)

    for k, kval in enumerate(waven):
        a = np.zeros((pressure_pa.size, pressure_pa.size), dtype=np.float64)
        b = np.zeros((pressure_pa.size, pressure_pa.size), dtype=np.float64)
        for n in range(pressure_pa.size):
            if n == pressure_pa.size - 1:
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n - 1])
                    - kval * (f**2) / dz2[n] * ((u[n - 1] - u[n]) / nz[n - 1])
                )
                a[n, n - 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n - 1])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (-1.0 / nz[n - 1])
                b[n, n - 1] = (f**2) / (dz2[n] * nz[n - 1])
            elif n == 0:
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n])
                    - kval * (f**2) / dz2[n] * ((u[n + 1] - u[n]) / nz[n])
                )
                a[n, n + 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (-1.0 / nz[n])
                b[n, n + 1] = (f**2) / (dz2[n] * nz[n])
            else:
                a[n, n - 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n - 1])
                a[n, n + 1] = (kval * u[n]) * (f**2) / (dz2[n] * nz[n])
                a[n, n] = (
                    kval * beta
                    - (kval**2) * (kval * u[n])
                    + (kval * u[n]) * (f**2) / dz2[n] * (-1.0 / nz[n - 1] - 1.0 / nz[n])
                    - kval
                    * (f**2)
                    / dz2[n]
                    * ((u[n - 1] - u[n]) / nz[n - 1] - (u[n] - u[n + 1]) / nz[n])
                )
                b[n, n - 1] = (f**2) / (dz2[n] * nz[n - 1])
                b[n, n + 1] = (f**2) / (dz2[n] * nz[n])
                b[n, n] = -(kval**2) + (f**2) / dz2[n] * (
                    -1.0 / nz[n - 1] - 1.0 / nz[n]
                )

        eigvals = np.linalg.eigvals(np.linalg.solve(b, a))
        growth[k] = np.max(np.imag(eigvals))

    if smooth_window > 1:
        half_window = smooth_window // 2
        growth_smoothed = np.full_like(growth, np.nan)
        for i in range(growth.size):
            start = max(0, i - half_window)
            stop = min(growth.size, i + half_window + 1)
            window_values = growth[start:stop]
            finite_mask = np.isfinite(window_values)
            if np.any(finite_mask):
                growth_smoothed[i] = np.mean(window_values[finite_mask])
        growth = growth_smoothed

    turning_points = np.where(growth[:-1] - growth[1:] > 0.0)[0]
    if turning_points.size == 0:
        return float("nan")

    return float(np.max(growth[: turning_points[-1] + 2]))


class TestGrowthRateHelpers:
    """Focused coverage for the growth-rate Python wrapper helpers."""

    def test_helper_missing_and_pressure_branches(self, monkeypatch):
        """Public normalization should accept hPa pressure and the pressure alias."""

        captured = {"methods": [], "source_pressure": []}

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            captured["methods"].append(method)
            captured["source_pressure"].append(
                np.asarray(p_in, dtype=np.float64).copy()
            )
            return np.asarray(x, dtype=np.float64)

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, f_cor, beta, smooth_window, *extra: (
                1.0e-6,
                0,
            ),
        )

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([230.0, 250.0, 280.0]),
            np.array([300.0, 650.0, 1000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
            method="pressure",
        )

        assert_allclose(result, 1.0e-6)
        assert_allclose(
            growth_rate_core._as_profile_vector(np.array([300.0]), "lat", 2),
            np.array([300.0, 300.0]),
        )
        assert captured["methods"] == ["linear", "linear"]
        assert_allclose(
            captured["source_pressure"][0],
            np.array([30000.0, 65000.0, 100000.0]),
        )
        matrix = growth_rate_core._coerce_profile_matrix(
            np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]),
            "values",
            pressure_size=2,
            pressure_dim=None,
        )
        assert matrix.shape == (3, 2)

    def test_helper_validation_errors(self):
        """Helper validation should reject malformed shapes, coordinates, and units."""

        with pytest.raises(ValueError, match="must be one-dimensional"):
            growth_rate_core._as_1d_float64(np.ones((2, 2)), "values")

        with pytest.raises(ValueError, match="one- or two-dimensional"):
            growth_rate_core._coerce_profile_matrix(
                np.ones((2, 2, 2)),
                "values",
                pressure_size=2,
                pressure_dim=None,
            )

        with pytest.raises(ValueError, match="dimension mismatch"):
            growth_rate_core._same_shape_or_raise(
                ("a", np.array([1.0, 2.0])),
                ("b", np.array([1.0])),
            )

        with pytest.raises(ValueError, match="strictly monotonic"):
            growth_rate_core._require_monotonic(
                np.array([1000.0, 900.0, 950.0]), "pressure"
            )

        with pytest.raises(ValueError, match="nonzero zonal distance"):
            growth_rate_core._prepare_wavenumber_inputs(
                "low",
                np.array([30.0], dtype=np.float64),
            )

        with pytest.raises(ValueError, match="at least two low-resolution"):
            growth_rate_core._prepare_wavenumber_inputs(
                "low",
                np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float64),
            )

        with pytest.raises(ValueError, match="match the pressure coordinate"):
            growth_rate_core._coerce_profile_matrix(
                np.array([1.0, 2.0]),
                "u",
                pressure_size=3,
                pressure_dim=None,
            )

        with pytest.raises(ValueError, match="exactly one axis matching"):
            growth_rate_core._coerce_profile_matrix(
                np.ones((3, 3)),
                "u",
                pressure_size=3,
                pressure_dim=None,
            )

        with pytest.raises(ValueError, match="length 1 or match the profile count"):
            growth_rate_core._as_profile_vector(np.array([1.0, 2.0]), "lat", 3)

        with pytest.raises(ValueError, match="scalar or one-dimensional"):
            growth_rate_core._as_profile_vector(np.ones((2, 2)), "lat", 2)

        with pytest.raises(ValueError, match="match the requested profile count"):
            growth_rate_core._broadcast_profile_matrix(np.ones((2, 3)), 3, "u")

        data = xr.DataArray(np.ones((2, 3), dtype=np.float64), dims=("time", "level"))
        matrix = growth_rate_core._coerce_profile_matrix(
            data,
            "u",
            pressure_size=3,
            pressure_dim="level",
        )
        assert matrix.shape == (2, 3)

        with pytest.raises(ValueError, match="vertical axis must match"):
            growth_rate_core._coerce_profile_matrix(
                xr.DataArray(
                    np.ones((2, 3), dtype=np.float64),
                    dims=("level", "profile"),
                ),
                "u",
                pressure_size=4,
                pressure_dim="level",
            )

        assert_allclose(
            growth_rate_core._as_profile_vector(np.array([45.0]), "lat", 3),
            np.array([45.0, 45.0, 45.0]),
        )

    def test_prepare_wavenumber_inputs_synonyms_and_invalid_mode(self):
        """Helper should normalize mode aliases and reject unsupported values."""
        high_mode = growth_rate_core._prepare_wavenumber_inputs("high-resolution", None)
        low_mode = growth_rate_core._prepare_wavenumber_inputs(
            "low_resolution",
            np.array([0.0, 60.0, 120.0, 180.0, 240.0, 300.0], dtype=np.float64),
        )

        assert high_mode[0] == 1
        assert low_mode[0] == 2
        assert low_mode[1] == 2
        assert low_mode[2] is not None

        with pytest.raises(ValueError, match="must be 'high' or 'low'"):
            growth_rate_core._prepare_wavenumber_inputs("medium", None)

    def test_latitude_band_helper_matches_chemke_reference_formula(self):
        """Latitude-band helper should reproduce the Chemke sin/cos averages."""

        f_cor, beta = growth_rate_core._f_beta_from_lat_bounds((30.0, 60.0))
        sin_avg = (np.sin(np.deg2rad(30.0)) + np.sin(np.deg2rad(60.0))) / 2.0
        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0

        assert_allclose(f_cor, 2.0 * OMEGA * sin_avg, rtol=0.0, atol=1e-15)
        assert_allclose(beta, 2.0 * OMEGA * cos_avg / RADIUS, rtol=0.0, atol=1e-21)

    def test_infer_tropopause_reorders_descending_pressure_for_wmo_solver(
        self, monkeypatch
    ):
        """The WMO tropopause helper should reorder descending pressure inputs."""

        captured = {}

        def fake_tropopause(temperature_profile, pressure_profile, pressure_unit="Pa"):
            captured["temperature"] = np.array(temperature_profile, copy=True)
            captured["pressure"] = np.array(pressure_profile, copy=True)
            captured["pressure_unit"] = pressure_unit
            return {
                "pressure": 250.0,
                "height": 11000.0,
                "level_index": 2,
                "lapse_rate": 1.8,
                "success": True,
            }

        monkeypatch.setattr(growth_rate_core, "_trop_wmo_profile", fake_tropopause)

        result = growth_rate_core._infer_tropopause_pressure_pa(
            np.array([285.0, 260.0, 230.0]),
            np.array([100000.0, 70000.0, 30000.0]),
        )

        assert_allclose(result, 25000.0)
        assert captured["pressure_unit"] == "Pa"
        assert_allclose(captured["pressure"], np.array([30000.0, 70000.0, 100000.0]))
        assert_allclose(captured["temperature"], np.array([230.0, 260.0, 285.0]))

    def test_build_solver_pressure_grid_rejects_non_tropospheric_bounds(
        self, monkeypatch
    ):
        """The public API should reject tropopauses at or below the lower bound."""

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 100000.0,
        )

        with pytest.raises(ValueError, match="diagnosed tropopause pressure"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([230.0, 260.0, 285.0]),
                np.array([30000.0, 70000.0, 100000.0]),
                lat=45.0,
            )

    def test_baroc_growth_rate_batch_rejects_conflicting_xarray_profile_dims(self):
        """Batched xarray inputs should share the same non-pressure dimension."""

        u = xr.DataArray(
            np.ones((2, 3), dtype=np.float64),
            dims=("time", "level"),
            coords={"time": [0, 1], "level": [300.0, 600.0, 1000.0]},
        )
        temperature = xr.DataArray(
            np.ones((2, 3), dtype=np.float64),
            dims=("year", "level"),
            coords={"year": [2000, 2001], "level": [300.0, 600.0, 1000.0]},
        )
        with pytest.raises(ValueError, match="share the same non-pressure dimension"):
            baroc_growth_rate(
                u,
                temperature,
                u["level"],
                lat=45.0,
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_xarray_latitude_band_reduces_lat_axis(self, monkeypatch):
        """Latitude-band xarray input should be cosine-weighted over the lat axis."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            captured["source_pressure"] = np.array(
                source_pressure, dtype=np.float64, copy=True
            )
            captured["f_cor"] = np.array(f_cor, dtype=np.float64, copy=True)
            captured["beta"] = np.array(beta, dtype=np.float64, copy=True)
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        lat_coord = np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64)
        u = xr.DataArray(
            np.array(
                [
                    [[1.0, 3.0, 5.0, 7.0], [2.0, 4.0, 6.0, 8.0], [3.0, 5.0, 7.0, 9.0]],
                    [[2.0, 4.0, 6.0, 8.0], [3.0, 5.0, 7.0, 9.0], [4.0, 6.0, 8.0, 10.0]],
                ],
                dtype=np.float64,
            ),
            dims=("time", "level", "lat"),
            coords={
                "time": [2000, 2001],
                "level": [300.0, 600.0, 1000.0],
                "lat": lat_coord,
            },
        )
        temperature = xr.DataArray(
            np.array(
                [
                    [
                        [220.0, 222.0, 224.0, 226.0],
                        [230.0, 232.0, 234.0, 236.0],
                        [240.0, 242.0, 244.0, 246.0],
                    ],
                    [
                        [221.0, 223.0, 225.0, 227.0],
                        [231.0, 233.0, 235.0, 237.0],
                        [241.0, 243.0, 245.0, 247.0],
                    ],
                ],
                dtype=np.float64,
            ),
            dims=("time", "level", "lat"),
            coords=u.coords,
        )

        result = baroc_growth_rate(
            u,
            temperature,
            u["level"],
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        selected_lat = xr.DataArray(
            lat_coord[1:3], dims=("lat",), coords={"lat": lat_coord[1:3]}
        )
        weights = xr.DataArray(
            np.cos(np.deg2rad(selected_lat.values)),
            dims=("lat",),
            coords={"lat": selected_lat.values},
        )
        expected_u = u.sel(lat=slice(30.0, 60.0)).weighted(weights).mean("lat")
        expected_temperature = (
            temperature.sel(lat=slice(30.0, 60.0)).weighted(weights).mean("lat")
        )
        sin_avg = (np.sin(np.deg2rad(30.0)) + np.sin(np.deg2rad(60.0))) / 2.0
        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("time",)
        assert_allclose(result.values, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(captured["u_input"], expected_u.values.T)
        assert_allclose(captured["temperature_input"], expected_temperature.values.T)
        assert_allclose(
            captured["source_pressure"],
            np.array([30000.0, 60000.0, 100000.0]),
        )
        assert_allclose(captured["f_cor"], np.full(2, 2.0 * OMEGA * sin_avg))
        assert_allclose(captured["beta"], np.full(2, 2.0 * OMEGA * cos_avg / RADIUS))

    def test_baroc_growth_rate_xarray_latitude_band_skips_inputs_without_lat_axis(
        self, monkeypatch
    ):
        """Latitude-band reduction should leave already-collapsed xarray inputs untouched."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        u = xr.DataArray(
            np.array(
                [[3.0, 4.0, 5.0], [4.0, 5.0, 6.0]],
                dtype=np.float64,
            ),
            dims=("time", "level"),
            coords={"time": [2000, 2001], "level": [300.0, 600.0, 1000.0]},
        )
        temperature = xr.DataArray(
            np.array(
                [
                    [
                        [220.0, 222.0, 224.0, 226.0],
                        [230.0, 232.0, 234.0, 236.0],
                        [240.0, 242.0, 244.0, 246.0],
                    ],
                    [
                        [221.0, 223.0, 225.0, 227.0],
                        [231.0, 233.0, 235.0, 237.0],
                        [241.0, 243.0, 245.0, 247.0],
                    ],
                ],
                dtype=np.float64,
            ),
            dims=("time", "level", "lat"),
            coords={
                "time": [2000, 2001],
                "level": [300.0, 600.0, 1000.0],
                "lat": [20.0, 40.0, 60.0, 80.0],
            },
        )

        result = baroc_growth_rate(
            u,
            temperature,
            u["level"],
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        weights = xr.DataArray(
            np.cos(np.deg2rad(np.array([40.0, 60.0], dtype=np.float64))),
            dims=("lat",),
            coords={"lat": [40.0, 60.0]},
        )
        expected_temperature = (
            temperature.sel(lat=slice(30.0, 60.0)).weighted(weights).mean("lat")
        )

        assert_allclose(result.values, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(captured["u_input"], u.values.T)
        assert_allclose(captured["temperature_input"], expected_temperature.values.T)

    def test_baroc_growth_rate_numpy_latitude_band_reduces_lat_axis(self, monkeypatch):
        """NumPy latitude-band input should use the explicit latitude coordinate."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            captured["f_cor"] = np.array(f_cor, dtype=np.float64, copy=True)
            captured["beta"] = np.array(beta, dtype=np.float64, copy=True)
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        lat_coord = np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64)
        u = np.array(
            [
                [[1.0, 3.0, 5.0, 7.0], [2.0, 4.0, 6.0, 8.0], [3.0, 5.0, 7.0, 9.0]],
                [[2.0, 4.0, 6.0, 8.0], [3.0, 5.0, 7.0, 9.0], [4.0, 6.0, 8.0, 10.0]],
            ],
            dtype=np.float64,
        )
        temperature = np.array(
            [
                [
                    [220.0, 222.0, 224.0, 226.0],
                    [230.0, 232.0, 234.0, 236.0],
                    [240.0, 242.0, 244.0, 246.0],
                ],
                [
                    [221.0, 223.0, 225.0, 227.0],
                    [231.0, 233.0, 235.0, 237.0],
                    [241.0, 243.0, 245.0, 247.0],
                ],
            ],
            dtype=np.float64,
        )

        result = baroc_growth_rate(
            u,
            temperature,
            np.array([300.0, 600.0, 1000.0], dtype=np.float64),
            lat=lat_coord,
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        weights = np.cos(np.deg2rad(lat_coord[1:3]))
        expected_u = np.average(u[:, :, 1:3], axis=2, weights=weights)
        expected_temperature = np.average(
            temperature[:, :, 1:3], axis=2, weights=weights
        )
        sin_avg = (np.sin(np.deg2rad(30.0)) + np.sin(np.deg2rad(60.0))) / 2.0
        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(captured["u_input"], expected_u.T)
        assert_allclose(captured["temperature_input"], expected_temperature.T)
        assert_allclose(captured["f_cor"], np.full(2, 2.0 * OMEGA * sin_avg))
        assert_allclose(captured["beta"], np.full(2, 2.0 * OMEGA * cos_avg / RADIUS))

    def test_baroc_growth_rate_numpy_latitude_band_requires_lat_coordinate(self):
        """Three-dimensional NumPy latitude-band input needs an explicit latitude coordinate."""

        with pytest.raises(ValueError, match="require `lat`"):
            baroc_growth_rate(
                np.ones((2, 3, 4), dtype=np.float64),
                np.full((2, 3, 4), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_xarray_latitude_band_rejects_multiple_lat_dims(self):
        """Latitude-band xarray input must not expose two latitude dimensions."""

        with pytest.raises(ValueError, match="at most one latitude dimension"):
            baroc_growth_rate(
                xr.DataArray(
                    np.ones((1, 3, 2, 2), dtype=np.float64),
                    dims=("time", "level", "lat", "latitude"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [30.0, 60.0],
                        "latitude": [30.0, 60.0],
                    },
                ),
                xr.DataArray(
                    np.ones((1, 3, 2, 2), dtype=np.float64) * 240.0,
                    dims=("time", "level", "lat", "latitude"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [30.0, 60.0],
                        "latitude": [30.0, 60.0],
                    },
                ),
                xr.DataArray([300.0, 600.0, 1000.0], dims=("level",)),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_xarray_latitude_band_requires_shared_lat_dim(self):
        """Latitude-band xarray inputs should agree on the latitude-dimension name."""

        with pytest.raises(ValueError, match="share the same latitude dimension"):
            baroc_growth_rate(
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64),
                    dims=("time", "level", "lat"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [30.0, 60.0],
                    },
                ),
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64) * 240.0,
                    dims=("time", "level", "latitude"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "latitude": [30.0, 60.0],
                    },
                ),
                xr.DataArray([300.0, 600.0, 1000.0], dims=("level",)),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_xarray_latitude_band_uses_embedded_coordinate(self):
        """Latitude-band xarray input should reject an extra public ``lat`` argument."""

        with pytest.raises(
            ValueError, match="latitude coordinate is taken from the DataArray"
        ):
            baroc_growth_rate(
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64),
                    dims=("time", "level", "lat"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [30.0, 60.0],
                    },
                ),
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64) * 240.0,
                    dims=("time", "level", "lat"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [30.0, 60.0],
                    },
                ),
                xr.DataArray([300.0, 600.0, 1000.0], dims=("level",)),
                lat=45.0,
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_xarray_latitude_band_ignores_nan_missing_values(
        self, monkeypatch
    ):
        """Latitude-band xarray reduction should ignore NaN samples."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        coords = {
            "time": [2000, 2001],
            "level": [300.0, 600.0, 1000.0],
            "lat": [20.0, 40.0, 60.0, 80.0],
        }
        u = xr.DataArray(
            np.array(
                [
                    [
                        [1.0, np.nan, 5.0, 7.0],
                        [2.0, np.nan, 6.0, 8.0],
                        [3.0, np.nan, 7.0, 9.0],
                    ],
                    [
                        [2.0, 4.0, np.nan, 8.0],
                        [3.0, 5.0, np.nan, 9.0],
                        [4.0, 6.0, np.nan, 10.0],
                    ],
                ],
                dtype=np.float64,
            ),
            dims=("time", "level", "lat"),
            coords=coords,
        )
        temperature = xr.DataArray(
            np.array(
                [
                    [
                        [220.0, np.nan, 224.0, 226.0],
                        [230.0, np.nan, 234.0, 236.0],
                        [240.0, np.nan, 244.0, 246.0],
                    ],
                    [
                        [221.0, 223.0, np.nan, 227.0],
                        [231.0, 233.0, np.nan, 237.0],
                        [241.0, 243.0, np.nan, 247.0],
                    ],
                ],
                dtype=np.float64,
            ),
            dims=("time", "level", "lat"),
            coords=coords,
        )

        result = baroc_growth_rate(
            u,
            temperature,
            u["level"],
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        assert_allclose(result.values, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(
            captured["u_input"],
            np.array([[5.0, 4.0], [6.0, 5.0], [7.0, 6.0]], dtype=np.float64),
        )
        assert_allclose(
            captured["temperature_input"],
            np.array(
                [[224.0, 223.0], [234.0, 233.0], [244.0, 243.0]], dtype=np.float64
            ),
        )

    def test_baroc_growth_rate_xarray_latitude_band_requires_selected_points(self):
        """Latitude-band xarray input should fail when the requested band is empty."""

        with pytest.raises(ValueError, match="does not select any latitude points"):
            baroc_growth_rate(
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64),
                    dims=("time", "level", "lat"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [10.0, 20.0],
                    },
                ),
                xr.DataArray(
                    np.ones((1, 3, 2), dtype=np.float64) * 240.0,
                    dims=("time", "level", "lat"),
                    coords={
                        "time": [0],
                        "level": [300.0, 600.0, 1000.0],
                        "lat": [10.0, 20.0],
                    },
                ),
                xr.DataArray([300.0, 600.0, 1000.0], dims=("level",)),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_numpy_latitude_band_requires_one_dimensional_lat(self):
        """NumPy latitude-band input should reject multi-dimensional latitude coordinates."""

        with pytest.raises(ValueError, match="must be one-dimensional"):
            baroc_growth_rate(
                np.ones((2, 3, 4), dtype=np.float64),
                np.full((2, 3, 4), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat=np.ones((2, 2), dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_numpy_latitude_band_requires_selected_points(self):
        """NumPy latitude-band input should fail when the requested band is empty."""

        with pytest.raises(ValueError, match="does not select any latitude points"):
            baroc_growth_rate(
                np.ones((2, 3, 4), dtype=np.float64),
                np.full((2, 3, 4), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat=np.array([70.0, 75.0, 80.0, 85.0], dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_numpy_latitude_band_rejects_invalid_array_rank(self):
        """NumPy latitude-band reduction only supports 2D or 3D raw arrays."""

        with pytest.raises(ValueError, match="two- or three-dimensional"):
            baroc_growth_rate(
                np.ones((2, 3, 4, 1), dtype=np.float64),
                np.full((2, 3, 4, 1), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat=np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_numpy_latitude_band_requires_one_pressure_axis(self):
        """NumPy latitude-band input should preserve the single-pressure-axis rule."""

        with pytest.raises(
            ValueError, match="exactly one axis matching the pressure coordinate"
        ):
            baroc_growth_rate(
                np.ones((2, 2, 4), dtype=np.float64),
                np.full((2, 2, 4), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat=np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_numpy_latitude_band_requires_one_latitude_axis(self):
        """NumPy latitude-band input should identify one explicit latitude axis."""

        with pytest.raises(
            ValueError, match="exactly one latitude axis matching `lat`"
        ):
            baroc_growth_rate(
                np.ones((2, 3, 5), dtype=np.float64),
                np.full((2, 3, 5), 240.0, dtype=np.float64),
                np.array([300.0, 600.0, 1000.0], dtype=np.float64),
                lat=np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64),
                lat_bounds=(30.0, 60.0),
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_numpy_latitude_band_ignores_nan_missing_values(
        self, monkeypatch
    ):
        """NumPy latitude-band reduction should ignore NaN samples."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        u = np.array(
            [
                [
                    [1.0, np.nan, 5.0, 7.0],
                    [2.0, np.nan, 6.0, 8.0],
                    [3.0, np.nan, 7.0, 9.0],
                ],
                [
                    [2.0, 4.0, np.nan, 8.0],
                    [3.0, 5.0, np.nan, 9.0],
                    [4.0, 6.0, np.nan, 10.0],
                ],
            ],
            dtype=np.float64,
        )
        temperature = np.array(
            [
                [
                    [220.0, np.nan, 224.0, 226.0],
                    [230.0, np.nan, 234.0, 236.0],
                    [240.0, np.nan, 244.0, 246.0],
                ],
                [
                    [221.0, 223.0, np.nan, 227.0],
                    [231.0, 233.0, np.nan, 237.0],
                    [241.0, 243.0, np.nan, 247.0],
                ],
            ],
            dtype=np.float64,
        )

        result = baroc_growth_rate(
            u,
            temperature,
            np.array([300.0, 600.0, 1000.0], dtype=np.float64),
            lat=np.array([20.0, 40.0, 60.0, 80.0], dtype=np.float64),
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(
            captured["u_input"],
            np.array([[5.0, 4.0], [6.0, 5.0], [7.0, 6.0]], dtype=np.float64),
        )
        assert_allclose(
            captured["temperature_input"],
            np.array(
                [[224.0, 223.0], [234.0, 233.0], [244.0, 243.0]], dtype=np.float64
            ),
        )


class TestGrowthRate:
    """Regression checks for the first-stage growth-rate API."""

    def test_barot_growth_rate_matches_reference(self):
        """Compiled barotropic growth should match the Chemke-style NumPy port."""

        lat = np.linspace(30.0, 60.0, 7)
        u = np.array([5.0, 12.0, 28.0, 42.0, 26.0, 11.0, 4.0], dtype=np.float64)

        expected = _barot_growth_reference(lat, u)
        result = barot_growth_rate(u, lat)

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_barot_growth_rate_constant_flow_is_nearly_neutral(self):
        """A constant barotropic flow should not produce meaningful growth."""

        lat = np.linspace(25.0, 55.0, 7)
        u = np.full(lat.shape, 20.0)

        result = barot_growth_rate(u, lat)

        assert np.isfinite(result)
        assert abs(result) < 1e-10

    def test_barot_growth_rate_rejects_custom_planetary_constants(self):
        """The compiled barotropic kernel should reject unsupported radius/omega values."""

        with pytest.raises(NotImplementedError, match="not implemented yet"):
            barot_growth_rate(
                np.array([10.0, 20.0, 15.0], dtype=np.float64),
                np.array([30.0, 45.0, 60.0], dtype=np.float64),
                radius=7_000_000.0,
            )

    def test_baroc_growth_rate_matches_reference_on_solver_grid_from_tropopause(self):
        """Compiled baroclinic growth should match the Chemke-style NumPy port."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        lat = 45.0

        expected = _baroc_growth_reference(u, temperature, pressure, lat)
        result = baroc_growth_rate(
            u,
            temperature,
            pressure,
            lat=lat,
            tropopause_pressure=200.0,
            solver_levels=5,
        )

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_baroc_growth_rate_requires_latitude(self):
        """Latitude is required for the Coriolis-dependent baroclinic solver."""

        with pytest.raises(ValueError, match="`lat` or `lat_bounds` is required"):
            baroc_growth_rate(
                np.array([10.0, 20.0, 30.0]),
                np.array([240.0, 260.0, 280.0]),
                np.array([30000.0, 60000.0, 100000.0]),
            )

    def test_baroc_growth_rate_rejects_conflicting_latitude_inputs(self):
        """Single-latitude and latitude-band controls cannot be mixed."""

        with pytest.raises(ValueError, match="cannot be provided together"):
            baroc_growth_rate(
                np.array([10.0, 20.0, 30.0]),
                np.array([240.0, 260.0, 280.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                lat_bounds=(30.0, 60.0),
            )

    def test_baroc_growth_rate_rejects_malformed_latitude_band(self):
        """Latitude-band control must contain exactly two endpoints."""

        with pytest.raises(ValueError, match="exactly two latitude values"):
            baroc_growth_rate(
                np.array([10.0, 20.0, 30.0]),
                np.array([240.0, 260.0, 280.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat_bounds=(30.0,),
            )

    def test_baroc_growth_rate_low_resolution_uses_default_chemke_lon_grid(
        self, monkeypatch
    ):
        """Low-resolution mode should fall back to the Chemke longitude grid."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["wavenumber_mode"] = int(extra[0])
            captured["wavenumber_count"] = int(extra[1])
            captured["zonal_length"] = float(extra[2])
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([10.0, 20.0, 30.0]),
            np.array([240.0, 260.0, 280.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            wavenumber_mode="low",
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        default_lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        expected_length = (
            abs(default_lon[0] - default_lon[-1])
            / 180.0
            * np.pi
            * RADIUS
            * np.cos(np.deg2rad(45.0))
        )

        assert_allclose(result, 1.0e-6)
        assert captured["wavenumber_mode"] == 2
        assert captured["wavenumber_count"] == 80
        assert_allclose(captured["zonal_length"], expected_length)

    def test_baroc_growth_rate_rejects_nonpositive_pressure(self):
        """Baroclinic pressure coordinates must remain strictly positive."""

        with pytest.raises(ValueError, match="strictly positive"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 0.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
            )

    def test_baroc_growth_rate_rejects_invalid_method(self):
        """The public interpolation method must be recognized."""

        with pytest.raises(ValueError, match="`method`"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                method="cubic",
            )

    def test_baroc_growth_rate_rejects_invalid_wavenumber_mode(self):
        """The public wavenumber mode must be recognized."""

        with pytest.raises(ValueError, match="`wavenumber_mode`"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                wavenumber_mode="medium",
            )

    def test_baroc_growth_rate_rejects_invalid_solver_levels_type(self):
        """The public solver-level count must be an integer."""

        with pytest.raises(TypeError, match="solver_levels"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3.5,
            )

    def test_baroc_growth_rate_rejects_too_small_solver_levels(self):
        """The public solver-level count must stay above the minimum."""

        with pytest.raises(ValueError, match="at least 2"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=1,
            )

    def test_baroc_growth_rate_rejects_invalid_smooth_window_type(self):
        """The public smoothing width must be an integer."""

        with pytest.raises(TypeError, match="smooth_window"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                smooth_window=2.5,
            )

    def test_baroc_growth_rate_rejects_too_small_smooth_window(self):
        """The public smoothing width must stay positive."""

        with pytest.raises(ValueError, match="at least 1"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                smooth_window=0,
            )

    def test_baroc_growth_rate_rejects_boolean_smooth_window(self):
        """The public smoothing width should not accept booleans."""

        with pytest.raises(TypeError, match="smooth_window"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                smooth_window=True,
            )

    def test_barot_growth_rate_returns_per_second_growth(self):
        """The public barotropic wrapper should return per-second growth."""

        lat = np.linspace(30.0, 60.0, 7)
        u = np.array([5.0, 12.0, 28.0, 42.0, 26.0, 11.0, 4.0], dtype=np.float64)

        growth = barot_growth_rate(u, lat)
        expected = _barot_growth_reference(lat, u)

        assert_allclose(growth, expected, rtol=1e-12, atol=1e-12)

    def test_barot_growth_rate_returns_nan_when_input_is_missing(self):
        """Missing barotropic inputs should short-circuit to NaN."""

        result = barot_growth_rate(
            np.array([10.0, np.nan, 20.0]),
            np.array([20.0, 30.0, 40.0]),
        )

        assert np.isnan(result)

    def test_barot_growth_rate_reorders_descending_latitudes_for_backend(
        self, monkeypatch
    ):
        """Descending latitude input should be reordered before the backend call."""

        captured = {}

        def fake_barot_backend(lat_values, u_values):
            captured["lat"] = np.array(lat_values, dtype=np.float64, copy=True)
            captured["u"] = np.array(u_values, dtype=np.float64, copy=True)
            return 2.5e-6, 0

        monkeypatch.setattr(
            growth_rate_core, "_dbarot_growth_rate_1d", fake_barot_backend
        )

        result = barot_growth_rate(
            np.array([30.0, 20.0, 10.0]),
            np.array([50.0, 40.0, 30.0]),
        )

        assert_allclose(result, 2.5e-6)
        assert_allclose(captured["lat"], np.array([30.0, 40.0, 50.0]))
        assert_allclose(captured["u"], np.array([10.0, 20.0, 30.0]))

    def test_barot_growth_rate_raises_when_backend_reports_error(self, monkeypatch):
        """Barotropic backend errors should raise a clear runtime error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_dbarot_growth_rate_1d",
            lambda lat_values, u_values: (0.0, 17),
        )

        with pytest.raises(RuntimeError, match="ier=17"):
            barot_growth_rate(
                np.array([5.0, 12.0, 28.0]),
                np.array([30.0, 45.0, 60.0]),
            )

    def test_baroc_growth_rate_uses_wmo_tropopause_when_target_grid_is_omitted(
        self, monkeypatch
    ):
        """When no explicit solver grid is given, WMO tropopause should be used."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        captured = {}

        monkeypatch.setattr(
            growth_rate_core,
            "_trop_wmo_profile",
            lambda temperature_profile, pressure_profile, pressure_unit="Pa": {
                "pressure": 230.0,
                "height": 11000.0,
                "level_index": 1,
                "lapse_rate": 1.8,
                "success": True,
            },
        )

        def fake_baroc_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["temperature_solver"] = np.array(
                temperature_solver, dtype=np.float64, copy=True
            )
            captured["f_cor"] = float(f_cor)
            captured["beta"] = float(beta)
            captured["smooth_window"] = int(smooth_window)
            return 1.25e-6, 0

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            fake_baroc_backend,
        )

        result = baroc_growth_rate(u, temperature, pressure, lat=45.0)

        assert_allclose(result, 1.25e-6, rtol=0.0, atol=0.0)
        assert_allclose(captured["f_cor"], 2.0 * OMEGA * np.sin(np.deg2rad(45.0)))
        assert_allclose(
            captured["beta"],
            2.0 * OMEGA * np.cos(np.deg2rad(45.0)) / RADIUS,
        )
        assert captured["smooth_window"] == 1
        assert_allclose(captured["pressure_solver"][0], 23000.0, rtol=0.0, atol=1e-12)
        assert_allclose(captured["pressure_solver"][-1], 100000.0, rtol=0.0, atol=1e-12)
        assert captured["temperature_solver"].shape == captured["pressure_solver"].shape

    def test_baroc_growth_rate_exposes_solver_levels_for_auto_grid(self, monkeypatch):
        """The public API should let callers control the automatic solver-grid size."""

        captured = {}

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 25000.0,
        )

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            captured["smooth_window"] = int(smooth_window)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0, 32.0, 38.0]),
            np.array([220.0, 235.0, 255.0, 272.0, 288.0]),
            np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0]),
            lat=45.0,
            solver_levels=11,
        )

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert captured["pressure_solver"].shape == (11,)
        assert captured["smooth_window"] == 1
        assert_allclose(captured["pressure_solver"][0], 25000.0, atol=1e-12)
        assert_allclose(captured["pressure_solver"][-1], 100000.0, atol=1e-12)

    def test_baroc_growth_rate_uses_explicit_tropopause_pressure_when_given(
        self, monkeypatch
    ):
        """Explicit tropopause pressure should override automatic WMO detection."""

        captured = {}

        def fake_wmo(*args, **kwargs):
            raise AssertionError("automatic WMO detection should not run")

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["pressure_solver"] = np.array(
                pressure_solver, dtype=np.float64, copy=True
            )
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_trop_wmo_profile", fake_wmo)
        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0, 32.0, 38.0]),
            np.array([220.0, 235.0, 255.0, 272.0, 288.0]),
            np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=250.0,
            solver_levels=5,
        )

        assert_allclose(result, 1.0e-6)
        assert_allclose(
            captured["pressure_solver"],
            np.linspace(25000.0, 100000.0, 5, dtype=np.float64),
        )

    def test_baroc_growth_rate_passes_smooth_window_to_fortran_backend(
        self, monkeypatch
    ):
        """The public API should forward smooth_window to the compiled backend."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["smooth_window"] = int(smooth_window)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
            smooth_window=5,
        )

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert captured["smooth_window"] == 5

    def test_baroc_growth_rate_defaults_to_high_resolution_wavenumbers(
        self, monkeypatch
    ):
        """The public API should default to the fixed high-resolution wavenumber grid."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["wavenumber_mode"] = int(extra[0])
            captured["wavenumber_count"] = int(extra[1])
            captured["zonal_length"] = float(extra[2])
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        assert_allclose(result, 1.0e-6)
        assert captured["wavenumber_mode"] == 1
        assert captured["wavenumber_count"] == 200
        assert_allclose(captured["zonal_length"], 1.0)

    def test_baroc_growth_rate_low_resolution_uses_lon_span(self, monkeypatch):
        """Low-resolution mode should forward the longitude-derived zonal length."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["wavenumber_mode"] = int(extra[0])
            captured["wavenumber_count"] = int(extra[1])
            captured["zonal_length"] = float(extra[2])
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat_bounds=(30.0, 60.0),
            lon=lon,
            wavenumber_mode="low",
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0
        expected_length = abs(lon[0] - lon[-1]) / 180.0 * np.pi * RADIUS * cos_avg

        assert_allclose(result, 1.0e-6)
        assert captured["wavenumber_mode"] == 2
        assert captured["wavenumber_count"] == 80
        assert_allclose(captured["zonal_length"], expected_length)

    def test_baroc_growth_rate_low_resolution_uses_single_latitude_geometry(
        self, monkeypatch
    ):
        """Low-resolution mode should use the supplied single latitude when no band is given."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["zonal_length"] = float(extra[2])
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            lon=lon,
            wavenumber_mode="low",
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        expected_length = (
            abs(lon[0] - lon[-1]) / 180.0 * np.pi * RADIUS * np.cos(np.deg2rad(45.0))
        )
        assert_allclose(result, 1.0e-6)
        assert_allclose(captured["zonal_length"], expected_length)

    def test_baroc_growth_rate_low_resolution_rejects_zero_zonal_length(self):
        """Low-resolution mode should reject unusable single-latitude geometry."""

        lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        with pytest.raises(ValueError, match="positive zonal length"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=np.nan,
                lon=lon,
                wavenumber_mode="low",
                tropopause_pressure=300.0,
                solver_levels=3,
            )

    def test_baroc_growth_rate_matches_reference_with_spectrum_smoothing(self):
        """Compiled spectrum smoothing should match the NumPy reference path."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0])
        u = np.array([8.0, 14.0, 24.0, 32.0, 38.0])
        lat = 45.0

        expected = _baroc_growth_reference(
            u,
            temperature,
            pressure,
            lat,
            smooth_window=5,
        )
        result = baroc_growth_rate(
            u,
            temperature,
            pressure,
            lat=lat,
            tropopause_pressure=200.0,
            solver_levels=5,
            smooth_window=5,
        )

        assert_allclose(result, expected, rtol=1e-8, atol=1e-11)

    def test_baroc_growth_rate_uses_chemke_latitude_band_dynamics(self, monkeypatch):
        """Latitude-band calls should forward Chemke-style band-mean ``f`` and ``beta``."""

        captured = {}

        def fake_backend(
            u_solver,
            theta_solver,
            pressure_solver,
            temperature_solver,
            f_cor,
            beta,
            smooth_window,
            *extra,
        ):
            captured["f_cor"] = float(f_cor)
            captured["beta"] = float(beta)
            return 1.0e-6, 0

        monkeypatch.setattr(growth_rate_core, "_dbaroc_growth_rate_1d", fake_backend)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        sin_avg = (np.sin(np.deg2rad(30.0)) + np.sin(np.deg2rad(60.0))) / 2.0
        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0

        assert_allclose(result, 1.0e-6, rtol=0.0, atol=0.0)
        assert_allclose(captured["f_cor"], 2.0 * OMEGA * sin_avg, rtol=0.0, atol=1e-15)
        assert_allclose(
            captured["beta"],
            2.0 * OMEGA * cos_avg / RADIUS,
            rtol=0.0,
            atol=1e-21,
        )

    def test_baroc_growth_rate_defaults_to_log_pressure_interpolation(
        self, monkeypatch
    ):
        """The default baroclinic interpolation mode should be log-pressure."""

        methods = []

        monkeypatch.setattr(
            growth_rate_core,
            "_infer_tropopause_pressure_pa",
            lambda temperature, pressure_pa: 30000.0,
        )

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            methods.append(method)
            return np.asarray(x, dtype=np.float64)

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, f_cor, beta, smooth_window, *extra: (
                1.0e-6,
                0,
            ),
        )

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([230.0, 250.0, 280.0]),
            np.array([30000.0, 65000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        assert_allclose(result, 1.0e-6)
        assert methods == ["log", "log"]

    def test_baroc_growth_rate_accepts_logp_alias(self, monkeypatch):
        """The legacy ``logp`` alias should still map to log-pressure interpolation."""

        methods = []

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            methods.append(method)
            return np.asarray(x, dtype=np.float64)

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, f_cor, beta, smooth_window, *extra: (
                1.0e-6,
                0,
            ),
        )

        baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([230.0, 250.0, 280.0]),
            np.array([30000.0, 65000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
            method="logp",
        )

        assert methods == ["log", "log"]

    def test_baroc_growth_rate_raises_when_auto_tropopause_detection_fails(
        self, monkeypatch
    ):
        """A failed WMO tropopause diagnosis should raise a clear error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_trop_wmo_profile",
            lambda temperature_profile, pressure_profile, pressure_unit="Pa": {
                "pressure": -999.0,
                "height": -999.0,
                "level_index": -1,
                "lapse_rate": -999.0,
                "success": False,
            },
        )

        with pytest.raises(RuntimeError, match="Automatic WMO tropopause detection"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0, 32.0, 38.0]),
                np.array([220.0, 235.0, 255.0, 272.0, 288.0]),
                np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0]),
                lat=45.0,
            )

    def test_baroc_growth_rate_returns_nan_when_input_is_missing(self):
        """Missing baroclinic inputs should short-circuit to NaN."""

        result = baroc_growth_rate(
            np.array([8.0, np.nan, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        assert np.isnan(result)

    def test_baroc_growth_rate_rejects_nonpositive_temperature(self):
        """Baroclinic temperature profiles must remain strictly positive."""

        with pytest.raises(ValueError, match="strictly positive"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 0.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
            )

    def test_baroc_growth_rate_returns_missing_when_interpolation_has_nan(
        self, monkeypatch
    ):
        """Interpolation failures should return NaN."""

        def fake_interp(x, p_in, p_out, *, method, extrapolate, missing_value):
            if float(np.asarray(x)[0]) == 8.0:
                return np.array([np.nan, np.nan, np.nan])
            return np.array([220.0, 235.0, 255.0])

        monkeypatch.setattr(growth_rate_core, "interp_pressure_1d", fake_interp)

        result = baroc_growth_rate(
            np.array([8.0, 14.0, 24.0]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=3,
        )

        assert np.isnan(result)

    def test_baroc_growth_rate_raises_when_backend_reports_error(self, monkeypatch):
        """Baroclinic backend errors should raise a clear runtime error."""

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_1d",
            lambda u_solver, theta_solver, pressure_solver, temperature_solver, f_cor, beta, smooth_window, *extra: (
                0.0,
                19,
            ),
        )

        with pytest.raises(RuntimeError, match="ier=19"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
            )

    def test_baroc_growth_rate_rejects_even_smooth_window(self):
        """Centered spectrum smoothing should require an odd window length."""

        with pytest.raises(ValueError, match="odd integer"):
            baroc_growth_rate(
                np.array([8.0, 14.0, 24.0]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
                solver_levels=3,
                smooth_window=4,
            )

    def test_baroc_growth_rate_batch_auto_diagnoses_tropopause(self, monkeypatch):
        """Batched profile input should use the compiled multi-profile WMO solver."""

        captured = {}

        def fake_trop_wmo(
            temperature,
            pressure,
            xdim,
            ydim,
            levdim,
            timedim=None,
            pressure_unit="hPa",
            lapse_criterion=2.0,
            missing_value=-999.0,
            check_pressure_order=True,
        ):
            captured["temperature"] = np.array(temperature, dtype=np.float64, copy=True)
            captured["pressure"] = np.array(pressure, dtype=np.float64, copy=True)
            captured["dims"] = (xdim, ydim, levdim, timedim)
            captured["pressure_unit"] = pressure_unit
            return {
                "pressure": np.array([300.0, 250.0], dtype=np.float64),
                "height": np.array([11000.0, 11500.0], dtype=np.float64),
                "level_index": np.array([1, 1], dtype=np.int32),
                "lapse_rate": np.array([1.8, 1.7], dtype=np.float64),
                "success": np.array([True, True]),
            }

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["target_pressure"] = np.array(
                target_pressure,
                dtype=np.float64,
                copy=True,
            )
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(growth_rate_core, "_trop_wmo", fake_trop_wmo)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([[220.0, 235.0, 255.0], [219.0, 234.0, 254.0]]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            solver_levels=4,
        )

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert captured["temperature"].shape == (2, 3)
        assert_allclose(captured["pressure"], np.array([30000.0, 60000.0, 100000.0]))
        assert captured["dims"] == (-1, 0, 1, None)
        assert captured["pressure_unit"] == "Pa"
        assert_allclose(
            captured["target_pressure"][:, 0],
            np.linspace(30000.0, 100000.0, 4, dtype=np.float64),
        )
        assert_allclose(
            captured["target_pressure"][:, 1],
            np.linspace(25000.0, 100000.0, 4, dtype=np.float64),
        )

    def test_baroc_growth_rate_batch_auto_tropopause_failure_returns_missing(
        self, monkeypatch
    ):
        """Profiles without a diagnosed WMO tropopause should return missing output."""

        captured = {}

        def fake_trop_wmo(
            temperature,
            pressure,
            xdim,
            ydim,
            levdim,
            timedim=None,
            pressure_unit="hPa",
            lapse_criterion=2.0,
            missing_value=-999.0,
            check_pressure_order=True,
        ):
            return {
                "pressure": np.array([300.0, -999.0], dtype=np.float64),
                "height": np.array([11000.0, -999.0], dtype=np.float64),
                "level_index": np.array([1, -999], dtype=np.int32),
                "lapse_rate": np.array([1.8, -999.0], dtype=np.float64),
                "success": np.array([True, False]),
            }

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["nprofile"] = u_input.shape[1]
            return np.array([1.0e-6]), np.array([0], dtype=np.int32)

        monkeypatch.setattr(growth_rate_core, "_trop_wmo", fake_trop_wmo)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([[220.0, 235.0, 255.0], [219.0, 234.0, 254.0]]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
        )

        assert captured["nprofile"] == 1
        assert_allclose(result[0], 1.0e-6)
        assert np.isnan(result[1])

    def test_baroc_growth_rate_batch_all_auto_tropopause_failures_skip_backend(
        self, monkeypatch
    ):
        """If no batch profile gets a valid WMO tropopause, the backend should not run."""

        def fake_trop_wmo(
            temperature,
            pressure,
            xdim,
            ydim,
            levdim,
            timedim=None,
            pressure_unit="hPa",
            lapse_criterion=2.0,
            missing_value=-999.0,
            check_pressure_order=True,
        ):
            return {
                "pressure": np.array([-999.0, -999.0], dtype=np.float64),
                "height": np.array([-999.0, -999.0], dtype=np.float64),
                "level_index": np.array([-999, -999], dtype=np.int32),
                "lapse_rate": np.array([-999.0, -999.0], dtype=np.float64),
                "success": np.array([False, False]),
            }

        def fail_backend(*args, **kwargs):
            raise AssertionError(
                "batch backend should not run when all tropopauses fail"
            )

        monkeypatch.setattr(growth_rate_core, "_trop_wmo", fake_trop_wmo)
        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fail_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([[220.0, 235.0, 255.0], [219.0, 234.0, 254.0]]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
        )

        assert np.all(np.isnan(result))

    def test_baroc_growth_rate_batch_broadcasts_climatology_to_backend(
        self, monkeypatch
    ):
        """The batched API should broadcast climatologies and call the batch backend once."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["u_input"] = np.array(u_input, dtype=np.float64, copy=True)
            captured["temperature_input"] = np.array(
                temperature_input, dtype=np.float64, copy=True
            )
            captured["source_pressure"] = np.array(
                source_pressure, dtype=np.float64, copy=True
            )
            captured["target_pressure"] = np.array(
                target_pressure, dtype=np.float64, copy=True
            )
            captured["f_cor"] = np.array(f_cor, dtype=np.float64, copy=True)
            captured["beta"] = np.array(beta, dtype=np.float64, copy=True)
            captured["interp_kind"] = int(interp_kind)
            captured["smooth_window"] = int(smooth_window)
            captured["wavenumber_mode"] = int(extra[0])
            captured["wavenumber_count"] = int(extra[1])
            captured["zonal_length"] = np.array(extra[2], dtype=np.float64, copy=True)
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert captured["u_input"].shape == (3, 2)
        assert captured["temperature_input"].shape == (3, 2)
        assert_allclose(
            captured["f_cor"],
            np.full(2, 2.0 * OMEGA * np.sin(np.deg2rad(45.0))),
        )
        assert_allclose(
            captured["beta"],
            np.full(2, 2.0 * OMEGA * np.cos(np.deg2rad(45.0)) / RADIUS),
        )
        assert captured["interp_kind"] == 2
        assert captured["smooth_window"] == 1
        assert captured["wavenumber_mode"] == 1
        assert captured["wavenumber_count"] == 200
        assert_allclose(captured["zonal_length"], np.ones(2, dtype=np.float64))
        assert_allclose(
            captured["target_pressure"][:, 0],
            np.linspace(30000.0, 100000.0, 4, dtype=np.float64),
        )

    def test_baroc_growth_rate_batch_low_resolution_uses_lon_span(self, monkeypatch):
        """Batch low-resolution mode should forward one zonal length per profile."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["wavenumber_mode"] = int(extra[0])
            captured["wavenumber_count"] = int(extra[1])
            captured["zonal_length"] = np.array(extra[2], dtype=np.float64, copy=True)
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        lat = np.array([40.0, 60.0], dtype=np.float64)
        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=lat,
            lon=lon,
            wavenumber_mode="low",
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        expected_lengths = (
            abs(lon[0] - lon[-1]) / 180.0 * np.pi * RADIUS * np.cos(np.deg2rad(lat))
        )

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert captured["wavenumber_mode"] == 2
        assert captured["wavenumber_count"] == 80
        assert_allclose(captured["zonal_length"], expected_lengths)

    def test_baroc_growth_rate_batch_low_resolution_rejects_zero_band_length(self):
        """Batch low-resolution mode should reject unusable latitude-band geometry."""

        lon = np.arange(0.0, 361.0, 1.5, dtype=np.float64)
        with pytest.raises(ValueError, match="positive zonal length"):
            baroc_growth_rate(
                np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat_bounds=(np.nan, 60.0),
                lon=lon,
                wavenumber_mode="low",
                tropopause_pressure=300.0,
                solver_levels=4,
            )

    def test_baroc_growth_rate_batch_uses_chemke_latitude_band_dynamics(
        self, monkeypatch
    ):
        """Batched latitude-band calls should pass band-mean ``f`` and ``beta`` arrays."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["f_cor"] = np.array(f_cor, dtype=np.float64, copy=True)
            captured["beta"] = np.array(beta, dtype=np.float64, copy=True)
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat_bounds=(30.0, 60.0),
            tropopause_pressure=300.0,
            solver_levels=4,
        )

        sin_avg = (np.sin(np.deg2rad(30.0)) + np.sin(np.deg2rad(60.0))) / 2.0
        cos_avg = (np.cos(np.deg2rad(30.0)) + np.cos(np.deg2rad(60.0))) / 2.0

        assert_allclose(result, np.array([1.0e-6, 2.0e-6]))
        assert_allclose(captured["f_cor"], np.full(2, 2.0 * OMEGA * sin_avg))
        assert_allclose(captured["beta"], np.full(2, 2.0 * OMEGA * cos_avg / RADIUS))

    def test_baroc_growth_rate_batch_returns_xarray_series(self, monkeypatch):
        """Batched xarray input should preserve the non-pressure coordinate."""

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        u = xr.DataArray(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]], dtype=np.float64),
            dims=("time", "level"),
            coords={"time": [2000, 2001], "level": [300.0, 600.0, 1000.0]},
        )
        temperature = xr.DataArray(
            np.array([220.0, 235.0, 255.0], dtype=np.float64),
            dims=("level",),
            coords={"level": [300.0, 600.0, 1000.0]},
        )

        result = baroc_growth_rate(
            u,
            temperature,
            temperature["level"],
            lat=45.0,
            tropopause_pressure=np.array([300.0, 250.0]),
            solver_levels=4,
        )

        assert isinstance(result, xr.DataArray)
        assert result.dims == ("time",)
        assert_allclose(result["time"], np.array([2000, 2001]))
        assert_allclose(result.values, np.array([1.0e-6, 2.0e-6]))

    def test_baroc_growth_rate_batch_masks_missing_profiles_with_nan_output(
        self, monkeypatch
    ):
        """Missing profiles should short-circuit to NaN output in batched mode."""

        captured = {}

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            captured["nprofile"] = u_input.shape[1]
            return np.array([1.0e-6]), np.array([0], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [np.nan, 13.0, 23.0]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
        )

        assert captured["nprofile"] == 1
        assert_allclose(result[0], 1.0e-6, rtol=1e-12, atol=1e-12)
        assert np.isnan(result[1])

    def test_baroc_growth_rate_batch_returns_missing_for_interp_failure_code(
        self, monkeypatch
    ):
        """Interpolation failures from the batch backend should map to missing output."""

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            return np.array([1.0e-6, -999.0]), np.array([0, 100], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        result = baroc_growth_rate(
            np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
        )

        assert_allclose(result[0], 1.0e-6)
        assert np.isnan(result[1])

    def test_baroc_growth_rate_batch_raises_when_profile_backend_reports_error(
        self, monkeypatch
    ):
        """Hard backend errors in batched mode should still raise."""

        def fake_backend(
            u_input,
            temperature_input,
            source_pressure,
            target_pressure,
            f_cor,
            beta,
            interp_kind,
            smooth_window,
            *extra,
        ):
            return np.array([1.0e-6, 2.0e-6]), np.array([0, 300], dtype=np.int32)

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fake_backend,
        )

        with pytest.raises(RuntimeError, match="profile indices \\[1\\]"):
            baroc_growth_rate(
                np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_batch_rejects_nonpositive_tropopause_pressure(self):
        """Explicit batch tropopause pressures must stay positive."""

        with pytest.raises(ValueError, match="strictly positive"):
            baroc_growth_rate(
                np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=np.array([300.0, 0.0]),
            )

    def test_baroc_growth_rate_batch_rejects_tropopause_at_lower_bound(self):
        """Explicit batch tropopause pressures must stay above the solver lower bound."""

        with pytest.raises(ValueError, match="lower-tropospheric pressure bound"):
            baroc_growth_rate(
                np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
                np.array([220.0, 235.0, 255.0]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=np.array([300.0, 1000.0]),
            )

    def test_baroc_growth_rate_batch_rejects_nonpositive_temperature(self):
        """Batched baroclinic temperature profiles must remain strictly positive."""

        with pytest.raises(ValueError, match="strictly positive"):
            baroc_growth_rate(
                np.array([[8.0, 14.0, 24.0], [7.0, 13.0, 23.0]]),
                np.array([[220.0, 235.0, 255.0], [219.0, 0.0, 254.0]]),
                np.array([30000.0, 60000.0, 100000.0]),
                lat=45.0,
                tropopause_pressure=300.0,
            )

    def test_baroc_growth_rate_batch_returns_nan_when_all_profiles_are_missing(
        self, monkeypatch
    ):
        """All-missing batched input should return NaN without backend work."""

        def fail_backend(*args, **kwargs):
            raise AssertionError("batch backend should not run for all-missing input")

        monkeypatch.setattr(
            growth_rate_core,
            "_dbaroc_growth_rate_profiles",
            fail_backend,
        )

        result = baroc_growth_rate(
            np.array([[np.nan, 14.0, 24.0], [7.0, 13.0, np.nan]]),
            np.array([220.0, 235.0, 255.0]),
            np.array([30000.0, 60000.0, 100000.0]),
            lat=45.0,
            tropopause_pressure=300.0,
        )

        assert np.all(np.isnan(result))

    def test_baroc_growth_rate_batch_matches_single_calls_with_explicit_tropopause(
        self,
    ):
        """The compiled batch path should match repeated single-profile calls."""

        pressure = np.array([20000.0, 40000.0, 60000.0, 80000.0, 100000.0])
        u = np.array(
            [[8.0, 14.0, 24.0, 32.0, 38.0], [9.0, 15.0, 22.0, 30.0, 36.0]],
            dtype=np.float64,
        )
        temperature = np.array([220.0, 235.0, 255.0, 272.0, 288.0], dtype=np.float64)
        tropopause = np.array([200.0, 250.0], dtype=np.float64)

        expected = np.array(
            [
                baroc_growth_rate(
                    u[idx],
                    temperature,
                    pressure,
                    lat=45.0,
                    tropopause_pressure=tropopause[idx],
                    solver_levels=5,
                )
                for idx in range(2)
            ]
        )
        result = baroc_growth_rate(
            u,
            temperature,
            pressure,
            lat=45.0,
            tropopause_pressure=tropopause,
            solver_levels=5,
        )

        assert_allclose(result, expected, rtol=1e-10, atol=1e-12)
