"""Ventilation index calculations for tropical cyclone research.

Implements the ventilated potential intensity (vPI) framework and the
ventilated genesis potential index (GPIv) following Chavas, Camargo &
Tippett (2025).

All functions accept and return ``xarray.DataArray`` objects for seamless
integration with gridded climate datasets (ERA5, CMIP6, etc.).
"""

from __future__ import annotations

import numpy as np
import xarray as xr

VI_MAX = 0.145  # Maximum ventilation index for nonzero vPI

# Thermodynamic constants used by the Chavas et al. (2025) entropy deficit.
DRY_AIR_CP = 1005.7  # J kg^-1 K^-1
DRY_AIR_GAS_CONSTANT = 287.04  # J kg^-1 K^-1
WATER_VAPOR_GAS_CONSTANT = 461.5  # J kg^-1 K^-1
LATENT_HEAT_VAPORIZATION = 2.501e6  # J kg^-1
TRIPLE_POINT_TEMPERATURE = 273.15  # K
REFERENCE_PRESSURE = 100000.0  # Pa
PRESSURE_600_HPA = 60000.0  # Pa


def _drop_scalar_coord(da: xr.DataArray, coord_name: str) -> xr.DataArray:
    """Remove scalar coordinate left by ``.sel(...)``."""
    if coord_name in da.coords and coord_name not in da.dims:
        return da.drop_vars(coord_name)
    return da


def _get_metpy_calc():
    """Lazy import of metpy.calc to avoid loading MetPy at module import time."""
    import metpy.calc as mpcalc

    return mpcalc


def _with_default_units(da: xr.DataArray, units: str) -> xr.DataArray:
    """Return a shallow copy with default units when metadata is missing."""
    if "units" in da.attrs:
        return da

    da_with_units = da.copy(deep=False)
    da_with_units.attrs = {**da.attrs, "units": units}
    return da_with_units


def _temperature_to_kelvin(da: xr.DataArray) -> xr.DataArray:
    """Return temperature in Kelvin, assuming Kelvin when units are missing."""
    units = str(da.attrs.get("units", "K")).strip().lower()
    celsius_units = {"c", "°c", "degc", "degree_celsius", "degrees_celsius"}
    if units in celsius_units or units in {"celsius", "degrees celsius"}:
        out = da + 273.15
        out.attrs = {**da.attrs, "units": "K"}
        return out
    return da


def _mixing_ratio_from_specific_humidity(q: xr.DataArray) -> xr.DataArray:
    """Convert specific humidity [kg/kg] to water-vapor mixing ratio [kg/kg]."""
    mpcalc = _get_metpy_calc()
    w = mpcalc.mixing_ratio_from_specific_humidity(q)
    if hasattr(w, "metpy"):
        w = w.metpy.dequantify()
    return w


def _saturation_vapor_pressure(temperature: xr.DataArray) -> xr.DataArray:
    """Calculate saturation vapor pressure over liquid water [Pa]."""
    mpcalc = _get_metpy_calc()
    t = _with_default_units(temperature, "K").metpy.quantify()
    return mpcalc.saturation_vapor_pressure(t).metpy.dequantify()


def _moist_entropy(
    pressure: float,
    temperature: xr.DataArray,
    mixing_ratio: xr.DataArray,
) -> xr.DataArray:
    """Calculate moist entropy [J kg^-1 K^-1]."""
    mpcalc = _get_metpy_calc()
    from metpy.units import units

    p = pressure * units.Pa
    t = _with_default_units(temperature, "K").metpy.quantify()
    w = _with_default_units(mixing_ratio, "kg/kg").metpy.quantify()

    vapor_pressure = mpcalc.vapor_pressure(p, w).metpy.dequantify()
    saturation_pressure = _saturation_vapor_pressure(temperature)
    relative_humidity = (vapor_pressure / saturation_pressure).clip(min=1e-10)
    dry_air_pressure = pressure - vapor_pressure

    return (
        DRY_AIR_CP * np.log(temperature / TRIPLE_POINT_TEMPERATURE)
        - DRY_AIR_GAS_CONSTANT * np.log(dry_air_pressure / REFERENCE_PRESSURE)
        + LATENT_HEAT_VAPORIZATION * mixing_ratio / temperature
        - WATER_VAPOR_GAS_CONSTANT * mixing_ratio * np.log(relative_humidity)
    )


def _saturation_entropy(
    pressure: float,
    temperature: xr.DataArray,
) -> xr.DataArray:
    """Calculate saturation moist entropy [J kg^-1 K^-1]."""
    mpcalc = _get_metpy_calc()
    from metpy.units import units

    p = pressure * units.Pa
    t = _with_default_units(temperature, "K").metpy.quantify()

    saturation_pressure = _saturation_vapor_pressure(temperature)
    saturation_mixing_ratio = mpcalc.saturation_mixing_ratio(p, t).metpy.dequantify()
    dry_air_pressure = pressure - saturation_pressure

    return (
        DRY_AIR_CP * np.log(temperature / TRIPLE_POINT_TEMPERATURE)
        - DRY_AIR_GAS_CONSTANT * np.log(dry_air_pressure / REFERENCE_PRESSURE)
        + LATENT_HEAT_VAPORIZATION * saturation_mixing_ratio / temperature
    )


def _ventilated_pi_ratio(vi_values: np.ndarray, vi_max: float) -> np.ndarray:
    """Calculate vPI/PI from VI values."""
    vi_values = np.asarray(vi_values, dtype=float)

    # Normalized ventilation ratio.
    r = vi_values / vi_max

    # Solve the cubic y³ - y + 2r/(3√3) = 0 for y = vPI / PI
    # using np.lib.scimath for type-safe complex intermediate values.
    with np.errstate(divide="ignore", invalid="ignore"):
        term1 = np.lib.scimath.sqrt(r**2 - 1.0)
        term2 = np.lib.scimath.power(term1 - r, 1.0 / 3.0)
        x = (1.0 / np.sqrt(3.0)) * term2
        ratio = np.real(x + 1.0 / (3.0 * x))

    mask = np.isfinite(vi_values) & (vi_values > 0) & (vi_values <= vi_max)
    return np.where(mask, ratio, np.nan).astype(float, copy=False)


def vertical_wind_shear(
    ds: xr.Dataset,
    u_var: str = "U",
    v_var: str = "V",
    level_coord: str = "level",
) -> xr.DataArray:
    """Calculate the 200–850 hPa vertical wind shear.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing wind components with a pressure-level coordinate.
    u_var : str
        Name of the zonal wind variable (default ``"U"``).
    v_var : str
        Name of the meridional wind variable (default ``"V"``).
    level_coord : str
        Name of the pressure-level coordinate (default ``"level"``).

    Returns
    -------
    xr.DataArray
        Wind shear magnitude in m/s with pressure-level dimension removed.
    """
    u200 = _drop_scalar_coord(ds[u_var].sel({level_coord: 200}), level_coord)
    u850 = _drop_scalar_coord(ds[u_var].sel({level_coord: 850}), level_coord)
    v200 = _drop_scalar_coord(ds[v_var].sel({level_coord: 200}), level_coord)
    v850 = _drop_scalar_coord(ds[v_var].sel({level_coord: 850}), level_coord)

    vws = np.hypot(u200 - u850, v200 - v850)
    vws.attrs = {
        "long_name": "Vertical Wind Shear (200-850 hPa)",
        "units": "m/s",
    }
    return vws


def entropy_deficit(
    ds: xr.Dataset,
    air_sea_disequilibrium: xr.DataArray,
    t_var: str = "T",
    q_var: str = "Q",
    level_coord: str = "level",
) -> xr.DataArray:
    """Calculate the entropy deficit parameter (Chi) at 600 hPa.

    Follows Chavas et al. (2025), Eq. (3)::

        χ = (s*_m(600) − s_m(600)) / (s*_SST − s_b)

    where the numerator is computed from 600-hPa temperature and
    specific humidity, and ``air_sea_disequilibrium`` supplies the
    denominator from the potential-intensity calculation.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with temperature and humidity variables on pressure levels.
    air_sea_disequilibrium : xr.DataArray
        Air-sea entropy disequilibrium term [J kg^-1 K^-1].
    t_var : str
        Name of the temperature variable in Kelvin (default ``"T"``).
    q_var : str
        Name of the specific humidity variable in kg/kg (default ``"Q"``).
    level_coord : str
        Name of the pressure-level coordinate (default ``"level"``).

    Returns
    -------
    xr.DataArray
        Dimensionless entropy deficit.
    """
    t600 = _temperature_to_kelvin(
        _drop_scalar_coord(ds[t_var].sel({level_coord: 600}), level_coord)
    )
    q600 = _drop_scalar_coord(ds[q_var].sel({level_coord: 600}), level_coord)
    mixing_ratio_600 = _mixing_ratio_from_specific_humidity(q600)

    moist_entropy_600 = _moist_entropy(PRESSURE_600_HPA, t600, mixing_ratio_600)
    saturation_entropy_600 = _saturation_entropy(PRESSURE_600_HPA, t600)
    chi = (saturation_entropy_600 - moist_entropy_600) / air_sea_disequilibrium

    chi = chi.where(np.isfinite(chi) & (chi > 0))
    chi = chi.clip(max=10.0)
    chi.attrs = {"long_name": "Entropy Deficit (Chi)", "units": "1"}
    return chi


def ventilation_index(
    vws: xr.DataArray,
    chi: xr.DataArray,
    pi: xr.DataArray,
) -> xr.DataArray:
    """Calculate the ventilation index (VI).

    VI = VWS * Chi / PI  (Tang & Emanuel 2012, Eq. (2))

    Parameters
    ----------
    vws : xr.DataArray
        Vertical wind shear [m/s].
    chi : xr.DataArray
        Entropy deficit [dimensionless].
    pi : xr.DataArray
        Potential intensity [m/s].

    Returns
    -------
    xr.DataArray
        Ventilation index (dimensionless), NaN where <= 0.
    """
    vi = vws * chi / pi
    vi = vi.where(np.isfinite(vi) & (vi > 0))
    vi.attrs = {"long_name": "Ventilation Index", "units": "1"}
    return vi


def air_sea_disequilibrium(
    pi: xr.DataArray,
    sst: xr.DataArray,
    outflow_temperature: xr.DataArray,
    ckcd: float = 0.9,
) -> xr.DataArray:
    """Calculate the air-sea entropy disequilibrium term used by ``Chi``.

    This inverts the Emanuel PI relationship using the PI diagnostics already
    produced by :mod:`skyborn.calc.GPI.xarray`.

    Parameters
    ----------
    pi : xr.DataArray
        Potential intensity [m/s].
    sst : xr.DataArray
        Sea-surface temperature [K] or [degC] if a units attribute is present.
    outflow_temperature : xr.DataArray
        Outflow temperature [K] or [degC] if a units attribute is present.
    ckcd : float
        Exchange coefficient ratio (default 0.9).

    Returns
    -------
    xr.DataArray
        Air-sea entropy disequilibrium [J kg^-1 K^-1].
    """
    sst_k = _temperature_to_kelvin(sst)
    outflow_temperature_k = _temperature_to_kelvin(outflow_temperature)
    asdeq = pi**2 * (1.0 / ckcd) * outflow_temperature_k
    asdeq = asdeq / (sst_k * (sst_k - outflow_temperature_k))
    asdeq = asdeq.where(np.isfinite(asdeq) & (asdeq > 0))
    asdeq.attrs = {
        "long_name": "Air-Sea Entropy Disequilibrium",
        "units": "J kg^-1 K^-1",
    }
    return asdeq


def ventilated_pi(
    pi: xr.DataArray,
    vi: xr.DataArray,
    vi_max: float = VI_MAX,
) -> xr.DataArray:
    """Calculate the ventilated potential intensity (vPI).

    Solves the cubic equation for normalized vPI using the analytic
    Cardano formula, matching Eqs. (5)-(6) of Chavas et al. (2025) and
    the reference implementation in tcvpigpiv v0.3.x.

    Boundary behaviour:

    - VI → 0⁺:  vPI → PI  (no reduction)
    - VI = VI_max:  vPI = PI / √3 ≈ 0.5774 × PI
    - VI ≤ 0 or VI > VI_max:  vPI = NaN  (genesis impossible)

    Parameters
    ----------
    pi : xr.DataArray
        Potential intensity [m/s].
    vi : xr.DataArray
        Ventilation index [dimensionless].
    vi_max : float
        Maximum VI for nonzero vPI (default 0.145).

    Returns
    -------
    xr.DataArray
        Ventilated potential intensity [m/s].
    """
    ratio_da = xr.apply_ufunc(
        _ventilated_pi_ratio,
        vi,
        kwargs={"vi_max": vi_max},
        dask="parallelized",
        output_dtypes=[float],
        keep_attrs=False,
    )
    vpi = pi * ratio_da
    vpi = vpi.where(np.isfinite(vpi) & (vpi > 0))
    vpi.attrs = {"long_name": "Ventilated Potential Intensity", "units": "m/s"}
    return vpi


def ventilation_components(
    ds: xr.Dataset,
    sst_var: str = "SSTK",
    psl_var: str = "SP",
    t_var: str = "T",
    q_var: str = "Q",
    u_var: str = "U",
    v_var: str = "V",
    level_coord: str = "level",
    lat_coord: str = "latitude",
    lon_coord: str = "longitude",
    ckcd: float = 0.9,
    vi_max: float = VI_MAX,
    vorticity_cap: float = 3.7e-5,
    dx: float = 2.0,
    dy: float = 2.0,
) -> xr.Dataset:
    """Calculate PI, ventilation metrics, and GPIv from one dataset.

    This high-level pipeline reuses the existing Skyborn GPI xarray backend for
    potential intensity diagnostics, then applies the ventilation framework.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing SST, surface pressure, temperature, humidity, and
        wind components.
    sst_var, psl_var, t_var, q_var, u_var, v_var : str
        Variable names in ``ds``.
    level_coord : str
        Name of the pressure-level coordinate.
    lat_coord : str
        Name of the latitude coordinate (default ``"latitude"``).
    lon_coord : str
        Name of the longitude coordinate (default ``"longitude"``).
    ckcd : float
        Exchange coefficient ratio passed to the GPI diagnostics.
    vi_max : float
        Maximum VI for nonzero vPI.
    vorticity_cap : float
        Upper bound for capped 850-hPa absolute vorticity.
    dx, dy : float
        Grid spacing in degrees for the final GPIv area factor.

    Returns
    -------
    xr.Dataset
        Dataset with PI diagnostics and ventilation components.
    """
    from skyborn.calc.GPI.xarray import pi_log_decomposition

    mixing_ratio = _mixing_ratio_from_specific_humidity(ds[q_var])
    mixing_ratio.attrs = {**ds[q_var].attrs, "units": "kg/kg"}

    pi_diagnostics = pi_log_decomposition(
        ds[sst_var],
        ds[psl_var],
        ds[level_coord],
        ds[t_var],
        mixing_ratio,
        CKCD=ckcd,
    )

    pi = pi_diagnostics["pi"]
    asdeq = air_sea_disequilibrium(pi, ds[sst_var], pi_diagnostics["t0"], ckcd=ckcd)
    vws = vertical_wind_shear(
        ds,
        u_var=u_var,
        v_var=v_var,
        level_coord=level_coord,
    )
    chi = entropy_deficit(
        ds,
        asdeq,
        t_var=t_var,
        q_var=q_var,
        level_coord=level_coord,
    )
    vi = ventilation_index(vws, chi, pi)
    vpi = ventilated_pi(pi, vi, vi_max=vi_max)
    eta_c = absolute_vorticity_850(
        ds,
        u_var=u_var,
        v_var=v_var,
        level_coord=level_coord,
        lat_coord=lat_coord,
        lon_coord=lon_coord,
        cap=vorticity_cap,
    )
    gpi_v = genesis_potential_index(vpi, eta_c, dx=dx, dy=dy, lat_coord=lat_coord)

    result = xr.Dataset(
        data_vars={
            "PI": pi,
            "vPI": vpi,
            "VWS": vws,
            "Chi": chi,
            "ventilation_index": vi,
            "eta_c": eta_c,
            "GPIv": gpi_v,
            "air_sea_disequilibrium": asdeq,
            "min_pressure": pi_diagnostics["min_pressure"],
            "error_flag": pi_diagnostics["error_flag"],
            "t0": pi_diagnostics["t0"],
        },
        attrs={
            "title": "Ventilated Genesis Potential Index Components",
            "method": "Chavas et al. (2025) ventilated PI/GPIv",
        },
    )
    return result


def absolute_vorticity_850(
    ds: xr.Dataset,
    u_var: str = "U",
    v_var: str = "V",
    level_coord: str = "level",
    lat_coord: str = "latitude",
    lon_coord: str = "longitude",
    cap: float = 3.7e-5,
) -> xr.DataArray:
    """Calculate the clipped 850-hPa absolute vorticity.

    Following Chavas et al. (2025), Eq. (7):

    .. math::
       \\eta_c =
       \\begin{cases}
         \\text{cap} & \\text{if } |f + \\zeta_{850}| > \\text{cap} \\\\
         f + \\zeta_{850} & \\text{otherwise}
       \\end{cases}

    where *f* is the Coriolis parameter, ζ₈₅₀ is the 850-hPa relative
    vorticity, and ``f + ζ₈₅₀`` is calculated with MetPy on the
    latitude-longitude grid.  Large-magnitude values are set to the positive
    cap, and non-positive values are masked out (NaN), matching the
    reference implementation used for Chavas et al. (2025).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with wind components on pressure levels.
    u_var, v_var : str
        Variable names for zonal and meridional wind.
    level_coord : str
        Name of the pressure-level coordinate (default ``"level"``).
    lat_coord : str
        Name of the latitude coordinate (default ``"latitude"``).
    lon_coord : str
        Name of the longitude coordinate (default ``"longitude"``).
    cap : float
        Upper bound on absolute vorticity (default 3.7e-5 s^-1).

    Returns
    -------
    xr.DataArray
        Clipped absolute vorticity [s^-1], NaN where non-positive.
    """
    u = _drop_scalar_coord(ds[u_var].sel({level_coord: 850}), level_coord)
    v = _drop_scalar_coord(ds[v_var].sel({level_coord: 850}), level_coord)

    mpcalc = _get_metpy_calc()
    u_for_vorticity = _with_default_units(u, "m/s").metpy.quantify()
    v_for_vorticity = _with_default_units(v, "m/s").metpy.quantify()
    eta = mpcalc.absolute_vorticity(
        u_for_vorticity,
        v_for_vorticity,
        latitude=ds[lat_coord],
        longitude=ds[lon_coord],
        x_dim=u.get_axis_num(lon_coord),
        y_dim=u.get_axis_num(lat_coord),
    ).metpy.dequantify()
    eta = eta.rename("absolute_vorticity")
    if "metpy_crs" in eta.coords:
        eta = eta.drop_vars("metpy_crs")

    eta = eta.assign_coords(
        {
            coord: u.coords[coord]
            for coord in u.coords
            if coord in eta.coords and coord != "metpy_crs"
        }
    )

    # Clip: set large-magnitude absolute vorticity to the positive cap.
    eta_c = xr.where(np.abs(eta) > cap, cap, eta)
    eta_c = eta_c.where(eta_c > 0)  # Remove non-positive values

    eta_c.attrs = {
        "long_name": "Capped 850-hPa Absolute Vorticity",
        "units": "s^-1",
    }
    return eta_c


def genesis_potential_index(
    vpi: xr.DataArray,
    eta_c: xr.DataArray,
    dx: float = 2.0,
    dy: float = 2.0,
    lat_coord: str = "latitude",
) -> xr.DataArray:
    """Calculate the ventilated genesis potential index (GPIv).

    Following Chavas et al. (2025), Eq. (8):

    .. math::
       \\mathrm{GPI}_v = (102.1 \\times \\mathrm{vPI} \\times \\eta_c)^{4.90}
       \\times \\cos\\varphi \\times \\Delta x \\times \\Delta y

    Parameters
    ----------
    vpi : xr.DataArray
        Ventilated potential intensity [m/s].
    eta_c : xr.DataArray
        Clipped absolute vorticity [s^-1].
    dx, dy : float
        Grid spacing in degrees (default 2.0, matching the paper).
    lat_coord : str
        Name of the latitude coordinate (default ``"latitude"``).

    Returns
    -------
    xr.DataArray
        Ventilated genesis potential index.
    """
    cos_lat = np.cos(np.deg2rad(vpi[lat_coord]))
    gpi_v = (102.1 * vpi * eta_c) ** 4.90 * cos_lat * dx * dy
    gpi_v = gpi_v.where(np.isfinite(gpi_v) & (gpi_v > 0))
    gpi_v.attrs = {
        "long_name": "Ventilated Genesis Potential Index",
        "units": "arbitrary",
    }
    return gpi_v
