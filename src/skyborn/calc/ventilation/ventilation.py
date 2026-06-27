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

# Physical constants
EARTH_RADIUS = 6_371_000.0  # m
EARTH_OMEGA = 7.2921159e-5  # s^-1
VI_MAX = 0.145  # Maximum ventilation index for nonzero vPI


def _drop_scalar_level(da: xr.DataArray) -> xr.DataArray:
    """Remove scalar level coordinate left by ``.sel(level=...)``."""
    if "level" in da.coords and "level" not in da.dims:
        return da.drop_vars("level")
    return da


def vertical_wind_shear(
    ds: xr.Dataset,
    u_var: str = "U",
    v_var: str = "V",
) -> xr.DataArray:
    """Calculate the 200–850 hPa vertical wind shear.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing wind components with a ``level`` dimension.
    u_var : str
        Name of the zonal wind variable (default ``"U"``).
    v_var : str
        Name of the meridional wind variable (default ``"V"``).

    Returns
    -------
    xr.DataArray
        Wind shear magnitude in m/s with ``level`` dimension removed.
    """
    u200 = _drop_scalar_level(ds[u_var].sel(level=200))
    u850 = _drop_scalar_level(ds[u_var].sel(level=850))
    v200 = _drop_scalar_level(ds[v_var].sel(level=200))
    v850 = _drop_scalar_level(ds[v_var].sel(level=850))

    vws = np.hypot(u200 - u850, v200 - v850)
    vws.attrs = {
        "long_name": "Vertical Wind Shear (200-850 hPa)",
        "units": "m/s",
    }
    return vws


def entropy_deficit(
    ds: xr.Dataset,
    rh_var: str = "R",
    q_var: str = "Q",
) -> xr.DataArray:
    """Calculate the entropy deficit parameter (Chi) at 600 hPa.

    Uses a robust relative-humidity-based proxy following the approach
    in Chavas et al. (2025), Eq. (3).  When 600-hPa relative humidity
    is available, the proxy is::

        χ = 1 − RH₆₀₀ / 100

    Otherwise, a specific-humidity-based normalization is used.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with humidity variables on pressure levels.
    rh_var : str
        Name of the relative humidity variable (default ``"R"``).
    q_var : str
        Name of the specific humidity variable (default ``"Q"``).

    Returns
    -------
    xr.DataArray
        Dimensionless entropy deficit, clipped to [0.02, 1.5].
    """
    if rh_var in ds:
        chi = 1.0 - _drop_scalar_level(ds[rh_var].sel(level=600)) / 100.0
    else:
        q600 = _drop_scalar_level(ds[q_var].sel(level=600))
        chi = 1.0 - q600 / q600.quantile(0.98)

    chi = chi.clip(min=0.02, max=1.5)
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


def ventilated_pi(
    pi: xr.DataArray,
    vi: xr.DataArray,
    vi_max: float = VI_MAX,
) -> xr.DataArray:
    """Calculate the ventilated potential intensity (vPI).

    Solves the cubic equation for normalized vPI using the analytic
    Cardano formula, matching Eqs. (5)-(6) of Chavas et al. (2025) and
    the reference implementation in tcvpigpiv v0.3.x.

    At VI = 0, vPI = PI (no reduction).
    At VI = VI_max, vPI = PI / sqrt(3) ≈ 0.5774 * PI.
    For VI > VI_max, vPI = 0 (genesis impossible).

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
    vi_limited = vi.where(vi <= vi_max)

    def _ratio(v):
        if not np.isfinite(v) or v <= 0 or v > vi_max:
            return np.nan
        ratio = v / vi_max
        term1 = (ratio**2 - 1.0) ** 0.5
        term2 = (term1 - ratio) ** (1.0 / 3.0)
        x = (1.0 / np.sqrt(3.0)) * term2
        return float(np.real(x + 1.0 / (3.0 * x)))

    ratio = xr.apply_ufunc(_ratio, vi_limited, vectorize=True, output_dtypes=[float])
    vpi = pi * ratio
    vpi = vpi.where(np.isfinite(vpi) & (vpi > 0))
    vpi.attrs = {"long_name": "Ventilated Potential Intensity", "units": "m/s"}
    return vpi


def absolute_vorticity_850(
    ds: xr.Dataset,
    u_var: str = "U",
    v_var: str = "V",
    cap: float = 3.7e-5,
) -> xr.DataArray:
    """Calculate the clipped 850-hPa absolute vorticity.

    Following Chavas et al. (2025), Eq. (7):

    .. math::
       \\eta_c = \\min(3.7 \\times 10^{-5}, |f + \\zeta_{850}|)

    where f is the Coriolis parameter and zeta_850 is the 850-hPa
    relative vorticity.  Values are clipped to ±cap and only positive
    values are retained.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with wind components on pressure levels.
    u_var, v_var : str
        Variable names for zonal and meridional wind.
    cap : float
        Upper bound on absolute vorticity (default 3.7e-5 s^-1).

    Returns
    -------
    xr.DataArray
        Clipped absolute vorticity [s^-1].
    """
    u = _drop_scalar_level(ds[u_var].sel(level=850))
    v = _drop_scalar_level(ds[v_var].sel(level=850))

    lat_rad = np.deg2rad(ds["latitude"])

    # Relative vorticity on the sphere
    dvdx = v.differentiate("longitude") / (
        np.deg2rad(1.0) * EARTH_RADIUS * np.cos(lat_rad)
    )
    dudy = u.differentiate("latitude") / (np.deg2rad(1.0) * EARTH_RADIUS)
    rel_vort = dvdx - dudy

    # Coriolis parameter
    coriolis = 2.0 * EARTH_OMEGA * np.sin(lat_rad)

    eta = rel_vort + coriolis

    # Clip: set |eta| > cap to sign(eta) * cap
    eta_c = xr.where(np.abs(eta) > cap, np.sign(eta) * cap, eta)
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

    Returns
    -------
    xr.DataArray
        Ventilated genesis potential index.
    """
    cos_lat = np.cos(np.deg2rad(vpi["latitude"]))
    gpi_v = (102.1 * vpi * eta_c) ** 4.90 * cos_lat * dx * dy
    gpi_v = gpi_v.where(np.isfinite(gpi_v) & (gpi_v > 0))
    gpi_v.attrs = {
        "long_name": "Ventilated Genesis Potential Index",
        "units": "arbitrary",
    }
    return gpi_v
