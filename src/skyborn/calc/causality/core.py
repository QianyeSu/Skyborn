#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for Liang and Granger causality analysis in atmospheric and climate data

This module provides causality analysis tools for time series data commonly
used in atmospheric and climate research.

Utilities for Liang and Granger causality analysis
https://github.com/LinkedEarth/Pyleoclim_util
"""

__all__ = [
    "ar1_fit_evenly",
    "sm_ar1_sim",
    "Surrogate",
    "granger_causality",
    "phaseran",
    "liang_causality",
    "liang",
    "signif_isopersist",
    "signif_isospec",
]
import sys
from importlib import import_module
from importlib import util as importlib_util
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
from types import ModuleType
from typing import Any, Dict, List, Optional, Union

import numpy as np
from scipy import fft as scipy_fft
from scipy.stats.mstats import mquantiles
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.arima_process import arma_generate_sample
from statsmodels.tsa.stattools import grangercausalitytests


def _load_liang_backend() -> Optional[ModuleType]:
    """Load the compiled Liang backend from source or an installed package."""
    backend_dir = Path(__file__).resolve().parent
    module_name = f"{__package__}._causality_core"

    for probe_dir in (backend_dir / "build", backend_dir):
        for suffix in EXTENSION_SUFFIXES:
            candidate = probe_dir / f"_causality_core{suffix}"
            if not candidate.exists():
                continue

            spec = importlib_util.spec_from_file_location(module_name, candidate)
            if spec is None or spec.loader is None:
                continue
            try:
                module = importlib_util.module_from_spec(spec)
                spec.loader.exec_module(module)
                sys.modules.setdefault(module_name, module)
                return module
            except Exception:
                continue

    try:
        return import_module(module_name)
    except Exception:
        return None


_LIANG_BACKEND = _load_liang_backend()


def _liang_backend_available() -> bool:
    """Return whether the compiled Liang backend is importable."""
    return _LIANG_BACKEND is not None


def _ar1_native_available() -> bool:
    """Return whether the compiled AR(1) filtering entry point is available."""
    return _LIANG_BACKEND is not None and hasattr(_LIANG_BACKEND, "ar1_filter_batch")


def _require_liang_backend() -> ModuleType:
    """Return the compiled Liang backend or raise a build error."""
    if _LIANG_BACKEND is None:
        raise ImportError(
            "The compiled Skyborn causality backend is unavailable. "
            "Build src/skyborn/calc/causality with Meson first."
        )
    return _LIANG_BACKEND


def _liang_single_backend(
    y1: np.ndarray, y2: np.ndarray, npt: int
) -> tuple[float, float, float, float, float]:
    """Run the compiled Liang kernel for one pair of time series."""
    backend = _require_liang_backend()
    return tuple(
        float(value)
        for value in backend.liang_single(
            np.asarray(y1, dtype=np.float64),
            np.asarray(y2, dtype=np.float64),
            int(npt),
        )
    )


def _liang_batch_backend(
    y1: np.ndarray, y2: np.ndarray, npt: int
) -> tuple[np.ndarray, np.ndarray]:
    """Run the compiled Liang kernel over surrogate columns."""
    backend = _require_liang_backend()
    t21, tau21 = backend.liang_batch(
        np.asfortranarray(y1, dtype=np.float64),
        np.asfortranarray(y2, dtype=np.float64),
        int(npt),
    )
    return np.asarray(t21, dtype=np.float64), np.asarray(tau21, dtype=np.float64)


def ar1_fit_evenly(y: np.ndarray) -> float:
    """Returns the lag-1 autocorrelation from AR(1) fit.

    Uses `statsmodels.tsa.arima.model.ARIMA <https://www.statsmodels.org/devel/generated/statsmodels.tsa.arima.model.ARIMA.html>`_. to
    calculate lag-1 autocorrelation

    MARK FOR DEPRECATION once uar1_fit is adopted

    Parameters
    ----------
    y : array
        Vector of (float) numbers as a time series

    Returns
    -------
    g : float
        Lag-1 autocorrelation coefficient

    """
    # syntax compatible with statsmodels v0.11.1
    # ar1_mod = sm.tsa.ARMA(y, (1, 0), missing='drop').fit(trend='nc', disp=0)
    # syntax compatible with statsmodels v0.12
    ar1_mod = ARIMA(y, order=(1, 0, 0), missing="drop", trend="ct").fit()
    g = ar1_mod.params[2]

    if g > 1:
        print(
            "Warning: AR(1) fitted autocorrelation greater than 1; setting to 1-eps^{1/4}"
        )
        eps = np.spacing(1.0)
        g = 1.0 - eps ** (1 / 4)

    return g


def _ar1_filter_native(innovations: np.ndarray, g: float, nout: int) -> np.ndarray:
    """Filter AR(1) innovations through the compiled Fortran backend."""
    backend = _require_liang_backend()
    if not hasattr(backend, "ar1_filter_batch"):
        raise ImportError("The compiled AR(1) filtering backend is unavailable")
    return np.asarray(
        backend.ar1_filter_batch(
            np.asfortranarray(innovations, dtype=np.float64),
            float(g),
            int(nout),
        ),
        dtype=np.float64,
    )


def sm_ar1_sim(
    n: int,
    p: int,
    g: float,
    sig: float,
    *,
    backend: str = "auto",
) -> np.ndarray:
    """Produce p realizations of an AR1 process of length n.

    Parameters
    ----------
    n : int
        row dimensions
    p : int
        column dimensions

    g : float
        lag-1 autocorrelation coefficient

    sig : float
        the standard deviation of the original time series
    backend : {"auto", "python", "fortran", "native"}, optional
        ``"auto"`` uses the Fortran/OpenMP filter for sufficiently large
        surrogate ensembles when it is available. ``"python"`` preserves the
        statsmodels implementation. ``"fortran"`` and ``"native"`` require
        the compiled filter.

    Returns
    -------
    red : numpy matrix
        n rows by p columns matrix of an AR1 process

    See also
    --------

    skyborn.causality.granger_causality : Granger causality analysis
    skyborn.causality.liang_causality : Liang information flow analysis

    """
    if backend not in {"auto", "python", "fortran", "native"}:
        raise ValueError(
            "backend must be one of: 'auto', 'python', 'fortran', 'native'"
        )
    if n < 1 or p < 0:
        raise ValueError("n must be positive and p must be non-negative")

    # Specify model parameters (statsmodels wants lag-0 coefficients as unity).
    ar = np.r_[1, -g]  # AR model parameter
    ma = np.r_[1, 0.0]  # MA model parameters
    # theoretical noise variance for red to achieve the same variance as X
    sig_n = sig * np.sqrt(1 - g**2)

    if p == 0:
        return np.empty(shape=(n, 0), order="F")

    burnin = 50
    native_available = _ar1_native_available()
    if backend in {"fortran", "native"}:
        if not native_available:
            raise ImportError(
                "The compiled AR(1) filtering backend is unavailable. "
                "Build src/skyborn/calc/causality with Meson first."
            )
        use_native = True
    else:
        use_native = backend == "auto" and p >= 32 and native_available

    if use_native:
        # Draw innovations in the same row-major order as the historical
        # per-column statsmodels loop, then let Fortran perform burn-in and
        # recurrence in parallel across surrogate columns.
        innovations = np.asfortranarray(
            np.random.normal(scale=sig_n, size=(p, n + burnin)).T
        )
        return _ar1_filter_native(innovations, g, n)

    # The Python path remains the exact compatibility implementation.
    return np.asfortranarray(
        arma_generate_sample(
            ar=ar,
            ma=ma,
            nsample=(p, n),
            axis=1,
            burnin=burnin,
            scale=sig_n,
        ).T
    )


# -------
# Main functions
# --------


def granger_causality(
    y1: np.ndarray,
    y2: np.ndarray,
    maxlag: Union[int, List[int]] = 1,
    addconst: bool = True,
    verbose: bool = True,
) -> Dict[Any, Any]:
    """Granger causality tests

    Four tests for the Granger non-causality of 2 time series.

    All four tests give similar results. params_ftest and ssr_ftest are equivalent based on F test which is identical to lmtest:grangertest in R.

    Wrapper for the functions described in statsmodels (https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.grangercausalitytests.html)

    Parameters
    ----------
    y1, y2: array
        vectors of (real) numbers with identical length, no NaNs allowed
    maxlag : int or int iterable, optional
        If an integer, computes the test for all lags up to maxlag. If an iterable, computes the tests only for the lags in maxlag.
    addconst : bool, optional
        Include a constant in the model.
    verbose : bool, optional
        Print results

    Returns
    -------
    dict
        All test results, dictionary keys are the number of lags. For each lag the values are a tuple, with the first element a dictionary with test statistic,
        pvalues, degrees of freedom, the second element are the OLS estimation results for the restricted model, the unrestricted model and the restriction (contrast)
        matrix for the parameter f_test.

    Notes
    -----

    The null hypothesis for Granger causality tests is that y2, does NOT Granger cause y1. Granger causality means that past values of y2 have a statistically significant effect on the current value of y1, taking past values of y1 into account as regressors. We reject the null hypothesis that y2 does not Granger cause y1 if the p-values are below a desired threshold (e.g. 0.05).

    The null hypothesis for all four test is that the coefficients corresponding to past values of the second time series are zero.

    ‘params_ftest’, ‘ssr_ftest’ are based on the F distribution

    ‘ssr_chi2test’, ‘lrtest’ are based on the chi-square distribution

    See also
    --------

    skyborn.causality.liang_causality : Information flow estimated using the Liang algorithm

    skyborn.causality.signif_isopersist : Significance test with AR(1) with same persistence

    skyborn.causality.signif_isospec : Significance test with surrogates with randomized phases

    References
    ----------

    Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. Econometrica, 37(3), 424-438.

    Granger, C. W. J. (1980). Testing for causality: A personal viewpoont. Journal of Economic Dynamics and Control, 2, 329-352.

    Granger, C. W. J. (1988). Some recent development in a concept of causality. Journal of Econometrics, 39(1-2), 199-211.

    """

    if len(y1) != len(y2):
        raise ValueError("Timeseries must be of same length")

    x = np.array([y1, y2]).T
    res = grangercausalitytests(x, maxlag=maxlag, addconst=addconst, verbose=verbose)
    return res


def phaseran(recblk: np.ndarray, nsurr: int) -> np.ndarray:
    """Simultaneous phase randomization of a set of time series

    It creates blocks of surrogate data with the same second order properties as the original
    time series dataset by transforming the original data into the frequency domain, randomizing the
    phases simultaneoulsy across the time series and converting the data back into the time domain.

    Written by Carlos Gias for MATLAB

    http://www.mathworks.nl/matlabcentral/fileexchange/32621-phase-randomization/content/phaseran.m

    Parameters
    ----------
    recblk : numpy array
        2D array , Row: time sample. Column: recording.
        An odd number of time samples (height) is expected.
        If that is not the case, recblock is reduced by 1 sample before the surrogate data is created.
        The class must be double and it must be nonsparse.

    nsurr : int
        is the number of image block surrogates that you want to generate.

    Returns
    -------
    surrblk : numpy array
        3D multidimensional array image block with the surrogate datasey along the third dimension

    See also
    --------

    skyborn.causality.liang_causality : Liang-Kleeman information flow analysis
    skyborn.causality.granger_causality : Granger causality analysis

    References
    ----------

    - Prichard, D., Theiler, J. Generating Surrogate Data for Time Series with Several Simultaneously Measured Variables (1994) Physical Review Letters, Vol 73, Number 7

    - Carlos Gias (2020). Phase randomization, MATLAB Central File Exchange
    """
    recblk = np.asarray(recblk, dtype=np.float64)
    if recblk.ndim not in (1, 2):
        raise ValueError("recblk must be a 1D or 2D array")
    if nsurr < 1:
        raise ValueError("nsurr must be a positive integer")

    was_1d = recblk.ndim == 1
    if was_1d:
        recblk = recblk[:, np.newaxis]

    # Get parameters
    nfrms = recblk.shape[0]

    if nfrms % 2 == 0:
        nfrms = nfrms - 1
        recblk = recblk[0:nfrms]

    len_ser = int((nfrms - 1) / 2)

    # Generate all requested phases in bounded batches. The transpose keeps
    # NumPy's random stream grouped exactly as in the historical per-surrogate
    # loop: each surrogate consumes one contiguous phase vector.
    #
    # Real-valued input only needs the positive-frequency half spectrum. The
    # frequency axis is kept contiguous so scipy.fft can process each batch
    # without repeated Python-level FFT calls.
    fft_recblk = scipy_fft.rfft(recblk, axis=0)
    nseries = recblk.shape[1]
    if was_1d:
        surrblk = np.empty((nfrms, nsurr), dtype=np.float64, order="F")
    else:
        surrblk = np.empty((nfrms, nseries, nsurr), dtype=np.float64)
    batch_size = min(int(nsurr), 256)

    for start in range(0, int(nsurr), batch_size):
        stop = min(start + batch_size, int(nsurr))
        batch_count = stop - start
        ph_rnd = np.random.rand(batch_count, len_ser).T
        ph_interv1 = np.exp(2 * np.pi * 1j * ph_rnd)
        fft_recblk_surr = np.empty(
            (batch_count, nseries, len_ser + 1), dtype=np.complex128
        )
        fft_recblk_surr[...] = fft_recblk.T[np.newaxis, :, :]
        fft_recblk_surr[:, :, 1:] *= ph_interv1.T[:, np.newaxis, :]

        irfft_batch = np.real(scipy_fft.irfft(fft_recblk_surr, n=nfrms, axis=-1))
        if was_1d:
            surrblk[:, start:stop] = irfft_batch[:, 0, :].T
        else:
            surrblk[:, :, start:stop] = irfft_batch.transpose(2, 1, 0)

    if was_1d:
        return surrblk
    return surrblk


class Surrogate:
    """Generate null-model surrogate ensembles for causality analysis.

    Parameters
    ----------
    data : array-like
        One-dimensional time series for AR(1) or phase-randomized surrogates,
        or a two-dimensional ``(time, series)`` block for phase
        randomization.

    Notes
    -----
    The class is a small state-free facade over the existing generation
    functions. It keeps surrogate construction separate from the Liang
    significance calculation while preserving the historical function APIs.
    """

    def __init__(self, data: np.ndarray):
        data = np.asarray(data, dtype=np.float64)
        if data.ndim not in (1, 2):
            raise ValueError("Surrogate data must be a 1D or 2D array")
        if data.shape[0] < 1:
            raise ValueError("Surrogate data must contain at least one sample")
        if not np.all(np.isfinite(data)):
            raise ValueError("Surrogate data must contain only finite values")
        self.data = data

    def phase_randomized(self, nsurr: int) -> np.ndarray:
        """Return phase-randomized surrogates for the stored data."""
        return phaseran(self.data, nsurr)

    def ar1(
        self,
        nsurr: int,
        *,
        g: Optional[float] = None,
        sig: Optional[float] = None,
        backend: str = "auto",
    ) -> np.ndarray:
        """Return AR(1) surrogates fitted to the stored one-dimensional data."""
        if self.data.ndim != 1:
            raise ValueError("AR(1) surrogates require one-dimensional data")
        if g is None:
            g = ar1_fit_evenly(self.data)
        if sig is None:
            sig = float(np.std(self.data))
        return sm_ar1_sim(
            self.data.size,
            nsurr,
            float(g),
            float(sig),
            backend=backend,
        )

    def generate(
        self,
        method: str,
        nsurr: int,
        *,
        backend: str = "auto",
        **kwargs: Any,
    ) -> np.ndarray:
        """Generate surrogates using a named null model.

        Accepted method names are ``"isospec"``, ``"phase_randomized"``,
        ``"isopersist"``, and ``"ar1"``.
        """
        normalized = method.lower()
        if normalized in {"isospec", "phase_randomized"}:
            if kwargs:
                unexpected = next(iter(kwargs))
                raise TypeError(f"Unexpected keyword argument: {unexpected}")
            return self.phase_randomized(nsurr)
        if normalized in {"isopersist", "ar1"}:
            return self.ar1(nsurr, backend=backend, **kwargs)
        raise KeyError(f"{method} is not a valid surrogate method")


def liang_causality(
    y1: np.ndarray,
    y2: np.ndarray,
    npt: int = 1,
    signif_test: str = "isospec",
    nsim: int = 1000,
    qs: List[float] = [0.005, 0.025, 0.05, 0.95, 0.975, 0.995],
    surrogate_backend: str = "auto",
) -> Dict[str, Any]:
    """Liang-Kleeman information flow

    Estimate the Liang information transfer from series y2 to series y1 with
    significance estimates using either an AR(1) tests with series with the same
    persistence or surrogates with randomized phases.

    Parameters
    ----------
    y1, y2 : array
        vectors of (real) numbers with identical length, no NaNs allowed
    npt : int >=1
        time advance in performing Euler forward differencing,
        e.g., 1, 2. Unless the series are generated with a highly chaotic deterministic system,
        npt=1 should be used
    signif_test : str; {'isopersist', 'isospec'}
        the method for significance test
        see signif_isospec and signif_isopersist for details.
    nsim : int
        the number of AR(1) surrogates for significance test
    qs : list
        the quantiles for significance test
    surrogate_backend : {"auto", "python", "fortran", "native"}, optional
        Backend used for AR(1) surrogate filtering when
        ``signif_test="isopersist"``.

    Returns
    -------
    res : dict
        A dictionary of results including:

        - T21 : float - information flow from y2 to y1 (Note: not y1 -> y2!)
        - tau21 : float - the standardized information flow from y2 to y1
        - Z : float - the total information flow from y2 to y1
        - dH1_star : float - dH*/dt (Liang, 2016)
        - dH1_noise : float
        - signif_qs : the quantiles for significance test
        - T21_noise : list - the quantiles of the information flow from noise2 to noise1 for significance testing
        - tau21_noise : list - the quantiles of the standardized information flow from noise2 to noise1 for significance testing

    See also
    --------
    skyborn.causality.liang : Information flow estimated using the Liang algorithm
    skyborn.causality.granger_causality : Information flow estimated using the Granger algorithm
    skyborn.causality.signif_isopersist : Significance test with AR(1) with same persistence
    skyborn.causality.signif_isospec : Significance test with surrogates with randomized phases

    References
    ----------
    Liang, X.S. (2013) The Liang-Kleeman Information Flow: Theory and Applications. Entropy, 15, 327-360, doi:10.3390/e15010327

    Liang, X.S. (2014) Unraveling the cause-effect relation between timeseries. Physical review, E 90, 052150

    Liang, X.S. (2015) Normalizing the causality between time series. Physical review, E 92, 022126

    Liang, X.S. (2016) Information flow and causality as rigorous notions ab initio. Physical review, E 94, 052201

    """

    y1 = np.asarray(y1, dtype=np.float64)
    y2 = np.asarray(y2, dtype=np.float64)
    result = liang(y1, y2, npt=npt)

    signif_test_func = {
        "isopersist": signif_isopersist,
        "isospec": signif_isospec,
    }

    if signif_test not in signif_test_func:
        raise KeyError(f"{signif_test} is not a valid significance test")

    # Preserve the historical surrogate length, which used the differenced
    # series after removing the final npt samples.
    y1_signif = y1[:-npt]
    y2_signif = y2[:-npt]
    signif_dict = signif_test_func[signif_test](
        y1_signif,
        y2_signif,
        method="liang",
        nsim=nsim,
        qs=qs,
        npt=npt,
        surrogate_backend=surrogate_backend,
    )
    T21_noise_qs = signif_dict["T21_noise_qs"]
    tau21_noise_qs = signif_dict["tau21_noise_qs"]

    res = {
        "T21": result["T21"],
        "tau21": result["tau21"],
        "Z": result["Z"],
        "dH1_star": result["dH1_star"],
        "dH1_noise": result["dH1_noise"],
        "signif_qs": qs,
        "T21_noise": T21_noise_qs,
        "tau21_noise": tau21_noise_qs,
    }

    return res


def _liang_python(y1: np.ndarray, y2: np.ndarray, npt: int = 1) -> Dict[str, float]:
    """
    Estimate the Liang information transfer from series y2 to series y1

    Parameters
    ----------
    y1, y2 : array
        Vectors of (real) numbers with identical length, no NaNs allowed

    npt : int  >=1
        Time advance in performing Euler forward differencing,
        e.g., 1, 2. Unless the series are generated with a highly chaotic deterministic system,
        npt=1 should be used

    Returns
    -------
    res : dict
        A dictionary of results including:

            - T21 (float): information flow from y2 to y1 (Note: not y1 -> y2!)
            - tau21 (float): the standardized information flow from y2 to y1
            - Z (float): the total information flow from y2 to y1
            - dH1_star (float): dH*/dt (Liang, 2016)
            - dH1_noise (float)

    See also
    --------

    skyborn.causality.liang_causality : Information flow estimated using the Liang algorithm
    skyborn.causality.granger_causality : Information flow estimated using the Granger algorithm
    skyborn.causality.signif_isopersist : Significance test with AR(1) with same persistence
    skyborn.causality.signif_isospec : Significance test with surrogates with randomized phases

    References
    ----------

    Liang, X.S. (2013) The Liang-Kleeman Information Flow: Theory and
            Applications. Entropy, 15, 327-360, doi:10.3390/e15010327

    Liang, X.S. (2014) Unraveling the cause-effect relation between timeseries.
        Physical review, E 90, 052150

    Liang, X.S. (2015) Normalizing the causality between time series.
        Physical review, E 92, 022126

    Liang, X.S. (2016) Information flow and causality as rigorous notions ab initio.
        Physical review, E 94, 052201

    """
    y1 = np.asarray(y1, dtype=np.float64)
    y2 = np.asarray(y2, dtype=np.float64)
    if y1.ndim != 1 or y2.ndim != 1:
        raise ValueError("Liang inputs must be 1D arrays")
    if y1.shape != y2.shape:
        raise ValueError("Liang inputs must share a length")
    if npt < 1 or y1.size - npt < 2:
        raise ValueError("npt must leave at least two differenced samples")
    if not np.all(np.isfinite(y1)) or not np.all(np.isfinite(y2)):
        raise ValueError("Liang inputs must contain only finite values")

    dt = 1
    nm = np.size(y1)

    grad1 = (y1[0 + npt :] - y1[0:-npt]) / (npt)
    grad2 = (y2[0 + npt :] - y2[0:-npt]) / (npt)

    y1 = y1[:-npt]
    y2 = y2[:-npt]

    N = nm - npt
    C = np.cov(y1, y2)
    detC = np.linalg.det(C)
    if not np.isfinite(detC) or C[0, 0] == 0 or detC == 0:
        raise ValueError("Liang information flow is undefined for singular covariance")

    dC = np.ndarray((2, 2))
    dC[0, 0] = np.sum((y1 - np.mean(y1)) * (grad1 - np.mean(grad1)))
    dC[0, 1] = np.sum((y1 - np.mean(y1)) * (grad2 - np.mean(grad2)))
    dC[1, 0] = np.sum((y2 - np.mean(y2)) * (grad1 - np.mean(grad1)))
    dC[1, 1] = np.sum((y2 - np.mean(y2)) * (grad2 - np.mean(grad2)))

    dC /= N - 1

    a11 = C[1, 1] * dC[0, 0] - C[0, 1] * dC[1, 0]
    a12 = -C[0, 1] * dC[0, 0] + C[0, 0] * dC[1, 0]

    a11 /= detC
    a12 /= detC

    f1 = np.mean(grad1) - a11 * np.mean(y1) - a12 * np.mean(y2)
    R1 = grad1 - (f1 + a11 * y1 + a12 * y2)
    Q1 = np.sum(R1 * R1)
    b1 = np.sqrt(Q1 * dt / N)

    NI = np.ndarray((4, 4))
    NI[0, 0] = N * dt / b1**2
    NI[1, 1] = dt / b1**2 * np.sum(y1 * y1)
    NI[2, 2] = dt / b1**2 * np.sum(y2 * y2)
    NI[3, 3] = 3 * dt / b1**4 * np.sum(R1 * R1) - N / b1**2
    NI[0, 1] = dt / b1**2 * np.sum(y1)
    NI[0, 2] = dt / b1**2 * np.sum(y2)
    NI[0, 3] = 2 * dt / b1**3 * np.sum(R1)
    NI[1, 2] = dt / b1**2 * np.sum(y1 * y2)
    NI[1, 3] = 2 * dt / b1**3 * np.sum(R1 * y1)
    NI[2, 3] = 2 * dt / b1**3 * np.sum(R1 * y2)

    NI[1, 0] = NI[0, 1]
    NI[2, 0] = NI[0, 2]
    NI[2, 1] = NI[1, 2]
    NI[3, 0] = NI[0, 3]
    NI[3, 1] = NI[1, 3]
    NI[3, 2] = NI[2, 3]

    invNI = np.linalg.pinv(NI)
    var_a12 = invNI[2, 2]
    T21 = C[0, 1] / C[0, 0] * (-C[1, 0] * dC[0, 0] + C[0, 0] * dC[1, 0]) / detC
    var_T21 = (C[0, 1] / C[0, 0]) ** 2 * var_a12

    dH1_star = a11
    dH1_noise = b1**2 / (2 * C[0, 0])

    Z = np.abs(T21) + np.abs(dH1_star) + np.abs(dH1_noise)

    tau21 = T21 / Z
    dH1_star = dH1_star / Z
    dH1_noise = dH1_noise / Z

    res = {
        "T21": T21,
        "tau21": tau21,
        "Z": Z,
        "dH1_star": dH1_star,
        "dH1_noise": dH1_noise,
    }

    return res


def liang(y1: np.ndarray, y2: np.ndarray, npt: int = 1) -> Dict[str, float]:
    """Estimate Liang information flow using the compiled backend when available."""
    y1 = np.asarray(y1, dtype=np.float64)
    y2 = np.asarray(y2, dtype=np.float64)
    if y1.ndim != 1 or y2.ndim != 1:
        raise ValueError("Liang inputs must be 1D arrays")
    if y1.shape != y2.shape:
        raise ValueError("Liang inputs must share a length")
    if npt < 1 or y1.size - npt < 2:
        raise ValueError("npt must leave at least two differenced samples")
    if not np.all(np.isfinite(y1)) or not np.all(np.isfinite(y2)):
        raise ValueError("Liang inputs must contain only finite values")

    if _liang_backend_available():
        values = _liang_single_backend(y1, y2, int(npt))
        keys = ("T21", "tau21", "Z", "dH1_star", "dH1_noise")
        return dict(zip(keys, values))

    return _liang_python(y1, y2, npt=npt)


def _liang_batch(
    y1: np.ndarray, y2: np.ndarray, npt: int
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate T21 and tau21 for a matrix of surrogate series."""
    y1 = np.asarray(y1, dtype=np.float64)
    y2 = np.asarray(y2, dtype=np.float64)

    if _liang_backend_available():
        return _liang_batch_backend(y1, y2, int(npt))

    if y1.shape != y2.shape or y1.ndim != 2:
        raise ValueError("Surrogate arrays must be 2D arrays with the same shape")

    t21 = np.empty(y1.shape[1], dtype=np.float64)
    tau21 = np.empty(y1.shape[1], dtype=np.float64)
    for index in range(y1.shape[1]):
        result = liang(y1[:, index], y2[:, index], npt=npt)
        t21[index] = result["T21"]
        tau21[index] = result["tau21"]
    return t21, tau21


def signif_isopersist(
    y1: np.ndarray,
    y2: np.ndarray,
    method: str,
    nsim: int = 1000,
    qs: list[float] = [0.005, 0.025, 0.05, 0.95, 0.975, 0.995],
    surrogate_backend: str = "auto",
    **kwargs,
) -> dict[str, np.ndarray]:
    """significance test with AR(1) with same persistence

    Parameters
    ----------
    y1, y2 : array
        vectors of (real) numbers with identical length, no NaNs allowed

    method : str; {'liang'}
        estimates for the Liang method

    nsim : int
        the number of AR(1) surrogates for significance test

    qs : list
        the quantiles for significance test
    surrogate_backend : {"auto", "python", "fortran", "native"}, optional
        Backend used for AR(1) surrogate filtering.

    Returns
    -------
    res_dict : dict

        A dictionary with the following information:

          - T21_noise_qs : list
            the quantiles of the information flow from noise2 to noise1 for significance testing
          - tau21_noise_qs : list
            the quantiles of the standardized information flow from noise2 to noise1 for significance testing

    See also
    --------

    skyborn.causality.liang_causality : Information flow estimated using the Liang algorithm
    skyborn.causality.granger_causality : Information flow estimated using the Granger algorithm
    skyborn.causality.signif_isospec : Significance test with surrogates with randomized phases

    """
    g1 = ar1_fit_evenly(y1)
    g2 = ar1_fit_evenly(y2)
    sig1 = np.std(y1)
    sig2 = np.std(y2)
    n = np.size(y1)
    noise1 = Surrogate(y1).generate(
        "isopersist",
        nsim,
        backend=surrogate_backend,
        g=g1,
        sig=sig1,
    )
    noise2 = Surrogate(y2).generate(
        "isopersist",
        nsim,
        backend=surrogate_backend,
        g=g2,
        sig=sig2,
    )

    if method == "liang":
        npt = kwargs["npt"] if "npt" in kwargs else 1
        T21_noise, tau21_noise = _liang_batch(noise1, noise2, npt=npt)
        tau21_noise_qs = mquantiles(tau21_noise, qs)
        T21_noise_qs = mquantiles(T21_noise, qs)

        res_dict = {
            "tau21_noise_qs": tau21_noise_qs,
            "T21_noise_qs": T21_noise_qs,
        }
    # TODO add granger method
    else:
        raise KeyError(f"{method} is not a valid method")

    return res_dict


def signif_isospec(
    y1: np.ndarray,
    y2: np.ndarray,
    method: str,
    nsim: int = 1000,
    qs: list[float] = [0.005, 0.025, 0.05, 0.95, 0.975, 0.995],
    **kwargs,
) -> dict[str, np.ndarray]:
    """significance test with surrogates with randomized phases

    Parameters
    ----------
    y1, y2 : array
        vectors of (real) numbers with identical length, no NaNs allowed
    method : str; {'liang'}
        estimates for the Liang method
    nsim : int
        the number of surrogates for significance test
    qs : list
        the quantiles for significance test
    kwargs : dict
        keyword arguments for the causality method (e.g. npt for Liang-Kleeman)

    Returns
    -------
    res_dict : dict
        A dictionary with the following information:
          - T21_noise_qs : list
                        the quantiles of the information flow from noise2 to noise1 for significance testing
          - tau21_noise_qs : list
                          the quantiles of the standardized information flow from noise2 to noise1 for significance testing

    See also
    --------

    skyborn.causality.liang_causality : Information flow estimated using the Liang algorithm
    skyborn.causality.granger_causality : Information flow estimated using the Granger algorithm
    skyborn.causality.signif_isopersist : Significance test with AR(1) with same persistence

    """

    noise1 = Surrogate(y1).generate("isospec", nsim)
    noise2 = Surrogate(y2).generate("isospec", nsim)

    if method == "liang":
        npt = kwargs["npt"] if "npt" in kwargs else 1
        T21_noise, tau21_noise = _liang_batch(noise1, noise2, npt=npt)
        tau21_noise_qs = mquantiles(tau21_noise, qs)
        T21_noise_qs = mquantiles(T21_noise, qs)

        res_dict = {
            "tau21_noise_qs": tau21_noise_qs,
            "T21_noise_qs": T21_noise_qs,
        }
    else:
        raise KeyError(f"{method} is not a valid method")

    return res_dict
