#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for Liang and Granger causality analysis
https://github.com/LinkedEarth/Pyleoclim_util
"""

__all__ = ["liang_causality", "granger_causality"]
from statsmodels.tsa.arima_process import arma_generate_sample
from statsmodels.tsa.arima.model import ARIMA
import numpy as np
from statsmodels.tsa.stattools import grangercausalitytests
from tqdm import tqdm
from scipy.stats.mstats import mquantiles


def ar1_fit_evenly(y):
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


def sm_ar1_sim(n, p, g, sig):
    """Produce p realizations of an AR1 process of length n with lag-1 autocorrelation g using statsmodels

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

    Returns
    -------
    red : numpy matrix
        n rows by p columns matrix of an AR1 process

    See also
    --------

    pyleoclim.utils.correlation.corr_sig : Estimates the Pearson's correlation and associated significance between two non IID time series
    pyleoclim.utils.correlation.fdr : Determine significance based on the false discovery rate

    """
    # specify model parameters (statsmodel wants lag0 coefficents as unity)
    ar = np.r_[1, -g]  # AR model parameter
    ma = np.r_[1, 0.0]  # MA model parameters
    # theoretical noise variance for red to achieve the same variance as X
    sig_n = sig * np.sqrt(1 - g**2)

    red = np.empty(shape=(n, p))  # declare array

    # simulate AR(1) model for each column
    for i in np.arange(p):
        red[:, i] = arma_generate_sample(
            ar=ar, ma=ma, nsample=n, burnin=50, scale=sig_n
        )

    return red


# -------
# Main functions
# --------


def granger_causality(y1, y2, maxlag=1, addconst=True, verbose=True):
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

    pyleoclim.utils.causality.liang_causality : information flow estimated using the Liang algorithm

    pyleoclim.utils.causality.signif_isopersist : significance test with AR(1) with same persistence

    pyleoclim.utils.causality.signif_isospec : significance test with surrogates with randomized phases

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


def phaseran(recblk, nsurr):
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

    pyleoclim.utils.correlation.corr_sig : Estimates the Pearson's correlation and associated significance between two non IID time series
    pyleoclim.utils.correlation.fdf : Determine significance based on the false discovery rate

    References
    ----------

    - Prichard, D., Theiler, J. Generating Surrogate Data for Time Series with Several Simultaneously Measured Variables (1994) Physical Review Letters, Vol 73, Number 7

    - Carlos Gias (2020). Phase randomization, MATLAB Central File Exchange
    """
    # Get parameters
    nfrms = recblk.shape[0]

    if nfrms % 2 == 0:
        nfrms = nfrms - 1
        recblk = recblk[0:nfrms]

    len_ser = int((nfrms - 1) / 2)
    interv1 = np.arange(1, len_ser + 1)
    interv2 = np.arange(len_ser + 1, nfrms)

    # Fourier transform of the original dataset
    fft_recblk = np.fft.fft(recblk)

    surrblk = np.zeros((nfrms, nsurr))

    #  for k in tqdm(np.arange(nsurr)):
    for k in np.arange(nsurr):
        ph_rnd = np.random.rand(len_ser)

        # Create the random phases for all the time series
        ph_interv1 = np.exp(2 * np.pi * 1j * ph_rnd)
        ph_interv2 = np.conj(np.flipud(ph_interv1))

        # Randomize all the time series simultaneously
        fft_recblk_surr = np.copy(fft_recblk)
        fft_recblk_surr[interv1] = fft_recblk[interv1] * ph_interv1
        fft_recblk_surr[interv2] = fft_recblk[interv2] * ph_interv2

        # Inverse transform
        surrblk[:, k] = np.real(np.fft.ifft(fft_recblk_surr))

    return surrblk


def liang_causality(
    y1,
    y2,
    npt=1,
    signif_test="isospec",
    nsim=1000,
    qs=[0.005, 0.025, 0.05, 0.95, 0.975, 0.995],
):
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
    pyleoclim.utils.causality.liang : information flow estimated using the Liang algorithm
    pyleoclim.utils.causality.granger_causality : information flow estimated using the Granger algorithm
    pyleoclim.utils.causality.signif_isopersist : significance test with AR(1) with same persistence
    pyleoclim.utils.causality.causality.signif_isospec : significance test with surrogates with randomized phases

    References
    ----------
    Liang, X.S. (2013) The Liang-Kleeman Information Flow: Theory and Applications. Entropy, 15, 327-360, doi:10.3390/e15010327

    Liang, X.S. (2014) Unraveling the cause-effect relation between timeseries. Physical review, E 90, 052150

    Liang, X.S. (2015) Normalizing the causality between time series. Physical review, E 92, 022126

    Liang, X.S. (2016) Information flow and causality as rigorous notions ab initio. Physical review, E 94, 052201

    """

    dt = 1
    nm = np.size(y1)

    grad1 = (y1[0 + npt :] - y1[0:-npt]) / (npt)
    grad2 = (y2[0 + npt :] - y2[0:-npt]) / (npt)

    y1 = y1[:-npt]
    y2 = y2[:-npt]

    N = nm - npt
    C = np.cov(y1, y2)
    detC = np.linalg.det(C)

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

    signif_test_func = {
        "isopersist": signif_isopersist,
        "isospec": signif_isospec,
    }

    signif_dict = signif_test_func[signif_test](
        y1, y2, method="liang", nsim=nsim, qs=qs, npt=npt
    )
    T21_noise_qs = signif_dict["T21_noise_qs"]
    tau21_noise_qs = signif_dict["tau21_noise_qs"]

    res = {
        "T21": T21,
        "tau21": tau21,
        "Z": Z,
        "dH1_star": dH1_star,
        "dH1_noise": dH1_noise,
        "signif_qs": qs,
        "T21_noise": T21_noise_qs,
        "tau21_noise": tau21_noise_qs,
    }

    return res


def liang(y1, y2, npt=1):
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

    pyleoclim.utils.causality.liang_causality : information flow estimated using the Liang algorithm
    pyleoclim.utils.causality.granger_causality : information flow estimated using the Granger algorithm
    pyleoclim.utils.causality.signif_isopersist : significance test with AR(1) with same persistence
    pyleoclim.utils.causality.signif_isospec : significance test with surrogates with randomized phases

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
    dt = 1
    nm = np.size(y1)

    grad1 = (y1[0 + npt :] - y1[0:-npt]) / (npt)
    grad2 = (y2[0 + npt :] - y2[0:-npt]) / (npt)

    y1 = y1[:-npt]
    y2 = y2[:-npt]

    N = nm - npt
    C = np.cov(y1, y2)
    detC = np.linalg.det(C)

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


def signif_isopersist(
    y1, y2, method, nsim=1000, qs=[0.005, 0.025, 0.05, 0.95, 0.975, 0.995], **kwargs
):
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

    pyleoclim.utils.causality.liang_causality : information flow estimated using the Liang algorithm
    pyleoclim.utils.causality.granger_causality : information flow estimated using the Granger algorithm
    pyleoclim.utils.causality.signif_isospec : significance test with surrogates with randomized phases

    """
    g1 = ar1_fit_evenly(y1)
    g2 = ar1_fit_evenly(y2)
    sig1 = np.std(y1)
    sig2 = np.std(y2)
    n = np.size(y1)
    noise1 = sm_ar1_sim(n, nsim, g1, sig1)
    noise2 = sm_ar1_sim(n, nsim, g2, sig2)

    if method == "liang":
        npt = kwargs["npt"] if "npt" in kwargs else 1
        T21_noise = []
        tau21_noise = []
        for i in tqdm(range(nsim), desc="Calculating causality between surrogates"):
            res_noise = liang(noise1[:, i], noise2[:, i], npt=npt)
            tau21_noise.append(res_noise["tau21"])
            T21_noise.append(res_noise["T21"])
        tau21_noise = np.array(tau21_noise)
        T21_noise = np.array(T21_noise)
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
    y1, y2, method, nsim=1000, qs=[0.005, 0.025, 0.05, 0.95, 0.975, 0.995], **kwargs
):
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

    pyleoclim.utils.causality.liang_causality : information flow estimated using the Liang algorithm
    pyleoclim.utils.causality.granger_causality : information flow estimated using the Granger algorithm
    pyleoclim.utils.causality.signif_isopersist : significance test with AR(1) with same persistence

    """

    noise1 = phaseran(y1, nsim)
    noise2 = phaseran(y2, nsim)

    if method == "liang":
        npt = kwargs["npt"] if "npt" in kwargs else 1
        T21_noise = []
        tau21_noise = []
        for i in tqdm(range(nsim), desc="Calculating causality between surrogates"):
            res_noise = liang(noise1[:, i], noise2[:, i], npt=npt)
            tau21_noise.append(res_noise["tau21"])
            T21_noise.append(res_noise["T21"])
        tau21_noise = np.array(tau21_noise)
        T21_noise = np.array(T21_noise)
        tau21_noise_qs = mquantiles(tau21_noise, qs)
        T21_noise_qs = mquantiles(T21_noise, qs)

        res_dict = {
            "tau21_noise_qs": tau21_noise_qs,
            "T21_noise_qs": T21_noise_qs,
        }
    # TODO Recode with Surrogate class
    else:
        raise KeyError(f"{method} is not a valid method")

    return res_dict
