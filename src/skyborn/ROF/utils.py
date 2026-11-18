# import os
from typing import Any, Dict, List, Tuple, Union

import numpy as np

# import pandas as pd
# from glob import glob
from scipy import stats

__all__ = [
    "speco",
    "chi2_test",
    "project_vectors",
    "unproject_vectors",
    "SSM",
    "get_nruns",
    "Cm_estimate",
    "Cv_estimate",
]


def speco(C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    This function computes eigenvalues and eigenvectors, in descending order
    :param C: numpy.ndarray
        A p x p symmetric real matrix
    :return:
    eigenvectors: numpy.ndarray
        The eigenvectors (eigenvectors[:, i] is the i-th eigenvector)
    eigenvalues_diag: numpy.ndarray
        The eigenvalues as a diagonal matrix
    """
    C_arr = np.asarray(C)
    if C_arr.ndim != 2 or C_arr.shape[0] != C_arr.shape[1]:
        raise ValueError("Matrix must be square")

    # Verify matrix symmetry before eigendecomposition
    if np.iscomplexobj(C_arr):
        imag_scale = max(float(np.max(np.abs(C_arr))), 1.0)
        if float(np.max(np.abs(np.imag(C_arr)))) > 1e-12 * imag_scale:
            raise ValueError("Matrix is not symmetric")
        C_real = np.real(C_arr).astype(np.float64, copy=False)
    else:
        C_real = C_arr.astype(np.float64, copy=False)

    if not np.allclose(C_real, C_real.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Matrix is not symmetric")

    # Symmetrize floating-point noise and use eigh for an orthonormal basis
    C_sym = 0.5 * (C_real + C_real.T)
    eigenvals_asc, eigenvecs_asc = np.linalg.eigh(C_sym)

    # Sort eigenvalues and eigenvectors in descending order
    eigenvectors = eigenvecs_asc[:, ::-1]
    eigenvalues_diag = np.diag(eigenvals_asc[::-1])

    return eigenvectors, eigenvalues_diag


def chi2_test(d_cons: float, df: int) -> float:
    """
    Check whether it is from a chi-squared distribution or not
    :param d_cons: float
        -2 log-likelihood
    :param df: int
        Degrees of freedom
    :return:
    pv_cons: float
        p-value for the test
    """
    rien = stats.chi2.cdf(d_cons, df=df)
    pv_cons = 1.0 - rien

    return pv_cons


def project_vectors(nt: int, X: np.ndarray) -> np.ndarray:
    """
    This function provides a projection matrix U that can be applied to X to ensure its covariance matrix to be
    full-ranked. Projects to a nt-1 subspace (ref: Ribes et al., 2013).
    :param nt: int
        number of time steps
    :param X: numpy.ndarray
        nt x nf array to be projected
    :return:
    projected_X: numpy.ndarray
        nt - 1 x nf array of projected timeseries
    """
    # Create centering matrix M = I - (1/nt) * ones_matrix
    # This removes the mean from each time series
    centering_matrix = np.eye(nt, nt) - np.ones((nt, nt)) / nt

    # Compute eigen-decomposition of centering matrix
    # Note: rank(M) = nt-1, so M has one eigenvalue equal to 0
    eigenvectors, _ = speco(centering_matrix)

    # Select first (nt-1) eigenvectors corresponding to non-zero eigenvalues
    # These form the projection matrix U
    projection_matrix = eigenvectors[:, : nt - 1].T

    # Apply projection to input data
    projected_X = np.dot(projection_matrix, X)

    return projected_X


def unproject_vectors(nt: int, Xc: np.ndarray) -> np.ndarray:
    """
    Unproject data from (nt-1) subspace back to nt subspace for trend computation.

    :param nt: int
        Number of time steps
    :param Xc: numpy.ndarray
        Projected data matrix (nt-1 x nf) to be unprojected
    :return:
    numpy.ndarray
        Unprojected data matrix (nt x nf)
    """
    # Create centering matrix: I - (1/nt)*ones
    centering_matrix = np.eye(nt, nt) - np.ones((nt, nt)) / nt

    # Get eigenvectors
    eigenvectors, _ = speco(centering_matrix)

    # Transpose of orthonormal eigenvectors is its inverse on the reduced subspace
    inverse_projection = eigenvectors[:, : nt - 1]

    # Apply inverse projection
    unprojected_X = np.dot(inverse_projection, Xc)

    return unprojected_X


def SSM(X_dict: Dict[str, np.ndarray], X_mm: np.ndarray) -> np.ndarray:
    """
    Calculates the squared difference between each models ensemble mean and the multi-model mean. Based on
    (Ribes et al., 2017)
    :param X_dict: dict
        Dictionary where keys are experiment names and values are arrays (n_members, n_time)
    :param X_mm: numpy.ndarray
        Multi-model ensemble mean, shape (n_time,)
    :return:
    D @ D.T: numpy.ndarray
        (nt - 1) x (nt - 1) full inter-model covariance matrix (sum of outer products)
    """
    # Make sure X_mm has right shape
    if X_mm.ndim == 1:
        X_mm = X_mm.reshape((len(X_mm), 1))

    # Calculate ensemble mean for each experiment
    exp_means = []
    for exp_name, exp_data in X_dict.items():
        # Get ensemble mean for this experiment
        ensemble_mean = np.mean(exp_data, axis=0)  # shape: (n_time,)
        exp_means.append(ensemble_mean)

    # Stack all experiment means: (n_time, n_experiments)
    X = np.column_stack(exp_means)

    # Apply projection (this is the default behavior in original SSM)
    Xc = project_vectors(X.shape[0], X)
    Xc_mm = project_vectors(X.shape[0], X_mm)

    # Full inter-model covariance (sum of outer products), as in Ribes et al. (2017).
    D = Xc - Xc_mm
    return D @ D.T

def get_nruns(X_dict: Dict[str, np.ndarray]) -> np.ndarray:
    """
    Gets the number of runs for each CESM2 experiment

    :param X_dict: dict
        Dictionary where keys are experiment names (e.g., 'CESM2-GHG', 'CESM2-EE')
        and values are arrays with shape (n_members, n_time)
    :return:
    nruns: numpy.ndarray
        Array with number of runs for each experiment
    """
    nruns = []
    for exp_name, exp_data in X_dict.items():
        nruns.append(exp_data.shape[0])  # number of members
    return np.array(nruns)


def Cm_estimate(
    X_dict: Dict[str, np.ndarray], Cv: np.ndarray, X_mm: np.ndarray
) -> np.ndarray:
    """
    Estimated covariance matrix for model error (Ribes et al., 2017)
    Modified for CESM2 experiments with ensemble member arrays

    :param X_dict: dict
        Dictionary where keys are experiment names (e.g., 'CESM2-GHG', 'CESM2-EE')
        and values are arrays with shape (n_members, n_time)
    :param Cv: numpy.ndarray
        Array with internal variability covariance matrix
    :param X_mm: numpy.ndarray
        Array with multi-model ensemble mean
    :return:
    Cm_pos_hat: numpy.ndarray
        Estimated covariance matrix for model error
    """

    # Calculate model differences using our modified SSM function
    _SSM = SSM(X_dict, X_mm)

    # Get number of runs and number of models using our modified function
    nruns = get_nruns(X_dict)
    nm = len(nruns)  # number of experiments

    # Calculate Cv_all based on number of runs for each experiment
    Cv_all = np.zeros(Cv.shape)
    for nr in nruns:
        Cv_all += Cv / nr

    if nm < 2:
        raise ValueError(
            "At least two models (nm >= 2) are required to estimate Cm."
        )

    # First estimation of Cm
    Cm_hat = (1.0 / (nm - 1.0)) * (_SSM - ((nm - 1.0) / nm) * Cv_all)
    Cm_sym = 0.5 * (Cm_hat + Cm_hat.T)

    # Set negative eigenvalues to zero and recompose the signal using eigh
    # so rank-deficient symmetric covariance matrices stay real and symmetric.
    S, X = np.linalg.eigh(Cm_sym)
    S = np.maximum(S, 0.0)
    Cm_pos_hat = (X * S) @ X.T  # spectral decomposition

    Cm_pos_hat = (1.0 + (1.0 / nm)) * Cm_pos_hat

    return Cm_pos_hat


def Cv_estimate(X_dict: Dict[str, np.ndarray], Cv: np.ndarray) -> np.ndarray:
    """
    Estimated covariance matrix for internal variability considering multiple models (Ribes et al., 2017)
    Modified for CESM2 experiments with ensemble member arrays

    :param X_dict: dict
        Dictionary where keys are experiment names and values are arrays (n_members, n_time)
    :param Cv: numpy.ndarray
        Array with internal variability covariance matrix
    :return:
    Cv_estimate: numpy.ndarray
        Estimated covariance matrix for internal variability considering multiple models
    """
    # Get number of runs and number of models
    nruns = get_nruns(X_dict)
    nm = len(nruns)

    Cv_all = np.zeros(Cv.shape)
    for nr in nruns:
        Cv_all += Cv / nr

    Cv_estimate = (1.0 / (nm**2.0)) * Cv_all

    return Cv_estimate


if __name__ == "__main__":
    T = 11
    M = np.eye(T, T) - np.ones((T, T)) / T
    speco(M)
