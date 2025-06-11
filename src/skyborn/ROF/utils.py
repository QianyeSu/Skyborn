import os
import numpy as np
import pandas as pd
from glob import glob
from scipy import stats
from typing import Dict, List, Tuple, Union, Any


def speco(C: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    This function computes eigenvalues and eigenvectors, in descending order
    :param C: numpy.ndarray
        A p x p symetric real matrix
    :return:
    P: numpy.ndarray
        The eigenvectors (P[:, i] is the ist eigenvector)
    D: numpy.ndarray
        The eigenvalues as a diagonal matrix
    """
    # Compute eigenvalues and eigenvectors (the eigenvectors are non unique so the values may change from one software
    # to another e.g. python, matlab, scilab)
    D0, P0 = np.linalg.eig(C)

    # Take real part (to avoid numeric noise, eg small complex numbers)
    if np.max(np.imag(D0)) / np.max(np.real(D0)) > 1e-12:
        raise ValueError("Matrix is not symmetric")

    # Check that C is symetric (<=> real eigen-values/-vectors)
    P1 = np.real(P0)
    D1 = np.real(D0)

    # sort eigenvalues in descending order and
    # get their indices to order the eigenvector
    Do = np.sort(D1)[::-1]
    o = np.argsort(D1)[::-1]

    P = P1[:, o]
    D = np.diag(Do)

    return P, D


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
    pv_cons = 1. - rien

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
    np.dot(U, X): numpy.ndarray
        nt - 1 x nf array of projected timeseries
    """
    M = np.eye(nt, nt) - np.ones((nt, nt)) / nt

    # Eigen-vectors/-values of M; note that rk(M)=nt-1, so M has one eigenvalue equal to 0.
    u, d = speco(M)

    # (nt-1) first eigenvectors (ie the ones corresponding to non-zero eigenvalues)
    U = u[:, :nt - 1].T

    return np.dot(U, X)


def unproject_vectors(nt: int, Xc: np.ndarray) -> np.ndarray:
    """
    This function provides unprojects a matrix nt subspace to we can compute the trends
    :param nt: int
        number of time steps
    :param Xc: numpy.ndarray
        nt x nf array to be unprojected
    :return:
    np.dot(U, X): numpy.ndarray
        nt - 1 x nf array of projected timeseries
    """
    M = np.eye(nt, nt) - np.ones((nt, nt)) / nt

    # Eigen-vectors/-values of M; note that rk(M)=nt-1, so M has one eigenvalue equal to 0.
    u, d = speco(M)

    # inverse of the projection matrix
    Ui = np.linalg.inv(u.T)[:, :nt - 1]

    return np.dot(Ui, Xc)


def SSM(X_dict: Dict[str, np.ndarray], X_mm: np.ndarray) -> np.ndarray:
    """
    Calculates the squared difference between each models ensemble mean and the multi-model mean. Based on
    (Ribes et al., 2017)
    :param X_dict: dict
        Dictionary where keys are experiment names and values are arrays (n_members, n_time)
    :param X_mm: numpy.ndarray
        Multi-model ensemble mean, shape (n_time,)
    :return:
    np.diag(((Xc - Xc_mm) ** 2.).sum(axis=1)): numpy.ndarray
        nt -1 x nt - 1 array of the difference between each model ensemble mean the multi-model mean
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

    return np.diag(((Xc - Xc_mm) ** 2.).sum(axis=1))


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


def Cm_estimate(X_dict: Dict[str, np.ndarray], Cv: np.ndarray, X_mm: np.ndarray) -> np.ndarray:
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

    # First estimation of Cm
    Cm_hat = (1. / (nm - 1.)) * (_SSM - ((nm - 1.) / nm) * Cv_all)

    # Set negative eigenvalues to zero and recompose the signal
    S, X = np.linalg.eig(Cm_hat)
    S[S < 0] = 0
    Cm_pos_hat = np.linalg.multi_dot(
        [X, np.diag(S), np.linalg.inv(X)])  # spectral decomposition

    Cm_pos_hat = (1. + (1. / nm)) * Cm_pos_hat

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

    Cv_estimate = (1. / (nm ** 2.)) * Cv_all

    return Cv_estimate


if __name__ == "__main__":
    T = 11
    M = np.eye(T, T) - np.ones((T, T)) / T
    speco(M)
