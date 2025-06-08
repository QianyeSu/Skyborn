import os
import numpy as np
import pandas as pd
from glob import glob
from scipy import stats


def speco(C):
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

#################################################################


def chi2_test(d_cons, df):
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

#################################################################


def project_vectors(nt, X):
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

#################################################################


def unproject_vectors(nt, Xc):
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

#################################################################


def SSM(exp, X_mm, init=1955, end=1995):
    """
    Calculates the squared difference between each models ensemble mean and the multi-model mean. Based on
    (Ribes et al., 2017)
    :param exp: str
        Experiment to calculate the difference (e.g., 'historical', 'historicalNat')
    :param X_mm: numpy.ndarray
        Array with multi-model ensemble mean
    :param init: int
        Correspondent year to start the analysis
    :param end: int
        Correspondent year to finish the analysis
    :return:
    np.diag(((Xc - Xc_mm) ** 2.).sum(axis=1)): numpy.ndarray
        nt -1 x nt - 1 array of the difference between each model ensemble mean the multi-model mean
    """
    # reads ensemble mean for each model
    ifiles = glob('data/model/%s/ensmean/*_%s_%s.csv' % (exp, init, end))

    df = pd.DataFrame()
    for ifile in ifiles:
        df_temp = pd.read_csv(ifile, index_col=0, parse_dates=True)
        df = pd.concat([df, df_temp['anomaly'].to_frame(
            os.path.basename(ifile)[:-4])], axis=1)

    # remove columns (ensemble members with nan)
    df.dropna(inplace=True, axis=1)
    # gets ensemble values and multi model (mm) ensemble
    X = df.values

    # project the data
    Xc = project_vectors(X.shape[0], X)
    Xc_mm = project_vectors(X.shape[0], X_mm.reshape((X.shape[0], 1)))

    return np.diag(((Xc - Xc_mm) ** 2.).sum(axis=1))

#################################################################


def get_nruns(exp, how='pandas', init=1955, end=1995):
    """
    Reads the number of runs for each model
    :param exp: str
        Experiment to calculate the difference (e.g., 'historical', 'historicalNat')
    :param how: str
        Used to see if the number of runs is calculated using the pandas dataframes or text file ('historicalOA' for
        example)
    :param init: int
        Correspondent year to start the analysis
    :param end: int
        Correspondent year to finish the analysis
    :return:
    nruns: numpy.ndarray
       Array with the number of runs for each model
    """
    if how == 'pandas':
        ifiles = glob('data/model/%s/ensemble/*_%s_%s.csv' % (exp, init, end))
        nruns = []

        for ifile in sorted(ifiles):
            df_temp = pd.read_csv(ifile, index_col=0, parse_dates=True)
            nruns.append(len(df_temp.columns))

        nruns = np.array(nruns)
    elif how == 'loadtxt':
        nruns = np.loadtxt('data/model/%s/ensemble/nruns_%s.txt' % (exp, exp))

    return nruns

#################################################################


def Cm_estimate(exp, Cv, X_mm, how_nr='pandas', init=1955, end=1995):
    """
    Estimated covariance matrix for model error (Ribes et al., 2017)
    :param exp: str
        Experiment to calculate the difference (e.g., 'historical', 'historicalNat')
    :param Cv: numpy.ndarray
        Array with internal variability covariance matrix
    :param X_mm: numpy.ndarray
        Array with multi-model ensemble mean
    :param how_nr:
        Used to see if the number of runs is calculated using the pandas dataframes or text file ('historicalOA' for
        example)
    :param init: int
        Correspondent year to start the analysis
    :param end: int
        Correspondent year to finish the analysis
    :return:
    Cm_pos_hat: numpy.ndarray
        Estimated covariance matrix for model error
    """

    # model difference
    _SSM = SSM(exp, X_mm, init=init, end=end)

    # nruns - number of runs / nm - number of models
    nruns = get_nruns(exp, how=how_nr, init=init, end=end)
    nm = len(nruns)

    Cv_all = np.zeros(Cv.shape)
    for nr in nruns:
        Cv_all += Cv / nr

    # first estimation of Cm
    Cm_hat = (1. / (nm - 1.)) * (_SSM - ((nm - 1.) / nm) * Cv_all)

    # set negative eigenvalues to zero and recompose the signal
    S, X = np.linalg.eig(Cm_hat)
    S[S < 0] = 0
    Cm_pos_hat = np.linalg.multi_dot(
        [X, np.diag(S), np.linalg.inv(X)])  # spectral decomposition

    Cm_pos_hat = (1. + (1. / nm)) * Cm_pos_hat

    return Cm_pos_hat

#################################################################


def Cv_estimate(exp, Cv, how_nr='pandas', init=1955, end=1995):
    """
    Estimated covariance matrix for internal variability considering multiple models (Ribes et al., 2017)
    :param exp: str
        Experiment to calculate the difference (e.g., 'historical', 'historicalNat')
    :param Cv: numpy.ndarray
        Array with internal variability covariance matrix
    :param how_nr:
        Used to see if the number of runs is calculated using the pandas dataframes or text file ('historicalOA' for
        example)
    :param init: int
        Correspondent year to start the analysis
    :param end: int
        Correspondent year to finish the analysis
    :return:
    Cv_estimate: numpy.ndarray
        Estimated covariance matrix for internal variability considering multiple models
    """
    # nruns - number of runs / nm - number of models
    nruns = get_nruns(exp, how=how_nr, init=init, end=end)
    nm = len(nruns)

    Cv_all = np.zeros(Cv.shape)
    for nr in nruns:
        Cv_all += Cv / nr

    Cv_estimate = (1. / (nm ** 2.) * Cv_all)

    return Cv_estimate

#################################################################


if __name__ == "__main__":
    T = 11
    M = np.eye(T, T) - np.ones((T, T)) / T
    speco(M)
