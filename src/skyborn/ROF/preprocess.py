import numpy as np
from typing import Tuple, Literal
from .utils import speco


class PreProcess:
    """
    A class to preprocess model data to be used in attribution methods like OLS, TLS and others. The script is
    heavily based on Aurelien Ribes (CNRM-GAME) scilab code, so for more information the user should consult the
    following reference:
        Ribes A., S. Planton, L. Terray (2012) Regularised optimal fingerprint for attribution.
        Part I: method, properties and idealised analysis. Climate Dynamics.

    :attribute y: numpy.ndarray
        Array of size nt with observations as a timeseries
    :attribute X: numpy.ndarray
        Array of size nt x nf, where nf is the number of forcings, with the model output as a timeseries
    :attribute Z: numpy.ndarray
        Array of ensembles with internal variability used to compute covariance matrix
    :attribute nt: int
        Number of timesteps for the timeseries

    :method extract_Z2(self, method='regular', frac=0.5):
        Split big sample Z in Z1 and Z2
    :method proj_fullrank(self, Z1, Z2):
        Provides a projection matrix to ensure its covariance matrix to be full-ranked
    :method creg(self, X, method='ledoit', alpha1=1e-10, alpha2=1):
        Computes regularized covariance matrix
    """

    def __init__(self, y: np.ndarray, X: np.ndarray, Z: np.ndarray):
        """
        :param y: numpy.ndarray
            Array of size nt with observations as a timeseries
        :param X: numpy.ndarray
            Array of size nt x nf, where nf is the number of forcings, with the model output as a timeseries
        :param Z: numpy.ndarray
            Array of ensembles with internal variability used to compute covariance matrix
        """
        self.y = y
        self.X = X
        self.Z = Z

        self.nt = y.shape[0]

    def extract_Z2(
        self, method: Literal["regular"] = "regular", frac: float = 0.5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        This function is used to split a big sample Z (dimension: nz x p, containing nz iid realisation of a random
        vector of size p) into two samples Z1 and Z2 (respectively of dimension nz1 x p and nz2 x p, with
        nz = nz1 + nz2). Further explanations in Ribes et al. (2012).
        :param method: str
            type of sampling used, for now may be only 'regular'
        :param frac: float
            fraction of realizations to put in Z2, the remaining is used in Z1
        :return:
        Z1: numpy.ndarray
            Array of size (nz1 x p)
        Z2: numpy.ndarray
            Array of size (nz2 x p)
        """
        nz = self.Z.shape[1]
        ind_z2 = np.zeros(nz)

        if method == "regular":
            # if frac = 0.5 ix would be [1, 3, 5, ..., ], so gets the index
            # for every two points
            # -1 is because python starts index at 0
            ix = np.arange(1 / frac - 1, nz, 1 / frac).astype(int)
            ind_z2[ix] = 1
            Z2 = self.Z[:, ind_z2 == 1]
            Z1 = self.Z[:, ind_z2 == 0]
        else:
            raise NotImplementedError("Method not implemented yet")

        return Z1, Z2

    def proj_fullrank(
        self, Z1: np.ndarray, Z2: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        This function provides a projection matrix U that can be applied to y, X, Z1 and Z2 to ensure its covariance
        matrix to be full-ranked. Uses variables defined in 'self', Z1 and Z2 computed in 'extract_Z2' method.
        :param Z1: numpy.ndarray
            Array of size (nz1 x p) of control simulation
        :param Z2: numpy.ndarray
            Array of size (nz1 x p) of control simulation
        :return:
        yc: numpy.ndarray
            y projected in U
        Xc: numpy.ndarray
            X projected in U
        Z1c: numpy.ndarray
            Z1 projected in U
        Z2c: numpy.ndarray
            Z2 projected in U
        """
        # M: the matrix corresponding to the temporal centering
        M = np.eye(self.nt, self.nt) - np.ones((self.nt, self.nt)) / self.nt

        # Eigen-vectors/-values of M; note that rk(M)=nt-1, so M has one eigenvalue equal to 0.
        u, d = speco(M)

        # (nt-1) first eigenvectors (ie the ones corresponding to non-zero eigenvalues)
        U = u[:, : self.nt - 1].T

        # Project all input data
        yc = np.dot(U, self.y)
        Xc = np.dot(U, self.X)
        Z1c = np.dot(U, Z1)
        Z2c = np.dot(U, Z2)

        return yc, Xc, Z1c, Z2c

    def creg(
        self,
        X: np.ndarray,
        method: Literal["ledoit", "specified"] = "ledoit",
        alpha1: float = 1e-10,
        alpha2: float = 1,
    ) -> np.ndarray:
        """
        This function compute the regularised covariance matrix estimate following the equation
        'Cr = alpha1 * Ip + alpha2 * CE' where alpha1 and alpha2 are parameters Ip is the p x p identity matrix and CE
        is the sample covariance matrix
        :param X: numpy.ndarray
            A n x p sample, meaning n iid realization of a random vector of size p.
        :param method: str
            method to compute the regularized covariance matrix
            - 'ledoit' uses Ledoit and Wolf (2003) estimate (default)
            - 'specified' uses specified values of alpha1 and alpha2
        :param alpha1: float
            Specified value for alpha1 (not used if method different than 'specified')
        :param alpha2: float
            Specified value for alpha1 (not used if method different than 'specified')
        :return:
        Cr: numpy.ndarray
            Regularized covariance matrix
        """
        n, p = X.shape

        CE = np.dot(X.T, X) / n  # sample covariance
        Ip = np.eye(p, p)

        # method for the regularised covariance matrix estimate as introduced by Ledoit & Wolf (2004) more specifically
        # on pages 379-380
        if method == "ledoit":
            m = np.trace(np.dot(CE, Ip)) / p  # first estimate in L&W
            XP = CE - m * Ip
            d2 = np.trace(np.dot(XP, XP.T)) / p  # second estimate in L&W

            bt = np.zeros(n)

            for i in range(n):
                Xi = X[i, :].reshape((1, p))
                Mi = np.dot(Xi.T, Xi)
                bt[i] = np.trace(np.dot((Mi - CE), (Mi - CE).T)) / p

            bb2 = (1.0 / n**2.0) * bt.sum()
            b2 = min(bb2, d2)  # third estimate in L&W
            a2 = d2 - b2  # fourth estimate in L&W

            alpha1 = b2 * m / d2
            alpha2 = a2 / d2

        elif method != "specified":
            raise NotImplementedError("Method not implemented yet")

        Cr = alpha1 * Ip + alpha2 * CE

        return Cr


if __name__ == "__main__":
    pass
    # observations
    # y = np.loadtxt('../example/HadCRUT4_1901-2010.asc')

    # # forcing patterns
    # X1 = np.loadtxt('../example/CNRM-CM5_ANT_1901-2010.asc')
    # X2 = np.loadtxt('../example/CNRM-CM5_NAT_1901-2010.asc')

    # X = np.stack([X1, X2], axis=1)

    # # internal variability ensemble
    # Z = np.loadtxt('../example/CTLruns.asc')
    # nt = len(y)
    # nz = int(len(Z) / nt)
    # Z = Z.reshape((nz, nt)).T  # every nt timesteps vector should be a column

    # p = PreProcess(y, X, Z)
    # Z1, Z2 = p.extract_Z2(frac=0.5)
    # yc, Xc, Z1c, Z2c = p.proj_fullrank(Z1, Z2)
    # Cr = p.creg(Z1c.T, method='ledoit')
