"""
This script is designed to support optimal fingerprinting attribution analysis by providing essential data processing and analysis functions.

Qianye Su
suqianye2000@gmail.com

Reference
 - ROF:
     https://github.com/pinplex/PyDnA
     https://github.com/rafaelcabreu/attribution, by Rafael Abreu
"""

import numpy as np
import scipy.stats as stats
from typing import Dict, Any, Tuple, Literal
from .preprocess import PreProcess
from .utils import chi2_test


class AttributionModel:
    """
    A class for attribution models. The OLS implementation is heavily based on Aurelien Ribes (CNRM-GAME) scilab code
    (see more in 'preprocess.py'). Also, Aurelien Ribes model proposed in 2017 is implemented following the reference:
        Ribes, Aurelien, et al. (2017) A new statistical approach to climate change detection and attribution.
        Climate Dynamics.

    :attribute X: numpy.ndarray
        Array of size nt x nf, where nf is the number of forcings, with the model output as a timeseries
    :attribute y: numpy.ndarray
        Array of size nt with observations as a timeseries

    :method ols(self, Cf, Proj, Z2, cons_test='AT99'):
        Ordinary Least Square (OLS) estimation of beta from the linear model y = beta * X + epsilon
    """

    def __init__(self, X: np.ndarray, y: np.ndarray):
        """
        :param X: numpy.ndarray
            Array of size nt x nf, where nf is the number of forcings, with the model output as a timeseries
        :param y: numpy.ndarray
            Array of size nt with observations as a timeseries
        """
        self.y = y
        self.X = X
        self.nt = y.shape[0]  # Number of time steps
        # 1 stands for the number of spatial patterns (dealing only with timeseries)
        self.n_reduced = self.nt - 1  # Number of reduced dimensions
        self.num_forcings = X.shape[1]  # Number of forcing factors

    def ols(
        self,
        Cf: np.ndarray,
        Proj: np.ndarray,
        Z2: np.ndarray,
        cons_test: Literal["AT99"] = "AT99",
    ) -> Dict[str, np.ndarray]:
        """
        Ordinary Least Square (OLS) estimation of beta from the linear model y = beta * X + epsilon as discussed in the
        following reference:
            Allen, Myles R., and Simon FB Tett (1999) Checking for model consistency in optimal fingerprinting.
            Climate Dynamics.
        :param Cf: numpy.ndarray
            Covariance matrix. Be sure that Cf is invertible to use this model (look at PreProcess class)
        :param Proj: numpy.ndarray
            Array of zeros and ones, indicating which forcings in each simulation
        :param Z2: numpy.ndarray
            Array of size (nz1 x p) of control simulation used to compute consistency test
        :param cons_test: str
            Which consistency test to be used
            - 'AT99' the formula provided by Allen & Tett (1999) (default)
        :return:
        Beta_hat: dict
            Dictionary with estimation of beta_hat and the upper and lower confidence intervals
        """

        # computes the covariance inverse
        Cf1 = np.linalg.inv(Cf)

        _Ft = np.linalg.multi_dot([self.X.T, Cf1, self.X])
        _Ft1 = np.linalg.inv(_Ft)
        Ft = np.linalg.multi_dot([_Ft1, self.X.T, Cf1]).T

        _y = self.y.reshape((self.nt, 1))
        beta_hat = np.linalg.multi_dot([_y.T, Ft, Proj.T])

        # 1-D confidence interval
        nz2 = Z2.shape[1]
        Z2t = Z2.T
        Var_valid = np.dot(Z2t.T, Z2t) / nz2
        Var_beta_hat = np.linalg.multi_dot([Proj, Ft.T, Var_valid, Ft, Proj.T])

        beta_hat_inf = beta_hat - 2.0 * stats.t.cdf(0.95, df=nz2) * np.sqrt(
            np.diag(Var_beta_hat).T
        )
        beta_hat_sup = beta_hat + 2.0 * stats.t.cdf(0.95, df=nz2) * np.sqrt(
            np.diag(Var_beta_hat).T
        )

        # consistency check
        epsilon = _y - np.linalg.multi_dot([self.X, np.linalg.inv(Proj), beta_hat.T])

        if cons_test == "AT99":  # formula provided by Allen & Tett (1999)
            d_cons = np.linalg.multi_dot(
                [epsilon.T, np.linalg.pinv(Var_valid), epsilon]
            ) / (self.n_reduced - self.num_forcings)
            rien = stats.f.cdf(d_cons, dfn=self.n_reduced - self.num_forcings, dfd=nz2)
            pv_cons = 1 - rien

        print("Consistency test: %s p-value: %.5f" % (cons_test, pv_cons))

        Beta_hat = {
            "beta_hat": beta_hat[0],
            "beta_hat_inf": beta_hat_inf[0],
            "beta_hat_sup": beta_hat_sup[0],
        }

        return Beta_hat

    def ribes(
        self, Cxi: np.ndarray, Cy: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Aurelien Ribes model proposed in 2017 is implemented following the reference:
        Ribes, Aurelien, et al. (2017) A new statistical approach to climate change detection and attribution.
        Climate Dynamics. It considers the following set of equations:

            Y_star = sum(X_star_i) for i from 1 to nf where nf is the number of forcings
            Y = Y_star + epsilon_y
            Xi = X_star_i + epsilon_xi

        Where epislon_y ~ N(0, Cy) and epislon_xi ~ N(0, Cxi)

        :param Cxi: numpy.ndarray
            Covariance matrix for each of the forcings Xi. Should be a 3D array (nt, nt, nf)
        :param Cy: numpy.ndarray
            Covariance matrix for the observations.
        :return:
        """
        X = self.X.sum(axis=1)
        Cx = Cxi.sum(axis=0)

        # Estimate the true state of variables (y) and (Xi) y_star and X_star_i using the MLE y_star_hat and
        # Xi_star_hat, respectively
        Xi_star_hat = np.zeros(self.X.shape)
        y_star_hat = self.y + np.linalg.multi_dot(
            [Cy, np.linalg.inv(Cy + Cx), (X - self.y)]
        )
        for i in range(Xi_star_hat.shape[1]):
            Xi_star_hat[:, i] = self.X[:, i] + np.linalg.multi_dot(
                [Cxi[i], np.linalg.inv(Cy + Cx), (self.y - X)]
            )

        # calculates variance for Y_star_hat
        Cy_star_hat = np.linalg.inv(np.linalg.inv(Cy) + np.linalg.inv(Cx))

        # calculates variance for Xi_star_hat
        Cxi_star_hat = np.zeros(Cxi.shape)
        for i in range(Cxi_star_hat.shape[0]):
            Cxi_temp = Cxi * 1.0
            # sum for every j different than i
            Cxi_temp[i] = 0.0
            Cxi_sum = Cxi_temp.sum(axis=0)

            Cxi_star_hat[i] = np.linalg.inv(
                np.linalg.inv(Cxi[i]) + np.linalg.inv(Cy + Cxi_sum)
            )

        # hypothesis test: compare with chi-square distribution
        print("#" * 60)
        print(
            "Hypothesis testing p-value for Chi-2 distribution and Maximum Likelihood ..."
        )

        # (internal variability only)
        d_cons = np.linalg.multi_dot([self.y.T, np.linalg.inv(Cy), self.y])
        print(
            "%30s: %.7f (%.7f)"
            % (
                "Internal variability only",
                chi2_test(d_cons, self.nt),
                np.exp(d_cons / -2.0),
            )
        )

        # (all forcings)
        d_cons = np.linalg.multi_dot(
            [(self.y - X).T, np.linalg.inv(Cy + Cx), (self.y - X)]
        )
        print(
            "%30s: %.7f (%.7f)"
            % ("All forcings", chi2_test(d_cons, self.nt), np.exp(d_cons / -2.0))
        )

        # (individual forcings)
        for i in range(self.X.shape[1]):
            d_cons = np.linalg.multi_dot(
                [
                    (self.y - self.X[:, i]).T,
                    np.linalg.inv(Cy + Cxi[i]),
                    (self.y - self.X[:, i]),
                ]
            )
            print(
                "%30s: %.7f (%.7f)"
                % (
                    "Forcing no %d only" % (i + 1),
                    chi2_test(d_cons, self.nt),
                    np.exp(d_cons / -2.0),
                )
            )

        return y_star_hat, Xi_star_hat, Cy_star_hat, Cxi_star_hat


if __name__ == "__main__":
    pass
