import numpy as np
import pandas as pd
from typing import List, Tuple, Union
from .utils import unproject_vectors


def calculate_trend(y: np.ndarray) -> float:
    """
    Calculate the trend by unprojecting a vector to the nt subspace and the using OLS estimation
    :param y: numpy.ndarray
        nt -1 array used to calculate the trend
    :return:
    beta_hat: float
        value of the scaling factor of the OLS adjustment
    """
    nt = len(y) + 1
    y_un = unproject_vectors(nt, y)  # unproject the data

    X = np.vstack([np.ones(nt), np.arange(nt)]).T

    # use OLS to calculate the trend
    beta_hat = np.linalg.multi_dot(
        [np.linalg.inv(np.linalg.multi_dot([X.T, X])), X.T, y_un]
    )

    return beta_hat[1]  # only return the trend


def calculate_uncertainty(
    y: np.ndarray, Cy: np.ndarray, alpha: float = 0.05, nsamples: int = 4000
) -> np.ndarray:
    """
    Calculate trend uncertainty by generating multiple series and the calculating the confidence interval
    :param y: numpy.ndarray
        nt -1 array used to calculate the trend
    :param Cy: numpy.ndarray
        nt -1 x nt -1 covariance matrix from the y vector
    :param alpha: float
        significance level
    :param nsamples: int
        number of repetitions
    :return:
    np.array([trend_min, trend_max]): np.ndarray
        array with the minimum and maximum values from the confidence interval
    """
    trends = np.zeros(nsamples)
    for i in range(nsamples):
        # generate random vector based on the mean and cov matrix
        y_random = np.random.multivariate_normal(y, Cy)

        # calculate the trend
        trends[i] = calculate_trend(y_random)

    trend_min = np.percentile(trends, (alpha * 100) / 2.0)
    trend_max = np.percentile(trends, 100 - (alpha * 100) / 2.0)

    return np.array([trend_min, trend_max])


def all_trends(
    y_star_hat: np.ndarray,
    Xi_star_hat: np.ndarray,
    Cy_star_hat: np.ndarray,
    Cxi_star_hat: np.ndarray,
) -> pd.DataFrame:
    """
    Calculate all trends (observations and each individual forcing) and save it to a csv file
    :param y_star_hat: np.ndarray
        vector of observations
    :param Xi_star_hat: np.ndarray
        nt -1 x nf matrix of forcings where nt is the number of time steps and nf is the number of forcings
    :param Cy_star_hat: np.ndarray
        nt -1 x nt -1 covariance matrix for the observations
    :param Cxi_star_hat: np.ndarray
        nf x nt -1 x nt -1 covariance matrix for each individual forcing
    :return:
    df: pandas.DataFrame
        dataframe with the trend for the observations and each of the forcings
    """
    trends_list: List[List[Union[str, float]]] = []

    trend = calculate_trend(y_star_hat)
    confidence_interval = calculate_uncertainty(y_star_hat, Cy_star_hat, alpha=0.05)

    trends_list.append(
        ["Observation", trend, confidence_interval[0], confidence_interval[1]]
    )

    print("-" * 60)
    print("Trends from the analysis ...")
    print(
        "%30s: %.3f (%.3f, %.3f)"
        % ("Observation", trend, confidence_interval[0], confidence_interval[1])
    )

    nf = Xi_star_hat.shape[1]
    for i in range(nf):
        trend = calculate_trend(Xi_star_hat[:, i])
        confidence_interval = calculate_uncertainty(
            Xi_star_hat[:, i], Cxi_star_hat[i], alpha=0.05
        )
        print(
            "%30s: %.3f (%.3f, %.3f)"
            % (
                "Forcing no %d only" % (i + 1),
                trend,
                confidence_interval[0],
                confidence_interval[1],
            )
        )

        trends_list.append(
            [
                "Forcing no %d only" % (i + 1),
                trend,
                confidence_interval[0],
                confidence_interval[1],
            ]
        )

    # save data as csv
    df = pd.DataFrame(
        trends_list, columns=["forcing", "trend", "trend_min", "trend_max"]
    )
    # df.to_csv('trends.csv', index=False)

    return df
