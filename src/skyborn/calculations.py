import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from typing import Tuple, Union
from sklearn.feature_selection import f_regression


def linearRegression(mat, sequence):
    '''
    回归函数 第一个是需要回归的变量，第二个为一个数组
    '''
    if len(mat) != len(sequence):
        raise ValueError('Data, array must be must be equal!!!!!!!'
                         ' %s and %s' % (mat.shape[0], len(sequence)))
    A = np.column_stack((sequence, np.ones(sequence.shape[0])))
    (a, b, c) = mat.shape
    mat_rshp = mat.reshape((a, b*c))
    mat_rshp_reg = np.linalg.lstsq(A, mat_rshp, rcond=None)[0][0]
    mat_reg = mat_rshp_reg.reshape((b, c))
    pvalue = f_regression(mat_rshp, sequence)[1].reshape(b, c)
    return mat_reg, pvalue


def convert_longitude_range(data: Union[xr.DataArray, xr.Dataset],
                            lon: str = 'lon',
                            center_on_180: bool = True
                            ) -> Union[xr.DataArray, xr.Dataset]:
    '''
    Wrap longitude coordinates of DataArray or Dataset to either -180..179 or 0..359.

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        An xarray DataArray or Dataset object containing longitude coordinates.
    lon : str, optional
        The name of the longitude coordinate, default is 'lon'.
    center_on_180 : bool, optional
        If True, wrap longitude from 0..359 to -180..179;
        If False, wrap longitude from -180..179 to 0..359.

    Returns
    -------
    xr.DataArray or xr.Dataset
        The DataArray or Dataset with wrapped longitude coordinates.
    '''
    # Wrap -180..179 to 0..359
    if center_on_180:
        data = data.assign_coords(**{lon: (lambda x: (x[lon] % 360))})
    # Wrap 0..359 to -180..179
    else:
        data = data.assign_coords(
            **{lon: (lambda x: ((x[lon] + 180) % 360) - 180)})
    return data.sortby(lon, ascending=True)
