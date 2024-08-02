import xarray as xr
def wrap_lon_to_180(data: xr.DataArray | xr.Dataset, 
                    lon: str = 'lon', 
                    center_on_180: bool = True
                    ) -> xr.DataArray | xr.Dataset:
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
        data = data.assign_coords(**{lon: (lambda x: ((x[lon] + 180) % 360) - 180)})
    return data.sortby(lon, ascending=True)

data = xr.open_dataset(r"E:/CumData/HadISST_sst.nc")
# print(data)
data = wrap_lon_to_180(data, lon = 'longitude', center_on_180 = True)
print(data)
data = wrap_lon_to_180(data, lon = 'longitude', center_on_180 = False)
print(data)