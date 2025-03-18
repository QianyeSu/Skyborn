import xarray as xr
import cfgrib

def read_netcdf(file_path: str) -> xr.Dataset:
    """读取NetCDF格式气象数据
    
    Args:
        file_path: 文件路径
        
    Returns:
        xarray Dataset对象
    """
    return xr.open_dataset(file_path)

def read_grib(file_path: str, **kwargs) -> xr.Dataset:
    """读取GRIB格式气象数据
    
    Args:
        file_path: 文件路径
        kwargs: 传递给cfgrib.open_dataset的额外参数
        
    Returns:
        xarray Dataset对象
    """
    return cfgrib.open_dataset(file_path, **kwargs)
