快速开始指南
============

本指南将帮助您在几分钟内开始使用 Skyborn。

第一步
------

导入 Skyborn
~~~~~~~~~~~~~

.. code-block:: python

   import skyborn
   import xarray as xr
   import numpy as np

基础数据转换
~~~~~~~~~~~~

将 GRIB 文件转换为 NetCDF 格式：

.. code-block:: python

   # 简单转换
   skyborn.grib2nc('input.grib', 'output.nc')

   # 带选项的转换
   skyborn.convert_grib_to_nc_simple(
       'input.grib',
       'output.nc',
       high_precision=True,
       compress=True
   )

数据处理
--------

加载和分析数据
~~~~~~~~~~~~~~

.. code-block:: python

   # 加载 NetCDF 数据
   ds = xr.open_dataset('output.nc')

   # 基本数据信息
   print(ds)
   print(ds.data_vars)

计算梯度
~~~~~~~~

.. code-block:: python

   # 经向梯度 (∂/∂lat)
   temp_grad_lat = skyborn.calculate_meridional_gradient(
       ds['temperature'],
       'latitude'
   )

   # 纬向梯度 (∂/∂lon)
   temp_grad_lon = skyborn.calculate_zonal_gradient(
       ds['temperature'],
       'longitude'
   )

   # 通用梯度
   pressure_grad = skyborn.calculate_gradient(
       ds['pressure'],
       'level'
   )

数据操作
~~~~~~~~

.. code-block:: python

   # 转换经度范围
   ds_converted = skyborn.convert_longitude_range(ds, target_range=(-180, 180))

   # 线性回归
   slope, intercept, r_value = skyborn.linear_regression(
       ds['temperature'].values.flatten(),
       ds['pressure'].values.flatten()
   )

高级功能
--------

插值和重网格化
~~~~~~~~~~~~~~

.. code-block:: python

   from skyborn.interp import interpolation, regridding

   # 数据插值
   interpolated = interpolation.interpolate_data(
       ds['temperature'],
       target_coords={'latitude': np.arange(-90, 91, 1)}
   )

   # 重网格化到新网格
   regridded = regridding.regrid_data(
       ds['temperature'],
       target_grid=(180, 360)  # 新分辨率
   )

因果关系分析
~~~~~~~~~~~~

.. code-block:: python

   # Liang 因果关系分析
   causality_result = skyborn.liang_causality(
       ds['temperature'].values,
       ds['precipitation'].values
   )

   # Granger 因果关系
   granger_result = skyborn.granger_causality(
       ds['temperature'].values,
       ds['pressure'].values,
       max_lag=5
   )

数据可视化
~~~~~~~~~~

.. code-block:: python

   from skyborn.plot import plotting, modplot

   # 基础绘图
   fig, ax = plotting.plot_contour(
       ds['temperature'].isel(time=0),
       levels=20,
       title='温度分布'
   )

   # 专业大气绘图
   fig = modplot.plot_wind_field(
       ds['u_wind'].isel(time=0),
       ds['v_wind'].isel(time=0)
   )

批量处理
--------

处理多个文件
~~~~~~~~~~~~

.. code-block:: python

   # 批量转换 GRIB 文件
   converted_files = skyborn.batch_convert_grib_to_nc(
       input_directory='./grib_data/',
       output_directory='./netcdf_data/',
       pattern='*.grb',
       high_precision=True,
       compress=True
   )

   print(f"已转换 {len(converted_files)} 个文件")

常用工作流程
------------

ERA5 数据处理
~~~~~~~~~~~~~

.. code-block:: python

   # 针对 ERA5 数据优化
   skyborn.convert_grib_to_nc(
       grib_files='era5_data.grib',
       output_file='era5_processed.nc',
       ignore_keys=['method', 'type', 'stream'],
       split_keys=['param', 'levtype'],
       data_type='NC_FLOAT',
       unlimited_dimension='time',
       file_kind=4,
       deflate_level=4
   )

气候分析流水线
~~~~~~~~~~~~~~

.. code-block:: python

   # 完整分析工作流程
   def analyze_climate_data(grib_file, output_dir):
       # 1. 转换数据
       nc_file = f"{output_dir}/converted_data.nc"
       skyborn.grib2nc(grib_file, nc_file, compress=True)

       # 2. 加载和处理
       ds = xr.open_dataset(nc_file)

       # 3. 计算梯度
       temp_grad = skyborn.calculate_meridional_gradient(
           ds['temperature'], 'latitude'
       )

       # 4. 因果关系分析
       causality = skyborn.liang_causality(
           ds['temperature'].values,
           ds['pressure'].values
       )

       # 5. 保存结果
       results = xr.Dataset({
           'temperature_gradient': temp_grad,
           'causality_strength': (['time'], causality)
       })
       results.to_netcdf(f"{output_dir}/analysis_results.nc")

       return results

最佳实践
--------

性能提示
~~~~~~~~

1. **使用压缩** 处理大型数据集
2. **适当指定数据类型** (NC_SHORT vs NC_FLOAT)
3. **分块处理** 非常大的文件
4. **使用适当的 ignore_keys** 针对您的数据类型

内存管理
~~~~~~~~

.. code-block:: python

   # 对于大型数据集，使用 dask
   import dask.array as da

   # 延迟加载数据
   ds = xr.open_dataset('large_file.nc', chunks={'time': 100})

   # 使用 dask 处理
   result = skyborn.calculate_gradient(ds['temperature'], 'latitude')
   result = result.compute()  # 执行计算

错误处理
~~~~~~~~

.. code-block:: python

   from skyborn.conversion import GribToNetCDFError

   try:
       skyborn.grib2nc('input.grib', 'output.nc')
   except GribToNetCDFError as e:
       print(f"转换失败: {e}")
   except FileNotFoundError as e:
       print(f"文件未找到: {e}")

下一步
------

* 浏览 :doc:`api/index` 获取详细的函数文档
* 查看 :doc:`examples/index` 获取更全面的示例
* 阅读 :doc:`modules/index` 了解具体模块
