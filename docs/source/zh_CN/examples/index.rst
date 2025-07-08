示例
====

使用示例和教程。

基本使用
--------

以下是一些 Skyborn 的基本使用示例：

GRIB to NetCDF 转换
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import skyborn

    # 简单转换
    skyborn.grib2nc('input.grib', 'output.nc')

    # 高精度转换
    skyborn.convert_grib_to_nc_simple('input.grib', 'output.nc', high_precision=True)

数据处理
~~~~~~~~

.. code-block:: python

    import skyborn
    import xarray as xr

    # 加载数据
    ds = xr.open_dataset('data.nc')

    # 转换经度范围
    ds_converted = skyborn.convert_longitude_range(ds, lon='longitude', center_on_180=True)

梯度计算
~~~~~~~~

.. code-block:: python

    import skyborn
    import numpy as np

    # 计算梯度
    data = np.random.rand(100)
    coords = np.linspace(-90, 90, 100)
    gradient = skyborn.calculate_gradient(data, coords)

因果分析
~~~~~~~~

.. code-block:: python

    import skyborn
    import numpy as np

    # 生成示例时间序列
    y1 = np.random.randn(1000)
    y2 = np.random.randn(1000)

    # Liang因果分析
    result = skyborn.liang_causality(y1, y2)
    print(f"信息流: {result['T21']}")

更多示例
--------

更多详细示例请参考：

- :doc:`../quickstart` - 快速开始指南
- :doc:`../api/index` - API 参考文档
