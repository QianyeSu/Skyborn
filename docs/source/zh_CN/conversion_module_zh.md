# GRIB 到 NetCDF 转换模块

`skyborn.conversion` 模块提供了一个 Python 接口来使用 eccodes 的 `grib_to_netcdf` 工具将 GRIB 文件转换为 NetCDF 格式。

## 安装要求

确保您已安装 eccodes 库：

```bash
# conda 安装
conda install -c conda-forge eccodes

# 或者使用 pip
pip install eccodes
```

## 快速开始

### 基础用法

```python
import skyborn

# 简单转换
skyborn.grib2nc('input.grib', 'output.nc')

# 或者使用完整函数名
skyborn.convert_grib_to_nc_simple('input.grib', 'output.nc')
```

### 高级用法

```python
from skyborn.conversion import convert_grib_to_nc

# 高精度转换
convert_grib_to_nc(
    grib_files='input.grib',
    output_file='output.nc',
    data_type='NC_FLOAT',           # 高精度浮点数
    ignore_keys=['type', 'step'],   # 忽略特定键
    unlimited_dimension='time',     # 设置时间为无限维度
    file_kind=4,                    # netCDF-4 格式
    deflate_level=6,                # 压缩级别
    shuffle=True                    # 启用洗牌压缩
)
```

### 批量转换

```python
from skyborn.conversion import batch_convert_grib_to_nc

# 转换目录中的所有 GRIB 文件
converted_files = batch_convert_grib_to_nc(
    input_directory='/path/to/grib_files',
    output_directory='/path/to/output',
    pattern='*.grib*',
    high_precision=True,
    compress=True
)
```

## 支持的功能

### 数据类型选项
- `NC_BYTE`: 8位整数
- `NC_SHORT`: 16位整数（默认）
- `NC_INT`: 32位整数
- `NC_FLOAT`: 32位浮点数
- `NC_DOUBLE`: 64位浮点数

### NetCDF 格式选项
- `1`: netCDF classic 格式
- `2`: netCDF 64位 classic 格式（默认）
- `3`: netCDF-4 格式
- `4`: netCDF-4 classic 模式格式

### 压缩选项
- 支持 netCDF-4 格式的压缩（0-9级别）
- 支持洗牌压缩以提高压缩率

## 函数参考

### convert_grib_to_nc()
完整功能的 GRIB 到 NetCDF 转换函数。

### convert_grib_to_nc_simple()
简化接口，适合快速转换。

### batch_convert_grib_to_nc()
批量转换目录中的多个文件。

### grib2nc()
`convert_grib_to_nc_simple()` 的简短别名。

### grib_to_netcdf()
`convert_grib_to_nc()` 的别名。

## 注意事项

1. 仅支持规则经纬度网格和规则高斯网格
2. 需要安装 eccodes 库
3. 对于大文件，建议使用压缩选项
4. 对于高精度需求，使用 `NC_FLOAT` 或 `NC_DOUBLE` 数据类型

## 示例脚本

查看 `examples/grib_to_netcdf_examples.py` 获取更多使用示例。

## 错误处理

模块会抛出 `GribToNetCDFError` 异常当转换失败时。确保：
- eccodes 已正确安装
- 输入文件存在且可读
- 输出目录有写权限
- GRIB 文件格式受支持
