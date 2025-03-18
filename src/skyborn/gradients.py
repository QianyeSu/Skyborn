import numpy as np
from typing import Union, Optional

EARTH_RADIUS = 6371e3  # 地球平均半径 (m)


def calculate_gradient(field: np.ndarray, coordinates: np.ndarray,
                       axis: int = -1, radius: float = 6371000.0) -> np.ndarray:
    """计算任意维度数组沿指定坐标的梯度

    Args:
        field: 待计算梯度的数据场，可以是任意维度的数组，如(time, level, lat, lon)
        coordinates: 沿着计算梯度的坐标数组，如纬度或经度值
        axis: 指定计算梯度的维度轴，默认为最后一个维度(-1)
        radius: 地球半径，默认为6371000.0米，用于纬度梯度计算

    Returns:
        与输入场相同形状的梯度场
    """
    # 检查输入数据维度是否匹配
    if coordinates.size != field.shape[axis]:
        raise ValueError(
            f"坐标数组大小({coordinates.size})与场数据在指定轴上的大小({field.shape[axis]})不匹配")

    # 确定是经度还是纬度坐标
    is_latitude = False
    if np.min(coordinates) >= -90 and np.max(coordinates) <= 90:
        is_latitude = True
        # 对于纬度，计算实际距离（单位：米）
        if is_latitude:
            # 将纬度转换为实际距离
            distances = coordinates * np.pi / 180.0 * radius
        else:
            # 对于经度，我们需要考虑纬度的影响，但这需要额外的纬度信息
            # 这里简单处理为直接使用经度差
            distances = coordinates

    # 创建与输入相同形状的输出数组
    gradient = np.zeros_like(field, dtype=float)

    # 为了使用numpy的高级索引，我们需要创建索引数组
    ndim = field.ndim
    idx_ranges = [slice(None)] * ndim

    # 对内部点使用中心差分
    inner_range = slice(1, field.shape[axis]-1)
    idx_forward = idx_ranges.copy()
    idx_forward[axis] = slice(2, field.shape[axis])

    idx_center = idx_ranges.copy()
    idx_center[axis] = inner_range

    idx_backward = idx_ranges.copy()
    idx_backward[axis] = slice(0, field.shape[axis]-2)

    # 使用矢量化操作计算内部点的梯度
    forward_dists = np.diff(distances[1:])
    backward_dists = np.diff(distances[:-1])
    total_dists = distances[2:] - distances[:-2]

    # 创建系数数组，形状适合广播
    shape = [1] * ndim
    shape[axis] = len(forward_dists)

    a0 = forward_dists.reshape(shape)
    b0 = backward_dists.reshape(shape)
    c0 = total_dists.reshape(shape)

    # 使用加权差分公式计算梯度
    gradient[tuple(idx_center)] = (
        b0 / a0 / c0 * field[tuple(idx_forward)] -
        a0 / b0 / c0 * field[tuple(idx_backward)] +
        (a0 - b0) / a0 / b0 * field[tuple(idx_center)]
    )

    # 处理边界点（前向和后向差分）
    # 左边界
    left_idx = idx_ranges.copy()
    left_idx[axis] = 0
    left_idx_plus = idx_ranges.copy()
    left_idx_plus[axis] = 1
    gradient[tuple(left_idx)] = (field[tuple(left_idx_plus)] -
                                 field[tuple(left_idx)]) / (distances[1] - distances[0])

    # 右边界
    right_idx = idx_ranges.copy()
    right_idx[axis] = -1
    right_idx_minus = idx_ranges.copy()
    right_idx_minus[axis] = -2
    gradient[tuple(right_idx)] = (field[tuple(right_idx)] -
                                  field[tuple(right_idx_minus)]) / (distances[-1] - distances[-2])

    return gradient


def calculate_meridional_gradient(field: np.ndarray, latitudes: np.ndarray,
                                  lat_axis: int = -1, radius: float = 6371000.0) -> np.ndarray:
    """计算经向梯度（沿纬度方向的梯度）

    Args:
        field: 待计算梯度的数据场，可以是任意维度的数组
        latitudes: 纬度数组（度）
        lat_axis: 指定纬度所在的轴，默认为最后一个维度(-1)
        radius: 地球半径，默认为6371000.0米

    Returns:
        经向梯度场
    """
    return calculate_gradient(field, latitudes, axis=lat_axis, radius=radius)


def calculate_vertical_gradient(field: np.ndarray,
                                pressure: np.ndarray,
                                pressure_axis: int = -3) -> np.ndarray:
    """计算垂直梯度（沿气压方向的梯度）

    Args:
        field: 待计算梯度的数据场
        pressure: 气压数组（Pa），必须为单调递减
        pressure_axis: 指定气压所在的轴，默认为倒数第三个维度(-3)

    Returns:
        垂直梯度场
    """
    return calculate_gradient(field, pressure, axis=pressure_axis, radius=None)


def calculate_zonal_gradient(field: np.ndarray, longitudes: np.ndarray, latitudes: np.ndarray,
                             lon_axis: int = -1, lat_axis: int = -2, radius: float = 6371000.0) -> np.ndarray:
    """计算纬向梯度（沿经度方向的梯度）

    Args:
        field: 待计算梯度的数据场，可以是任意维度的数组
        longitudes: 经度数组（度）
        latitudes: 纬度数组（度），用于计算不同纬度下经度的实际距离
        lon_axis: 指定经度所在的轴，默认为最后一个维度(-1)
        lat_axis: 指定纬度所在的轴，默认为倒数第二个维度(-2)
        radius: 地球半径，默认为6371000.0米

    Returns:
        纬向梯度场
    """
    # 获取纬度因子，用于调整不同纬度下经度间的实际距离
    cos_lat = np.cos(np.radians(latitudes))

    # 如果场是4D (time, level, lat, lon)
    if field.ndim == 4 and lon_axis == -1 and lat_axis == -2:
        # 创建一个广播形状的纬度因子数组
        cos_lat_expanded = cos_lat.reshape(1, 1, -1, 1)

        # 将经度转换为考虑纬度的实际距离
        effective_distances = np.radians(
            longitudes) * radius * cos_lat_expanded

        # 现在计算梯度
        return calculate_gradient(field, effective_distances, axis=lon_axis, radius=1.0)

    # 如果场是3D (time, lat, lon)
    elif field.ndim == 3 and lon_axis == -1 and lat_axis == -2:
        cos_lat_expanded = cos_lat.reshape(1, -1, 1)
        effective_distances = np.radians(
            longitudes) * radius * cos_lat_expanded
        return calculate_gradient(field, effective_distances, axis=lon_axis, radius=1.0)

    else:
        # 对于其他维度组合，需要创建适当的广播形状
        broadcast_shape = [1] * field.ndim
        broadcast_shape[lat_axis] = len(latitudes)
        cos_lat_expanded = cos_lat.reshape(broadcast_shape)

        # 创建有效距离数组
        effective_longitudes = np.radians(longitudes) * radius

        # 对于每个纬度，计算梯度
        result = np.zeros_like(field)

        # 循环处理每个纬度（这部分实现取决于具体数据结构，可能需要调整）
        for i in range(len(latitudes)):
            idx = [slice(None)] * field.ndim
            idx[lat_axis] = i

            # 调整当前纬度的经度距离
            current_effective_dist = effective_longitudes * cos_lat[i]

            # 计算当前纬度的梯度
            result[tuple(idx)] = calculate_gradient(
                field[tuple(idx)], current_effective_dist, axis=lon_axis, radius=1.0)

        return result
