"""Miller projection and grid indexing."""
import numpy as np

def lon_lat_2dikaer(lon, lat):
    """
    米勒圆柱投影 (Miller Cylindrical Projection)
    将经纬度 (lon, lat) 转换为笛卡尔平面坐标 (X, Y) [单位: km]
    """
    M_PI = 3.14159265358979323846
    L = 6381372.0 * M_PI * 2.0  # 地球周长 (米)
    W = L                       # 展开 x 轴总长
    H = L / 2.0                 # y 轴总长
    mill = 2.3                  # 米勒投影常数

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    x = lon * M_PI / 180.0
    y = lat * M_PI / 180.0
    y1 = 1.25 * np.log(np.tan(0.25 * M_PI + 0.4 * y))

    dikaerX = (W / 2.0) + (W / (2.0 * M_PI)) * x
    dikaerY = (H / 2.0) - (H / (2.0 * mill)) * y1

    return dikaerX / 1000.0, dikaerY / 1000.0


def xy_to_coor(x_km, y_km):
    """
    米勒圆柱投影逆变换
    将笛卡尔平面坐标 (X, Y) [单位: km] 还原为经纬度 (lon, lat)
    """
    M_PI = 3.14159265358979323846
    L = 6381372.0 * M_PI * 2.0
    W = L
    H = L / 2.0
    mill = 2.3

    x = np.asarray(x_km, dtype=float) * 1000.0
    y = np.asarray(y_km, dtype=float) * 1000.0

    lat_rad = ((H / 2.0 - y) * 2.0 * mill) / (1.25 * H)
    lat = ((np.arctan(np.exp(lat_rad)) - 0.25 * M_PI) * 180.0) / (0.4 * M_PI)
    lon = (x - W / 2.0) * 360.0 / W

    return lon, lat


def get_matlab_grid_idx(grid_ticks, values):
    """
    严格对齐 MATLAB 的网格索引分配逻辑:
    xgrididx = find((lengrid - xlon) > 0); grididx = xgrididx(1);
    返回 1-based 网格索引 (1 <= idx <= len(grid_ticks))
    """
    # np.searchsorted(side='right') + 1 严格等于 MATLAB find((ticks - val) > 0, 1)
    indices = np.searchsorted(grid_ticks, values, side='right') + 1
    # 上限保护不超过 ticks 长度
    indices = np.clip(indices, 1, len(grid_ticks))
    return indices


