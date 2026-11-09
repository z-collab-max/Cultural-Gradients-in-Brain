from pathlib import Path
import json,re,struct,math
import numpy as np

def read_poly_shapefile(path: Path) -> list[dict]:
    """Read Polygon/Polyline records needed by this figure.

    Returns one item per SHP record so it stays aligned with the DBF table.
    Coordinates are numpy arrays with shape (n, 2), one per part/ring.
    """
    records: list[dict] = []
    with path.open('rb') as f:
        header = f.read(100)
        if len(header) != 100 or struct.unpack('>i', header[:4])[0] != 9994:
            raise ValueError(f'Not a valid ESRI shapefile: {path}')
        while True:
            rh = f.read(8)
            if not rh:
                break
            if len(rh) != 8:
                raise ValueError(f'Truncated shapefile record header: {path}')
            (_, words) = struct.unpack('>2i', rh)
            content = f.read(words * 2)
            if len(content) != words * 2:
                raise ValueError(f'Truncated shapefile record: {path}')
            shape_type = struct.unpack_from('<i', content, 0)[0]
            if shape_type == 0:
                records.append({'bbox': None, 'parts': []})
                continue
            if shape_type not in (3, 5):
                raise ValueError(f'Unsupported shape type {shape_type} in {path}')
            bbox = struct.unpack_from('<4d', content, 4)
            (n_parts, n_points) = struct.unpack_from('<2i', content, 36)
            part_starts = list(struct.unpack_from(f'<{n_parts}i', content, 44))
            point_offset = 44 + 4 * n_parts
            pts = np.frombuffer(content, dtype='<f8', count=2 * n_points, offset=point_offset).reshape(-1, 2).copy()
            bounds = part_starts + [n_points]
            parts = [pts[bounds[i]:bounds[i + 1]] for i in range(n_parts)]
            records.append({'bbox': bbox, 'parts': parts})
    return records

def decode_echarts_geometry(g: dict) -> list[list[np.ndarray]]:

    def zig(v: int) -> int:
        return v >> 1 ^ -(v & 1)

    def decode_ring(s: str, offset: list[int]) -> np.ndarray:
        (px, py) = offset
        pts = []
        for i in range(0, len(s), 2):
            px += zig(ord(s[i]) - 64)
            py += zig(ord(s[i + 1]) - 64)
            pts.append((px / 1024.0, py / 1024.0))
        return np.asarray(pts, dtype=float)
    if 'encodeOffsets' not in g:
        if g['type'] == 'Polygon':
            return [[np.asarray(r, float) for r in g['coordinates']]]
        return [[np.asarray(r, float) for r in p] for p in g['coordinates']]
    if g['type'] == 'Polygon':
        return [[decode_ring(s, g['encodeOffsets'][i]) for (i, s) in enumerate(g['coordinates'])]]
    return [[decode_ring(s, g['encodeOffsets'][pi][i]) for (i, s) in enumerate(poly)] for (pi, poly) in enumerate(g['coordinates'])]

def load_usa_geojson(js_path: Path) -> dict:
    text = js_path.read_text(encoding='utf-8')
    call = re.search('echarts\\.registerMap\\([^,]+,\\s*', text)
    if call is None:
        raise ValueError(f'Could not find registerMap data in {js_path}')
    start = text.index('{', call.end())
    (data, _) = json.JSONDecoder(strict=False).raw_decode(text, start)
    for feature in data['features']:
        decoded = decode_echarts_geometry(feature['geometry'])
        feature['geometry'] = {'type': 'MultiPolygon', 'coordinates': [[ring.tolist() for ring in poly] for poly in decoded]}
    return data

def albers(lon, lat, lon0, lat0, phi1, phi2):
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    lam = np.deg2rad(lon)
    phi = np.deg2rad(lat)
    (lam0, phi0) = (math.radians(lon0), math.radians(lat0))
    (p1, p2) = (math.radians(phi1), math.radians(phi2))
    n = 0.5 * (math.sin(p1) + math.sin(p2))
    c = math.cos(p1) ** 2 + 2 * n * math.sin(p1)
    rho = np.sqrt(np.maximum(0.0, c - 2 * n * np.sin(phi))) / n
    rho0 = math.sqrt(c - 2 * n * math.sin(phi0)) / n
    theta = n * (lam - lam0)
    return (rho * np.sin(theta), rho0 - rho * np.cos(theta))

ROBIN_X = np.asarray([1.0, 0.9986, 0.9954, 0.99, 0.9822, 0.973, 0.96, 0.9427, 0.9216, 0.8962, 0.8679, 0.835, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322])

ROBIN_Y = np.asarray([0.0, 0.062, 0.124, 0.186, 0.248, 0.31, 0.372, 0.434, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0])

def robinson(lon, lat):
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    a = np.clip(np.abs(lat), 0, 90) / 5.0
    i = np.minimum(np.floor(a).astype(int), 17)
    f = a - i
    xx = ROBIN_X[i] * (1 - f) + ROBIN_X[i + 1] * f
    yy = ROBIN_Y[i] * (1 - f) + ROBIN_Y[i + 1] * f
    return (0.8487 * np.deg2rad(lon) * xx, 1.3523 * np.sign(lat) * yy)