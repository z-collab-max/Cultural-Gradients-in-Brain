"""Assign origin coordinates to provinces and aggregate scanning centers."""
from pathlib import Path
from collections import Counter
import csv,math,struct
import numpy as np
import h5py

def read_poly_shapefile(path: Path) -> list[dict]:
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

def read_dbf(path: Path, encoding: str='gbk') -> list[dict[str, str]]:
    with path.open('rb') as f:
        h = f.read(32)
        n_records = struct.unpack_from('<I', h, 4)[0]
        header_len = struct.unpack_from('<H', h, 8)[0]
        record_len = struct.unpack_from('<H', h, 10)[0]
        fields: list[tuple[str, int]] = []
        while f.tell() < header_len - 1:
            d = f.read(32)
            if not d or d[0] == 13:
                break
            name = d[:11].split(b'\x00', 1)[0].decode('ascii', 'ignore')
            fields.append((name, d[16]))
        f.seek(header_len)
        out: list[dict[str, str]] = []
        for _ in range(n_records):
            raw = f.read(record_len)
            if len(raw) != record_len:
                break
            if raw[:1] == b'*':
                out.append({})
                continue
            row: dict[str, str] = {}
            pos = 1
            for (name, width) in fields:
                value = raw[pos:pos + width]
                pos += width
                row[name] = value.rstrip(b' \x00').decode(encoding, 'replace')
                if name in ['ADCODE99','ADCODE93'] and row[name]:row[name]=str(int(float(row[name])))
            out.append(row)
    return out

def points_in_ring(x: np.ndarray, y: np.ndarray, ring: np.ndarray) -> np.ndarray:
    if len(ring) < 3:
        return np.zeros(len(x), dtype=bool)
    inside = np.zeros(len(x), dtype=bool)
    (x0, y0) = ring[-1]
    for (x1, y1) in ring:
        crossing = (y1 > y) != (y0 > y)
        denom = y0 - y1
        if denom != 0:
            x_at_y = (x0 - x1) * (y - y1) / denom + x1
            inside ^= crossing & (x < x_at_y)
        (x0, y0) = (x1, y1)
    return inside

def assign_province_codes(lon: np.ndarray, lat: np.ndarray, polygons: list[dict], dbf_rows: list[dict]) -> np.ndarray:
    if len(polygons) != len(dbf_rows):
        raise ValueError('Province SHP and DBF record counts do not match')
    codes = np.full(len(lon), '', dtype='U12')
    for (shape, attrs) in zip(polygons, dbf_rows):
        if not shape['parts']:
            continue
        code = attrs.get('ADCODE99') or attrs.get('ADCODE93') or ''
        if not code:
            continue
        (xmin, ymin, xmax, ymax) = shape['bbox']
        idx = np.flatnonzero((codes == '') & (lon >= xmin) & (lon <= xmax) & (lat >= ymin) & (lat <= ymax))
        if not len(idx):
            continue
        hit = np.zeros(len(idx), dtype=bool)
        for ring in shape['parts']:
            hit ^= points_in_ring(lon[idx], lat[idx], ring)
        codes[idx[hit]] = code
    valid_shapes = []
    for (shape, attrs) in zip(polygons, dbf_rows):
        code = attrs.get('ADCODE99') or attrs.get('ADCODE93') or ''
        if code and shape['parts']:
            valid_shapes.append((shape, code))
    fallback_distances = []
    far_outliers = []
    for idx in np.flatnonzero(codes == ''):
        (x, y) = (float(lon[idx]), float(lat[idx]))
        bbox_d2 = []
        for (shape, _) in valid_shapes:
            (xmin, ymin, xmax, ymax) = shape['bbox']
            dx = max(xmin - x, 0.0, x - xmax)
            dy = max(ymin - y, 0.0, y - ymax)
            bbox_d2.append(dx * dx + dy * dy)
        take = min(16, len(valid_shapes))
        candidates = np.argpartition(np.asarray(bbox_d2), take - 1)[:take]
        (best_code, best_d2) = ('', float('inf'))
        lon_scale = max(math.cos(math.radians(y)), 0.35)
        for j in candidates:
            (shape, code) = valid_shapes[int(j)]
            for ring in shape['parts']:
                d2 = ((ring[:, 0] - x) * lon_scale) ** 2 + (ring[:, 1] - y) ** 2
                local = float(d2.min())
                if local < best_d2:
                    (best_code, best_d2) = (code, local)
        best_distance = math.sqrt(best_d2)
        if best_distance <= 2.0:
            codes[idx] = best_code
            fallback_distances.append(best_distance)
        else:
            far_outliers.append((idx, best_distance))
    if fallback_distances:
        print(f'China: nearest-boundary assignment used for {len(fallback_distances)} coordinates (max angular distance {max(fallback_distances):.3f} deg)')
    if far_outliers:
        print(f'China: {len(far_outliers)} far coordinate outlier(s) left for center-modal assignment')
    return codes

def load_china_data(mat_path: Path, province_poly_path: Path) -> dict:
    with h5py.File(mat_path, 'r') as f:
        g = f['brain_data']
        lat = np.ravel(g['latitude_18y'][...]).astype(float)
        lon = np.ravel(g['longitude_18y'][...]).astype(float)
        center = np.ravel(g['centerID_rerange'][...]).astype(float)
    keep = np.isfinite(lat) & np.isfinite(lon) & np.isfinite(center)
    (lat, lon, center) = (lat[keep], lon[keep], center[keep].astype(int))
    polygons = read_poly_shapefile(province_poly_path)
    attrs = read_dbf(province_poly_path.with_suffix('.dbf'))
    province = assign_province_codes(lon, lat, polygons, attrs)
    far = np.flatnonzero(province == '')
    for idx in far:
        peer = (center == center[idx]) & (province != '')
        if np.any(peer):
            province[idx] = Counter(province[peer].tolist()).most_common(1)[0][0]
    if len(far):
        print(f'China: center-modal province used for {len(far)} far outlier(s)')
    assigned = province != ''
    if not assigned.all():
        raise RuntimeError(f'Only {assigned.sum()}/{len(assigned)} Chinese coordinates were assigned to province polygons; check the CRS/map file.')
    prov_codes = sorted(set(province))
    center_ids = sorted(set(center.tolist()))
    prov_index = {v: i for (i, v) in enumerate(prov_codes)}
    center_index = {v: i for (i, v) in enumerate(center_ids)}
    counts = np.zeros((len(prov_codes), len(center_ids)), dtype=int)
    prov_lat = np.zeros(len(prov_codes))
    prov_lon = np.zeros(len(prov_codes))
    for (p, code) in enumerate(prov_codes):
        m = province == code
        prov_lat[p] = float(lat[m].mean())
        prov_lon[p] = float(lon[m].mean())
    for (pcode, cid) in zip(province, center):
        counts[prov_index[pcode], center_index[int(cid)]] += 1
    center_lat = np.zeros(len(center_ids))
    center_lon = np.zeros(len(center_ids))
    for c in range(len(center_ids)):
        p = int(np.argmax(counts[:, c]))
        (center_lat[c], center_lon[c]) = (prov_lat[p], prov_lon[p])
    totals = counts.sum(axis=1)
    visible = totals > 0
    return {'n': int(len(lat)), 'center_ids': center_ids, 'center_lat': center_lat, 'center_lon': center_lon, 'prov_codes': np.asarray(prov_codes)[visible], 'prov_lat': prov_lat[visible], 'prov_lon': prov_lon[visible], 'counts': counts[visible], 'totals': totals[visible], 'n_excluded_small': int(totals[~visible].sum()), 'province_shapes': polygons}

def load_usa(path: Path) -> list[dict]:
    groups: dict[str, dict] = {}
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            state = row['state']
            g = groups.setdefault(state, {'name': state.replace('_', ' '), 'n': 0, 'lat': 0.0, 'lon': 0.0, 'datasets': set()})
            g['n'] += 1
            g['lat'] += float(row['latitude'])
            g['lon'] += float(row['longitude'])
            g['datasets'].add(row['dataset'])
    out = []
    for g in groups.values():
        out.append({'name': g['name'], 'n': g['n'], 'lat': g['lat'] / g['n'], 'lon': g['lon'] / g['n'], 'sites': len(g['datasets'])})
    return out
