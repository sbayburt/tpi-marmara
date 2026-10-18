"""Minimal geospatial readers: GeoPackage (SQLite+WKB) and GeoTIFF, numpy only."""
import sqlite3, struct
import numpy as np
import tifffile


# ---------- GeoPackage ----------
def gpkg_layers(path):
    con = sqlite3.connect(path)
    rows = con.execute(
        "SELECT table_name, data_type, srs_id FROM gpkg_contents").fetchall()
    con.close()
    return rows


def _wkb_points(blob):
    """Parse a GPKG geometry blob; return list of (x, y) for Point/MultiPoint."""
    if blob is None:
        return []
    # GPKG binary header: magic 'GP', version, flags, srs_id, envelope
    assert blob[:2] == b"GP", "not a GPKG geometry blob"
    flags = blob[3]
    env_code = (flags >> 1) & 0x07
    env_sizes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}
    hdr = 8 + env_sizes[env_code]
    wkb = blob[hdr:]
    return _parse_wkb(wkb)


def _parse_wkb(wkb):
    byte_order = wkb[0]
    e = "<" if byte_order == 1 else ">"
    gtype = struct.unpack(e + "I", wkb[1:5])[0]
    base = gtype % 1000
    has_z = (gtype // 1000) in (1, 3)
    has_m = (gtype // 1000) in (2, 3)
    ncoord = 2 + int(has_z) + int(has_m)
    pos = 5
    out = []
    if base == 1:  # Point
        vals = struct.unpack(e + "d" * ncoord, wkb[pos:pos + 8 * ncoord])
        out.append((vals[0], vals[1]))
    elif base == 4:  # MultiPoint
        n = struct.unpack(e + "I", wkb[pos:pos + 4])[0]
        pos += 4
        for _ in range(n):
            sub = wkb[pos:]
            out.extend(_parse_wkb(sub))
            pos += 1 + 4 + 8 * ncoord
    else:
        raise ValueError("geometry type %d not a point type" % base)
    return out


def read_points(path, table=None):
    """Return Nx2 float array of point coordinates from a GPKG point layer."""
    layers = gpkg_layers(path)
    if table is None:
        table = layers[0][0]
    con = sqlite3.connect(path)
    gcol = con.execute(
        "SELECT column_name FROM gpkg_geometry_columns WHERE table_name=?",
        (table,)).fetchone()[0]
    blobs = con.execute('SELECT "%s" FROM "%s"' % (gcol, table)).fetchall()
    con.close()
    pts = []
    for (b,) in blobs:
        pts.extend(_wkb_points(b))
    return np.asarray(pts, dtype=float)


def _parse_wkb_polys(wkb):
    """Return list of exterior rings (Nx2 arrays) for Polygon / MultiPolygon."""
    byte_order = wkb[0]
    e = "<" if byte_order == 1 else ">"
    gtype = struct.unpack(e + "I", wkb[1:5])[0]
    base = gtype % 1000
    has_z = (gtype // 1000) in (1, 3)
    has_m = (gtype // 1000) in (2, 3)
    nc = 2 + int(has_z) + int(has_m)
    pos = 5
    rings = []
    if base == 3:  # Polygon
        nring = struct.unpack(e + "I", wkb[pos:pos + 4])[0]
        pos += 4
        for r in range(nring):
            npt = struct.unpack(e + "I", wkb[pos:pos + 4])[0]
            pos += 4
            buf = np.frombuffer(wkb[pos:pos + 8 * nc * npt],
                                dtype=(e + "f8")).reshape(npt, nc)
            pos += 8 * nc * npt
            if r == 0:                      # exterior ring only
                rings.append(buf[:, :2].copy())
        return rings, pos
    elif base == 6:  # MultiPolygon
        npoly = struct.unpack(e + "I", wkb[pos:pos + 4])[0]
        pos += 4
        for _ in range(npoly):
            sub, used = _parse_wkb_polys(wkb[pos:])
            rings.extend(sub)
            pos += used
        return rings, pos
    raise ValueError("geometry type %d not a polygon type" % base)


def read_polygons(path, table=None, fields=()):
    """Return (list_of_ring_lists, dict_of_field_arrays)."""
    layers = gpkg_layers(path)
    if table is None:
        table = layers[0][0]
    con = sqlite3.connect(path)
    gcol = con.execute(
        "SELECT column_name FROM gpkg_geometry_columns WHERE table_name=?",
        (table,)).fetchone()[0]
    cols = ['"%s"' % gcol] + ['"%s"' % f for f in fields]
    rows = con.execute('SELECT %s FROM "%s"' % (",".join(cols), table)).fetchall()
    con.close()
    geoms, attrs = [], {f: [] for f in fields}
    for row in rows:
        blob = row[0]
        if blob is None:
            geoms.append([])
        else:
            assert blob[:2] == b"GP"
            flags = blob[3]
            env = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(flags >> 1) & 0x07]
            rings, _ = _parse_wkb_polys(blob[8 + env:])
            geoms.append(rings)
        for i, f in enumerate(fields):
            attrs[f].append(row[1 + i])
    return geoms, {k: np.array(v, dtype=object) for k, v in attrs.items()}


def table_columns(path, table):
    con = sqlite3.connect(path)
    cols = [r[1] for r in con.execute('PRAGMA table_info("%s")' % table)]
    con.close()
    return cols


def read_attrs(path, table, fields):
    con = sqlite3.connect(path)
    q = 'SELECT %s FROM "%s"' % (",".join('"%s"' % f for f in fields), table)
    rows = con.execute(q).fetchall()
    con.close()
    return rows


# ---------- GeoTIFF ----------
def read_geotiff(path):
    """Return (array, transform_dict) with origin / pixel size from GeoTIFF tags."""
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        arr = page.asarray()
        tags = {t.name: t.value for t in page.tags.values()}
    scale = tags.get("ModelPixelScaleTag") or tags.get("ModelPixelScale")
    tie = tags.get("ModelTiepointTag") or tags.get("ModelTiepoint")
    nodata = tags.get("GDAL_NODATA")
    meta = {"height": arr.shape[0], "width": arr.shape[1]}
    if scale and tie:
        meta["px"] = float(scale[0])
        meta["py"] = float(scale[1])
        meta["x0"] = float(tie[3])
        meta["y0"] = float(tie[4])
    if nodata is not None:
        try:
            meta["nodata"] = float(nodata)
        except (TypeError, ValueError):
            pass
    return arr, meta


def write_geotiff(path, arr, meta):
    """Write a single-band float32 GeoTIFF carrying the same grid definition."""
    extratags = [
        (33550, "d", 3, (meta["px"], meta["py"], 0.0)),
        (33922, "d", 6, (0.0, 0.0, 0.0, meta["x0"], meta["y0"], 0.0)),
    ]
    tifffile.imwrite(path, arr.astype(np.float32), extratags=extratags)


def sample_raster(arr, meta, xy):
    """Nearest-cell sample of a raster at Nx2 map coordinates. NaN when outside."""
    x = np.asarray(xy[:, 0]); y = np.asarray(xy[:, 1])
    col = np.floor((x - meta["x0"]) / meta["px"]).astype(int)
    row = np.floor((meta["y0"] - y) / meta["py"]).astype(int)
    ok = (col >= 0) & (col < meta["width"]) & (row >= 0) & (row < meta["height"])
    out = np.full(x.shape, np.nan)
    out[ok] = arr[row[ok], col[ok]]
    return out


def rasterize_rings(ring_lists, meta, fill=-1):
    """Burn polygon indices onto the grid using matplotlib point-in-polygon.

    Even-odd rule across a feature's exterior rings; later features win overlaps.
    """
    from matplotlib.path import Path
    h, w = meta["height"], meta["width"]
    out = np.full((h, w), fill, dtype=np.int32)
    for idx, rings in enumerate(ring_lists):
        if not rings:
            continue
        xs = np.concatenate([r[:, 0] for r in rings])
        ys = np.concatenate([r[:, 1] for r in rings])
        c0 = max(0, int(np.floor((xs.min() - meta["x0"]) / meta["px"])))
        c1 = min(w - 1, int(np.ceil((xs.max() - meta["x0"]) / meta["px"])))
        r0 = max(0, int(np.floor((meta["y0"] - ys.max()) / meta["py"])))
        r1 = min(h - 1, int(np.ceil((meta["y0"] - ys.min()) / meta["py"])))
        if c1 < c0 or r1 < r0:
            continue
        cc, rr = np.meshgrid(np.arange(c0, c1 + 1), np.arange(r0, r1 + 1))
        px = meta["x0"] + (cc + 0.5) * meta["px"]
        py = meta["y0"] - (rr + 0.5) * meta["py"]
        pts = np.column_stack([px.ravel(), py.ravel()])
        inside = np.zeros(pts.shape[0], dtype=bool)
        for ring in rings:
            inside ^= Path(ring).contains_points(pts)
        if inside.any():
            sub = out[r0:r1 + 1, c0:c1 + 1]
            sub[inside.reshape(cc.shape)] = idx
            out[r0:r1 + 1, c0:c1 + 1] = sub
    return out
