"""Minimal ESRI Shapefile reader (polygon geometry + dBASE attributes), numpy only."""
import struct
import numpy as np


def read_shp_polygons(path):
    """Return a list of features; each is a list of rings (Nx2 float arrays)."""
    with open(path, "rb") as fh:
        buf = fh.read()
    pos = 100  # main header
    feats = []
    n = len(buf)
    while pos < n:
        _num, clen = struct.unpack(">ii", buf[pos:pos + 8])
        pos += 8
        rec = buf[pos:pos + clen * 2]
        pos += clen * 2
        shp_type = struct.unpack("<i", rec[0:4])[0]
        if shp_type == 0:            # null shape
            feats.append([])
            continue
        if shp_type not in (5, 15, 25):   # Polygon, PolygonZ, PolygonM
            feats.append([])
            continue
        nparts, npoints = struct.unpack("<ii", rec[36:44])
        parts = np.frombuffer(rec[44:44 + 4 * nparts], dtype="<i4")
        off = 44 + 4 * nparts
        pts = np.frombuffer(rec[off:off + 16 * npoints],
                            dtype="<f8").reshape(npoints, 2)
        rings = []
        for i in range(nparts):
            a = parts[i]
            b = parts[i + 1] if i + 1 < nparts else npoints
            rings.append(pts[a:b].copy())
        feats.append(rings)
    return feats


def read_dbf(path, encoding="utf-8"):
    """Return (field_names, list_of_row_dicts) from a dBASE III/IV table."""
    with open(path, "rb") as fh:
        buf = fh.read()
    nrec, hlen, rlen = struct.unpack("<IHH", buf[4:12])
    fields, pos = [], 32
    while buf[pos] != 0x0D:
        raw = buf[pos:pos + 32]
        name = raw[:11].split(b"\x00")[0].decode(encoding, "replace")
        ftype = chr(raw[11])
        flen = raw[16]
        fdec = raw[17]
        fields.append((name, ftype, flen, fdec))
        pos += 32
    start = hlen
    rows = []
    for r in range(nrec):
        base = start + r * rlen
        if buf[base:base + 1] == b"*":        # deleted
            continue
        off = base + 1
        row = {}
        for (name, ftype, flen, fdec) in fields:
            raw = buf[off:off + flen]
            off += flen
            s = raw.decode(encoding, "replace").strip()
            if ftype in ("N", "F"):
                try:
                    row[name] = float(s) if (fdec or "." in s) else int(s)
                except ValueError:
                    row[name] = np.nan
            elif ftype == "L":
                row[name] = s.upper() in ("Y", "T")
            else:
                row[name] = s
            if row[name] == "":
                row[name] = None
        rows.append(row)
    return [f[0] for f in fields], rows


def read_shapefile(base, encoding="utf-8"):
    """base without extension -> (geoms, field_names, rows)."""
    geoms = read_shp_polygons(base + ".shp")
    names, rows = read_dbf(base + ".dbf", encoding)
    return geoms, names, rows
