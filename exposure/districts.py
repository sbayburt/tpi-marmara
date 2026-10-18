"""Rasterise the 128 district polygons of the published analysis onto the TPI
grid and validate against the district means reported in the paper."""
import sys
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules

from shp import read_shp_polygons, read_dbf
from config import *
from geo import read_geotiff, rasterize_rings
from proj import albers_to_utm

TPI = TPI_RASTER
SHP = DISTRICT_SHP

tpi, meta = read_geotiff(TPI)
valid = tpi > 0
feats = read_shp_polygons(SHP + ".shp")
names, rows = read_dbf(SHP + ".dbf")
utm_feats = []
for rings in feats:
    out = []
    for r in rings:
        x, y = albers_to_utm(r[:, 0], r[:, 1])
        out.append(np.column_stack([x, y]))
    utm_feats.append(out)
dist = rasterize_rings(utm_feats, meta, fill=-1)
print("district raster: cells assigned", int((dist >= 0).sum()), "of", dist.size)

# validation against the published zonal statistics (mean over cells > 0)
diffs, cnt_ok = [], 0
for i, r in enumerate(rows):
    m = (dist == i) & valid
    if m.sum() == 0:
        if r["tpi_hucre"] and r["tpi_hucre"] > 0:
            print("  no cells for", r["Adı"], "paper n=", r["tpi_hucre"])
        continue
    mean = float(tpi[m].mean())
    diffs.append(abs(mean - float(r["tpi_ort"])))
    cnt_ok += int(m.sum() == int(r["tpi_hucre"]))
diffs = np.array(diffs)
print("districts validated: %d; |mean diff| max %.4f, median %.4f; cell counts identical for %d districts"
      % (len(diffs), diffs.max(), np.median(diffs), cnt_ok))
np.savez_compressed(os.path.join(OUT_DIR, "districts.npz"), dist=dist.astype(np.int16),
                    names=np.array([r["Adı"] for r in rows], dtype=object),
                    tpi_ort=np.array([r["tpi_ort"] for r in rows], dtype=float),
                    tpi_hucre=np.array([r["tpi_hucre"] for r in rows], dtype=float),
                    gi_sinif=np.array([r["gi_sinif"] for r in rows], dtype=object),
                    alan_km2=np.array([r["alan_km2"] for r in rows], dtype=float))
