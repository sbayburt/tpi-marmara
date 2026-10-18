"""Derived layers from the WorldCover counts: class fractions, sea mask,
distance to coast.  Cross-checked against the published valid-cell mask."""
import sys
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
from scipy import ndimage

from config import *
from geo import read_geotiff, write_geotiff

d = np.load(os.path.join(OUT_DIR, "wc_counts.npz"), allow_pickle=True)
counts = d["counts"].astype(np.int64)            # H x W x K
classes = d["classes"]
H, W, K = counts.shape
tot = counts.sum(axis=2)
frac = np.where(tot[..., None] > 0, counts / np.maximum(tot, 1)[..., None], 0.0).astype(np.float32)
k80 = int(np.where(classes == 80)[0][0])
water = frac[..., k80]

# ---- sea mask: water-dominated cells connected to the grid edge or > 1,000 km2
wet = water > 0.5
lab, n = ndimage.label(wet, structure=[[0, 1, 0], [1, 1, 1], [0, 1, 0]])
sizes = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
edge_ids = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])))
edge_ids.discard(0)
sea = np.zeros_like(wet)
big = 0
for i in range(1, n + 1):
    if i in edge_ids or sizes[i - 1] * 0.0625 > 1000:
        sea |= lab == i
        if sizes[i - 1] * 0.0625 > 1000 and i not in edge_ids:
            big += 1
print("water components:", n, " sea components (edge-connected):", len(edge_ids),
      " large non-edge components counted as sea:", big)
print("sea cells: %d (%.0f km2)  lake/river water cells not sea: %d"
      % (sea.sum(), sea.sum() * 0.0625, (wet & ~sea).sum()))

# ---- distance to coast (cell centre to nearest sea cell, minus half a cell)
edt = ndimage.distance_transform_edt(~sea)
dist_m = np.clip((edt - 0.5) * 250.0, 0, None).astype(np.float32)
dist_m[sea] = -1

# ---- cross-check with the published valid cells
tpi, meta = read_geotiff(TPI_RASTER)
valid = tpi > 0
print("valid cells:", valid.sum(), " valid cells with no WorldCover pixel:", int((valid & (tot == 0)).sum()),
      " valid cells classified as sea:", int((valid & sea).sum()))
nz = valid & sea
if nz.any():
    print("  water fraction of those:", np.round(np.sort(water[nz])[-10:], 2))
print("valid cells within 1 km of coast: %d (%.1f%%)" % ((valid & (dist_m <= 1000)).sum(), 100 * (valid & (dist_m <= 1000)).sum() / valid.sum()))

np.savez_compressed(os.path.join(OUT_DIR, "derived.npz"), frac=frac, tot=tot.astype(np.int32),
                    sea=sea, dist_m=dist_m, classes=classes)
write_geotiff(os.path.join(OUT_DIR, "out_dist_to_coast_m.tif"), dist_m, meta)
write_geotiff(os.path.join(OUT_DIR, "out_builtup_fraction.tif"), frac[..., int(np.where(classes == 50)[0][0])], meta)
print("saved derived.npz")
