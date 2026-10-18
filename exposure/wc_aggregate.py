"""Aggregate ESA WorldCover 2021 (10 m, EPSG:4326) onto the published 250 m
TPI grid (EPSG:32635, 777 x 1542, origin 414390.994 / 4600180.695).

Output: wc_counts.npz with counts[row, col, k] = number of 10 m pixels of
class CLASSES[k] whose centre falls in the cell.  Pure numpy + tifffile.
"""
import sys, time
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import tifffile

from config import *
from proj import tm_forward

SRC = RAW_DIR + os.sep
TILES = {"N39E024": (24.0, 42.0), "N39E027": (27.0, 42.0), "N39E030": (30.0, 42.0)}
PX = 1.0 / 12000.0                      # degrees per WorldCover pixel
X0, Y0, CELL, H, W = 414390.99421868543, 4600180.6948133595, 250.0, 777, 1542
CLASSES = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100], dtype=np.uint8)
K = len(CLASSES)
LUT = np.full(256, 255, dtype=np.uint8)
for k, c in enumerate(CLASSES):
    LUT[c] = k

# geographic window generously covering the UTM grid rectangle
LAT_N, LAT_S, LON_W, LON_E = 41.65, 39.65, 25.85, 30.65
BR, BC = 1200, 3000                      # block size (rows, cols) in pixels

counts = np.zeros(H * W * K, dtype=np.int64)
t0 = time.time()
for name, (tx0, ty0) in TILES.items():
    path = WC_TILE % name
    r0 = max(0, int((ty0 - LAT_N) / PX)); r1 = min(36000, int(np.ceil((ty0 - LAT_S) / PX)))
    c0 = max(0, int((LON_W - tx0) / PX)); c1 = min(36000, int(np.ceil((LON_E - tx0) / PX)))
    if r1 <= r0 or c1 <= c0:
        continue
    with tifffile.TiffFile(path) as tf:
        page = tf.pages[0]
        arr = page.asarray(maxworkers=2)
    print(name, "decoded", arr.shape, "window rows %d-%d cols %d-%d  %.0fs" % (r0, r1, c0, c1, time.time() - t0), flush=True)
    npx = 0
    for rb in range(r0, r1, BR):
        rows = np.arange(rb, min(rb + BR, r1))
        lat = ty0 - (rows + 0.5) * PX
        for cb in range(c0, c1, BC):
            cols = np.arange(cb, min(cb + BC, c1))
            lon = tx0 + (cols + 0.5) * PX
            blk = arr[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
            kk = LUT[blk]
            valid = kk != 255
            if not valid.any():
                continue
            LON, LAT = np.meshgrid(lon, lat)
            x, y = tm_forward(LON[valid], LAT[valid])
            c = np.floor((x - X0) / CELL).astype(np.int64)
            r = np.floor((Y0 - y) / CELL).astype(np.int64)
            ok = (c >= 0) & (c < W) & (r >= 0) & (r < H)
            if not ok.any():
                continue
            flat = (r[ok] * W + c[ok]) * K + kk[valid][ok].astype(np.int64)
            counts += np.bincount(flat, minlength=H * W * K)
            npx += int(ok.sum())
    print(name, "pixels assigned: %d   %.0fs" % (npx, time.time() - t0), flush=True)
    del arr

counts = counts.reshape(H, W, K)
np.savez_compressed(os.path.join(OUT_DIR, "wc_counts.npz"), counts=counts.astype(np.uint32),
                    classes=CLASSES, x0=X0, y0=Y0, cell=CELL)
tot = counts.sum(axis=2)
print("done %.0fs; cells with any pixel: %d; pixels per cell median %d, min %d, max %d"
      % (time.time() - t0, int((tot > 0).sum()), int(np.median(tot[tot > 0])), int(tot[tot > 0].min()), int(tot.max())))
