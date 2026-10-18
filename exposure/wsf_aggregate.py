"""WSF Tracker (MindEarth/DLR, 10 m, first-detection half-year 1..20 = 2016-07 .. 2026-01)
aggregated onto the 250 m TPI grid: settlement area per cell per period, in m²."""
import sys, time
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules

from config import *
from proj import tm_forward

SRC = WSF_NPZ
X0, Y0, CELL, H, W = 414390.99421868543, 4600180.6948133595, 250.0, 777, 1542
NP = 20
t0 = time.time()
z = np.load(SRC)
data, xs, ys = z["data"], z["x"].astype(np.float64), z["y"].astype(np.float64)
print("array", data.shape, data.dtype, "x step %.6f° y step %.6f°" % (xs[1] - xs[0], ys[0] - ys[1]))
dlon = abs(float(xs[1] - xs[0])); dlat = abs(float(ys[0] - ys[1]))
area = np.zeros(H * W * (NP + 1), dtype=np.float64)        # m² of settlement by (cell, period)
BR = 1000
npx = 0
for r0 in range(0, data.shape[0], BR):
    blk = data[r0:r0 + BR]
    rr, cc = np.nonzero(blk)
    if rr.size == 0:
        continue
    per = blk[rr, cc].astype(np.int64)
    lat = ys[r0 + rr]; lon = xs[cc]
    x, y = tm_forward(lon, lat)
    c = np.floor((x - X0) / CELL).astype(np.int64); r = np.floor((Y0 - y) / CELL).astype(np.int64)
    ok = (c >= 0) & (c < W) & (r >= 0) & (r < H) & (per >= 1) & (per <= NP)
    pix_area = (dlat * 111320.0) * (dlon * 111320.0 * np.cos(np.deg2rad(lat[ok])))
    flat = (r[ok] * W + c[ok]) * (NP + 1) + per[ok]
    area += np.bincount(flat, weights=pix_area, minlength=area.size)
    npx += int(ok.sum())
    if (r0 // BR) % 5 == 0:
        print("  rows %d/%d  %.0fs" % (r0 + BR, data.shape[0], time.time() - t0), flush=True)
area = area.reshape(H, W, NP + 1)[..., 1:]                   # drop index 0
np.savez_compressed(os.path.join(OUT_DIR, "wsf_area.npz"), area=area.astype(np.float32))
tot = area.sum(axis=2)
print("done %.0fs; settlement pixels assigned %d; settlement 2016-07: %.0f km², 2026-01: %.0f km² (+%.1f%%)"
      % (time.time() - t0, npx, area[..., 0].sum() / 1e6, tot.sum() / 1e6, 100 * (tot.sum() / area[..., 0].sum() - 1)))
