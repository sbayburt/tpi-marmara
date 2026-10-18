"""Table 8: the ten districts carrying the most pressure on natural or semi-natural cover."""
import os, json
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
from config import *
from geo import read_geotiff
from scipy import stats

tpi, meta = read_geotiff(TPI_RASTER); valid = tpi > 0; high = valid & (tpi >= 0.6)
DS = np.load(os.path.join(OUT_DIR, "districts.npz"), allow_pickle=True); dist_id = DS["dist"]; names = list(DS["names"])
C = np.load(os.path.join(OUT_DIR, "cells.npz")); dom, terr, f_nat = C["dom"], C["terr"], C["f_nat"]
rows = []
for i, nm in enumerate(names):
    m = valid & (dist_id == i)
    if m.sum() == 0:
        continue
    NP = float((tpi[m] * f_nat[m]).sum()); share = NP / float(tpi[m].sum())
    hm = m & high; hn = hm & (terr > 0) & (dom >= 2) & (dom <= 5)        # dominant cover natural
    rows.append(dict(d=nm, n=int(m.sum()), mean=float(tpi[m].mean()), NP=NP, share=share, high=int(hm.sum()), hnat=int(hn.sum())))
rk_mean = {r["d"]: k + 1 for k, r in enumerate(sorted(rows, key=lambda r: -r["mean"]))}
rk_np = {r["d"]: k + 1 for k, r in enumerate(sorted(rows, key=lambda r: -r["NP"]))}
top = sorted(rows, key=lambda r: -r["NP"])[:10]
for r in top:
    print("%-14s n %5d mean %.3f (rk %3d)  NP %6.1f (rk %2d) share %.2f  high nat/all %3d/%3d" % (r["d"], r["n"], r["mean"], rk_mean[r["d"]], r["NP"], rk_np[r["d"]], r["share"], r["hnat"], r["high"]))
print("Spearman(mean TPI, NP) over %d districts: %.3f" % (len(rows), stats.spearmanr([r["mean"] for r in rows], [r["NP"] for r in rows]).statistic))
json.dump({"top10": [dict(r, rk_mean=rk_mean[r["d"]], rk_np=rk_np[r["d"]]) for r in top]}, open(os.path.join(OUT_DIR, "table8.json"), "w"), ensure_ascii=False, indent=1)
