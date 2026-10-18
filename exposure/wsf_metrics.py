"""Settlement change 2016-07 -> 2026-01 (WSF Tracker) across the pressure gradient."""
import json, sys
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules

from config import *
from geo import read_geotiff

OUT = OUT_DIR + os.sep
tpi, meta = read_geotiff(TPI_RASTER)
valid = tpi > 0; high = valid & (tpi >= 0.6)
D = np.load(OUT + "derived.npz"); sea, tot = D["sea"], D["tot"]
dist_m = D["dist_m"].copy(); dist_m[valid & (dist_m < 0)] = 0
DS = np.load(OUT + "districts.npz", allow_pickle=True); dist_id = DS["dist"]; names = list(DS["names"])
M = np.load(OUT + "pa_masks.npz"); marine = M["marine"]
A = np.load(OUT + "wsf_area.npz")["area"].astype(np.float64)       # H x W x 20, m²
s16 = A[..., 0]; s26 = A.sum(axis=2); s21 = A[..., :11].sum(axis=2); new = s26 - s16
region_land = (dist_id >= 0) & ~sea & (tot > 0)
land0 = region_land & ~valid                                        # land cells with TPI = 0

R = {"strata": [], "bands": []}
lines = []
P = lines.append
tot_new_region = new[region_land].sum()
P("Region land: settlement 2016 %.0f km², 2021 %.0f km², 2026 %.0f km²; new 2016-2026 %.0f km² (+%.1f%%)" %
  (s16[region_land].sum() / 1e6, s21[region_land].sum() / 1e6, s26[region_land].sum() / 1e6, tot_new_region / 1e6, 100 * tot_new_region / s16[region_land].sum()))
P("Valid cells (TPI > 0): new %.0f km² = %.1f%% of regional new settlement; valid cells are %.1f%% of region land" %
  (new[valid & region_land].sum() / 1e6, 100 * new[valid & region_land].sum() / tot_new_region, 100 * (valid & region_land).sum() / region_land.sum()))

CLS = [("Low (0.00–0.20)", 0.0, 0.2), ("Moderate (0.20–0.40)", 0.2, 0.4), ("Moderate-high (0.40–0.60)", 0.4, 0.6),
       ("High (0.60–0.80)", 0.6, 0.8), ("Very high (0.80–1.00)", 0.8, 1.0001)]
P("\n%-28s %7s %8s %8s %8s %7s %8s %9s %9s" % ("stratum", "cells", "s2016km2", "s2026km2", "new km2", "growth%", "%newReg", "built16%", "built26%"))
def row(label, m):
    a16, a26 = s16[m].sum() / 1e6, s26[m].sum() / 1e6; nw = a26 - a16; cell_area = m.sum() * 0.0625
    d = dict(label=label, cells=int(m.sum()), s2016=round(a16, 1), s2026=round(a26, 1), new=round(nw, 1),
             growth_pct=round(100 * nw / a16, 1) if a16 > 0 else None, share_new=round(100 * (new[m].sum()) / tot_new_region, 1),
             built16_pct=round(100 * a16 / cell_area, 1), built26_pct=round(100 * a26 / cell_area, 1), land_pct=round(100 * (m & region_land).sum() / region_land.sum(), 1))
    P("%-28s %7d %8.1f %8.1f %8.1f %7.1f %8.1f %9.1f %9.1f" % (label, d["cells"], a16, a26, nw, d["growth_pct"] or 0, d["share_new"], d["built16_pct"], d["built26_pct"]))
    return d
for name, lo, hi in CLS:
    R["strata"].append(row(name, valid & (tpi >= lo) & (tpi < hi)))
R["strata"].append(row("All valid cells", valid))
R["strata"].append(row("Land cells with TPI = 0", land0))
R["strata"].append(row("Regional land, total", region_land))
P("")
for name, lo, hi in [("0–1 km", 0, 1000), ("1–3 km", 1000, 3000), ("3–10 km", 3000, 10000), ("> 10 km", 10000, 1e12)]:
    R["bands"].append(row(name, region_land & (dist_m >= lo) & (dist_m < hi)))
# high cells only, and high cells in band 0-1 km
hb = row("High cells (TPI >= 0.60)", high); R["high"] = hb
R["high_coastal"] = row("High cells within 1 km", high & (dist_m < 1000))
R["sepa_islands"] = row("Adalar + Erdek (SEPA land)", valid & marine & np.isin(dist_id, [names.index("Adalar"), names.index("Erdek")]))
R["sile_kandira"] = row("Şile + Kandıra", valid & np.isin(dist_id, [names.index("Şile"), names.index("Kandıra")]))
# new-settlement density per km² of cell area by class, and per high cell
P("\nnew settlement per km² of cell area: " + ", ".join("%s %.1f ha/km²" % (r["label"][:9], 100 * r["new"] / (r["cells"] * 0.0625)) for r in R["strata"][:6]))
# half-year series for region and for high cells (for text: is growth slowing?)
ser_r = A[region_land].sum(axis=0) / 1e6; ser_h = A[high].sum(axis=0) / 1e6
P("region new settlement per half-year (km²), periods 2..20: " + " ".join("%.1f" % v for v in ser_r[1:]))
P("high cells new per half-year (km²): " + " ".join("%.2f" % v for v in ser_h[1:]))
# districts: top 10 by new settlement in valid cells, and the ten Table-8 districts
dn = []
for i, nm in enumerate(names):
    m = region_land & (dist_id == i)
    if m.sum() == 0: continue
    dn.append((float(new[m].sum() / 1e6), nm, float(100 * new[m].sum() / max(s16[m].sum(), 1)), float(new[m & valid].sum() / 1e6)))
dn.sort(reverse=True)
P("\nTop 10 districts by new settlement 2016-2026 (km², growth %, of which in valid cells):")
for v, nm, g, vv in dn[:10]: P("  %-16s %6.1f km²  +%5.1f%%  in valid cells %5.1f" % (nm, v, g, vv))
R["districts_top"] = [dict(district=nm, new=round(v, 1), growth=round(g, 1), new_valid=round(vv, 1)) for v, nm, g, vv in dn[:10]]
open(OUT + "wsf_metrics.txt", "w").write("\n".join(lines)); json.dump(R, open(OUT + "wsf_metrics.json", "w"), ensure_ascii=False, indent=1)
print("\n".join(lines))
