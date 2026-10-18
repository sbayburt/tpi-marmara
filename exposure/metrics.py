"""Exposure metrics for the EMA revision (section 4.5).

A  land-cover composition of each TPI class (terrestrial pixels)
B  coastal distance bands
C  Ramsar sites in the study area (the only WDPA polygons present)
D  district 'pressure on natural cover'
E  relation between TPI and built-up fraction
"""
import sys, json
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
from scipy import stats, ndimage

from config import *
from geo import read_geotiff, rasterize_rings
from shp import read_shp_polygons, read_dbf
from proj import tm_forward

OUT = OUT_DIR + os.sep
tpi, meta = read_geotiff(TPI_RASTER)
valid = tpi > 0
D = np.load(OUT + "derived.npz", allow_pickle=True)
frac, tot, sea, dist_m, classes = D["frac"], D["tot"], D["sea"], D["dist_m"].copy(), list(D["classes"])
dist_m[valid & (dist_m < 0)] = 0.0        # valid cells lying on water at 10 m are coastal by definition
DS = np.load(OUT + "districts.npz", allow_pickle=True)
dist_id, dnames = DS["dist"], DS["names"]

kidx = {c: classes.index(c) for c in classes}
GROUPS = [("Built-up", [50]), ("Cropland", [40]), ("Tree cover", [10]),
          ("Grassland / shrubland", [20, 30]), ("Bare / sparse (beach, dune, rock)", [60]),
          ("Herbaceous wetland", [90]), ("Inland water", [80])]
NATURAL = [10, 20, 30, 60, 90]
counts = (frac * tot[..., None]).round().astype(np.int64)          # back to pixel counts
terr = tot - counts[..., kidx[80]]                                   # terrestrial pixels per cell
f_nat = np.where(terr > 0, sum(counts[..., kidx[c]] for c in NATURAL) / np.maximum(terr, 1), 0.0)
f_built = np.where(terr > 0, counts[..., kidx[50]] / np.maximum(terr, 1), 0.0)

CLS = [("Low (0.00–0.20)", 0.0, 0.2), ("Moderate (0.20–0.40)", 0.2, 0.4),
       ("Moderate-high (0.40–0.60)", 0.4, 0.6), ("High (0.60–0.80)", 0.6, 0.8),
       ("Very high (0.80–1.00)", 0.8, 1.0001)]
high = valid & (tpi >= 0.6)
R = {}
lines = []
P = lines.append

# ---------------------------------------------------------------- water-only cells
wonly = valid & (terr == 0)
P("Valid cells %d; entirely water at 10 m: %d (%.2f%%), of which TPI>=0.6: %d; max TPI among them %.3f"
  % (valid.sum(), wonly.sum(), 100 * wonly.sum() / valid.sum(), (wonly & high).sum(), tpi[wonly].max()))
R["water_only_cells"] = int(wonly.sum()); R["water_only_high"] = int((wonly & high).sum())

# ---------------------------------------------------------------- A composition
P("\nA. Land-cover composition of terrestrial pixels, by TPI class (per cent)")
hdr = "%-28s %8s" % ("class", "cells") + "".join("%10s" % g[0][:9] for g in GROUPS if g[0] != "Inland water") + "%10s" % "natural"
P(hdr)
R["A"] = []
for name, lo, hi in CLS + [("All valid cells", 0.0, 1.0001)]:
    m = valid & (tpi >= lo) & (tpi < hi)
    T = terr[m].sum()
    row = {"class": name, "cells": int(m.sum())}
    s = "%-28s %8d" % (name, m.sum())
    for g, codes in GROUPS:
        if g == "Inland water":
            continue
        share = 100 * sum(counts[..., kidx[c]][m].sum() for c in codes) / T
        row[g] = round(float(share), 1); s += "%10.1f" % share
    nat = 100 * sum(counts[..., kidx[c]][m].sum() for c in NATURAL) / T
    row["natural"] = round(float(nat), 1); s += "%10.1f" % nat
    R["A"].append(row); P(s)

# cell-level view: dominant cover of high-pressure cells, and low-built-up high cells
dom = np.argmax(np.stack([sum(counts[..., kidx[c]] for c in codes) for g, codes in GROUPS[:-1]], axis=-1), axis=-1)
P("\nDominant terrestrial cover of the %d high + very-high cells (cells):" % high.sum())
R["A_dom_high"] = {}
for i, (g, _) in enumerate(GROUPS[:-1]):
    n = int((high & (terr > 0) & (dom == i)).sum()); R["A_dom_high"][g] = n
    P("   %-36s %5d  (%.1f%%)" % (g, n, 100 * n / (high & (terr > 0)).sum()))
lowb = high & (terr > 0) & (f_built < 0.25)
P("High-pressure cells with built-up fraction < 0.25: %d (%.1f%% of terrestrial high cells)" % (lowb.sum(), 100 * lowb.sum() / (high & (terr > 0)).sum()))
R["A_high_lowbuilt"] = int(lowb.sum())

# ---------------------------------------------------------------- B coastal bands
P("\nB. Coastal distance bands (distance from cell centre to nearest sea cell)")
BANDS = [("0–1 km", 0, 1000), ("1–3 km", 1000, 3000), ("3–10 km", 3000, 10000), ("> 10 km", 10000, 1e12)]
mass = tpi[valid].sum()
region_land = (dist_id >= 0) & ~sea & (tot > 0)          # land inside the 128 district polygons
P("%-10s %8s %8s %8s %9s %9s %10s %9s %9s %9s" % ("band", "%land", "cells", "%cells", "meanTPI", "%mass", "high", "%high", "%natural", "%natdom"))
R["B"] = []
for name, lo, hi in BANDS:
    m = valid & (dist_m >= lo) & (dist_m < hi)
    hm = m & high
    T = terr[hm].sum()
    nat = 100 * sum(counts[..., kidx[c]][hm].sum() for c in NATURAL) / T if T > 0 else float("nan")
    hmt = hm & (terr > 0)
    natdom = 100 * (hmt & (dom >= 2) & (dom <= 5)).sum() / max(hmt.sum(), 1)
    row = dict(band=name, pct_land=round(100 * (region_land & (dist_m >= lo) & (dist_m < hi)).sum() / region_land.sum(), 1),
               cells=int(m.sum()), pct_cells=round(100 * m.sum() / valid.sum(), 1),
               mean_tpi=round(float(tpi[m].mean()), 3), pct_mass=round(100 * tpi[m].sum() / mass, 1),
               high=int(hm.sum()), pct_high=round(100 * hm.sum() / high.sum(), 1), nat_high=round(float(nat), 1),
               natdom_high=round(float(natdom), 1), mean_built_high=round(float(f_built[hmt].mean()), 2))
    R["B"].append(row)
    P("%-10s %8.1f %8d %8.1f %9.3f %9.1f %10d %9.1f %9.1f %9.1f" % (name, row["pct_land"], row["cells"], row["pct_cells"], row["mean_tpi"], row["pct_mass"], row["high"], row["pct_high"], row["nat_high"], row["natdom_high"]))
R["region_land_km2"] = round(float(region_land.sum() * 0.0625), 0)
# land share within 1 km for context: all land cells (not only valid)
land = (~sea) & (tot > 0)
P("Context: of all land cells in the grid, %.1f%% lie within 1 km of the coast; of valid cells %.1f%%; of high cells %.1f%%"
  % (100 * (land & (dist_m <= 1000)).sum() / land.sum(), 100 * (valid & (dist_m < 1000)).sum() / valid.sum(), 100 * (high & (dist_m < 1000)).sum() / high.sum()))
R["B_land_within1km"] = round(float(100 * (land & (dist_m <= 1000)).sum() / land.sum()), 1)

# ---------------------------------------------------------------- C Ramsar sites
P("\nC. Ramsar sites inside the study area (WDPA)")
import glob
ram = []
for d in sorted(glob.glob(os.path.join(WDPA_DIR, "*", ""))):
    g = read_shp_polygons(d + "WDPA_WDOECM_Oct2026_Public_TUR_shp-polygons.shp")
    names, rows = read_dbf(d + "WDPA_WDOECM_Oct2026_Public_TUR_shp-polygons.dbf")
    for feat, r in zip(g, rows):
        xs = np.concatenate([q[:, 0] for q in feat])
        if 25.9 < xs.mean() < 30.6 and r["NAME"] in ("Manyas", "Lake Uluabat"):
            utm = []
            for q in feat:
                x, y = tm_forward(q[:, 0], q[:, 1]); utm.append(np.column_stack([x, y]))
            ram.append((r["NAME"], utm))
pa = rasterize_rings([u for _, u in ram], meta, fill=-1)
R["C"] = []
for i, (nm, _) in enumerate(ram):
    inside = pa == i
    ring3 = (ndimage.distance_transform_edt(~inside) * 250 <= 3000) & ~inside
    vi, v3 = valid & inside, valid & ring3
    row = dict(name=nm, cells=int(inside.sum()), valid_in=int(vi.sum()), mean_in=round(float(tpi[vi].mean()), 3) if vi.any() else None,
               max_in=round(float(tpi[vi].max()), 3) if vi.any() else None, valid_3km=int(v3.sum()),
               mean_3km=round(float(tpi[v3].mean()), 3) if v3.any() else None, max_3km=round(float(tpi[v3].max()), 3) if v3.any() else None,
               high_in_3km=int((v3 & high).sum() + (vi & high).sum()))
    R["C"].append(row); P("   " + json.dumps(row, ensure_ascii=False))

# ---------------------------------------------------------------- D districts
P("\nD. District pressure on natural cover  (NP = sum TPI x natural fraction; share = NP / sum TPI)")
R["D"] = []
for i, nm in enumerate(dnames):
    m = valid & (dist_id == i)
    if m.sum() == 0:
        continue
    npres = float((tpi[m] * f_nat[m]).sum()); tsum = float(tpi[m].sum())
    R["D"].append(dict(district=str(nm), n=int(m.sum()), mean_tpi=round(float(tpi[m].mean()), 4), NP=round(npres, 2),
                       share_nat=round(npres / tsum, 3), high=int((m & high).sum()), high_nat=int((m & high & (f_nat >= 0.5)).sum()),
                       mean_fnat=round(float(f_nat[m].mean()), 3)))
Dd = sorted(R["D"], key=lambda r: -r["NP"])
rank_mean = {r["district"]: k + 1 for k, r in enumerate(sorted(R["D"], key=lambda r: -r["mean_tpi"]))}
P("%-18s %6s %8s %8s %8s %6s %8s %9s" % ("district", "n", "meanTPI", "rk_mean", "NP", "rk_NP", "shareNat", "highNat/high"))
for k, r in enumerate(Dd[:20]):
    P("%-18s %6d %8.3f %8d %8.1f %6d %8.2f %5d/%d" % (r["district"], r["n"], r["mean_tpi"], rank_mean[r["district"]], r["NP"], k + 1, r["share_nat"], r["high_nat"], r["high"]))
P("Districts with the highest share of pressure on natural cover (n >= 50 cells):")
for r in sorted([r for r in R["D"] if r["n"] >= 50], key=lambda r: -r["share_nat"])[:12]:
    P("   %-18s share %.2f  meanTPI %.3f (rank %d)  high %d" % (r["district"], r["share_nat"], r["mean_tpi"], rank_mean[r["district"]], r["high"]))
rs = stats.spearmanr([r["mean_tpi"] for r in R["D"]], [r["NP"] for r in R["D"]])
P("Spearman(mean TPI, NP) across %d districts = %.3f" % (len(R["D"]), rs.statistic))
R["D_rho"] = round(float(rs.statistic), 3)

# ---------------------------------------------------------------- E built-up relation
P("\nE. TPI vs built-up fraction over valid terrestrial cells")
mt = valid & (terr > 0)
rho = stats.spearmanr(tpi[mt], f_built[mt]).statistic
P("Spearman rho = %.3f (n = %d)" % (rho, mt.sum()))
for name, lo, hi in CLS:
    m = mt & (tpi >= lo) & (tpi < hi)
    P("   %-28s mean built-up %.2f  median %.2f" % (name, f_built[m].mean(), np.median(f_built[m])))
R["E_rho"] = round(float(rho), 3)
R["E_built_by_class"] = [round(float(f_built[mt & (tpi >= lo) & (tpi < hi)].mean()), 3) for _, lo, hi in CLS]

open(OUT + "metrics.txt", "w").write("\n".join(lines))
json.dump(R, open(OUT + "metrics.json", "w"), ensure_ascii=False, indent=1,
          default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
np.savez_compressed(OUT + "cells.npz", f_nat=f_nat.astype(np.float32), f_built=f_built.astype(np.float32), dom=dom.astype(np.int8), terr=terr.astype(np.int32), pa=pa.astype(np.int8))
print("\n".join(lines))
