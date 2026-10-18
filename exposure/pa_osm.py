"""Protected areas from the OSM export (Overpass, 18 Oct 2026): rasterise onto the
TPI grid and compute exposure metric C.  WDPA (2 Ramsar polygons) is merged in."""
import json, re, sys, glob
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
from scipy import ndimage
from matplotlib.path import Path

from config import *
from geo import read_geotiff
from shp import read_shp_polygons, read_dbf
from proj import tm_forward

OUT = OUT_DIR + os.sep
tpi, meta = read_geotiff(TPI_RASTER)
valid = tpi > 0; high = valid & (tpi >= 0.6)
H, W, X0, Y0 = meta["height"], meta["width"], meta["x0"], meta["y0"]
D = np.load(OUT + "derived.npz"); sea = D["sea"]; tot = D["tot"]
DS = np.load(OUT + "districts.npz", allow_pickle=True); dist_id = DS["dist"]
region_land = (dist_id >= 0) & ~sea & (tot > 0)

g = json.load(open(OSM_PA_GEOJSON))
feats = []
for f in g["features"]:
    geom = f["geometry"]; p = f["properties"]
    if geom["type"] == "Polygon":
        polys = [geom["coordinates"]]
    elif geom["type"] == "MultiPolygon":
        polys = geom["coordinates"]
    else:
        continue
    name = p.get("name") or p.get("name:tr") or p.get("name:en") or "(unnamed)"
    title = p.get("protection_title") or ""
    desig = name + " " + title + " " + str(p.get("boundary")) + " " + str(p.get("leisure"))
    low = name.lower()
    if re.search(r"mesire|kent orman|araştırma orman|forest fakültesi|üniversitesi", low):
        cat = "recreation"            # recreation / urban-forest designations: not conservation
    elif re.search(r"özel çevre koruma", low) or re.search(r"deniz mill", low):
        cat = "marine"                # marine special environmental protection areas / marine NP
    elif re.search(r"[а-яА-Я]|Δ|ρ|Εθν", name):
        cat = "foreign"               # Bulgarian / Greek sites inside the bounding box
    else:
        cat = "terrestrial"
    rings = []
    for poly in polys:
        for k, ring in enumerate(poly):
            r = np.asarray(ring, dtype=float)
            x, y = tm_forward(r[:, 0], r[:, 1]); rings.append(np.column_stack([x, y]))
    feats.append(dict(name=name, cat=cat, title=title, boundary=p.get("boundary"), rings=rings))

# add the two WDPA Ramsar polygons (Manyas, Uluabat) as terrestrial PAs if not already present
for d_ in sorted(glob.glob(os.path.join(WDPA_DIR, "*", ""))):
    gg = read_shp_polygons(d_ + "WDPA_WDOECM_Oct2026_Public_TUR_shp-polygons.shp")
    names, rows = read_dbf(d_ + "WDPA_WDOECM_Oct2026_Public_TUR_shp-polygons.dbf")
    for feat, r in zip(gg, rows):
        if r["NAME"] in ("Manyas", "Lake Uluabat"):
            rings = []
            for q in feat:
                x, y = tm_forward(q[:, 0], q[:, 1]); rings.append(np.column_stack([x, y]))
            feats.append(dict(name=r["NAME"] + " (Ramsar, WDPA)", cat="terrestrial", title="Ramsar", boundary="wdpa", rings=rings))


def mask_of(rings):
    xs = np.concatenate([r[:, 0] for r in rings]); ys = np.concatenate([r[:, 1] for r in rings])
    c0 = max(0, int(np.floor((xs.min() - X0) / 250))); c1 = min(W - 1, int(np.ceil((xs.max() - X0) / 250)))
    r0 = max(0, int(np.floor((Y0 - ys.max()) / 250))); r1 = min(H - 1, int(np.ceil((Y0 - ys.min()) / 250)))
    m = np.zeros((H, W), dtype=bool)
    if c1 < c0 or r1 < r0:
        return m
    cc, rr = np.meshgrid(np.arange(c0, c1 + 1), np.arange(r0, r1 + 1))
    pts = np.column_stack([X0 + (cc.ravel() + 0.5) * 250, Y0 - (rr.ravel() + 0.5) * 250])
    inside = np.zeros(pts.shape[0], dtype=bool)
    for ring in rings:
        if len(ring) >= 3:
            inside ^= Path(ring).contains_points(pts)
    m[r0:r1 + 1, c0:c1 + 1] = inside.reshape(cc.shape)
    return m


terr_mask = np.zeros((H, W), dtype=bool); marine_mask = np.zeros((H, W), dtype=bool)
per = []
for f in feats:
    m = mask_of(f["rings"]); f["mask"] = m
    if f["cat"] == "terrestrial":
        terr_mask |= m
    elif f["cat"] == "marine":
        marine_mask |= m
    v = valid & m
    if f["cat"] in ("terrestrial", "marine") and v.sum() > 0:
        ring1 = (ndimage.distance_transform_edt(~m) * 250 <= 1000) & ~m
        per.append(dict(name=f["name"], cat=f["cat"], title=f["title"], cells=int(m.sum()), area_km2=round(m.sum() * 0.0625, 1),
                        land_pct=round(100 * (m & region_land).sum() / region_land.sum(), 2), valid_pct=round(100 * v.sum() / valid.sum(), 2),
                        mass_pct=round(100 * float(tpi[v].sum()) / float(tpi[valid].sum()), 2), high_pct=round(100 * (v & high).sum() / high.sum(), 2),
                        valid_in=int(v.sum()), mean_in=round(float(tpi[v].mean()), 3), max_in=round(float(tpi[v].max()), 3),
                        high_in=int((v & high).sum()), high_1km=int((valid & ring1 & high).sum()), valid_1km=int((valid & ring1).sum())))

cats = {c: sum(1 for f in feats if f["cat"] == c) for c in ("terrestrial", "marine", "recreation", "foreign")}
print("features by category:", cats)
print("terrestrial PA cells in region land: %d (%.1f%% of region land, %.0f km2); marine SEPA cells: %d" %
      ((terr_mask & region_land).sum(), 100 * (terr_mask & region_land).sum() / region_land.sum(), (terr_mask & region_land).sum() * 0.0625, marine_mask.sum()))
dist_t = ndimage.distance_transform_edt(~terr_mask) * 250
R = {}
for label, m in [("inside", terr_mask), ("within 1 km", (dist_t <= 1000)), ("within 3 km", (dist_t <= 3000))]:
    v = valid & m; hv = high & m
    R[label] = dict(land_pct=round(100 * (region_land & m).sum() / region_land.sum(), 1), valid=int(v.sum()), valid_pct=round(100 * v.sum() / valid.sum(), 1),
                    mean=round(float(tpi[v].mean()), 3), high=int(hv.sum()), high_pct=round(100 * hv.sum() / high.sum(), 1), mass_pct=round(100 * tpi[v].sum() / tpi[valid].sum(), 1))
    print("terrestrial PAs %-12s land %5.1f%%  valid %6d (%4.1f%%)  mean %.3f  high %4d (%4.1f%%)  mass %4.1f%%" %
          (label, R[label]["land_pct"], v.sum(), R[label]["valid_pct"], R[label]["mean"], hv.sum(), R[label]["high_pct"], R[label]["mass_pct"]))
v = valid & marine_mask; hv = high & marine_mask
R["overall_mean"] = round(float(tpi[valid].mean()), 3)
R["high_outside_any"] = int((high & ~(terr_mask | marine_mask)).sum())
R["n_terrestrial_polygons"] = sum(1 for f in feats if f["cat"] == "terrestrial")
R["marine"] = dict(valid=int(v.sum()), mean=round(float(tpi[v].mean()), 3), high=int(hv.sum()), high_pct=round(100 * hv.sum() / high.sum(), 1),
                   land_valid=int((v & ~sea).sum()), high_land=int((hv & ~sea).sum()))
print("marine SEPA/NP: valid cells inside %d (on land %d), mean %.3f, high %d (%.1f%%), high on land %d" %
      (v.sum(), (v & ~sea).sum(), R["marine"]["mean"], hv.sum(), R["marine"]["high_pct"], (hv & ~sea).sum()))
# which districts hold the marine-SEPA land cells
names = list(DS["names"])
ids, cnt = np.unique(dist_id[v & ~sea & (dist_id >= 0)], return_counts=True)
print("  by district:", sorted([(names[i], int(c)) for i, c in zip(ids, cnt)], key=lambda t: -t[1])[:10])
print("\nper-PA (terrestrial + marine) sorted by high cells inside, then mean:")
per.sort(key=lambda r: (-r["high_in"], -r["mean_in"]))
for r in per[:25]:
    print("  %-52s %-11s area %7.1f km2  valid_in %5d  mean %.3f  max %.3f  high_in %4d  high_1km %4d" %
          (r["name"][:52], r["cat"], r["area_km2"], r["valid_in"], r["mean_in"], r["max_in"], r["high_in"], r["high_1km"]))
json.dump(dict(summary=R, per=per, cats=cats), open(OUT + "pa_metrics.json", "w"), ensure_ascii=False, indent=1, default=lambda o: float(o))
np.savez_compressed(OUT + "pa_masks.npz", terr=terr_mask, marine=marine_mask)
import pickle
pickle.dump([dict(name=f["name"], cat=f["cat"], rings=f["rings"]) for f in feats], open(OUT + "pa_rings.pkl", "wb"))
# location of the maximum inside the two largest terrestrial designations
for nm in ("Uludağ", "Çanakkale Savaşları Gelibolu"):
    f = next(f for f in feats if f["name"].startswith(nm)); m = f["mask"] & valid
    r, c = np.unravel_index(np.argmax(np.where(m, tpi, -1)), tpi.shape)
    print("max in %s: TPI %.3f at district %s, %.1f km from coast" % (nm, tpi[r, c], DS["names"][dist_id[r, c]] if dist_id[r, c] >= 0 else "?", D["dist_m"][r, c] / 1000))
