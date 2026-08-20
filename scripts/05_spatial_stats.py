# -*- coding: utf-8 -*-
"""
05_spatial_stats.py — global Moran's I, LISA and Getis-Ord Gi* (Section 4.2, Figure 5).

Aggregates the TPI surface to district level by zonal statistics over non-zero
cells, then computes global and local spatial statistics. Insular districts are
connected by a nearest-neighbour link so that island units are not left without
neighbours.

Run:  python scripts/05_spatial_stats.py
"""

import numpy as np
import geopandas as gpd
from rasterstats import zonal_stats
from libpysal.weights import Queen, KNN, w_union
from esda.moran import Moran, Moran_Local
from esda.getisord import G_Local
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import (TPI_MASKED, TPI_RASTER, DISTRICTS_SHP, DISTRICT_RESULT,
                    RES, FIG, CRS_METRIC, CELL_SIZE_M, MIN_DISTRICT_AREA_KM2,
                    MIN_TPI_CELLS, N_PERMUTATIONS, HOTSPOT_COLOURS, FIG_DPI, banner)

CRS_CANDIDATES = {
    "UTM_35N": "EPSG:32635",
    "Albers_TR": ("+proj=aea +lat_1=36.5 +lat_2=41.0 +lat_0=0 +lon_0=35 "
                  "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),
    "WGS84": "EPSG:4326",
    "WebMercator": "EPSG:3857",
}


def read_districts() -> gpd.GeoDataFrame:
    """Read district polygons, detecting the CRS if the .prj is unusable."""
    if not DISTRICTS_SHP.exists():
        raise SystemExit(f"Missing input: {DISTRICTS_SHP}\nSee data/README.md.")

    prj = DISTRICTS_SHP.with_suffix(".prj")
    backup = DISTRICTS_SHP.with_suffix(".prj.bak")
    moved = False
    try:
        if prj.exists():
            prj.rename(backup)
            moved = True
        raw = gpd.read_file(DISTRICTS_SHP, encoding="utf-8")
    finally:
        if moved and backup.exists():
            backup.rename(prj)

    for label, crs in CRS_CANDIDATES.items():
        try:
            test = raw.set_crs(crs, allow_override=True).to_crs("EPSG:4326")
            minx, miny, maxx, maxy = test.total_bounds
            if 24 < minx < 46 and 35 < miny < 43 and 24 < maxx < 46 and 35 < maxy < 43:
                print(f"  CRS detected: {label}")
                return raw.set_crs(crs, allow_override=True).to_crs(CRS_METRIC)
        except Exception:
            continue
    raise SystemExit("Could not determine the CRS of the district layer.")


def classify(z: float, p: float) -> str:
    if p >= 0.10:
        return "Not sig."
    if z > 0:
        return "Hot 99%" if p < 0.01 else ("Hot 95%" if p < 0.05 else "Hot 90%")
    return "Cold 99%" if p < 0.01 else ("Cold 95%" if p < 0.05 else "Cold 90%")


def main() -> None:
    banner("05_spatial_stats.py", "Moran's I, LISA and Getis-Ord Gi* at district level")

    src = TPI_MASKED if TPI_MASKED.exists() else TPI_RASTER
    g = read_districts()
    name_col = next((c for c in g.columns
                     if c.lower() in ("name", "adi", "adı", "district", "ilce", "ilçe")),
                    g.columns[0])
    print(f"  Districts read: {len(g)} (name field: '{name_col}')")

    # dissolve multipart districts sharing a name, drop slivers
    g = g.dissolve(by=name_col, as_index=False)
    g["area_km2"] = g.geometry.area / 1e6
    g = g[g["area_km2"] >= MIN_DISTRICT_AREA_KM2].reset_index(drop=True)

    # zonal statistics over non-zero cells only (nodata=0)
    print(f"  Zonal statistics against {src.name}")
    zs = zonal_stats(g, str(src), stats=["mean", "max", "count"], nodata=0)
    g["tpi_mean"] = [round(z["mean"], 6) if z["mean"] else 0.0 for z in zs]
    g["tpi_max"] = [round(z["max"], 6) if z["max"] else 0.0 for z in zs]
    g["tpi_cells"] = [z["count"] or 0 for z in zs]

    before = len(g)
    g = g[g["tpi_cells"] >= MIN_TPI_CELLS].reset_index(drop=True)
    print(f"  Districts with measurable TPI: {before} -> {len(g)}")

    # spatial weights
    w = Queen.from_dataframe(g, use_index=False)
    if len(w.islands) > 0:
        print(f"  {len(w.islands)} insular district(s); adding nearest-neighbour links")
        w = w_union(w, KNN.from_dataframe(g, k=1))
    w.transform = "r"
    y = g["tpi_mean"].values

    mi = Moran(y, w, permutations=N_PERMUTATIONS)
    print(f"\n  GLOBAL MORAN'S I")
    print(f"    I = {mi.I:.3f}   z = {mi.z_sim:.2f}   p = {mi.p_sim:.5f} "
          f"({N_PERMUTATIONS} permutations)")

    lm = Moran_Local(y, w, permutations=N_PERMUTATIONS)
    quadrant = {1: "High-High", 2: "Low-High", 3: "Low-Low", 4: "High-Low"}
    g["LISA"] = [quadrant[q] if p < 0.05 else "Not sig."
                 for q, p in zip(lm.q, lm.p_sim)]

    gi = G_Local(y, w, star=True, permutations=N_PERMUTATIONS)
    g["GiZ"] = np.round(gi.Zs, 3)
    g["GiP"] = np.round(gi.p_sim, 4)
    g["GiClass"] = [classify(z, p) for z, p in zip(gi.Zs, gi.p_sim)]

    hot = int(g["GiClass"].str.startswith("Hot").sum())
    cold = int(g["GiClass"].str.startswith("Cold").sum())
    print(f"\n  GETIS-ORD Gi*")
    print(g["GiClass"].value_counts().to_string().replace("\n", "\n    "))
    print(f"    Summary: {hot} hot spots, {cold} cold spots, "
          f"{len(g)-hot-cold} not significant")

    top = g.nlargest(20, "GiZ")[[name_col, "tpi_mean", "GiZ", "GiClass"]]
    print(f"\n  Twenty highest Gi* z-scores:\n{top.to_string(index=False)}")

    hh = g[g["LISA"] == "High-High"][name_col].astype(str).tolist()
    hl = g[g["LISA"] == "High-Low"][name_col].astype(str).tolist()

    g.to_file(DISTRICT_RESULT, driver="GPKG")
    with open(RES / "spatial_statistics.txt", "w", encoding="utf-8") as f:
        f.write(f"SPATIAL STATISTICS ({len(g)} districts)\n{'='*60}\n\n")
        f.write(f"Global Moran's I = {mi.I:.3f}, z = {mi.z_sim:.2f}, "
                f"p = {mi.p_sim:.5f} ({N_PERMUTATIONS} permutations)\n")
        f.write("Queen contiguity + nearest-neighbour link for islands, row-standardised\n\n")
        f.write(f"Getis-Ord Gi*\n{g['GiClass'].value_counts().to_string()}\n")
        f.write(f"\nSummary: {hot} hot spots, {cold} cold spots\n\n")
        f.write(f"Twenty highest Gi* z-scores\n{top.to_string(index=False)}\n\n")
        f.write(f"High-High cluster ({len(hh)}): {', '.join(hh)}\n")
        f.write(f"High-Low outliers ({len(hl)}): {', '.join(hl) if hl else 'none'}\n")

    fig, ax = plt.subplots(figsize=(11, 7))
    for label, colour in HOTSPOT_COLOURS.items():
        sub = g[g["GiClass"] == label]
        if len(sub):
            sub.plot(ax=ax, color=colour, edgecolor="#999999", lw=0.3, label=label)
    ax.legend(fontsize=8, loc="lower left", title="Gi* classification")
    ax.set_axis_off()
    ax.set_title("Getis-Ord Gi* hot-spot analysis of the Tourism Pressure Index")
    fig.tight_layout()
    fig.savefig(FIG / "figure05_hotspots.png", dpi=FIG_DPI)

    print(f"\n  Districts with statistics -> {DISTRICT_RESULT.name}")
    print(f"  Figure 5 (quick-look) -> figure05_hotspots.png")
    print("  For the published figure, symbolise GiClass in a GIS using the")
    print("  colours defined in config.HOTSPOT_COLOURS.")


if __name__ == "__main__":
    main()
