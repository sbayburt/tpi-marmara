# -*- coding: utf-8 -*-
"""
07_validate.py — validation against three independent open sources
(Section 4.4, Table 5, Figures 7-9).

  1. Official district-level accommodation statistics (arrivals, overnight stays, occupancy)
  2. User review density from a places API, assigned to districts by spatial containment
  3. Officially designated Tourism Centres and Conservation/Development Regions (ROC, AUC)

Run:  export GOOGLE_PLACES_API_KEY=...   &&   python scripts/07_validate.py
"""

import re
import json
import time
import unicodedata

import numpy as np
import pandas as pd
import geopandas as gpd
import requests
from scipy import stats
from sklearn.metrics import roc_curve, auc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import (DISTRICT_RESULT, MINISTRY_XLSX, DESIGNATED_XLSX, RES, FIG,
                    GOOGLE_CACHE, PLACES_API_KEY, PLACES_SEARCH_URL,
                    PLACES_QUERIES, PLACES_RADIUS_M, MARMARA_PROVINCES,
                    FIG_DPI, banner)

FIELD_MASK = "places.id,places.userRatingCount,places.location"
REPORT = []


def log(line: str = "") -> None:
    print(line)
    REPORT.append(line)


def norm(s) -> str:
    """Normalise Turkish district names for matching."""
    if pd.isna(s):
        return ""
    s = str(s).strip().lower()
    for a, b in {"ı": "i", "ş": "s", "ğ": "g", "ü": "u",
                 "ö": "o", "ç": "c", "â": "a", "î": "i", "û": "u"}.items():
        s = s.replace(a, b)
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]", "", s)


def correlate(x: np.ndarray, y: np.ndarray) -> dict:
    """Spearman on raw values, Pearson on log-transformed, plus RMSE and MAE."""
    ylog = np.log1p(y)
    sp_r, sp_p = stats.spearmanr(x, y)
    pe_r, pe_p = stats.pearsonr(x, ylog)
    xn = (x - x.min()) / (x.max() - x.min())
    yn = (ylog - ylog.min()) / (ylog.max() - ylog.min())
    return {"n": len(x), "spearman": sp_r, "spearman_p": sp_p,
            "pearson": pe_r, "pearson_p": pe_p,
            "rmse": float(np.sqrt(np.mean((xn - yn) ** 2))),
            "mae": float(np.mean(np.abs(xn - yn)))}


# ------------------------------------------------------------- source 1
def validate_ministry(g: gpd.GeoDataFrame, name_col: str) -> pd.DataFrame:
    log("\n" + "=" * 60)
    log("SOURCE 1 — Official accommodation statistics")
    log("=" * 60)
    if not MINISTRY_XLSX.exists():
        log(f"  Missing {MINISTRY_XLSX.name}; skipped.")
        return pd.DataFrame()

    df = pd.read_excel(MINISTRY_XLSX, sheet_name=0, header=None, skiprows=3)
    df.columns = ["province", "district",
                  "arr_foreign", "arr_domestic", "arr_total",
                  "night_foreign", "night_domestic", "night_total",
                  "stay_foreign", "stay_domestic", "stay_total",
                  "occ_foreign", "occ_domestic", "occ_total"]
    df["province"] = df["province"].ffill()
    df = df[df["district"].notna()]
    df = df[df["district"].astype(str).str.strip().str.lower() != "toplam"]
    df = df[df["province"].isin(MARMARA_PROVINCES)].copy()
    df["key"] = df["district"].apply(norm)
    df.loc[df["key"] == "merkez", "key"] = (
        df.loc[df["key"] == "merkez", "province"].apply(norm) + "merkez")
    for col in ("arr_total", "night_total", "occ_total"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    merged = g[[name_col, "key", "tpi_mean"]].merge(
        df[["province", "district", "key", "arr_total", "night_total", "occ_total"]],
        on="key", how="inner")
    log(f"  Matched districts: {len(merged)}")

    rows = []
    for col, label in [("arr_total", "Arrivals"),
                       ("night_total", "Overnight stays"),
                       ("occ_total", "Occupancy rate")]:
        sub = merged[(merged["tpi_mean"] > 0) & (merged[col] > 0)]
        if len(sub) < 10:
            continue
        r = correlate(sub["tpi_mean"].values, sub[col].values)
        rows.append({"Source": "Ministry statistics", "Measure": label, **r})
        log(f"\n  {label} (n = {r['n']})")
        log(f"    Spearman rho = {r['spearman']:.3f} (p = {r['spearman_p']:.2e})")
        log(f"    Pearson r    = {r['pearson']:.3f} (p = {r['pearson_p']:.2e})")
        log(f"    RMSE = {r['rmse']:.3f}   MAE = {r['mae']:.3f}")
        if col == "night_total":
            plot_scatter(sub, name_col, "night_total", r,
                         "log(1 + overnight stays)",
                         "TPI vs. official overnight stays",
                         FIG / "figure07_overnight_stays.png", "#C0504D")
    return pd.DataFrame(rows)


# ------------------------------------------------------------- source 2
def fetch_places(district: str, lat: float, lon: float, cache: dict) -> tuple:
    total_reviews, count = 0, 0
    for query in PLACES_QUERIES:
        key = f"{district}|{query}"
        if key in cache:
            data = cache[key]
        else:
            headers = {"Content-Type": "application/json",
                       "X-Goog-Api-Key": PLACES_API_KEY,
                       "X-Goog-FieldMask": FIELD_MASK}
            payload = {"textQuery": f"{query} {district}",
                       "languageCode": "tr", "regionCode": "TR",
                       "maxResultCount": 20,
                       "locationBias": {"circle": {
                           "center": {"latitude": lat, "longitude": lon},
                           "radius": PLACES_RADIUS_M}}}
            data = {}
            for attempt in range(3):
                try:
                    resp = requests.post(PLACES_SEARCH_URL, headers=headers,
                                         json=payload, timeout=30)
                    if resp.status_code == 200:
                        data = resp.json()
                        break
                    if resp.status_code == 429:
                        time.sleep(2 ** attempt)
                        continue
                    break
                except requests.RequestException:
                    time.sleep(2 ** attempt)
            cache[key] = data
            GOOGLE_CACHE.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")
            time.sleep(0.12)
        for place in (data.get("places") or []):
            total_reviews += place.get("userRatingCount", 0) or 0
            count += 1
    return total_reviews, count


def validate_reviews(g: gpd.GeoDataFrame, name_col: str) -> pd.DataFrame:
    log("\n" + "=" * 60)
    log("SOURCE 2 — Review density")
    log("=" * 60)

    cache = json.loads(GOOGLE_CACHE.read_text(encoding="utf-8")) if GOOGLE_CACHE.exists() else {}

    if not cache and not PLACES_API_KEY:
        log("  No cache and no API key; skipped.")
        return pd.DataFrame()

    if not cache:
        log("  Building cache (this consumes API quota)")
        g84 = g.to_crs("EPSG:4326")
        for i, row in g.iterrows():
            centroid = g84.geometry.iloc[i].centroid
            fetch_places(str(row[name_col]), centroid.y, centroid.x, cache)

    # assign each establishment to the district that contains it
    points = []
    for _, data in cache.items():
        for place in (data.get("places") or []):
            loc = place.get("location") or {}
            if place.get("id") and loc.get("latitude"):
                points.append({"place_id": place["id"],
                               "lat": loc["latitude"], "lon": loc["longitude"],
                               "reviews": place.get("userRatingCount", 0) or 0})
    pdf = pd.DataFrame(points).drop_duplicates(subset="place_id")
    gp = gpd.GeoDataFrame(pdf, geometry=gpd.points_from_xy(pdf.lon, pdf.lat),
                          crs="EPSG:4326")
    g84 = g.to_crs("EPSG:4326")
    joined = gpd.sjoin(gp, g84[[name_col, "geometry"]], how="inner", predicate="within")
    log(f"  Unique establishments: {len(pdf):,}; assigned to a district: {len(joined):,} "
        f"({100*len(joined)/max(len(pdf),1):.1f}%)")

    agg = joined.groupby(name_col).agg(reviews=("reviews", "sum"),
                                       places=("place_id", "count")).reset_index()
    merged = g[[name_col, "tpi_mean", "area_km2"]].merge(agg, on=name_col, how="left")
    merged[["reviews", "places"]] = merged[["reviews", "places"]].fillna(0)
    merged["reviews_km2"] = merged["reviews"] / merged["area_km2"].clip(lower=1)

    rows = []
    sub = merged[(merged["tpi_mean"] > 0) & (merged["reviews"] > 0)]
    for col, label in [("reviews", "Total reviews"), ("reviews_km2", "Reviews per km2")]:
        r = correlate(sub["tpi_mean"].values, sub[col].values)
        rows.append({"Source": "Review data", "Measure": label, **r})
        log(f"\n  {label} (n = {r['n']})")
        log(f"    Spearman rho = {r['spearman']:.3f} (p = {r['spearman_p']:.2e})")
        log(f"    Pearson r    = {r['pearson']:.3f} (p = {r['pearson_p']:.2e})")
        if col == "reviews_km2":
            plot_scatter(sub, name_col, "reviews_km2", r,
                         "log(1 + reviews per km2)",
                         "TPI vs. review density",
                         FIG / "figure08_review_density.png", "#4A7C9E")
    return pd.DataFrame(rows)


# ------------------------------------------------------------- source 3
def validate_designated(g: gpd.GeoDataFrame, name_col: str) -> pd.DataFrame:
    log("\n" + "=" * 60)
    log("SOURCE 3 — Officially designated tourism areas")
    log("=" * 60)
    if not DESIGNATED_XLSX.exists():
        log(f"  Missing {DESIGNATED_XLSX.name}; skipped.")
        return pd.DataFrame()

    df = pd.read_excel(DESIGNATED_XLSX, sheet_name=0, header=None, skiprows=3).iloc[:, :8]
    df.columns = ["no", "code", "name", "area_ha", "type", "province", "district", "theme"]
    df = df[df["district"].notna()]

    records = []
    for _, row in df.iterrows():
        for part in str(row["district"]).split(","):
            if part.strip():
                records.append({"district": part.strip(), "theme": row["theme"]})
    designated = pd.DataFrame(records)
    designated["key"] = designated["district"].apply(norm)

    g = g.copy()
    g["designated"] = g["key"].isin(set(designated["key"])).astype(int)
    positives = int(g["designated"].sum())
    log(f"  Districts containing a designated area: {positives} / {len(g)}")

    if positives < 3:
        log("  Too few positives for a ROC analysis; skipped.")
        return pd.DataFrame()

    fpr, tpr, _ = roc_curve(g["designated"].values, g["tpi_mean"].values)
    area = auc(fpr, tpr)
    u, p = stats.mannwhitneyu(g[g["designated"] == 1]["tpi_mean"],
                              g[g["designated"] == 0]["tpi_mean"],
                              alternative="greater")
    log(f"\n  ROC AUC = {area:.3f}")
    log(f"  Mann-Whitney U = {u:.0f}, p = {p:.4f} (one-tailed)")
    log(f"  Mean TPI: designated {g[g['designated']==1]['tpi_mean'].mean():.3f} "
        f"vs other {g[g['designated']==0]['tpi_mean'].mean():.3f}")

    fig, ax = plt.subplots(figsize=(6, 5.5))
    ax.plot(fpr, tpr, lw=2.2, color="#C0504D", label=f"TPI (AUC = {area:.3f})")
    ax.plot([0, 1], [0, 1], lw=1, ls="--", color="#999999", label="Random (0.50)")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC: TPI predicting designated tourism areas")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
    fig.tight_layout()
    fig.savefig(FIG / "figure09_roc.png", dpi=FIG_DPI)
    log(f"  Figure 9 -> figure09_roc.png")

    return pd.DataFrame([{"Source": "Designated areas", "Measure": f"ROC AUC = {area:.3f}",
                          "n": len(g), "spearman": np.nan, "spearman_p": np.nan,
                          "pearson": np.nan, "pearson_p": p, "rmse": np.nan, "mae": np.nan}])


def plot_scatter(sub, name_col, col, r, ylabel, title, path, colour) -> None:
    fig, ax = plt.subplots(figsize=(7, 5.5))
    x = sub["tpi_mean"].values
    y = np.log1p(sub[col].values)
    ax.scatter(x, y, s=32, c=colour, edgecolors="#333333", linewidths=0.4, alpha=0.85)
    for _, row in sub.nlargest(10, col).iterrows():
        ax.annotate(str(row[name_col]), (row["tpi_mean"], np.log1p(row[col])),
                    fontsize=7, alpha=0.75, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("District mean TPI")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nSpearman rho = {r['spearman']:.3f}, n = {r['n']}")
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)


def main() -> None:
    banner("07_validate.py", "Validation against three independent open sources")
    if not DISTRICT_RESULT.exists():
        raise SystemExit(f"Missing {DISTRICT_RESULT.name}. Run 05_spatial_stats.py first.")

    g = gpd.read_file(DISTRICT_RESULT)
    name_col = next((c for c in g.columns
                     if c.lower() in ("name", "adi", "adı", "district", "ilce", "ilçe")),
                    g.columns[0])
    g["key"] = g[name_col].apply(norm)
    log(f"  Districts: {len(g)}")

    tables = [validate_ministry(g, name_col),
              validate_reviews(g, name_col),
              validate_designated(g, name_col)]
    table5 = pd.concat([t for t in tables if len(t)], ignore_index=True)

    if len(table5):
        table5.to_csv(RES / "table5_validation.csv", index=False, encoding="utf-8-sig")
        log(f"\n  Table 5 -> table5_validation.csv")

    (RES / "validation_report.txt").write_text("\n".join(REPORT), encoding="utf-8")
    log(f"  Full report -> validation_report.txt")


if __name__ == "__main__":
    main()
