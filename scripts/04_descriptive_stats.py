# -*- coding: utf-8 -*-
"""
04_descriptive_stats.py — Table 2, Table 3 and Figure 4.

Computes descriptive statistics of the TPI surface (including the Gini
coefficient of spatial concentration), the distribution across the five
pressure classes, and the histogram used as Figure 4.

Run:  python scripts/04_descriptive_stats.py
"""

import numpy as np
import pandas as pd
import rasterio
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import (TPI_MASKED, TPI_RASTER, RES, FIG, CELL_SIZE_M,
                    CLASS_BREAKS, CLASS_LABELS, FIG_DPI, banner)


def gini(x: np.ndarray) -> float:
    """Gini coefficient of a non-negative array."""
    s = np.sort(x)
    n = s.size
    return float((2 * np.sum(np.arange(1, n + 1) * s) / (n * s.sum())) - (n + 1) / n)


def main() -> None:
    banner("04_descriptive_stats.py", "Descriptive statistics and the class histogram")

    src_path = TPI_MASKED if TPI_MASKED.exists() else TPI_RASTER
    print(f"  Reading {src_path.name}")
    with rasterio.open(src_path) as src:
        arr = src.read(1).astype("float64")
        nodata = src.nodata

    valid = np.isfinite(arr) & (arr > 0)
    if nodata is not None:
        valid &= (arr != nodata)
    v = arr[valid]

    cell_km2 = (CELL_SIZE_M / 1000) ** 2
    q1, med, q3 = np.percentile(v, [25, 50, 75])

    table2 = {
        "Valid cells (n)": f"{v.size:,}",
        "Area covered (km2)": f"{v.size * cell_km2:,.0f}",
        "Minimum": f"{v.min():.3f}",
        "Maximum": f"{v.max():.3f}",
        "Mean": f"{v.mean():.3f}",
        "Median": f"{med:.3f}",
        "First quartile (Q1)": f"{q1:.3f}",
        "Third quartile (Q3)": f"{q3:.3f}",
        "Standard deviation": f"{v.std():.3f}",
        "Skewness": f"{stats.skew(v):.2f}",
        "Kurtosis": f"{stats.kurtosis(v):.2f}",
        "Gini coefficient": f"{gini(v):.3f}",
    }
    print("\n  TABLE 2 — Descriptive statistics")
    for k, val in table2.items():
        print(f"    {k:24s} {val}")
    pd.DataFrame(list(table2.items()), columns=["Statistic", "Value"]).to_csv(
        RES / "table2_descriptive.csv", index=False, encoding="utf-8-sig")

    rows = []
    for i in range(5):
        lo, hi = CLASS_BREAKS[i], CLASS_BREAKS[i + 1]
        m = (v > lo) & (v <= hi)
        rows.append({
            "Class": CLASS_LABELS[i],
            "Range": f"{lo:.2f}-{hi:.2f}",
            "Cells": int(m.sum()),
            "Area_km2": round(m.sum() * cell_km2, 1),
            "Share_pct": round(100 * m.sum() / v.size, 1),
        })
    t3 = pd.DataFrame(rows)
    print("\n  TABLE 3 — Class distribution")
    print(t3.to_string(index=False))
    t3.to_csv(RES / "table3_classes.csv", index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(v, bins=100, color="#4A7C9E", edgecolor="none")
    for b in CLASS_BREAKS[1:-1]:
        ax.axvline(b, color="#C00000", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("TPI value")
    ax.set_ylabel("Number of cells (log scale)")
    ax.set_title("Distribution of TPI values with class boundaries")
    fig.tight_layout()
    out = FIG / "figure04_histogram.png"
    fig.savefig(out, dpi=FIG_DPI)
    print(f"\n  Figure 4 -> {out.name}")


if __name__ == "__main__":
    main()
