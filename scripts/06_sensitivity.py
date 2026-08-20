# -*- coding: utf-8 -*-
"""
06_sensitivity.py — three-scenario weighting sensitivity (Section 4.3, Table 4, Figure 6).

Recomputes the index under equal, tourism-oriented and infrastructure-oriented
weighting schemes and reports pairwise agreement.

Run:  python scripts/06_sensitivity.py
"""

import numpy as np
import pandas as pd
import rasterio
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import (PROC, RES, FIG, INDICATORS, SCENARIOS, CLASS_BREAKS,
                    FIG_DPI, banner)


def main() -> None:
    banner("06_sensitivity.py", "Weighting sensitivity across three scenarios")

    surfaces, profile = {}, None
    for name in INDICATORS:
        path = PROC / f"kde_{name}.tif"
        if not path.exists():
            raise SystemExit(f"Missing {path.name}. Run 03_build_tpi.py first.")
        with rasterio.open(path) as src:
            surfaces[name] = src.read(1).astype("float64")
            if profile is None:
                profile = src.profile.copy()

    def combine(weights: dict) -> np.ndarray:
        out = np.zeros_like(next(iter(surfaces.values())))
        for name, w in weights.items():
            out += w * surfaces[name]
        return out / out.max() if out.max() > 0 else out

    scenario_surfaces = {}
    for label, weights in SCENARIOS.items():
        surface = combine(weights)
        scenario_surfaces[label] = surface
        out = PROC / f"tpi_{label}.tif"
        profile.update(dtype="float32", count=1, compress="lzw")
        with rasterio.open(out, "w", **profile) as dst:
            dst.write(surface.astype("float32"), 1)
        print(f"  {label:20s} -> {out.name}")

    valid = np.ones_like(next(iter(scenario_surfaces.values())), dtype=bool)
    for surface in scenario_surfaces.values():
        valid &= np.isfinite(surface) & (surface > 0)
    print(f"\n  Cells compared: {valid.sum():,}")

    def classes(a): return np.digitize(a, CLASS_BREAKS[1:-1])

    rows, labels = [], list(SCENARIOS)
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a = scenario_surfaces[labels[i]][valid]
            b = scenario_surfaces[labels[j]][valid]
            rows.append({
                "Scenario pair": f"{labels[i]} vs {labels[j]}",
                "Pearson r": round(stats.pearsonr(a, b)[0], 3),
                "Spearman rho": round(stats.spearmanr(a, b)[0], 3),
                "Class change %": round(100 * np.mean(classes(a) != classes(b)), 1),
            })

    table = pd.DataFrame(rows)
    print(f"\n  TABLE 4 — Scenario agreement\n{table.to_string(index=False)}")
    table.to_csv(RES / "table4_sensitivity.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, label in zip(axes, labels):
        img = ax.imshow(np.where(valid, scenario_surfaces[label], np.nan),
                        cmap="RdYlGn_r", vmin=0, vmax=1)
        ax.set_title(label.replace("_", " "))
        ax.set_axis_off()
    fig.colorbar(img, ax=axes, shrink=0.7, label="TPI")
    fig.savefig(FIG / "figure06_scenarios.png", dpi=FIG_DPI, bbox_inches="tight")
    print(f"\n  Figure 6 -> figure06_scenarios.png")

    if (table["Pearson r"] > 0.95).all():
        print("\n  All pairings exceed r = 0.95: the index is robust to moderate reweighting.")


if __name__ == "__main__":
    main()
