"""Figure 10 — environmental exposure of the TPI surface.
(a) high and very-high pressure cells by dominant terrestrial land cover
(b) land-cover composition of each pressure class
(c) coastal distance bands: share of regional land, of pressure mass and of high-pressure cells
174 mm wide, 600 dpi PNG (combination art) + PDF.
"""
import sys, json
import os
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from config import *
from geo import read_geotiff
from shp import read_shp_polygons
from proj import albers_to_utm, tm_forward

OUT = OUT_DIR + os.sep
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7.5, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "pdf.fonttype": 42})
COL = {"Built-up": "#6250d6", "Cropland": "#e87ba4", "Tree cover": "#008300",
       "Grassland / shrubland": "#eda100", "Bare / sparse": "#9a9a9a", "Herbaceous wetland": "#1baf7a"}
GREY = ["#c9c9c9", "#8c8c8c", "#3a3a3a"]
INK, INK2 = "#1a1a1a", "#555555"

tpi, meta = read_geotiff(TPI_RASTER)
valid = tpi > 0; high = valid & (tpi >= 0.6)
D = np.load(OUT + "derived.npz"); sea = D["sea"]
C = np.load(OUT + "cells.npz"); dom, terr = C["dom"], C["terr"]
R = json.load(open(OUT + "metrics.json"))
X0, Y0, H, W = meta["x0"], meta["y0"], meta["height"], meta["width"]
ext = [X0, X0 + W * 250, Y0 - H * 250, Y0]

fig = plt.figure(figsize=(174 / 25.4, 176 / 25.4))
gs_top = fig.add_gridspec(1, 1, left=0.012, right=0.988, top=0.992, bottom=0.495)
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1], wspace=0.36, left=0.175, right=0.975, top=0.415, bottom=0.145)

# ------------------------------------------------------------------ (a) map
ax = fig.add_subplot(gs_top[0, 0])
bg = np.where(sea, 0, 1).astype(float)
ax.imshow(bg, extent=ext, cmap=matplotlib.colors.ListedColormap(["#dfe7ee", "#f8f8f6"]), vmin=0, vmax=1,
          interpolation="nearest", zorder=0)
feats = read_shp_polygons(DISTRICT_SHP + ".shp")
for rings in feats:
    for r in rings:
        x, y = albers_to_utm(r[:, 0], r[:, 1])
        ax.plot(x, y, color="#b9b9b9", lw=0.25, zorder=1)
import pickle
PAS = pickle.load(open(OUT + "pa_rings.pkl", "rb"))
from matplotlib.patches import Polygon as MplPolygon
for f in PAS:
    if f["cat"] == "terrestrial":
        for k, ring in enumerate(f["rings"]):
            if len(ring) >= 3:
                ax.add_patch(MplPolygon(ring, closed=True, fc="#cfe8c9" if k == 0 else "#f8f8f6", ec="#6aa84f", lw=0.3, zorder=1.5))
    elif f["cat"] == "marine" and ("Marmara" in f["name"] or "Saros" in f["name"]):
        for ring in f["rings"][:1]:
            ax.plot(ring[:, 0], ring[:, 1], color="#1f6f78", lw=0.55, ls=(0, (3, 2)), zorder=1.8)
rr, cc = np.where(high)
xs = X0 + (cc + 0.5) * 250; ys = Y0 - (rr + 0.5) * 250
dcls = dom[rr, cc]; tr = terr[rr, cc]
groups = [("Built-up", (dcls == 0) & (tr > 0)), ("Tree cover", (dcls == 2) & (tr > 0)),
          ("Grassland / shrubland", (dcls == 3) & (tr > 0)),
          ("Other (cropland, bare, wetland, water)", ~(((dcls == 0) | (dcls == 2) | (dcls == 3)) & (tr > 0)))]
for name, m in [groups[0], groups[3], groups[2], groups[1]]:          # draw the small classes last
    col = COL.get(name, "#9a9a9a")
    ax.scatter(xs[m], ys[m], s=3.2, marker="s", c=col, linewidths=0, zorder=3 if name != "Built-up" else 2)
ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.set_aspect("equal")
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_edgecolor("#999999"); s.set_linewidth(0.5)

places = [("İstanbul", 28.90, 41.07, "left", 0, 9000), ("İzmit", 29.92, 40.77, "left", 4000, -9000),
          ("Bursa", 29.06, 40.19, "left", 4000, -4000), ("Uludağ", 29.22, 40.07, "left", 4000, -8000),
          ("Şile", 29.61, 41.18, "left", 3000, 4000), ("Kandıra", 30.15, 41.07, "left", 3000, 4000),
          ("Erdek", 27.79, 40.40, "right", -3500, 2500), ("Bandırma", 27.97, 40.35, "left", 3500, -6500),
          ("Bozcaada", 26.07, 39.83, "left", 5000, 5000), ("Gelibolu", 26.67, 40.41, "right", -3000, 1000),
          ("Şarköy", 27.11, 40.61, "left", -9000, 7000), ("Adalar", 29.09, 40.87, "left", 6000, -6000),
          ("Çanakkale", 26.41, 40.15, "left", 4000, -5000), ("Marmara Is.", 27.60, 40.60, "center", 0, -9000),
          ("Tekirdağ", 27.51, 40.98, "center", 0, 6500), ("Yalova", 29.28, 40.66, "left", 3000, -7000)]
for name, lon, lat, ha, dx, dy in places:
    x, y = tm_forward(lon, lat)
    ax.plot(x, y, marker="o", ms=1.6, mfc="white", mec=INK, mew=0.4, zorder=5, ls="none")
    ax.annotate(name, (x, y), xytext=(x + dx, y + dy), fontsize=6.2, color=INK, ha=ha, va="center", zorder=6)
# scale bar 50 km
sx, sy = ext[0] + 12000, ext[2] + 14000
ax.plot([sx, sx + 50000], [sy, sy], color=INK, lw=1.2, solid_capstyle="butt", zorder=6)
ax.text(sx + 25000, sy + 3500, "50 km", ha="center", va="bottom", fontsize=6.2, color=INK, zorder=6)
ax.text(ext[0] + 5000, ext[3] - 6000, "(a)", fontsize=8.5, fontweight="bold", va="top", color=INK, zorder=6)
handles = [Line2D([], [], marker="s", ls="none", ms=4.2, mfc=COL["Built-up"], mec="none", label="Built-up"),
           Line2D([], [], marker="s", ls="none", ms=4.2, mfc=COL["Tree cover"], mec="none", label="Tree cover"),
           Line2D([], [], marker="s", ls="none", ms=4.2, mfc=COL["Grassland / shrubland"], mec="none", label="Grassland / shrubland"),
           Line2D([], [], marker="s", ls="none", ms=4.2, mfc="#9a9a9a", mec="none", label="Other"),
           Patch(fc="#cfe8c9", ec="#6aa84f", lw=0.4, label="Terrestrial protected area"),
           Line2D([], [], color="#1f6f78", lw=0.9, ls=(0, (3, 2)), label="Special Environmental Protection Area")]
leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.004, 0.935), frameon=True, fontsize=6.4, framealpha=0.95, edgecolor="#cccccc",
                title="High and very-high pressure cells (TPI ≥ 0.60)\nby dominant terrestrial land cover", title_fontsize=6.4,
                borderpad=0.6, handletextpad=0.5, labelspacing=0.35)
leg.get_frame().set_linewidth(0.4); leg._legend_box.align = "left"

# ------------------------------------------------------------------ (b) composition
ax = fig.add_subplot(gs[0, 0])
rows = R["A"]
labels = [r["class"].replace(" (", "\n(") for r in rows]
order = ["Built-up", "Cropland", "Tree cover", "Grassland / shrubland", "Bare / sparse (beach, dune, rock)", "Herbaceous wetland"]
short = {"Bare / sparse (beach, dune, rock)": "Bare / sparse"}
ypos = np.arange(len(rows))[::-1]
left = np.zeros(len(rows))
for g in order:
    vals = np.array([r[g] for r in rows])
    key = short.get(g, g)
    ax.barh(ypos, vals, left=left, height=0.62, color=COL[key], edgecolor="white", linewidth=0.8, label=key)
    for yv, lv, v in zip(ypos, left, vals):
        if v >= 8:
            ax.text(lv + v / 2, yv, "%.0f" % v, ha="center", va="center", fontsize=6.2,
                    color="white" if key in ("Built-up", "Tree cover") else INK)
    left += vals
ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=6.4)
ax.set_xlim(0, 100); ax.set_xlabel("Share of terrestrial surface (%)", fontsize=6.8)
ax.tick_params(axis="x", labelsize=6.4, length=2); ax.tick_params(axis="y", length=0)
for s in ["top", "right", "left"]:
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color("#999999")
ax.legend(loc="upper center", bbox_to_anchor=(0.38, -0.21), ncol=3, frameon=False, fontsize=6.2, handlelength=1.0,
          columnspacing=0.9, handletextpad=0.4)
ax.text(-0.42, 1.03, "(b)", transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=INK)

# ------------------------------------------------------------------ (c) coastal bands
ax = fig.add_subplot(gs[0, 1])
B = R["B"]
bands = [b["band"] for b in B]
series = [("Regional land area", [b["pct_land"] for b in B]), ("Pressure mass (Σ TPI)", [b["pct_mass"] for b in B]),
          ("High-pressure cells", [b["pct_high"] for b in B])]
x = np.arange(len(bands)); w = 0.26
for i, (name, vals) in enumerate(series):
    bars = ax.bar(x + (i - 1) * w, vals, width=w - 0.03, color=GREY[i], label=name, linewidth=0)
    for xb, v in zip(x + (i - 1) * w, vals):
        ax.text(xb, v + 1.2, "%.1f" % v if v < 10 else "%.0f" % v, ha="center", va="bottom", fontsize=5.9, color=INK)
ax.set_xticks(x); ax.set_xticklabels(bands, fontsize=6.4)
ax.set_xlabel("Distance from the coast", fontsize=6.8)
ax.set_ylabel("Share of regional total (%)", fontsize=6.8)
ax.set_ylim(0, 82); ax.tick_params(axis="y", labelsize=6.4, length=2); ax.tick_params(axis="x", length=0)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color("#999999"); ax.spines["bottom"].set_color("#999999")
ax.yaxis.grid(True, color="#e6e6e6", lw=0.4); ax.set_axisbelow(True)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=6.2, handlelength=1.0,
          columnspacing=1.0, handletextpad=0.4)
ax.text(-0.2, 1.03, "(c)", transform=ax.transAxes, fontsize=8.5, fontweight="bold", color=INK)

fig.savefig(os.path.join(OUT_DIR, "Fig10_exposure.png"), dpi=600)
fig.savefig(os.path.join(OUT_DIR, "Fig10_exposure.pdf"))
fig.savefig(os.path.join(OUT_DIR, "Fig10_exposure_preview.png"), dpi=150)
print("saved")
