"""Figure 2 redrawn: the methodological workflow as implemented (no log1p),
88,102 valid cells, correct comparison sources, and the exposure analysis
added to the analytical band.  Same aspect ratio as the original (1.26:1)."""
import os
import matplotlib
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
from config import OUT_DIR
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42})
W, H = 12.6, 11.2
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
INK, INK2 = "#1a1a1a", "#4d4d4d"
STYLE = {"data": ("#e8f1fb", "#9dbbe0"), "calc": ("#e6f4ee", "#93c7ae"), "agg": ("#ece9f8", "#b1a7e3"),
         "out": ("#fbeae6", "#e0a697"), "ana": ("#f1efe9", "#c9c4b4")}


def box(x, y, w, h, title, lines, kind, title_size=14, size=11.5, pad=0.02):
    fc, ec = STYLE[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.18", fc=fc, ec=ec, lw=1.2))
    ax.text(x + w / 2, y + h - 0.38, title, ha="center", va="center", fontsize=title_size, fontweight="bold", color=INK)
    for k, ln in enumerate(lines):
        ax.text(x + 0.25, y + h - 0.85 - k * 0.37, ln, ha="left", va="center", fontsize=size, color=INK2)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16, lw=1.3, color="#8a8a8a", shrinkA=2, shrinkB=2))


# ---- five processing stages
bw, bh, y0 = 2.75, 2.55, 8.15
xs = [0.35, 3.45, 6.55, 9.65]
box(xs[0], y0, bw, bh, "1 · Acquisition", ["OSM Overpass API", "7 POI layers", "Blue Flag list (TÜRÇEV)", "Places geocoding", "n = 10,890 features"], "data")
box(xs[1], y0, bw, bh, "2 · Pre-process", ["Centroid extraction", "Reproject to EPSG:32635", "De-duplication", "Land mask"], "data")
box(xs[2], y0, bw, bh, "3 · Density (KDE)", ["250 m grid", "Gaussian kernel", "σ = 750 m (3 cells)", "One surface per indicator"], "calc")
box(xs[3], y0, bw, bh, "4 · Normalise", ["Clip at 95th percentile", "Divide by 95th percentile", "Floor: values < 0.05 → 0", "Common 0–1 scale"], "calc")
for a, b in zip(xs[:-1], xs[1:]):
    arrow(a + bw, y0 + bh / 2, b, y0 + bh / 2)
# stage 5 below stage 4
y5 = 5.15
box(xs[3], y5, bw, 2.3, "5 · Aggregate", ["Weighted sum, Σ wᵢ = 1", "Min–max rescaling to 0–1", "Five pressure classes"], "agg")
arrow(xs[3] + bw / 2, y0, xs[3] + bw / 2, y5 + 2.3)
# output
ox, ow, oh = 4.55, 3.9, 1.5
oy = y5 + (2.3 - oh) / 2
ax.add_patch(FancyBboxPatch((ox, oy), ow, oh, boxstyle="round,pad=0.02,rounding_size=0.18", fc=STYLE["out"][0], ec=STYLE["out"][1], lw=1.2))
ax.text(ox + ow / 2, oy + oh - 0.5, "TPI surface (0–1)", ha="center", va="center", fontsize=15, fontweight="bold", color=INK)
ax.text(ox + ow / 2, oy + 0.45, "250 m · 88,102 valid cells", ha="center", va="center", fontsize=11.5, color=INK2)
arrow(xs[3], y5 + 2.3 / 2, ox + ow, oy + oh / 2)
# analysis band
ay, ah, ax0, aw = 0.35, 2.75, 0.35, 12.05
# context layers that enter the evaluation but not the index
cx, cy, cw, ch = 0.35, y5 - 0.1, 3.7, 2.5
ax.add_patch(FancyBboxPatch((cx, cy), cw, ch, boxstyle="round,pad=0.02,rounding_size=0.18", fc="#fbf3e2", ec="#d9b36a", lw=1.2))
ax.text(cx + cw / 2, cy + ch - 0.38, "Context layers (not in the index)", ha="center", va="center", fontsize=12.5, fontweight="bold", color=INK)
for k, ln in enumerate(["Ministry statistics · review density", "Designated tourism areas (ROC)", "ESA WorldCover 2021 land cover (10 m)",
                        "Protected areas: WDPA + OpenStreetMap", "WSF Tracker settlement 2016–2026 (10 m)"]):
    ax.text(cx + 0.22, cy + ch - 0.82 - k * 0.36, ln, ha="left", va="center", fontsize=10.3, color=INK2)
arrow(cx + cw / 2, cy, cx + cw / 2, ay + ah)
ax.add_patch(FancyBboxPatch((ax0, ay), aw, ah, boxstyle="round,pad=0.02,rounding_size=0.18", fc=STYLE["ana"][0], ec=STYLE["ana"][1], lw=1.2))
ax.text(ax0 + aw / 2, ay + ah - 0.38, "Analysis and evaluation (Sections 4.1–4.5)", ha="center", va="center", fontsize=14, fontweight="bold", color=INK)
arrow(ox + ow / 2, oy, ox + ow / 2, ay + ah)
subs = [("Descriptive statistics", ["Gini coefficient", "class distribution"]),
        ("Spatial structure", ["Moran's I", "Getis-Ord Gi* hot spots"]),
        ("Weighting sensitivity", ["three scenarios", "class change"]),
        ("Independent indicators", ["ministry statistics", "review density · ROC"]),
        ("Environmental exposure", ["land cover · coast · protected", "areas · settlement change"])]
sw, sh, gap = 2.22, 1.65, 0.2
sx = ax0 + (aw - (5 * sw + 4 * gap)) / 2
for k, (t, ls) in enumerate(subs):
    x = sx + k * (sw + gap); y = ay + 0.3
    ax.add_patch(FancyBboxPatch((x, y), sw, sh, boxstyle="round,pad=0.02,rounding_size=0.14", fc="#f8f7f3", ec="#c9c4b4", lw=1.0))
    ax.text(x + sw / 2, y + sh - 0.38, t, ha="center", va="center", fontsize=11.5, fontweight="bold", color=INK)
    for j, ln in enumerate(ls):
        ax.text(x + sw / 2, y + sh - 0.8 - j * 0.36, ln, ha="center", va="center", fontsize=10.5, color=INK2)

fig.savefig(os.path.join(OUT_DIR, "Fig2_workflow.png"), dpi=600, facecolor="white")
fig.savefig(os.path.join(OUT_DIR, "Fig2_workflow.pdf"), facecolor="white")
fig.savefig(os.path.join(OUT_DIR, "Fig2_workflow_preview.png"), dpi=110, facecolor="white")
print("saved")
