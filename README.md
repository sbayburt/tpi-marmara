# Tourism Pressure Index (TPI) — Marmara Region

Reproducible workflow for the paper:

> **A Reproducible Tourism Pressure Index (TPI) from Open Geospatial Data: Methodology and Application to the Marmara Region, Türkiye**
> Serdar Bayburt, Gürcan Büyüksalih, Cem Gazioğlu (2026)

This repository contains the complete, script-based production chain used to build the
Tourism Pressure Index (TPI) for the Marmara Region of Türkiye at 250 m resolution, together
with the spatial-statistical analysis and the three-source validation reported in the paper.
Every processing decision — indicator selection, kernel bandwidth, outlier treatment,
weighting, land masking — is expressed in code and can be re-executed or varied by a third party.

---

## What this repository produces

Running the scripts in order reproduces:

1. Eight tourism-supply density surfaces and the combined **TPI raster** (250 m, 0–1).
2. **Descriptive statistics** of the TPI surface (mean, median, Gini, class distribution) and a histogram.
3. **Global Moran's I**, **Local Moran's I (LISA)** and **Getis-Ord Gi\*** hot-spot analysis at district level.
4. A **three-scenario weighting sensitivity** analysis.
5. **Validation** against three independent open sources (official accommodation statistics,
   review density, officially designated tourism areas), including scatter plots and a ROC curve.

---

## Repository structure

```
tpi-marmara/
├── README.md                     ← this file
├── LICENSE                       ← MIT (code) — see note on data licensing below
├── CITATION.cff                  ← how to cite this repository
├── requirements.txt              ← Python dependencies
├── environment.yml               ← optional conda environment
├── .gitignore
│
├── scripts/                      ← all processing and analysis code
│   ├── 01_download_osm.py        ← retrieve 7 POI layers from OpenStreetMap (Overpass API)
│   ├── 02_geocode_blueflag.py    ← geocode the official Blue Flag beach list (Places API)
│   ├── 03_build_tpi.py           ← KDE, outlier treatment, weighting → TPI raster
│   ├── 04_descriptive_stats.py   ← Table 2, Table 3, Figure 4 (histogram)
│   ├── 05_spatial_stats.py       ← Moran's I, LISA, Getis-Ord Gi*, Figure 5
│   ├── 06_sensitivity.py         ← three weighting scenarios, Table 4, Figure 6
│   ├── 07_validate.py            ← three-source validation, Table 5, Figures 7–9
│   └── config.py                 ← shared paths, CRS, weights, constants
│
├── data/
│   ├── raw/                      ← inputs you download (not tracked; see data/README.md)
│   └── processed/                ← derived layers written by the scripts
│
├── figures/                      ← figures written by the scripts (PNG, 300 dpi)
├── results/                      ← tables and text summaries written by the scripts
└── docs/
    └── data_sources.md           ← exact provenance of every input (Section 3.1)
```

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/<your-username>/tpi-marmara.git
cd tpi-marmara

# 2. Create environment (choose one)
pip install -r requirements.txt
#   or
conda env create -f environment.yml && conda activate tpi

# 3. Obtain input data (see data/README.md for the exact sources)
#    Place downloaded files in data/raw/

# 4. Set your Google Places API key (only needed for scripts 02 and 07)
export GOOGLE_PLACES_API_KEY="your-key"      # Windows: set GOOGLE_PLACES_API_KEY=...

# 5. Run the pipeline in order
python scripts/01_download_osm.py
python scripts/02_geocode_blueflag.py
python scripts/03_build_tpi.py
python scripts/04_descriptive_stats.py
python scripts/05_spatial_stats.py
python scripts/06_sensitivity.py
python scripts/07_validate.py
```

Each script prints what it reads and writes, and can be run independently once its
inputs exist. Paths, the coordinate reference system and the indicator weights are
centralised in `scripts/config.py` — edit there, not inside individual scripts.

---

## Method summary

| Stage | Operation | Key parameters |
|-------|-----------|----------------|
| Acquisition | OSM Overpass API + official Blue Flag list | 8 layers, n = 10,890 |
| Pre-process | centroid extraction, reprojection, land mask | EPSG:32635 |
| Density | Gaussian kernel density estimation | 250 m grid, σ = 750 m |
| Normalise | 95th-percentile clip → log1p → rescale 0–1 | per indicator |
| Aggregate | weighted linear sum | Σ wᵢ = 1 |

Indicator weights: hotels 0.25, beaches 0.20, restaurants/cafés 0.15, marinas 0.10,
camping 0.10, museums 0.10, ferry piers 0.05, Blue Flag 0.05.

---

## Data availability and licensing

- **Code** in this repository is released under the MIT License (see `LICENSE`).
- **Input data** are openly available but carry their own licences and are **not
  redistributed** here. `data/README.md` lists each source and how to obtain it.
  OpenStreetMap data are © OpenStreetMap contributors, available under the Open
  Database License (ODbL).
- **Derived indicator layers** produced by the scripts may be deposited in an archival
  repository (e.g. Zenodo) and linked from the paper's Data Availability statement.

---

## Citation

If you use this code or the derived layers, please cite the paper (see `CITATION.cff`).

---

## Contact

Serdar Bayburt — Denge Altyapı Proje Müşavirlik ve İnşaat A.Ş., İstanbul, Türkiye.
Issues and questions: please use the GitHub Issues tab.
