# -*- coding: utf-8 -*-
"""
config.py — shared configuration for the TPI workflow.

Every path, parameter and weight used by the pipeline is defined here.
Edit this file rather than individual scripts.
"""

from pathlib import Path
import os

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
FIG = ROOT / "figures"
RES = ROOT / "results"

for _d in (RAW, PROC, FIG, RES):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------- input files
# Place these in data/raw/ — see data/README.md for how to obtain them.
DISTRICTS_SHP = RAW / "districts.shp"          # district polygons (with a name field)
LAND_MASK_SHP = RAW / "provinces.shp"          # provincial polygons, dissolved for masking
BLUEFLAG_LIST = RAW / "blueflag_2026.csv"      # official Blue Flag list: district, beach_name
MINISTRY_XLSX = RAW / "ministry_accommodation_2025.xlsx"
DESIGNATED_XLSX = RAW / "designated_tourism_areas.xlsx"

# --------------------------------------------------------- output files
TPI_RASTER = PROC / "tpi_marmara_250m.tif"
TPI_MASKED = PROC / "tpi_marmara_250m_masked.tif"
DISTRICT_RESULT = PROC / "tpi_districts.gpkg"
GOOGLE_CACHE = PROC / "places_cache.json"

# ------------------------------------------------------------------ CRS
CRS_GEOGRAPHIC = "EPSG:4326"      # WGS 84, as delivered by the APIs
CRS_METRIC = "EPSG:32635"         # UTM Zone 35N — all distance operations

# ------------------------------------------------------ grid and kernel
CELL_SIZE_M = 250                 # raster cell size in metres
SIGMA_CELLS = 3                   # Gaussian kernel sd in cells (3 x 250 m = 750 m)
SIGMA_M = SIGMA_CELLS * CELL_SIZE_M

# --------------------------------------------------- outlier treatment
CLIP_PERCENTILE = 95              # clip densities above this percentile
# Transformation is log1p, applied after clipping; see 03_build_tpi.py

# ------------------------------------------------------------- classes
CLASS_BREAKS = [0.0, 0.20, 0.40, 0.60, 0.80, 1.0]
CLASS_LABELS = ["Low", "Moderate", "Moderate-high", "High", "Very high"]

# ------------------------------------------------------------- weights
# S2 = the scheme reported in Table 1 of the paper (the published index)
INDICATORS = [
    "hotels", "beaches", "restaurants", "marinas",
    "camping", "museums", "ferry_piers", "blueflag",
]

WEIGHTS_S2_TOURISM = {
    "hotels": 0.25, "beaches": 0.20, "restaurants": 0.15, "marinas": 0.10,
    "camping": 0.10, "museums": 0.10, "ferry_piers": 0.05, "blueflag": 0.05,
}
WEIGHTS_S1_EQUAL = {k: 0.125 for k in INDICATORS}
WEIGHTS_S3_INFRA = {
    "hotels": 0.15, "beaches": 0.15, "restaurants": 0.10, "marinas": 0.20,
    "camping": 0.15, "museums": 0.05, "ferry_piers": 0.15, "blueflag": 0.05,
}
SCENARIOS = {
    "S1_equal": WEIGHTS_S1_EQUAL,
    "S2_tourism": WEIGHTS_S2_TOURISM,
    "S3_infrastructure": WEIGHTS_S3_INFRA,
}

# ------------------------------------------------- Overpass API queries
# Bounding box for the Marmara Region: south, west, north, east
BBOX = (39.5, 25.6, 42.2, 31.5)

OVERPASS_URL = "https://overpass-api.de/api/interpreter"
OVERPASS_MIRRORS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
]

OSM_QUERIES = {
    "hotels":      '["tourism"~"^(hotel|motel|guest_house|hostel|apartment)$"]',
    "beaches":     '["natural"="beach"]',
    "restaurants": '["amenity"~"^(restaurant|cafe)$"]',
    "marinas":     '["leisure"="marina"]',
    "camping":     '["tourism"~"^(camp_site|caravan_site)$"]',
    "museums":     '["tourism"~"^(museum|artwork)$"]|["historic"]',
    "ferry_piers": '["amenity"="ferry_terminal"]',
}

# --------------------------------------------------- Marmara provinces
MARMARA_PROVINCES = [
    "Balıkesir", "Bilecik", "Bursa", "Çanakkale", "Edirne", "İstanbul",
    "Kırklareli", "Kocaeli", "Sakarya", "Tekirdağ", "Yalova",
]

# ---------------------------------------------------- district filters
MIN_DISTRICT_AREA_KM2 = 1.0       # discard sliver polygons below this
MIN_TPI_CELLS = 20                # district must contain at least this many non-zero cells

# ------------------------------------------------ spatial statistics
N_PERMUTATIONS = 9999
SIGNIFICANCE_LEVELS = [0.10, 0.05, 0.01]   # 90%, 95%, 99%

# ------------------------------------------------------- Places API
PLACES_API_KEY = os.environ.get("GOOGLE_PLACES_API_KEY", "")
PLACES_SEARCH_URL = "https://places.googleapis.com/v1/places:searchText"
PLACES_QUERIES = ["otel", "turistik yer", "plaj", "müze"]
PLACES_RADIUS_M = 8000

# ------------------------------------------------------ figure style
FIG_DPI = 300
HOTSPOT_COLOURS = {
    "Hot 99%": "#B2182B", "Hot 95%": "#EF8A62", "Hot 90%": "#FDDBC7",
    "Not sig.": "#F7F7F7",
    "Cold 90%": "#D1E5F0", "Cold 95%": "#67A9CF", "Cold 99%": "#2166AC",
}


def banner(script_name: str, description: str) -> None:
    """Print a consistent header so pipeline runs are easy to follow."""
    line = "=" * 70
    print(f"\n{line}\n{script_name}\n{description}\n{line}")
