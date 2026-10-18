"""Paths for the environmental-exposure analysis (Section 4.5 of the paper).

Everything is relative to the repository root unless overridden by the
environment variables TPI_RAW_DIR, TPI_OUT_DIR, TPI_RASTER, TPI_DISTRICTS.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = os.environ.get("TPI_RAW_DIR", os.path.join(ROOT, "data", "exposure_raw"))
OUT_DIR = os.environ.get("TPI_OUT_DIR", os.path.join(ROOT, "outputs", "exposure"))
# land-masked TPI surface used throughout the paper (88,102 valid cells, EPSG:32635, 250 m)
TPI_RASTER = os.environ.get("TPI_RASTER", os.path.join(ROOT, "outputs", "TPI_kara2.tif"))
# 128 district polygons of the published analysis (custom Albers CRS, see proj.py), path without extension
DISTRICT_SHP = os.environ.get("TPI_DISTRICTS", os.path.join(ROOT, "data", "districts", "TPI_ANALIZ_SONUC"))

WC_TILE = os.path.join(RAW_DIR, "ESA_WorldCover_10m_2021_v200_%s_Map.tif")   # % tile name
WDPA_DIR = os.path.join(RAW_DIR, "wdpa")                                     # unzipped WDPA_*_TUR_shp_0/1/2
OSM_PA_GEOJSON = os.path.join(RAW_DIR, "osm_protected_areas.geojson")
WSF_NPZ = os.path.join(RAW_DIR, "wsf_tracker_marmara.npz")

os.makedirs(OUT_DIR, exist_ok=True)
