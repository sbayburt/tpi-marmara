# Input data

Input data are **not redistributed** in this repository: each source carries its own
licence and terms of use. This file lists every input, where to obtain it, and the
format the scripts expect. Place all files in `data/raw/`.

---

## 1. Point-of-interest data — retrieved automatically

| File | Produced by | Source |
|------|-------------|--------|
| `osm_hotels.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_beaches.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_restaurants.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_marinas.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_camping.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_museums.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |
| `osm_ferry_piers.gpkg` | `01_download_osm.py` | OpenStreetMap via Overpass API |

No manual download is required — run `python scripts/01_download_osm.py`.
OpenStreetMap data are © OpenStreetMap contributors, available under the
[Open Database License (ODbL)](https://www.openstreetmap.org/copyright).

---

## 2. Blue Flag beaches — manual list, automatic geocoding

**File:** `blueflag_2026.csv`
**Source:** official national Blue Flag award list published by the Foundation for
Environmental Education in Türkiye (TÜRÇEV) for the relevant bathing season.

Expected format:

```csv
province,district,beach_name
Balıkesir,Erdek,Örnek Halk Plajı
Çanakkale,Bozcaada,Ayazma Plajı
...
```

Then run `python scripts/02_geocode_blueflag.py`, which geocodes each entry and writes
`osm_blueflag.gpkg`. Rows that cannot be matched confidently are listed in
`results/blueflag_review.csv` for manual checking.

---

## 3. Administrative boundaries

| File | Description | Required fields |
|------|-------------|-----------------|
| `districts.shp` | District polygons covering the study region | a name field (`name`, `adi`, `district` or similar) |
| `provinces.shp` | Provincial polygons, used to build the land mask | none |

Obtain from the national administrative boundary dataset published by the relevant
Turkish authority. The scripts detect the coordinate reference system automatically
and tolerate a missing or malformed `.prj`.

---

## 4. Comparison data (validation only)

These are used by `07_validate.py` and do **not** enter the index.

| File | Description | Source |
|------|-------------|--------|
| `ministry_accommodation_2025.xlsx` | District-level arrivals, overnight stays and occupancy | Ministry of Culture and Tourism, annual district accommodation statistics |
| `designated_tourism_areas.xlsx` | Register of Tourism Centres and Culture and Tourism Conservation and Development Regions | Ministry of Culture and Tourism |

The accommodation workbook is read with three header rows skipped and is expected to
have the standard published column order (province, district, then arrivals, overnight
stays, average stay and occupancy, each split foreign / domestic / total).

---

## 5. Cartographic base

The Esri World Terrain basemap is used for display only in the published figures.
It is not part of the computation and is not required to reproduce the results.

---

## Reproducibility note

Derived layers (the KDE surfaces, the TPI raster and the district-level results) are
written to `data/processed/` by the scripts. Large rasters are excluded from version
control by `.gitignore`; deposit them in an archival repository such as Zenodo and cite
the DOI in the paper's Data Availability statement.
