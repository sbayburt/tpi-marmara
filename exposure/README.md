# Environmental exposure of the TPI surface (paper Section 4.5)

Scripts that reproduce Tables 6–10 and Figure 10 of the paper: the land-cover
composition of the pressure classes, the coastal distance bands, the district
"pressure on natural cover" ranking, the protected-area comparison and the
2016–2026 settlement change. They consume the published TPI surface produced
by the main pipeline (`scripts/`) and four open context layers that do **not**
enter the index.

Everything is plain NumPy/SciPy — no GDAL, rasterio or GeoPandas — so the
chain runs on a bare Python installation. Map projections (UTM 35N, and the
custom Albers CRS of the district layer) are implemented in `proj.py`;
GeoTIFF, GeoPackage and Shapefile readers in `geo.py` / `shp.py`.

## Inputs

| Layer | Source | How to obtain | Put in `data/exposure_raw/` as |
|---|---|---|---|
| TPI surface | main pipeline | `outputs/TPI_kara2.tif` (land-masked, 88,102 valid cells) | — (see `config.py`) |
| Districts | main pipeline | 128 district polygons with results, custom Albers CRS | — (see `config.py`) |
| Land cover | ESA WorldCover 2021 v200, 10 m (Zanaga et al., 2022; CC BY 4.0) | tiles N39E024, N39E027, N39E030 from `https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_<TILE>_Map.tif` | `ESA_WorldCover_10m_2021_v200_<TILE>_Map.tif` |
| Protected areas (international) | WDPA, October 2026 (UNEP-WCMC & IUCN) | country download for Türkiye (shapefile) from protectedplanet.net, unzip the three inner archives | `wdpa/WDPA_*_TUR_shp_0/`, `_1/`, `_2/` |
| Protected areas (national) | OpenStreetMap, Overpass API, October 2026 (ODbL) | run the query below in overpass-turbo.eu and export as GeoJSON | `osm_protected_areas.geojson` |
| Settlement 2016–2026 | WSF Tracker (Marconcini & Metz-Marconcini, 2026; CC BY 4.0) | `python wsf_extract.py` (needs `xarray zarr obstore`) | `wsf_tracker_marmara.npz` |

Overpass query used for the protected areas (the eleven provinces of the region):

```
[out:json][timeout:180];
area["ISO3166-2"~"^TR-(10|11|16|17|22|34|39|41|54|59|77)$"]["admin_level"="4"]->.marmara;
(
  way["boundary"~"^(protected_area|national_park)$"](area.marmara);
  relation["boundary"~"^(protected_area|national_park)$"](area.marmara);
  way["leisure"="nature_reserve"](area.marmara);
  relation["leisure"="nature_reserve"](area.marmara);
  way["protection_title"](area.marmara);
  relation["protection_title"](area.marmara);
);
out body; >; out skel qt;
```

`pa_osm.py` excludes recreation areas and urban forests (name matching
"mesire", "kent ormanı", research forests) and sites outside Türkiye, and
treats the two Special Environmental Protection Areas and the North Aegean
Marine National Park separately from the terrestrial designations.

## Running

```
pip install numpy scipy tifffile matplotlib          # analysis
pip install xarray zarr obstore                      # only for wsf_extract.py
python wsf_extract.py        # ~1 GB download, writes wsf_tracker_marmara.npz
python districts.py          # rasterise the district layer; validates against the published district means
python wc_aggregate.py       # WorldCover 10 m -> 250 m class counts (~6 min, 7 GB RAM)
python derived.py            # class fractions, sea mask, distance to coast
python pa_osm.py             # protected areas: masks, per-area statistics (Table 9)
python metrics.py            # Tables 6–7 and Figure 10 inputs
python table8.py             # Table 8
python wsf_aggregate.py      # settlement area per cell and half-year
python wsf_metrics.py        # Table 10
python fig10.py              # Figure 10
python fig2.py               # Figure 2 (workflow diagram)
```

Paths are set in `config.py` and can be overridden with the environment
variables `TPI_RAW_DIR`, `TPI_OUT_DIR`, `TPI_RASTER`, `TPI_DISTRICTS`.

## Outputs (in `outputs/exposure/`)

`metrics.txt` / `metrics.json` (Tables 6–7, Section 4.5 figures), `table8.json`
(Table 8), `pa_metrics.json` (Table 9), `wsf_metrics.txt` / `.json` (Table 10),
`out_dist_to_coast_m.tif` and `out_builtup_fraction.tif` (250 m, EPSG:32635),
`Fig10_exposure.png/.pdf`, `Fig2_workflow.png/.pdf`.

## Notes

* `districts.py` reproduces the published district mean TPI and cell counts
  for all 128 units exactly; this validates the Albers → UTM transformation
  and the rasterisation against the original GIS workflow.
* 982 valid cells (1.1 %) lie entirely on water at 10 m resolution because the
  250 m land mask is coarser than WorldCover; they are kept in the first
  coastal band and excluded from the land-cover composition.
* The WSF Tracker pixel value is the half-year of first detection
  (1 = 2016-07 … 20 = 2026-01); the record is monotonic, so change measures
  expansion of the built footprint, not redevelopment within it.
* OpenStreetMap coverage of protected areas is incomplete and boundary
  accuracy varies; the protected-area statistics are lower bounds.
