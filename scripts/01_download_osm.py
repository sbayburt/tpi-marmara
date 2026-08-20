# -*- coding: utf-8 -*-
"""
01_download_osm.py — retrieve seven POI layers from OpenStreetMap.

Queries the Overpass API layer by layer for the Marmara bounding box,
reduces polygon/way features to centroids, and writes one GeoPackage
per indicator to data/raw/.

Run:  python scripts/01_download_osm.py
"""

import time
import json
import requests
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

from config import (RAW, BBOX, OVERPASS_MIRRORS, OSM_QUERIES,
                    CRS_GEOGRAPHIC, banner)

HEADERS = {"User-Agent": "TPI-Marmara/1.0 (academic research; contact via GitHub)"}


def build_query(selector: str) -> str:
    """Assemble an Overpass QL query for nodes, ways and relations."""
    s, w, n, e = BBOX
    parts = []
    for sel in selector.split("|"):
        for kind in ("node", "way", "relation"):
            parts.append(f"{kind}{sel}({s},{w},{n},{e});")
    return f"[out:json][timeout:180];({''.join(parts)});out center;"


def fetch(query: str, retries: int = 3):
    """POST the query, rotating mirrors and backing off on rate limits."""
    for attempt in range(retries):
        url = OVERPASS_MIRRORS[attempt % len(OVERPASS_MIRRORS)]
        try:
            r = requests.post(url, data={"data": query}, headers=HEADERS, timeout=200)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 504):
                wait = 10 * (attempt + 1)
                print(f"    HTTP {r.status_code}; waiting {wait}s")
                time.sleep(wait)
                continue
            print(f"    HTTP {r.status_code} from {url}")
        except requests.RequestException as exc:
            print(f"    {type(exc).__name__}; retrying")
            time.sleep(5 * (attempt + 1))
    return None


def to_points(payload) -> gpd.GeoDataFrame:
    """Convert an Overpass response to a point GeoDataFrame (centroids for ways)."""
    rows = []
    for el in (payload or {}).get("elements", []):
        if el.get("type") == "node":
            lat, lon = el.get("lat"), el.get("lon")
        else:
            centre = el.get("center") or {}
            lat, lon = centre.get("lat"), centre.get("lon")
        if lat is None or lon is None:
            continue
        tags = el.get("tags", {}) or {}
        rows.append({
            "osm_id": el.get("id"),
            "osm_type": el.get("type"),
            "name": tags.get("name", ""),
            "geometry": Point(lon, lat),
        })
    if not rows:
        return gpd.GeoDataFrame(columns=["osm_id", "osm_type", "name", "geometry"],
                                geometry="geometry", crs=CRS_GEOGRAPHIC)
    return gpd.GeoDataFrame(pd.DataFrame(rows), geometry="geometry", crs=CRS_GEOGRAPHIC)


def main() -> None:
    banner("01_download_osm.py", "Retrieving POI layers from the Overpass API")
    summary = {}
    for layer, selector in OSM_QUERIES.items():
        out = RAW / f"osm_{layer}.gpkg"
        if out.exists():
            existing = gpd.read_file(out)
            print(f"  {layer:12s} already present ({len(existing)} features) — skipping")
            summary[layer] = len(existing)
            continue

        print(f"  {layer:12s} querying...")
        payload = fetch(build_query(selector))
        gdf = to_points(payload)
        gdf = gdf.drop_duplicates(subset=["osm_id", "osm_type"]).reset_index(drop=True)
        gdf.to_file(out, driver="GPKG")
        summary[layer] = len(gdf)
        print(f"  {layer:12s} {len(gdf):>6} features -> {out.name}")
        time.sleep(3)   # be polite to the public API

    (RAW / "osm_counts.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\n  Total OSM features: {sum(summary.values()):,}")
    print(f"  Counts written to {RAW / 'osm_counts.json'}")


if __name__ == "__main__":
    main()
