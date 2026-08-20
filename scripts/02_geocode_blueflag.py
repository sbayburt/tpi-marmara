# -*- coding: utf-8 -*-
"""
02_geocode_blueflag.py — geocode the official Blue Flag beach list.

The Blue Flag award is only sporadically tagged in OpenStreetMap, so this
layer is compiled from the official national list and geocoded through the
Places API with district-centred location bias and bounding-box validation.
Each result receives a confidence label; low-confidence rows are written to a
separate file for manual review.

Input:  data/raw/blueflag_2026.csv  with columns: province, district, beach_name
Output: data/raw/osm_blueflag.gpkg  (name kept consistent with the OSM layers)
        results/blueflag_review.csv (rows needing manual verification)

Run:  export GOOGLE_PLACES_API_KEY=...   &&   python scripts/02_geocode_blueflag.py
"""

import json
import time
import requests
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from config import (RAW, RES, BBOX, BLUEFLAG_LIST, PLACES_API_KEY,
                    PLACES_SEARCH_URL, GOOGLE_CACHE, CRS_GEOGRAPHIC, banner)

FIELD_MASK = "places.id,places.displayName,places.formattedAddress,places.location"


def load_cache() -> dict:
    if GOOGLE_CACHE.exists():
        return json.loads(GOOGLE_CACHE.read_text(encoding="utf-8"))
    return {}


def save_cache(cache: dict) -> None:
    GOOGLE_CACHE.write_text(json.dumps(cache, ensure_ascii=False), encoding="utf-8")


def geocode(query: str, cache: dict) -> dict:
    """Text search with caching; returns the raw API payload."""
    if query in cache:
        return cache[query]
    headers = {"Content-Type": "application/json",
               "X-Goog-Api-Key": PLACES_API_KEY,
               "X-Goog-FieldMask": FIELD_MASK}
    payload = {"textQuery": query, "languageCode": "tr",
               "regionCode": "TR", "maxResultCount": 5}
    for attempt in range(3):
        try:
            r = requests.post(PLACES_SEARCH_URL, headers=headers, json=payload, timeout=30)
            if r.status_code == 200:
                data = r.json()
                cache[query] = data
                save_cache(cache)
                return data
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            print(f"    HTTP {r.status_code} for '{query}'")
            return {}
        except requests.RequestException:
            time.sleep(2 ** attempt)
    return {}


def in_bbox(lat: float, lon: float) -> bool:
    s, w, n, e = BBOX
    return (s <= lat <= n) and (w <= lon <= e)


def main() -> None:
    banner("02_geocode_blueflag.py", "Geocoding the official Blue Flag beach list")
    if not PLACES_API_KEY:
        raise SystemExit("Set GOOGLE_PLACES_API_KEY before running this script.")
    if not BLUEFLAG_LIST.exists():
        raise SystemExit(f"Missing input: {BLUEFLAG_LIST}\nSee data/README.md.")

    beaches = pd.read_csv(BLUEFLAG_LIST)
    required = {"province", "district", "beach_name"}
    if not required.issubset(beaches.columns):
        raise SystemExit(f"{BLUEFLAG_LIST.name} must contain columns: {sorted(required)}")

    cache = load_cache()
    rows, review = [], []

    for i, rec in beaches.iterrows():
        query = f"{rec['beach_name']} plaj {rec['district']} {rec['province']}"
        data = geocode(query, cache)
        places = data.get("places") or []

        if not places:
            review.append({**rec.to_dict(), "reason": "no result", "confidence": "NONE"})
            print(f"  [{i+1:>3}] {rec['beach_name'][:35]:35s} NO RESULT")
            continue

        top = places[0]
        loc = top.get("location", {})
        lat, lon = loc.get("latitude"), loc.get("longitude")

        if lat is None or not in_bbox(lat, lon):
            review.append({**rec.to_dict(), "reason": "outside bbox", "confidence": "LOW"})
            print(f"  [{i+1:>3}] {rec['beach_name'][:35]:35s} OUTSIDE BBOX")
            continue

        address = top.get("formattedAddress", "")
        confidence = "HIGH" if rec["district"].lower() in address.lower() else "MEDIUM"
        if confidence == "MEDIUM":
            review.append({**rec.to_dict(), "reason": "district not in address",
                           "confidence": "MEDIUM"})

        rows.append({
            "beach_name": rec["beach_name"],
            "district": rec["district"],
            "province": rec["province"],
            "confidence": confidence,
            "address": address,
            "geometry": Point(lon, lat),
        })
        print(f"  [{i+1:>3}] {rec['beach_name'][:35]:35s} {confidence}")
        time.sleep(0.12)

    gdf = gpd.GeoDataFrame(pd.DataFrame(rows), geometry="geometry", crs=CRS_GEOGRAPHIC)
    out = RAW / "osm_blueflag.gpkg"
    gdf.to_file(out, driver="GPKG")
    print(f"\n  Geocoded: {len(gdf)} / {len(beaches)} -> {out.name}")

    if review:
        rev = RES / "blueflag_review.csv"
        pd.DataFrame(review).to_csv(rev, index=False, encoding="utf-8-sig")
        print(f"  {len(review)} rows need manual review -> {rev.name}")


if __name__ == "__main__":
    main()
