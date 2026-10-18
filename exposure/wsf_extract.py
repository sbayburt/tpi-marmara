import os
"""Download the WSF Tracker subset for the study area from the public GeoZarr store.

Requires Python >= 3.11 and:  pip install numpy xarray zarr obstore
Output: <RAW_DIR>/wsf_tracker_marmara.npz with arrays data (int8), x (lon), y (lat).
Pixel value = half-year of first detection (1 = 2016-07 ... 20 = 2026-01), 0 = no settlement.
Source: Marconcini & Metz-Marconcini (2026), https://doi.org/10.5281/zenodo.20424537 (CC BY 4.0).
"""
import numpy as np
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling modules
import xarray as xr
from config import WSF_NPZ

MIN_LON, MIN_LAT, MAX_LON, MAX_LAT = 25.85, 39.65, 30.65, 41.65


def open_store():
    try:
        import obstore
        from zarr.storage import ObjectStore
        store = obstore.store.S3Store(bucket="us-west-2.opendata.source.coop",
                                      prefix="mindearth/wsf/World_WSF_20160701-20260101.zarr",
                                      config={"aws_skip_signature": "true", "aws_region": "us-west-2"})
        return xr.open_zarr(ObjectStore(store, read_only=True), group="0", decode_coords="all")
    except Exception as e:                       # fallback: plain HTTPS (pip install fsspec aiohttp)
        print("obstore route failed (%s); trying HTTPS" % e)
        return xr.open_zarr("https://data.source.coop/mindearth/wsf/World_WSF_20160701-20260101.zarr", group="0", decode_coords="all")


ds = open_store()
da = ds["wsf_tracker"].sel(x=slice(MIN_LON, MAX_LON), y=slice(MAX_LAT, MIN_LAT))
out = np.zeros(da.shape, dtype=np.int8)
for i in range(0, da.shape[0], 1500):
    out[i:i + 1500] = da.isel(y=slice(i, i + 1500)).values.astype(np.int8)
    print("  %d / %d rows" % (min(i + 1500, da.shape[0]), da.shape[0]), flush=True)
np.savez_compressed(WSF_NPZ, data=out, x=da["x"].values, y=da["y"].values)
print("written:", WSF_NPZ, "settlement pixels:", int((out > 0).sum()))
