# -*- coding: utf-8 -*-
"""
03_build_tpi.py — build the Tourism Pressure Index raster.

Stages (Sections 3.2-3.5 of the paper):
  1. read the eight indicator point layers, reproject to EPSG:32635
  2. rasterise each to a 250 m grid
  3. Gaussian kernel density estimation, sigma = 3 cells (750 m)
  4. 95th-percentile clipping -> log1p -> rescale to 0-1
  5. weighted linear sum, land mask, write GeoTIFF

Run:  python scripts/03_build_tpi.py
"""

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from rasterio.features import geometry_mask
from scipy import ndimage

from config import (RAW, PROC, TPI_RASTER, TPI_MASKED, LAND_MASK_SHP,
                    CRS_METRIC, CELL_SIZE_M, SIGMA_CELLS, CLIP_PERCENTILE,
                    INDICATORS, WEIGHTS_S2_TOURISM, banner)


def read_layer(name: str) -> gpd.GeoDataFrame:
    path = RAW / f"osm_{name}.gpkg"
    if not path.exists():
        raise SystemExit(f"Missing layer: {path}\nRun 01_download_osm.py / 02_geocode_blueflag.py first.")
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf.to_crs(CRS_METRIC)


def build_grid(layers: dict, pad_m: float = 5000.0):
    """Common grid covering all layers, padded so kernels are not truncated."""
    bounds = np.array([g.total_bounds for g in layers.values() if len(g)])
    minx, miny = bounds[:, 0].min() - pad_m, bounds[:, 1].min() - pad_m
    maxx, maxy = bounds[:, 2].max() + pad_m, bounds[:, 3].max() + pad_m
    width = int(np.ceil((maxx - minx) / CELL_SIZE_M))
    height = int(np.ceil((maxy - miny) / CELL_SIZE_M))
    transform = from_origin(minx, maxy, CELL_SIZE_M, CELL_SIZE_M)
    return transform, width, height


def rasterise_counts(gdf, transform, width, height) -> np.ndarray:
    """Count points falling in each cell."""
    counts = np.zeros((height, width), dtype="float64")
    if len(gdf) == 0:
        return counts
    inv = ~transform
    for geom in gdf.geometry:
        if geom is None or geom.is_empty:
            continue
        col, row = inv * (geom.x, geom.y)
        r, c = int(row), int(col)
        if 0 <= r < height and 0 <= c < width:
            counts[r, c] += 1.0
    return counts


def kde(counts: np.ndarray) -> np.ndarray:
    """Gaussian kernel density estimation by convolution."""
    return ndimage.gaussian_filter(counts, sigma=SIGMA_CELLS, mode="constant", cval=0.0)


def normalise(surface: np.ndarray) -> np.ndarray:
    """95th-percentile clip -> log1p -> rescale to [0, 1]."""
    positive = surface[surface > 0]
    if positive.size == 0:
        return surface
    q = np.percentile(positive, CLIP_PERCENTILE)
    clipped = np.minimum(surface, q)
    stretched = np.log1p(clipped)
    lo, hi = stretched.min(), stretched.max()
    if hi <= lo:
        return np.zeros_like(stretched)
    return (stretched - lo) / (hi - lo)


def main() -> None:
    banner("03_build_tpi.py", "Building the TPI raster (KDE -> normalise -> weighted sum)")

    print("  Reading indicator layers")
    layers = {}
    for name in INDICATORS:
        gdf = read_layer(name)
        layers[name] = gdf
        print(f"    {name:12s} {len(gdf):>6} features")

    transform, width, height = build_grid(layers)
    print(f"\n  Grid: {width} x {height} cells at {CELL_SIZE_M} m "
          f"({width*height/1e6:.2f} M cells)")

    tpi = np.zeros((height, width), dtype="float64")
    profile = {
        "driver": "GTiff", "height": height, "width": width, "count": 1,
        "dtype": "float32", "crs": CRS_METRIC, "transform": transform,
        "nodata": 0.0, "compress": "lzw",
    }

    print(f"\n  Kernel density estimation (sigma = {SIGMA_CELLS} cells "
          f"= {SIGMA_CELLS*CELL_SIZE_M} m)")
    for name in INDICATORS:
        counts = rasterise_counts(layers[name], transform, width, height)
        surface = normalise(kde(counts))
        weight = WEIGHTS_S2_TOURISM[name]
        tpi += weight * surface

        out = PROC / f"kde_{name}.tif"
        with rasterio.open(out, "w", **profile) as dst:
            dst.write(surface.astype("float32"), 1)
        print(f"    {name:12s} w={weight:.3f}  max={surface.max():.4f}  -> {out.name}")

    # final rescale so the index spans the full unit interval
    if tpi.max() > 0:
        tpi = tpi / tpi.max()

    with rasterio.open(TPI_RASTER, "w", **profile) as dst:
        dst.write(tpi.astype("float32"), 1)
    print(f"\n  Unmasked index -> {TPI_RASTER.name}")

    # ---- land mask ----
    if LAND_MASK_SHP.exists():
        land = gpd.read_file(LAND_MASK_SHP).to_crs(CRS_METRIC)
        mask = geometry_mask(land.geometry, out_shape=(height, width),
                             transform=transform, invert=True)
        masked = np.where(mask, tpi, 0.0)
        with rasterio.open(TPI_MASKED, "w", **profile) as dst:
            dst.write(masked.astype("float32"), 1)
        valid = int((masked > 0).sum())
        print(f"  Land-masked index -> {TPI_MASKED.name}")
        print(f"  Valid cells: {valid:,} "
              f"({valid * (CELL_SIZE_M/1000)**2:,.0f} km2)")
    else:
        print(f"  NOTE: {LAND_MASK_SHP.name} not found; land mask skipped.")


if __name__ == "__main__":
    main()
