# Data provenance (Section 3.1 of the paper)

This document records the exact provenance of every input, corresponding to
Table 1 and Section 3.1 of the manuscript.

## Indicator layers

| Indicator | Count | Weight | Source | Retrieved |
|-----------|-------|--------|--------|-----------|
| Hotels | 2,418 | 0.25 | OpenStreetMap / Overpass API | July 2026 |
| Beaches | 467 | 0.20 | OpenStreetMap / Overpass API | July 2026 |
| Restaurants and cafés | 7,209 | 0.15 | OpenStreetMap / Overpass API | July 2026 |
| Marinas | 60 | 0.10 | OpenStreetMap / Overpass API | July 2026 |
| Camping areas | 184 | 0.10 | OpenStreetMap / Overpass API | July 2026 |
| Museums and heritage sites | 293 | 0.10 | OpenStreetMap / Overpass API | July 2026 |
| Ferry piers | 214 | 0.05 | OpenStreetMap / Overpass API | July 2026 |
| Blue Flag beaches | 45 | 0.05 | TÜRÇEV 2026 award list, geocoded | July 2026 |
| **Total** | **10,890** | **1.00** | | |

## Processing parameters

| Parameter | Value |
|-----------|-------|
| Coordinate reference system | WGS 84 / UTM Zone 35N (EPSG:32635) |
| Cell size | 250 m |
| Kernel | Gaussian, σ = 3 cells = 750 m |
| Outlier treatment | 95th-percentile clip, then log1p |
| Normalisation | linear rescale to [0, 1] per indicator |
| Aggregation | weighted linear sum, Σ wᵢ = 1 |
| Land mask | dissolved provincial polygons |

## Analytical parameters

| Parameter | Value |
|-----------|-------|
| Spatial weights | Queen contiguity + nearest-neighbour link for islands |
| Permutations | 9,999 |
| Significance levels | 90%, 95%, 99% |
| District filter | area ≥ 1 km², at least 20 non-zero TPI cells |

## Validation sources

1. **Official accommodation statistics** — Ministry of Culture and Tourism,
   district-level arrivals, overnight stays and occupancy rate for 2025.
2. **Review density** — user review counts attached to tourism establishments,
   retrieved through the Google Places API and assigned to districts by spatial
   containment.
3. **Designated tourism areas** — official register of Tourism Centres and Culture
   and Tourism Conservation and Development Regions.

None of the three validation sources enters the construction of the index.
