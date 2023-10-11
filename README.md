# BioSCape_across_scales

**Goal:** Process airborne hyperspectral imagery to produce maps of functional richness, divergence and evenness.

**Scripts implemented in Python**

**S01:** Module containing functions that pre-process and compute diversity metrics.

**S02:** Full step-by-step workflow for single image.

**S03:** Workflow for looping through many images; user can choose site, NEON product and year in global parameters and produce FRic and CV values for a range of defined window sizes.

**S04:** Testing workflow on single image with AWS - parallelization in progress.

Processing scripts should be implemented as:
1. Correct (BRDF/Topo) and export flightlines as tiffs (mosaic_flightlines.py)
2. Clip flightlines to larger field site boundary (clip_flightlines.py)
3. Mosaic clipped flightlines together (merge_clipped.py)
4. Clip mosaic down to individual rasters for each field site (~1x1km) (clip_mosaic.py)
5. Implement workflow for functional richness (S03) using the raster around each field site.

**Authors:** M. Hayden, C. Amaral, M. Rossi, N. Stavros

Updated: Oct 11. 2023
