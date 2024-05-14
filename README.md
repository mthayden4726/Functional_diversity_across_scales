# BioSCape_across_scales

**Goal:** Process airborne hyperspectral imagery to produce maps of functional richness, divergence and evenness.

**Scripts implemented in Python**

See [Documentation.md](https://github.com/mthayden4726/BioSCape_across_scales/blob/25c5cc9100f1ee14097ed48f8af2e4acf129dd67/Documentation.md) for detailed documentation on the workflow, scripts, and required user inputs. 

**S01:** Module containing functions that pre-process and compute diversity metrics.

**S02:** Image Correction (BRDF & Topo)

**S03:** Image Clipping (Clip flightlines to shapefiles around field sites)

**S04:** Mosaic Images (Mosaic clipped rasters for each field site)

**S05:** Compute Functional Metrics.

**S04:** Model Scaling Relationships (Functional Metrics ~ Area)

Processing scripts should be implemented as:
1. Correct (BRDF/Topo), mask (NDVI threshold) and export flightlines as tiffs (**S02** scripts)
2. Clip flightlines to larger field site boundary (**S03** scripts)
3. Mosaic clipped flightlines together (**S04** scripts)
4. Compute functional richness and divergence metrics for each field site's raster (**S05** scripts)
5. Analyze the output to assess scaling of functional metrics with area (**S06** scripts)

**Authors:** M. Hayden, C. Amaral, M. Rossi, N. Stavros

Updated: May 14th, 2024

Next: Combining scripts into workflow based on FUNCTIONS (currently requires input for each site) --> script in progress: Workflow.py

Currently: Change site code and variables for each run; see scripts with various site code suffices.
