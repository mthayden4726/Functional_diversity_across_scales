#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that loops through individual files to correct, clip, mosaic and compute FRic & FDiv for plots within NEON site.

User must define parameters that feed into calling each script. 

Requires:
- correction coefficients
- shapefiles (plot boundaries)


Author: M. Hayden
Last Updated: 1/23/24
"""

# Load required libraries
import os

# Define variables
SITECODE = 'BART'
DOMAIN_ID = 'D01'
ID_NO = '5'
YEAR = '201908'
DATE = '20190825'
DATE_ID = '2019082513'
EPSG = '36219'
NDVI_THRESHOLD = 0.25


# Workflow to get Functional Richness and Functional Divergence for set of NEON plots (within a site)

# 1. Topo/BRDF corrections
script_string = "python 02_scripts/S02_TopoBRDF_Corrections.py --SITECODE " + SITECODE 
                + "--DOMAIN " + DOMAIN_ID 
                + "--ID_NO " + ID_NO  
                + "--DATE " + DATE 
                + "--DATE_ID " + DATE_ID 
                + "--EPSG " + EPSG 
                + "--NDVI " + NDVI_THRESHOLD
os.system(script_string)
print("Corrections complete.")

# 2. Clip to regions of interest
script_string = "python 02_scripts/S03_Clip_Corrected.py --SITECODE " + SITECODE 
                + "--YEAR " + YEAR 
os.system(script_string)
print("Flightlines clipped.")

# 3. Mosaic flightlines
script_string = "python 02_scripts/S04_Mosaic_Clipped_Raster.py --SITECODE " + SITECODE 
                + "--EPSG " + EPSG
os.system(script_string)
print("Flightlines mosaicked.")

# 4. Compute Functional Richness
script_string = "python 02_scripts/S05_Compute_FRic_FDiv.py --SITECODE " + SITECODE 
os.system(script_string)
print("FRic & FDiv computed for all mosaics.")

