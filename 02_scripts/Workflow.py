# Import functions
from S01_Functions import * 
import csv

# Set working directories
Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# Set global parameters
window_sizes = [60, 120, 240, 480, 700, 960, 1200, 1500, 2000, 2200]
ndvi_threshold = 0.4 # ndvi threshold for radiometric filtering
comps = 3 # default component numbers for PCA

# Define NEON sites of interest
NEON_sites = ['TEAK', 'SERC', 'SRER']

# Pull flightlines for each site
# download from S3
# Load list into array
with open('NEON_flightlines.csv', newline='') as csvfile:
    all_flights = list(csv.reader(csvfile))
my_flights = all_flights[all_flights["Site"].isin(NEON_sites)]

## PROCESSING WORKFLOW ##

# 1. Topo/BRDF corrections
topo_brdf_correct(my_flights)
print("Corrections complete")

# 2. Clip to regions of interest
clip_flightline(corrected_flightlines)
print("Flightlines Clipped")

# 3. Mosaic flightlines
mosaic_flightlines(clipped_flightlines)
print("Flightlines mosaicked")

# 4. Compute Functional Richness
compute_fric(all_mosaics)
print("FRic computed for all mosaics")

# 5. Compute Functional Evenness
compute_feve(all_mosaics)
print("FEve computed for all mosaics")

# 6. Compute Functional Divergence
compute_fdiv(all_mosaics)
print("FDiv computed for all mosaics")
