# Import functions
from S01_Functions import * 

# Universal variables

# Define NEON sites of interest
NEON_sites = ['TEAK', 'SERC', 'SRER']

# Pull flightlines for each site
# Write code here to pull in a table and subset to sites of interest - called 'all_flights'
# all_flights = 

## PROCESSING WORKFLOW ##

# 1. Topo/BRDF corrections
topo_brdf_correct(all_flights)
print("Corrections complete")

# 2. Clip to regions of interest

# 3. Mosaic flightlines

# 4. Compute Functional Richness

# 5. Compute Functional Evenness

# 6. Compute Functional Divergence
