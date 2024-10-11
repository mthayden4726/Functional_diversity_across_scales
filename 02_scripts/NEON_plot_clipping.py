## This script clips existing mosaics to the coordinates of the plot and creates a clip for the plot.

import pandas as pd
import rasterio
from rasterio.windows import from_bounds
from shapely.geometry import box
import pyproj
from rasterio.plot import show
import matplotlib.pyplot as plt

# Load the mosaic (raster) file
raster_path = "/Users/meha3816/Downloads/MOSAIC_BART_013.tif"
mosaic = rasterio.open(raster_path)

# Load the CSV with site data (site name, latitude, longitude)
csv_path = "/NOEN_sr_summaries/NEON_species_summaries.csv"
df = pd.read_csv(csv_path)

# Filter by name of mosaic to get the corresponding lat and long
#mosaicID = "BART_013"  # Replace this with the site name you're interested in
mosaicID = mosaic_filename.split("_")[1] + "_" + mosaic_filename.split("_")[2].split(".")[0]
filtered_site = df[df['plotID'] == mosaicID]

# Check if the site exists in the data
if filtered_site.empty:
    raise ValueError(f"Site {site_name} not found in the CSV.")

# Get the latitude and longitude for the selected site
lat = filtered_site.iloc[10]['latitude']
lon = filtered_site.iloc[10]['longitude']
coordinate = (lon, lat)

# Load the mosaic (raster) file
raster_path = "/Users/meha3816/Downloads/MOSAIC_BART_013.tif"
mosaic = rasterio.open(raster_path)

# Transform the coordinate to the CRS of the raster (usually a projected CRS)
project = pyproj.Transformer.from_crs("EPSG:4326", mosaic.crs, always_xy=True)  # Convert from WGS84 to the raster CRS
x_proj, y_proj = project.transform(coordinate[0], coordinate[1])

# Create a 20x20m bounding box (square) around the transformed coordinate
size = 10  # 10 meters in each direction from the center
square = box(x_proj - size, y_proj - size, x_proj + size, y_proj + size)

# Clip the mosaic using the bounding box
window = from_bounds(square.bounds[0], square.bounds[1], square.bounds[2], square.bounds[3], mosaic.transform)

# Read the clipped data
clipped_data = mosaic.read(window=window)

# Update metadata for the clipped raster
clipped_transform = mosaic.window_transform(window)
clipped_meta = mosaic.meta.copy()
clipped_meta.update({
    "driver": "GTiff",
    "height": clipped_data.shape[1],
    "width": clipped_data.shape[2],
    "transform": clipped_transform
})

# Save the clipped mosaic as a new file
output_path = f"/Users/meha3816/Downloads/clipped_mosaic_{site_name}.tif"
with rasterio.open(output_path, "w", **clipped_meta) as dest:
    dest.write(clipped_data)
