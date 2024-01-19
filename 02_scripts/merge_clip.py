import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
import glob
import os
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
from osgeo import gdal

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

def upload_to_s3(bucket_name, file_path, s3_key):
    """
    Upload a file from an EC2 instance to Amazon S3.

    :param bucket_name: Name of the S3 bucket
    :param file_path: Local path to the file on the EC2 instance
    :param s3_key: Destination key in the S3 bucket (e.g., folder/file_name.ext)
    """
    # Initialize the S3 client
    s3 = boto3.client('s3')

    try:
    # Upload the file
        s3.upload_file(file_path, bucket_name, s3_key)
        print(f'Successfully uploaded {file_path} to {bucket_name}/{s3_key}')
    except Exception as e:
        print(f"Error uploading file: {e}")

#search_criteria = "*.tif"
#dirpath = "TEAK_flightlines/"
#q = os.path.join(dirpath, search_criteria)

files = ['TEAK_flightlines/20190615_171251_output_0.tif', 
         'TEAK_flightlines/20190615_171251_output_1.tif',
         'TEAK_flightlines/20190615_171251_output_2.tif',
         'TEAK_flightlines/20190615_171251_output_3.tif']
        # 'TEAK_flightlines/20190615_171251_output_4.tif',
        # 'TEAK_flightlines/20190615_171251_output_5.tif',
        # 'TEAK_flightlines/20190615_171251_output_6.tif',
        # 'output_fullarray_170625.tif']

# Load the polygon for clipping ()
#clip_file = 'TEAK_clip.shp'
gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
gdal.SetConfigOption('CHECK_DISK_FREE_SPACE', 'FALSE')
#s3.download_file(bucket_name, clip_file, Out_Dir + '/clip_polygon.shp')
#clip_polygon = gpd.read_file(Out_Dir + '/clip_polygon.shp')

# Empty list with files to mosaic
#datasets = []

#for i, file in enumerate(files):
#    print('Loading file from S3')
#    s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
#    flight = Out_Dir + '/file_' + str(i) + '.tif'
#    print(flight)
#    with rasterio.open(flight) as src:
#        clipped_data, transform = rasterio.mask.mask(src, clip_polygon.geometry, crop=True)
#        clipped_dataset = rasterio.open('clipped_file_' + str(i) + '.tif', 'w', **src.profile)
#        print(clipped_dataset)
#        clipped_dataset.write(clipped_data)
#        datasets.append(clipped_dataset)
#print(datasets)
# Open in read mode and add to file list
src_files_to_mosaic = []
for i,file in enumerate(files):
    s3.download_file(bucket_name, file, Out_Dir + '/file_' + str(i) + '.tif')
    flight  = Out_Dir + '/file_' + str(i) + '.tif'
    src = rasterio.open(flight)
    print(src.meta)
    src_files_to_mosaic.append(src)

print(src_files_to_mosaic)

mosaic, out_trans = merge(src_files_to_mosaic, nodata = -9999)
print('Merge complete')

# Merge files
#bounds = [-119.9968, 36.925, -118.9968, 37.0104]
#mosaic, out_trans = merge(datasets)
#print('Merge complete')
#os.remove(clipped_dataset)
# Clip to smaller extent
# Create bounding box
# minx, miny = 36.9625, -118.9968
# maxx, maxy = 37.0104, -119.0708
# bbox = box(minx, miny, maxx, maxy)
# Insert box into GeoDataFrame
# geo = gpd.GeoDataFrame({'geometry':bbox},index = [0], crs = from_epsg(4326))
# Re-project into same as raster data
# geo = geo.to_crs(crs = mosaic.crs.data)
# Convert to rasterio ready coordinates
# def getFeatures(gdf):
#    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
#    import json
#    return [json.loads(gdf.to_json())['features'][0]['geometry']]
# coords = getFeatures(geo)
# print(coords)
# Clip the raster
# out_img, out_transform = mask(raster = mosaic, shapes = coords, crop = True)
# Update metadata
out_meta = src.meta.copy()
print(out_meta)
print(mosaic.shape)
out_meta.update({"driver": "GTiff",
   "height": mosaic.shape[1],
   "width": mosaic.shape[2],
   "transform": out_trans,
                 "crs": "+init=epsg:32611 +units=m +no_defs "})
print(out_meta)
# Write to computer, send to S3
mosaic_name = Out_Dir + "/mosaic.tif"
with rasterio.open(mosaic_name, "w", **out_meta) as dest:
    dest.write(mosaic)
print("File written")

# Push to S3 bucket
destination_s3_key = 'TEAK_flightlines/Mosaic_0-3.tif'
local_file_path = Out_Dir + '/mosaic.tif'
upload_to_s3(bucket_name, local_file_path, destination_s3_key)
print("File uploaded to S3")
os.remove(local_file_path)


