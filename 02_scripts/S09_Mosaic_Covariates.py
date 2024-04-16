import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio.mask import mask
from rasterio.warp import transform_bounds
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject
import glob
import os
import boto3
from botocore.exceptions import NoCredentialsError, PartialCredentialsError, ClientError
from shapely.geometry import box
import geopandas as gpd
import fiona
from fiona.crs import from_epsg
import pycrs
from osgeo import gdal
import numpy as np
from S01_Functions import * # add scripts folder to python path manager

Data_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/01_rawdata'
Out_Dir = '/home/ec2-user/BioSCape_across_scales/01_data/02_processed'
bucket_name = 'bioscape.gra'
s3 = boto3.client('s3')

# Set config options
gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
gdal.SetConfigOption('CHECK_DISK_FREE_SPACE', 'FALSE')

src_files_to_mosaic = []

file_ID = [
              '005',
  '008',
  '010',
  '011',
  '018',
  '019',
  '021',
  '024',
  '030',
  '043',
  '073'
]

for i,ID in enumerate(file_ID):
    src_files_to_mosaic = []
    # List files associated with a single buffer shape
    search_criteria1 = str(ID)
    search_criteria2 = "DTM"
    search_criteria3 = "Clipped"
    dirpath = "Environmental_Covariates/ONAQ"

    # List objects in the S3 bucket in the matching directory
    objects = s3.list_objects_v2(Bucket=bucket_name, Prefix=dirpath)['Contents']
    # Filter objects based on the search criteria
    files = [obj['Key'] for obj in objects if obj['Key'].endswith('.tif') and (search_criteria1 in obj['Key'])
            and (search_criteria2 in obj['Key']) and (search_criteria3 in obj['Key'])]
    print(files)
    for j,file in enumerate(files):
        flight  = Out_Dir + '/file_' + str(j) + '.tif'
        try:
            s3.download_file(bucket_name, file, flight)
            print(f"The file '{file}' exists.")
        except Exception as e:
            print(f"Error: {e}")
        src = rasterio.open(flight)
        src_files_to_mosaic.append(src)

    # Mosaic files
    print(src_files_to_mosaic)
    mosaic, out_trans = merge(src_files_to_mosaic, method = 'max')
    print('Merge complete')
    # Update metadata
    out_meta = src.meta.copy()
    print(out_meta)
    print(mosaic.shape)
    out_meta.update({"driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans})
    print(out_meta)

    # Write to computer, send to S3
    local_file_path = Out_Dir + "/mosaic_ONAQ.tif"
    with rasterio.open(local_file_path, "w", **out_meta) as dest:
        dest.write(mosaic)
    print("File written")

    # Push to S3 bucket
    destination_s3_key = 'Environmental_Covariates/ONAQ/DTM_Mosaic_'+str(ID)+'.tif'
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    print("File uploaded to S3")

    # Get summary metrics
    with rasterio.open(local_file_path) as DTM_src:
      DTM_data = DTM_src.read(1, masked=True)

    # initialize summaries
    data = [['Mean:', DEM_data.mean()], 
            ['Max:', DEM_data.mean()],
            ['Min:', DEM_data.mean()],
            ['Std:', DEM_data.mean()],
            ['Var:', DEM_data.mean()]
           ]
 
    # Create the pandas DataFrame
    df = pd.DataFrame(data, columns=['Variables', 'Values'])
    df.to_csv('summary.csv')
    local_csv_path = 'summary.csv'
    
    # Push to S3 bucket
    destination_s3_key = 'Environmental_Covariates/ONAQ/DTM_Mosaic_Summary'+str(ID)+'.csv'
    upload_to_s3(bucket_name, local_csv_path, destination_s3_key)
    print("File uploaded to S3")
    
    # Remove unneeded files (mosaic and shapefile)
    os.remove(local_file_path)
    os.remove(local_csv_path)
  
    mosaic = None
    src_files_to_mosaic = None 
