# Push to S3 bucket
    destination_s3_key = 'Environmental_Covariates/' + SITECODE + '/' + ENV + '_Mosaic_' +str(ID)+ '.tif'
    upload_to_s3(bucket_name, local_file_path, destination_s3_key)
    print("File uploaded to S3")

    # Get summary metrics
    with rasterio.open(local_file_path) as data_src:
      env_data = data_src.read(1, masked=True)

    # initialize summaries
    data = [{'Site': SITECODE,
             'Plot': str(ID), 
             'Env': ENV,
             'Mean': env_data.mean(),
             'Median': np.median(env_data),
             'Max': env_data.max(),
             'Min': env_data.min(),
             'Std': env_data.std(),
             'Var': env_data.var()
                            }]
    # append to dataframe
    summary_data = summary_data.append(data, ignore_index=True)
    print(summary_data)
    
    # Remove unneeded files (mosaic and shapefile)
    os.remove(local_file_path)
  
    mosaic = None
    src_files_to_mosaic = None 

# Create the pandas DataFrame
print(summary_data)
summary_data.to_csv('summary.csv')
local_csv_path = 'summary.csv'
    
# Push to S3 bucket
destination_s3_key = 'Environmental_Covariates/Summary_files/' + SITECODE + '_' + ENV + '_Summary.csv'
upload_to_s3(bucket_name, local_csv_path, destination_s3_key)
print("File uploaded to S3")

os.remove(local_csv_path)
