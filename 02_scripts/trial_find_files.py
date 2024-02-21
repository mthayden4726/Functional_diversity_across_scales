from S01_Functions import * # add scripts folder to python path manager
SITECODE = 'HEAL'
PRODUCTCODE = 'DP1.30006.001'
YEAR = '2019-08'

SERVER = 'http://data.neonscience.org/api/v0/'
# Build so that we can loop through reading files in
url = SERVER+'data/'+PRODUCTCODE+'/'+SITECODE+'/'+YEAR
# Request the url
data_request = requests.get(url)
# Convert the request to Python JSON object
data_json = data_request.json()
print(data_json)
# Create list of file paths of interest for given site, product, year
file_paths = []
for file in data_json['data']['files'][:20]:
  if 'reflectance.h5' in file['name']:
    file_paths.insert(1, file['url'])
    print(file['url'])
    
print(file_paths)
