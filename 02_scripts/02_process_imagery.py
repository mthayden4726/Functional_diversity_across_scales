#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for processing NEON imagery.

Goal: compute alpha diversity and functional richness for an array of
window sizes.


Created on Fri Mar 24 2023
Last updated on Mon Apr 17 2023

Authors: M. Hayden
"""

# 1. Set-up - load packages, data, and functions

# Import packages
import hytools as ht
import matplotlib.pyplot as plt
import numpy as np
import requests
import sklearn
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import kneed
from kneed import KneeLocator
from scipy.spatial import ConvexHull

# Read in functions from other scripts? 

# Set global parameters
window_sizes = [10, 15, 20, 30, 40, 50, 75, 100, 150]  # list of window sizes to test

# Set radiometric thresholds
ndvi_threshold = 0.4
# shade_threshold = 'NA' # should be low threshold across NIR region (15%)
# cloud_threshold = 'NA' # should be high threshold in blue band

# Set default component numbers for A) PCA and B) K-means clustering
comps = 4
nclusters = 15

# Import data
url = 'https://storage.googleapis.com/neon-aop-products/2019/FullSite/D14/2019_SRER_3/L1/Spectrometer/ReflectanceH5/2019090414/NEON_D14_SRER_DP1_20190904_165832_reflectance.h5'    
file = 'NEON_D14_SRER_DP1_20190904_165832_reflectance.h5'
# Request the url
data_request = requests.get(url)
print(data_request.status_code) #200 means success
open(file, 'wb').write(data_request.content)

# Read in file
neon_image = '/Users/meha3816/Desktop/bioscape_proj/bioscape_specdiv/01_data/01_rawdata/NEON_D14_SRER_DP1_20190904_165832_reflectance.h5'
neon = ht.HyTools()
neon.read_file(neon_image,'neon')

# Review file
neon.map_info

# Show image
show_rgb(neon, correct = [])

# 2. Image Preprocessing - remove non-vegetated pixels

# Set bad bands for removal
neon.create_bad_bands([[300,400],[1300,1450],[1780,2000],[2450,2600]])

# sensitivity analysis - this will need to vary depending on the ecosystem and pixel size
# include options for various ecosystems and pixel size

# Calculate NDVI
ndvi = neon.ndi()

# Set a mask from which to sample for PCA where NDVI is greater than threshold
# I.e., select for only live vegetation pixels
neon.mask['sample'] = (neon.mask['no_data']) & (ndvi > ndvi_threshold)

# View mask
plt.matshow(neon.mask['sample'])


# 3a. Dimensionality Reduction - Select subset of image to initiate PCA

# Choose percent of unmasked pixels to sample
sample_size = 0.1

# Subsample a set of pixels for PCA
sub_samples = np.zeros((neon.lines,neon.columns)).astype(bool)
idx = np.array(np.where(neon.mask['sample'])).T
idxRand= idx[np.random.choice(range(len(idx)),int(len(idx)*sample_size), replace = False)].T
sub_samples[idxRand[0],idxRand[1]] = True
neon.mask['samples'] = sub_samples

X = []

# Create subset of data that excludes bad bands
# should be a list of 336 arrays (or however many bands)
for band_num,band in enumerate(neon.bad_bands):
    if ~band:
        X.append(neon.get_band(band_num,mask='samples'))

X = np.array(X).T

# Attempting to speed that up by eliminating for loop

mask = ~neon.bad_bands #make mask the inverse, TRUE is now good bands
X = neon.get_band(mask, mask='samples') #sample pixels from the good bands


# 3b. Dimensionality Reduction - Apply PCA to subset - could be own function

# Center, scale and fit PCA transform - scales based on mean reflectance at each band
x_mean = X.mean(axis=0)[np.newaxis,:]
X = X.astype('float32') # necessary to manually convert to float for next function to work
X -=x_mean
x_std = X.std(axis=0,ddof=1)[np.newaxis,:]
X /=x_std
X = X[~np.isnan(X.sum(axis=1)) & ~np.isinf(X.sum(axis=1)),:]

# Perform initial PCA fit
pca = PCA(n_components=15) # set max number of components
pca.fit(X)

fig = plt.figure(figsize = (6,4))
ax = fig.add_subplot(111)
ax.plot(np.cumsum(pca.explained_variance_ratio_), c= 'r',lw=2)
ax.set_xlabel("Principal components")
_ = ax.set_ylabel("Cumulative explained variance ratio")

# Choose number of components that explain much variance
comps = int(input("Number of components to use: "))
# Using 2 for now

# Rerun with above number of components
pca = PCA(n_components=comps)
pca.fit(X)
pca_out = pca.fit(X)

# 3c. Dimensionality Reduction - Apply PCA to entire image - could be own function

pca_transform = np.zeros((neon.lines,neon.columns,comps))
iterator = neon.iterate(by = 'chunk',chunk_size = (500,500))

while not iterator.complete:
    chunk = iterator.read_next()
    X_chunk = chunk[:,:,~neon.bad_bands].astype(np.float32)
    X_chunk = X_chunk.reshape((X_chunk.shape[0]*X_chunk.shape[1],X_chunk.shape[2]))
    X_chunk -=x_mean
    X_chunk /=x_std
    X_chunk[np.isnan(X_chunk) | np.isinf(X_chunk)] = 0
    pca_chunk=  pca.transform(X_chunk)
    pca_chunk = pca_chunk.reshape((chunk.shape[0],chunk.shape[1],comps))
    pca_chunk[chunk[:,:,0] == neon.no_data] =0
    pca_transform[iterator.current_line:iterator.current_line+pca_chunk.shape[0],
                  iterator.current_column:iterator.current_column+pca_chunk.shape[1]] = pca_chunk

# See how much variation each component explains
exp_var_pca = pca.explained_variance_ratio_
exp_var_pca

# Plot PCA bands

red = pca_transform[:,:,0]
green = pca_transform[:,:,1]
blue = pca_transform[:,:,2]
prgb=  np.stack([red,green,blue])
prgb = np.moveaxis(prgb,0,-1).astype(float)
prgb[prgb ==neon.no_data] = np.nan
bottom = np.nanpercentile(prgb,5,axis = (0,1))
top = np.nanpercentile(prgb,95,axis = (0,1))
prgb = np.clip(prgb,bottom,top)
prgb = (prgb-np.nanmin(prgb,axis=(0,1)))/(np.nanmax(prgb,axis= (0,1))-np.nanmin(prgb,axis= (0,1)))
prgb[~neon.mask['no_data']] = 1

fig = plt.figure(figsize = (6,6),facecolor='w')
fig.tight_layout()
ax = fig.add_subplot(111)
ax.matshow(prgb)
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# 4a. Unsupervised clustering - Choose number of clusters, K
# Initiate K-means clustering
kmeans_kwargs = {
   "init": "random", # chooses n clusters non-randomly to speed convergence
   "n_init": 10,
   "max_iter": 300,
   "random_state": 42,
   }

pca_sample= pca.transform(X.astype(float))

# Loop through possible K values of 2-50 to determine optimal cluster number
# Record intertia for each value
SSE = []
for cluster in range(5, 50, 5):
    clusters = KMeans(n_clusters=cluster, **kmeans_kwargs)
    clusters.fit(pca_sample)
    SSE.append(clusters.inertia_)

# Plot inertia against number of clusters (k) and determine elbow

plt.style.use("fivethirtyeight")
plt.plot(range(5, 50, 5), SSE)
plt.xticks(range(5, 50, 5))
plt.xlabel("Number of Clusters")
plt.ylabel("SSE")
plt.show()

kl = KneeLocator(
   range(5, 50, 5), SSE, curve="convex", direction="decreasing"
   )
print("Best cluster number:", kl.elbow)

# reshape to 2d array to fit with kmeans algorithm below
nx, ny, nsamples = pca_transform.shape
reshape = pca_transform.reshape((nsamples,nx*ny)).T

# check that it matches original data
neon.lines*neon.columns

# 4b. Unsupervised Clustering - Run K-means algorithm on PCA

# Select number of clusters
nclusters = kl.elbow
# nclusters = int(input("Number of clusters: ")) #this was from previous without elbow plots

pca_trans= reshape
clusters = KMeans(n_clusters=nclusters)
clusters.fit(pca_trans)
classes = clusters.predict(pca_trans.reshape(neon.lines*neon.columns,comps)).reshape(neon.lines,neon.columns)
classes[~neon.mask['no_data']] = nclusters

fig = plt.figure(figsize = (6,6),facecolor='w')
fig.tight_layout()
ax = fig.add_subplot(111)
ax.matshow(classes,cmap = plt.get_cmap('tab20c'))
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')


# 5.  Alpha Diversity - Run algorithm to compute alpha diversity, plot results

# For now, not using the calc_alpha from sister_biodiv because of how slow it is
# Calculate alpha diversity for a range of window sizes

results = calc_alpha(arr = classes, windows = window_sizes)

# to access individual elements in dict: results_tab["window_25"]
# or take means of all items 
Shann_list = {k:np.nanmean(v) for k, v in results.items()}

# Plot results
names = list(Shann_list.keys())
values = list(Shann_list.values())

plt.bar(range(len(Shann_list)), values, tick_label=names)
plt.show()

# 6. Functional Richness
# Use similar structure but compute complex hull volume based on PCs

# If each row is a pixel and each column is a trait
# For now use PCs as traits
# First do with a small sample to see if it works
# Then use iter to iterate through chunks of the image (will be 200 chunks)

# points need to be 2-D

chv_results = calc_fun_rich(neon, window_sizes)
        

##########

