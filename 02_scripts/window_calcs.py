#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 10:25:17 2023

@author: meha3816
"""
## Load necessary packages ##
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
import subprocess
from urllib.request import urlretrieve
import parmap
import os
import tqdm 

def window_calcs(window_size, pca_chunk, comps):
    
    """ Calculate convex hull volume for a single PCA chunk and window size.
    FOR USE IN PARALLEL PROCESSING OF FUNCTIONAL RICHNESS.
    
    Parameters:
    -----------
    pca_chunk: PCA chunk from NEON image.
    window_sizes: list/array of integers
    comps: Number of PCs. Here, set to 4. 
    
    Returns:
    -----------
    volume_mean: functional richness for given window size and image.
    
    """
    half_window = window_size//2
    results_FR = []
    fric = np.zeros(pca_chunk.shape)
    for i in range(half_window, pca_chunk.shape[0]-half_window):
        for j in range(half_window, pca_chunk.shape[1]-half_window):
            sub_arr = pca_chunk[i-half_window:i+half_window+1, j-half_window:j+half_window+1, :]
            sub_arr = sub_arr.reshape((sub_arr.shape[0]*sub_arr.shape[1],comps))
            if np.nanmean(sub_arr) == 0.0:
                continue
            hull = ConvexHull(sub_arr) 
            fric[i,j]= hull.volume
    results_FR.append(np.nanmean(fric))        
    return results_FR