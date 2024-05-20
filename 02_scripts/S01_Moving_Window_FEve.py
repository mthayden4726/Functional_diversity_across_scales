## This script is in progress - as of May 2024, I did not get functional evenness to work ##


import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
import csv
from tqdm import tqdm
from scipy.spatial import distance


def calculate_FEve(mstvect, dist_matrix, nbSpecies):
    EW = np.zeros(nbSpecies - 1)
    flag = 0
    for i in range(1, nbSpecies):  # Skip the first node (index 0) as it has no incoming edges
        # Find the parent node in the minimum spanning tree
        parent = (mstvect[:i] == i).nonzero()[0]
        if len(parent) > 0:  # Ensure there's a parent node
            parent = parent[0]
            # Calculate the edge weight as the distance between the node and its parent
            EW[flag] = dist_matrix[i, parent]
            flag += 1
    return EW

def primms_mst(raster):
    # Step 1: Calculate pairwise distances
    distances = distance.pdist(raster.reshape(-1, raster.shape[-1]))
    dist_matrix = distance.squareform(distances)

    # Step 2: Construct the minimum spanning tree
    mst = minimum_spanning_tree(dist_matrix)

    # Step 3: Compute the Euclidean distance of each branch (EW)
    mst = mst.toarray()
    EW = np.sqrt(np.sum((raster[mst > 0][:, None] - raster[mst > 0][None, :])**2, axis=-1))

    # Step 4: Calculate partial weighted evenness (PEW)
    S = np.count_nonzero(mst)
    PEW = EW/np.sum(EW)

    return PEW, S

def calculate_FEve_villager(PEW, S):
    threshold = 1 / (S - 1)  # Threshold value
    FEve = (np.sum(np.minimum(PEW, threshold)) - threshold) / (1 - threshold)
    return FEve
  
def window_calcs_feve(args):
    
    """ Calculate functional evenness using MST for a single PCA chunk and window size.
    FOR USE IN PARALLEL PROCESSING OF FUNCTIONAL EVENNESS.
    
    Parameters:
    -----------
    pca_chunk: PCA chunk from NEON image.
    window_sizes: list/array of integers
    
    Returns:
    -----------
    FEve: functional evenness for given window size and image.
    
    """
    windows, pca_x, results_FE, local_file_path  = args
    for window in tqdm(windows, desc='Processing window for batch'):
        print(window)
        half_window = window // 2
        FEve_values = []
        for i in tqdm(range(half_window, pca_x.shape[0] - half_window, 30), desc='Processing window index'):
            for j in range(half_window, pca_x.shape[1] - half_window, 30):
                sub_arr = pca_x[i - half_window:i + half_window + 1, j - half_window:j + half_window + 1, :]
                #if sub_arr is None:
                #    print(f"sub_arr is None at position ({i}, {j})")
                #    continue  # Skip this iteration if sub_arr is None
                print(f"sub_arr shape: {sub_arr.shape}")
                PEW, S = primms_mst(sub_arr)
                print("PEW:", PEW)
                print("S:", S)
                FEve = calculate_FEve_villager(PEW, S)
                print("FEve (Villéger et al.):", FEve_villager)
                FEve_values.extend(FEve)
                
        with open(local_file_path, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['Window_Size', 'FEve'])  # Write header
            for FEve_value in FEve_values:
                csvwriter.writerow([window, FEve_value])

    
