import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import pdist, squareform
import csv
from tqdm import tqdm

def calculate_FEve(mstvect, dist_matrix, nb):
    EW = np.zeros(nb - 1)
    flag = 0
    for m in range((nb - 1) * nb // 2):
        if mstvect[m] != 0:
            EW[flag] = dist_matrix[m]
            flag += 1
    return EW
  
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
    windows, pca_chunk, results_FE, local_file_path  = args
    for window in tqdm(windows, desc='Processing window for batch'):
        print(window)
        half_window = window // 2
        FEve_values = []
        for i in tqdm(range(half_window, pca_chunk.shape[0] - half_window, 15), desc='Processing window index'):
            for j in range(half_window, pca_chunk.shape[1] - half_window, 15):
                sub_arr = pca_chunk[i - half_window:i + half_window + 1, j - half_window:j + half_window + 1, :]
                nbSpecies = sub_arr.shape[0] * sub_arr.shape[1]
                distances = pdist(sub_arr.reshape(nbSpecies, -1))
                dist_matrix = squareform(distances)
                mst = minimum_spanning_tree(dist_matrix)
                mstvect = mst.toarray().flatten()
                FEve = calculate_FEve(mstvect, dist_matrix, nbSpecies)
                print(FEve)
                FEve_values.extend(FEve)
                
        with open(local_file_path, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            if csvfile.tell() == 0:
                csvwriter.writerow(['Window_Size', 'FEve'])  # Write header
            for FEve_value in FEve_values:
                csvwriter.writerow([window, FEve_value])

    
