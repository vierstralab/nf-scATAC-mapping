from scipy.sparse import hstack, load_npz, save_npz
import sys
import pandas as pd
import numpy as np

def main(matrices, maps):
    return hstack(matrices), np.concatenate(maps)

if __name__ == '__main__':
    samples_meta = pd.read_table(sys.argv[1])
    matrices = [load_npz(f"{x}.barcodes.npz") for x in samples_meta['fragment_file_id']]
    maps = [np.loadtxt(f"{x}.unique_barcodes.map", dtype=str) for x in samples_meta['fragment_file_id']]
    concat_matrix, concat_map = main(matrices, maps)
    save_npz(f"{sys.argv[2]}.sparse_matrix.npz", concat_matrix)
    save_npz(f"{sys.argv[2]}.fragments.map", concat_map)
