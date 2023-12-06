import pandas as pd
from scipy.sparse import coo_matrix, save_npz
import numpy as np
import sys


def main(melt_df, dhs_mapping, fragments_mapping):
    melt_df = melt_df.groupby(['dhs_id', 'fragment_id']).size().reset_index(name='value')
    melt_df['dhs_index'] = melt_df['dhs_id'].map(dhs_mapping)
    melt_df['fragment_index'] = melt_df['fragment_id'].map(fragments_mapping)
    return coo_matrix((melt_df['value'], (melt_df['dhs_index'], melt_df['fragment_index'])),
                    shape=(len(dhs_mapping), len(fragments_mapping)), dtype=np.uint16)


def create_mapping(filename):
    mapping = np.loadtxt(filename, dtype=str)
    return {value: idx for idx, value in enumerate(mapping)}


if __name__ == '__main__':
    print("Reading files")
    data = pd.read_table(sys.argv[1], header=None, names=['dhs_id', 'fragment_id'])

    fragments_mapping = create_mapping(sys.argv[2])
    dhs_mapping = create_mapping(sys.argv[3])
    print("Converting to sparse matrix")
    sparse_matrix = main(data, fragments_mapping=fragments_mapping, dhs_mapping=dhs_mapping)
    print("Saving sparse matrix")
    save_npz(sys.argv[4], sparse_matrix)
