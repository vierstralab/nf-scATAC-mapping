import pandas as pd
from scipy.sparse import coo_matrix, save_npz
import numpy as np
import sys


def main(melt_df, dhs_mapping, fragments_mapping):
    melt_df['value'] = melt_df.groupby(['dhs_id', 'fragment_id']).transform('count')
    melt_df['dhs_index'] = melt_df['dhs_id'].map(dhs_mapping)
    melt_df['fragment_index'] = melt_df['fragment_id'].map(fragments_mapping)
    return coo_matrix((melt_df['value'], (melt_df['dhs_index'], melt_df['fragment_index'])))


def create_mapping(filename):
    mapping = np.loadtxt(filename, dtype=str)
    return {value: idx for idx, value in enumerate(mapping)}


if __name__ == '__main__':
    data = pd.read_table(sys.argv[1], header=None, names=['dhs_id', 'fragment_id'])

    fragments_mapping = create_mapping(sys.argv[32])
    dhs_mapping = create_mapping(sys.argv[3])

    sparse_matrix = main(data, fragments_mapping=fragments_mapping, dhs_mapping=dhs_mapping)
    save_npz(sys.argv[4], sparse_matrix)
