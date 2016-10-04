# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

import qiime
import biom
import skbio
import skbio.diversity
import numpy
import pandas as pd


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'unweighted_unifrac', 'weighted_unifrac'}


def non_phylogenetic_metrics():
    return {'cityblock', 'euclidean', 'seuclidean', 'sqeuclidean', 'cosine',
            'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
            'braycurtis', 'mahalanobis', 'yule', 'matching', 'dice',
            'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
            'sokalsneath', 'wminkowski'}


def beta_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                      metric: str)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        otu_ids=feature_ids,
        tree=phylogeny
    )


def beta(table: biom.Table, metric: str)-> skbio.DistanceMatrix:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids
    )


def bioenv(output_dir: str, distance_matrix: skbio.DistanceMatrix,
           metadata: qiime.Metadata) -> None:
    # convert metadata to numeric values where applicable, drop the non-numeric
    # values, and then drop samples that contain NaNs
    df = metadata.to_dataframe()
    df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    df = df.select_dtypes([numpy.number]).dropna()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(df.index)
    filtered_dm_length = distance_matrix.shape[0]

    result = skbio.stats.distance.bioenv(distance_matrix, df)
    index_fp = os.path.join(output_dir, 'index.html')
    with open(index_fp, 'w') as fh:
        fh.write('<html><body>')
        if initial_dm_length != filtered_dm_length:
            fh.write("<b>Warning</b>: Some samples were filtered from the "
                     "input distance matrix because they were missing "
                     "metadata values.<br><b>The input distance matrix "
                     "contained %d samples but bioenv was computed on "
                     "only %d samples.</b><p>"
                     % (initial_dm_length, filtered_dm_length))
        fh.write(result.to_html())
        fh.write('</body></html>')
