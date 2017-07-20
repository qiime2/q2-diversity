# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import biom
import skbio
import skbio.diversity
import skbio.tree
import sklearn.metrics


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
                      metric: str, n_jobs: int=1)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    try:
        results = skbio.diversity.beta_diversity(
            metric=metric,
            counts=counts,
            ids=sample_ids,
            otu_ids=feature_ids,
            tree=phylogeny,
            pairwise_func=sklearn.metrics.pairwise_distances,
            n_jobs=n_jobs
        )
    except skbio.tree.MissingNodeError as e:
        message = str(e).replace('otu_ids', 'feature_ids')
        message = message.replace('tree', 'phylogeny')
        raise skbio.tree.MissingNodeError(message)

    return results


def beta(table: biom.Table, metric: str, n_jobs: int=1)-> skbio.DistanceMatrix:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        pairwise_func=sklearn.metrics.pairwise_distances,
        n_jobs=n_jobs
    )
