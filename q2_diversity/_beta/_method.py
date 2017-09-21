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
import unifrac

from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat

from functools import partial

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


def all_metrics():
    return phylogenetic_metrics() | non_phylogenetic_metrics()


def beta_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                      metric: str, n_jobs: int=1)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")
    if n_jobs != 1 and metric == 'weighted_unifrac':
        raise ValueError("Weighted UniFrac is not parallelizable")

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


def beta_phylogenetic_hpc(table: BIOMV210Format, phylogeny: NewickFormat,
                          metric: str, n_jobs: int=1,
                          variance_adjusted: bool=False,
                          alpha=1.0,
                          bypass_tips: bool=False) -> skbio.DistanceMatrix:
    if metric == 'unweighted_unifrac':
        f = unifrac.unweighted
    elif metric == 'weighted_unnormalized_unifrac':
        f = unifrac.weighted_unnormalized
    elif metric == 'weighted_normalized_unifrac':
        f = unifrac.weighted_normalized
    elif metric == 'generalized_unifrac':
        f = partial(unifrac.generalized, alpha=alpha)
    else:
        raise ValueError("Unknown metric: %s" % metric)

    # unifrac processes tables and trees should be filenames
    return f(str(table), str(phylogeny), threads=n_jobs,
             variance_adjusted=variance_adjusted, bypass_tips=bypass_tips)


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
