# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
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
import psutil

from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat

from functools import partial


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'unweighted_unifrac', 'weighted_unifrac'}


def phylogenetic_metrics_alt_dict():
    return {'unweighted_unifrac': unifrac.unweighted,
            'weighted_unifrac': unifrac.weighted_unnormalized,
            'weighted_normalized_unifrac': unifrac.weighted_normalized,
            'generalized_unifrac': unifrac.generalized}


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


def beta_phylogenetic_alt(table: BIOMV210Format, phylogeny: NewickFormat,
                          metric: str, n_jobs: int=1,
                          variance_adjusted: bool=False,
                          alpha: float=None,
                          bypass_tips: bool=False) -> skbio.DistanceMatrix:

    metrics = phylogenetic_metrics_alt_dict()
    generalized_unifrac = 'generalized_unifrac'

    if metric not in metrics:
        raise ValueError("Unknown metric: %s" % metric)

    if alpha is not None and metric != generalized_unifrac:
        raise ValueError('The alpha parameter is only allowed when the choice'
                         ' of metric is generalized_unifrac')

    # this behaviour is undefined, so let's avoid a seg fault
    cpus = psutil.cpu_count(logical=False)
    if n_jobs > cpus:
        raise ValueError('The value of n_jobs cannot exceed the number of '
                         'processors (%d) available in this system.' % cpus)

    if metric == generalized_unifrac:
        alpha = 1.0 if alpha is None else alpha
        f = partial(metrics[metric], alpha=alpha)
    else:
        f = metrics[metric]

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
