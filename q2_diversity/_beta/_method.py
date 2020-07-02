# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
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
import numpy as np

from skbio.stats.composition import clr
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import jensenshannon


def phylogenetic_metrics_dict():
    return {'unweighted_unifrac': unifrac.unweighted,
            'weighted_unnormalized_unifrac': unifrac.weighted_unnormalized,
            'weighted_normalized_unifrac': unifrac.weighted_normalized,
            'generalized_unifrac': unifrac.generalized}


def phylogenetic_metrics():
    return set(phylogenetic_metrics_dict())


def non_phylogenetic_metrics():
    return {'cityblock', 'euclidean', 'seuclidean', 'sqeuclidean', 'cosine',
            'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
            'braycurtis', 'mahalanobis', 'yule', 'matching', 'dice',
            'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
            'sokalsneath', 'wminkowski', 'aitchison', 'canberra_adkins',
            'jensenshannon'}


def all_metrics():
    return phylogenetic_metrics() | non_phylogenetic_metrics()


def beta_phylogenetic(ctx, table, phylogeny,
                      metric, threads=1,
                      variance_adjusted=False,
                      alpha=None,
                      bypass_tips=False):

    func = ctx.get_action('diversity_lib', 'beta_phylogenetic_dispatch')
    distance_matrix = func(table, phylogeny, metric=metric, threads=threads,
                           variance_adjusted=variance_adjusted, alpha=alpha,
                           bypass_tips=bypass_tips)
    return tuple(distance_matrix)


def beta(table: biom.Table, metric: str,
         pseudocount: int = 1, n_jobs: int = 1) -> skbio.DistanceMatrix:

    if not (metric in non_phylogenetic_metrics()):
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().T

    def aitchison(x, y, **kwds):
        return euclidean(clr(x), clr(y))

    def canberra_adkins(x, y, **kwds):
        if (x < 0).any() or (y < 0).any():
            raise ValueError("Canberra-Adkins is only defined over positive "
                             "values.")

        nz = ((x > 0) | (y > 0))
        x_ = x[nz]
        y_ = y[nz]
        nnz = nz.sum()

        return (1. / nnz) * np.sum(np.abs(x_ - y_) / (x_ + y_))

    def jensen_shannon(x, y, **kwds):
        return jensenshannon(x, y)

    if metric == 'aitchison':
        counts += pseudocount
        metric = aitchison
    elif metric == 'canberra_adkins':
        metric = canberra_adkins
    elif metric == 'jensenshannon':
        metric = jensen_shannon

    if table.is_empty():
        raise ValueError("The provided table object is empty")

    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        validate=True,
        pairwise_func=sklearn.metrics.pairwise_distances,
        n_jobs=n_jobs
    )
