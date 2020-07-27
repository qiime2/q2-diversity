# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unifrac


# TODO: remove these collections ASAP
def phylogenetic_metrics_dict():
    return {'unweighted_unifrac': unifrac.unweighted,
            'weighted_unifrac': unifrac.weighted_unnormalized,
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


def beta(ctx, table, metric, pseudocount=1, n_jobs=1):

    func = ctx.get_action('diversity_lib', 'beta_dispatch')
    distance_matrix = func(table, metric, pseudocount, n_jobs)
    return tuple(distance_matrix)
