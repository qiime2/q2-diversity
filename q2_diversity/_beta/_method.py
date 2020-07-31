# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_diversity_lib.beta import METRICS

all_phylo_metrics = METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL']
all_nonphylo_metrics = METRICS['NONPHYLO']['IMPL'] \
                       | METRICS['NONPHYLO']['UNIMPL']


def all_metrics():
    return all_phylo_metrics | all_nonphylo_metrics


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
