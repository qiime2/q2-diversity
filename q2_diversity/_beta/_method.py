# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

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
    if metric not in all_phylo_metrics:
        raise ValueError("Unknown metric: %s" % metric)

    if alpha is not None and metric != 'generalized_unifrac':
        raise ValueError('The alpha parameter is only allowed when the choice'
                         ' of metric is generalized_unifrac')

    # HACK: this logic will be simpler once the remaining unifracs are done
    if metric in ('unweighted_unifrac', 'weighted_unifrac') \
            and not variance_adjusted:
        func = ctx.get_action('diversity_lib', metric)
    else:
        # handle unimplemented unifracs
        func = ctx.get_action('diversity_lib', 'beta_phylogenetic_passthrough')
        func = partial(func, metric=metric, alpha=alpha,
                       variance_adjusted=variance_adjusted)

    result = func(table, phylogeny, threads=threads,
                  bypass_tips=bypass_tips)
    return tuple(result)


def beta(ctx, table, metric, pseudocount=1, n_jobs=1):

    func = ctx.get_action('diversity_lib', 'beta_dispatch')
    distance_matrix = func(table, metric, pseudocount, n_jobs)
    return tuple(distance_matrix)
