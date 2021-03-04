# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from . import METRICS


def beta_phylogenetic(ctx,
                      table,
                      phylogeny,
                      metric,
                      threads=1,
                      variance_adjusted=False,
                      alpha=None,
                      bypass_tips=False):
    # TODO: remove when we can handle optional type-mapped parameters
    if alpha is not None and metric != 'generalized_unifrac':
        raise ValueError('The alpha parameter is only allowed when the choice'
                         ' of metric is generalized_unifrac')

    # TODO: this logic will be simpler once the remaining unifracs are
    # implemented in q2-diversity-lib
    if metric in ('unweighted_unifrac', 'weighted_unifrac') \
            and not variance_adjusted:
        metric = METRICS['NAME_TRANSLATIONS'][metric]
        action = ctx.get_action('diversity_lib', metric)
        dm, = action(table, phylogeny, threads=threads,
                     bypass_tips=bypass_tips)
    else:
        # handle unimplemented unifracs
        action = ctx.get_action('diversity_lib',
                                'beta_phylogenetic_passthrough')
        dm, = action(table, phylogeny, metric=metric, threads=threads,
                     variance_adjusted=variance_adjusted, alpha=alpha,
                     bypass_tips=bypass_tips)

    return dm


def beta(ctx, table, metric, pseudocount=1, n_jobs=1):
    if metric in METRICS['NONPHYLO']['IMPL']:
        metric = METRICS['NAME_TRANSLATIONS'][metric]
        action = ctx.get_action('diversity_lib', metric)
        dm, = action(table=table, n_jobs=n_jobs)
    else:
        action = ctx.get_action('diversity_lib', 'beta_passthrough')
        dm, = action(table=table, metric=metric, pseudocount=pseudocount,
                     n_jobs=n_jobs)

    return dm
