# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_diversity_lib.alpha import METRICS

all_phylo_metrics = METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL']
all_nonphylo_metrics = METRICS['NONPHYLO']['IMPL'] \
                       | METRICS['NONPHYLO']['UNIMPL']


def alpha_phylogenetic(ctx, table, phylogeny, metric):
    metric_tr = METRICS['NAME_TRANSLATIONS'][metric]

    f = ctx.get_action('diversity_lib', metric_tr)
    result = f(table, phylogeny)
    return tuple(result)


def alpha(ctx, table, metric):
    metric_tr = METRICS['NAME_TRANSLATIONS'][metric]

    if metric in METRICS['NONPHYLO']['IMPL']:
        func = ctx.get_action('diversity_lib', metric_tr)
        result = func(table=table)
    else:
        func = ctx.get_action('diversity_lib', 'alpha_passthrough')
        result = func(table=table, metric=metric_tr)

    return tuple(result)
