# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from . import METRICS


def alpha_phylogenetic(ctx, table, phylogeny, metric):
    metric = METRICS['NAME_TRANSLATIONS'][metric]

    action = ctx.get_action('diversity_lib', metric)
    vector, = action(table, phylogeny)

    return vector


def alpha(ctx, table, metric):
    if metric in METRICS['NONPHYLO']['IMPL']:
        metric = METRICS['NAME_TRANSLATIONS'][metric]
        action = ctx.get_action('diversity_lib', metric)
        vector, = action(table=table)
    else:
        action = ctx.get_action('diversity_lib', 'alpha_passthrough')
        vector, = action(table=table, metric=metric)

    return vector
