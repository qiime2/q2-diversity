# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_diversity_lib.alpha import METRICS
from q2_diversity_lib import translate_metric_name

all_phylo_metrics = METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL']
all_nonphylo_metrics = METRICS['NONPHYLO']['IMPL'] \
                       | METRICS['NONPHYLO']['UNIMPL']
metric_name_translations = METRICS['METRIC_NAME_TRANSLATIONS']


def alpha_phylogenetic(ctx, table, phylogeny, metric):
    metric_tr = translate_metric_name(metric, metric_name_translations)

    f = ctx.get_action('diversity_lib', metric_tr)
    result = f(table, phylogeny)
    return tuple(result)


def alpha(ctx, table, metric):
    implemented_metrics = METRICS['NONPHYLO']['IMPL']
    metric_tr = translate_metric_name(metric, metric_name_translations)

    if metric in implemented_metrics:
        func = ctx.get_action('diversity_lib', metric_tr)
        result = func(table=table)
    else:
        func = ctx.get_action('diversity_lib', 'alpha_passthrough')
        result = func(table=table, metric=metric_tr)

    return tuple(result)
