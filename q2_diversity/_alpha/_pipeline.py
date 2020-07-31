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


METRIC_NAME_TRANSLATIONS = {'shannon': 'shannon_entropy',
                            'pielou_e': 'pielou_evenness'}


def _translate_metric_name_for_div_lib(m: str) -> str:
    return METRIC_NAME_TRANSLATIONS[m] if m in METRIC_NAME_TRANSLATIONS else m


def alpha_phylogenetic(ctx, table, phylogeny, metric):
    metrics = all_phylo_metrics
    if metric not in metrics:
        raise ValueError("Unknown metric: %s" % metric)

    metric = _translate_metric_name_for_div_lib(metric)

    f = ctx.get_action('diversity_lib', metric)
    result = f(table, phylogeny)
    return tuple(result)


def alpha(ctx, table, metric):
    implemented_metrics = METRICS['NONPHYLO']['IMPL']
    if metric not in all_nonphylo_metrics:
        raise ValueError("Unknown metric: %s" % metric)

    metric = _translate_metric_name_for_div_lib(metric)

    if metric in implemented_metrics:
        func = ctx.get_action('diversity_lib', metric)
        result = func(table=table)
    else:
        func = ctx.get_action('diversity_lib', 'alpha_passthrough')
        result = func(table=table, metric=metric)

    return tuple(result)
