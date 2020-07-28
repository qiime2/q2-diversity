# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unifrac


# TODO: remove these collections ASAP
def phylogenetic_metrics():
    return {'faith_pd'}


# must contain an entry for every metric in phylogenetic_metrics
def _phylogenetic_functions():
    return {'faith_pd': unifrac.faith_pd}


def non_phylogenetic_metrics():
    return {'ace', 'chao1', 'chao1_ci', 'berger_parker_d', 'brillouin_d',
            'dominance', 'doubles', 'enspie', 'esty_ci', 'fisher_alpha',
            'goods_coverage', 'heip_e', 'kempton_taylor_q', 'margalef',
            'mcintosh_d', 'mcintosh_e', 'menhinick', 'michaelis_menten_fit',
            'observed_features', 'osd', 'pielou_e', 'robbins', 'shannon',
            'simpson', 'simpson_e', 'singles', 'strong', 'gini_index',
            'lladser_pe', 'lladser_ci'}


def alpha_phylogenetic(ctx, table, phylogeny, metric):
    f = ctx.get_action('diversity_lib', 'alpha_phylogenetic_dispatch')
    result = f(table, phylogeny, metric)
    return tuple(result)


def alpha(ctx, table, metric, drop_undefined_samples=False):
    f = ctx.get_action('diversity_lib', 'alpha_dispatch')
    result = f(table, metric, drop_undefined_samples)
    return tuple(result)
