# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import skbio.diversity


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'faith_pd'}


def non_phylogenetic_metrics():
    return {'ace', 'chao1', 'chao1_ci', 'berger_parker_d', 'brillouin_d',
            'dominance', 'doubles', 'enspie', 'esty_ci', 'fisher_alpha',
            'goods_coverage', 'heip_e', 'kempton_taylor_q', 'margalef',
            'mcintosh_d', 'mcintosh_e', 'menhinick', 'michaelis_menten_fit',
            'feature_counts', 'osd', 'pielou_e', 'robbins', 'shannon',
            'simpson', 'simpson_e', 'singles', 'strong', 'gini_index',
            'lladser_pe', 'lladser_ci'}


def alpha_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                       metric: str) -> pd.Series:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    try:
        result = skbio.diversity.alpha_diversity(metric=metric,
                                                 counts=counts,
                                                 ids=sample_ids,
                                                 otu_ids=feature_ids,
                                                 tree=phylogeny)
    except skbio.tree.MissingNodeError:
        raise skbio.tree.MissingNodeError("All ``feature_ids`` must be "
                                          "present as tip names in ``tree``.")
    result.name = metric
    return result


def alpha(table: biom.Table, metric: str) -> pd.Series:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)

    new_metric = metric if metric != 'feature_counts' else 'observed_otus'

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    result = skbio.diversity.alpha_diversity(metric=new_metric, counts=counts,
                                             ids=sample_ids)
    result.name = metric
    return result
