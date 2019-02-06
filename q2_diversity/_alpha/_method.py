# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import numpy as np
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
            'observed_otus', 'osd', 'pielou_e', 'robbins', 'shannon',
            'simpson', 'simpson_e', 'singles', 'strong', 'gini_index',
            'lladser_pe', 'lladser_ci'}


def _batch_table(table, batchsize=1000):
    # always have at least 1 partition
    n_partitions = (len(table.ids()) / batchsize) + 1
    for id_split in np.array_split(table.ids(), n_partitions):
        subset = table.filter(set(id_split), inplace=False).remove_empty()

        counts = subset.matrix_data.toarray().astype(int).T
        feature_ids = subset.ids(axis='observation')

        yield (counts, id_split, feature_ids)


def alpha_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                       metric: str) -> pd.Series:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    results = []
    for counts, sample_ids, feature_ids in _batch_table(table):
        try:
            result = skbio.diversity.alpha_diversity(metric=metric,
                                                     counts=counts,
                                                     ids=sample_ids,
                                                     otu_ids=feature_ids,
                                                     tree=phylogeny)
        except skbio.tree.MissingNodeError as e:
            message = str(e).replace('otu_ids', 'feature_ids')
            message = message.replace('tree', 'phylogeny')
            raise skbio.tree.MissingNodeError(message)

        result.name = metric
        results.append(result)

    return pd.concat(results)


def alpha(table: biom.Table, metric: str) -> pd.Series:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)
    if table.is_empty():
        raise ValueError("The provided table object is empty")

    results = []
    for counts, sample_ids, _ in _batch_table(table):
        result = skbio.diversity.alpha_diversity(metric=metric, counts=counts,
                                                 ids=sample_ids)
        result.name = metric
        results.append(result)

    return pd.concat(results)
