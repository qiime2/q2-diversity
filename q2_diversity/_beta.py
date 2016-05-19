# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import skbio
import skbio.diversity


def beta_diversity(metric: str, feature_table: biom.Table,
                   phylogeny: skbio.TreeNode=None) -> skbio.DistanceMatrix:
    counts = feature_table.matrix_data.toarray().astype(int).T
    sample_ids = feature_table.ids(axis='sample')

    if metric in ('unweighted_unifrac', 'weighted_unifrac'):
        feature_ids = feature_table.ids(axis='observation')
        return skbio.diversity.beta_diversity(
            metric=metric,
            counts=counts,
            ids=sample_ids,
            otu_ids=feature_ids,
            tree=phylogeny
        )
    else:
        if phylogeny is not None:
            raise TypeError(
                "Phylogeny was provided for non-phylogenetic metric %r"
                % metric)
        return skbio.diversity.beta_diversity(
            metric=metric,
            counts=counts,
            ids=sample_ids
        )
