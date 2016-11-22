# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import skbio
import pandas as pd

from q2_diversity import (alpha, alpha_phylogenetic, beta, beta_phylogenetic,
                          pcoa)
from q2_feature_table import rarefy


def core_metrics(table: biom.Table, phylogeny: skbio.TreeNode,
                 sampling_depth: int) -> (pd.Series,
                                          pd.Series,
                                          pd.Series,
                                          pd.Series,
                                          skbio.DistanceMatrix,
                                          skbio.DistanceMatrix,
                                          skbio.DistanceMatrix,
                                          skbio.DistanceMatrix,
                                          skbio.OrdinationResults,
                                          skbio.OrdinationResults,
                                          skbio.OrdinationResults,
                                          skbio.OrdinationResults):
    rarefied_table = rarefy(table=table, sampling_depth=sampling_depth)

    faith_pd_vector = alpha_phylogenetic(
        table=rarefied_table, phylogeny=phylogeny, metric='faith_pd')
    observed_otus_vector = alpha(table=rarefied_table, metric='observed_otus')
    shannon_vector = alpha(table=rarefied_table, metric='shannon')
    evenness_vector = alpha(table=rarefied_table, metric='pielou_e')

    unweighted_unifrac_distance_matrix = beta_phylogenetic(
        table=rarefied_table, phylogeny=phylogeny, metric='unweighted_unifrac')
    weighted_unifrac_distance_matrix = beta_phylogenetic(
        table=rarefied_table, phylogeny=phylogeny, metric='weighted_unifrac')
    jaccard_distance_matrix = beta(table=rarefied_table, metric='jaccard')
    bray_curtis_distance_matrix = beta(
        table=rarefied_table, metric='braycurtis')

    unweighted_unifrac_pcoa_results = pcoa(
        distance_matrix=unweighted_unifrac_distance_matrix)
    weighted_unifrac_pcoa_results = pcoa(
        distance_matrix=weighted_unifrac_distance_matrix)
    jaccard_pcoa_results = pcoa(distance_matrix=jaccard_distance_matrix)
    bray_curtis_pcoa_results = pcoa(
        distance_matrix=bray_curtis_distance_matrix)

    return (
        faith_pd_vector, observed_otus_vector, shannon_vector, evenness_vector,
        unweighted_unifrac_distance_matrix, weighted_unifrac_distance_matrix,
        jaccard_distance_matrix, bray_curtis_distance_matrix,
        unweighted_unifrac_pcoa_results, weighted_unifrac_pcoa_results,
        jaccard_pcoa_results, bray_curtis_pcoa_results
    )
