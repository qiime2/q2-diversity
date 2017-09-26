# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
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


_core_phylogenetic_output = (biom.Table,
                             pd.Series,
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
                             skbio.OrdinationResults)

_core_output = (biom.Table,
                pd.Series,
                pd.Series,
                pd.Series,
                skbio.DistanceMatrix,
                skbio.DistanceMatrix,
                skbio.OrdinationResults,
                skbio.OrdinationResults)


def _core_metrics(table, sampling_depth, n_jobs, phylogeny=None):
    rarefied_table = rarefy(table=table, sampling_depth=sampling_depth)

    # Non-Phylogenetic Metrics
    observed_otus_vector = alpha(table=rarefied_table, metric='observed_otus')
    shannon_vector = alpha(table=rarefied_table, metric='shannon')
    evenness_vector = alpha(table=rarefied_table, metric='pielou_e')

    jaccard_distance_matrix = beta(table=rarefied_table, metric='jaccard',
                                   n_jobs=n_jobs)
    bray_curtis_distance_matrix = beta(
        table=rarefied_table, metric='braycurtis', n_jobs=n_jobs)

    jaccard_pcoa_results = pcoa(distance_matrix=jaccard_distance_matrix)
    bray_curtis_pcoa_results = pcoa(
        distance_matrix=bray_curtis_distance_matrix)

    # Phylogenetic Metrics
    if phylogeny is not None:
        faith_pd_vector = alpha_phylogenetic(
            table=rarefied_table, phylogeny=phylogeny, metric='faith_pd')

        unweighted_unifrac_distance_matrix = beta_phylogenetic(
            table=rarefied_table, phylogeny=phylogeny,
            metric='unweighted_unifrac', n_jobs=n_jobs)
        weighted_unifrac_distance_matrix = beta_phylogenetic(
            table=rarefied_table, phylogeny=phylogeny,
            metric='weighted_unifrac')

        unweighted_unifrac_pcoa_results = pcoa(
            distance_matrix=unweighted_unifrac_distance_matrix)
        weighted_unifrac_pcoa_results = pcoa(
            distance_matrix=weighted_unifrac_distance_matrix)

        metrics = (rarefied_table, faith_pd_vector, observed_otus_vector,
                   shannon_vector, evenness_vector,
                   unweighted_unifrac_distance_matrix,
                   weighted_unifrac_distance_matrix, jaccard_distance_matrix,
                   bray_curtis_distance_matrix,
                   unweighted_unifrac_pcoa_results,
                   weighted_unifrac_pcoa_results, jaccard_pcoa_results,
                   bray_curtis_pcoa_results)
    else:
        metrics = (rarefied_table, observed_otus_vector, shannon_vector,
                   evenness_vector, jaccard_distance_matrix,
                   bray_curtis_distance_matrix, jaccard_pcoa_results,
                   bray_curtis_pcoa_results)

    return metrics


def core_metrics_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                              sampling_depth: int,
                              n_jobs: int=1) -> _core_phylogenetic_output:
    return _core_metrics(table, sampling_depth, n_jobs, phylogeny)


def core_metrics(table: biom.Table, sampling_depth: int,
                 n_jobs: int=1) -> _core_output:
    return _core_metrics(table, sampling_depth, n_jobs)
