# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def core_metrics(ctx, table, sampling_depth, metadata, with_replacement=False,
                 n_jobs=1):
    rarefy = ctx.get_action('feature_table', 'rarefy')
    observed_features = ctx.get_action('diversity_lib', 'observed_features')
    pielou_e = ctx.get_action('diversity_lib', 'pielou_evenness')
    shannon = ctx.get_action('diversity_lib', 'shannon_entropy')
    braycurtis = ctx.get_action('diversity_lib', 'bray_curtis')
    jaccard = ctx.get_action('diversity_lib', 'jaccard')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')

    results = []
    rarefied_table, = rarefy(table=table, sampling_depth=sampling_depth,
                             with_replacement=with_replacement)
    results.append(rarefied_table)

    for metric in (observed_features, shannon, pielou_e):
        results += metric(table=rarefied_table)

    dms = []
    for metric in (jaccard, braycurtis):
        beta_results = metric(table=rarefied_table, n_jobs=n_jobs)
        results += beta_results
        dms += beta_results

    pcoas = []
    for dm in dms:
        pcoa_results = pcoa(distance_matrix=dm)
        results += pcoa_results
        pcoas += pcoa_results

    for pcoa in pcoas:
        results += emperor_plot(pcoa=pcoa, metadata=metadata)

    return tuple(results)


def core_metrics_phylogenetic(ctx, table, phylogeny, sampling_depth, metadata,
                              n_jobs_or_threads=1):
    faith_pd = ctx.get_action('diversity_lib', 'faith_pd')
    unweighted_unifrac = ctx.get_action('diversity_lib', 'unweighted_unifrac')
    weighted_unifrac = ctx.get_action(
            'diversity_lib',
            'weighted_unifrac')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')
    core_metrics = ctx.get_action('diversity', 'core_metrics')

    cr = core_metrics(table=table, sampling_depth=sampling_depth,
                      metadata=metadata, n_jobs=n_jobs_or_threads)

    faith_pd_vector, = faith_pd(table=cr.rarefied_table,
                                phylogeny=phylogeny)

    dms = []
    dms += unweighted_unifrac(table=cr.rarefied_table, phylogeny=phylogeny,
                              threads=n_jobs_or_threads)
    dms += weighted_unifrac(table=cr.rarefied_table,
                            phylogeny=phylogeny,
                            threads=n_jobs_or_threads)

    pcoas = []
    for dm in dms:
        pcoas += pcoa(distance_matrix=dm)

    plots = []
    for pcoa in pcoas:
        plots += emperor_plot(pcoa=pcoa, metadata=metadata)

    return (
        cr.rarefied_table, faith_pd_vector, cr.observed_features_vector,
        cr.shannon_vector, cr.evenness_vector, *dms,
        cr.jaccard_distance_matrix, cr.bray_curtis_distance_matrix,
        *pcoas, cr.jaccard_pcoa_results, cr.bray_curtis_pcoa_results,
        *plots, cr.jaccard_emperor, cr.bray_curtis_emperor)
