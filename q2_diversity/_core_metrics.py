# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def core_metrics(ctx, table, sampling_depth, metadata, n_jobs=1):
    rarefy = ctx.get_action('feature_table', 'rarefy')
    alpha = ctx.get_action('diversity', 'alpha')
    beta = ctx.get_action('diversity', 'beta')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')

    results = []
    rarefied_table, = rarefy(table=table, sampling_depth=sampling_depth)
    results.append(rarefied_table)

    for metric in 'observed_otus', 'shannon', 'pielou_e':
        results += alpha(table=rarefied_table, metric=metric)

    dms = []
    for metric in 'jaccard', 'braycurtis':
        beta_results = beta(table=rarefied_table, metric=metric, n_jobs=n_jobs)
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
                              n_jobs=1):
    alpha_phylogenetic = ctx.get_action('diversity', 'alpha_phylogenetic')
    beta_phylogenetic = ctx.get_action('diversity', 'beta_phylogenetic')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')
    core_metrics = ctx.get_action('diversity', 'core_metrics')

    cr = core_metrics(table=table, sampling_depth=sampling_depth,
                      metadata=metadata, n_jobs=n_jobs)

    faith_pd_vector, = alpha_phylogenetic(table=cr.rarefied_table,
                                          phylogeny=phylogeny,
                                          metric='faith_pd')

    dms = []
    for metric in 'unweighted_unifrac', 'weighted_unifrac':
        dms += beta_phylogenetic(table=cr.rarefied_table, phylogeny=phylogeny,
                                 metric=metric, n_jobs=n_jobs)

    pcoas = []
    for dm in dms:
        pcoas += pcoa(distance_matrix=dm)

    plots = []
    for pcoa in pcoas:
        plots += emperor_plot(pcoa=pcoa, metadata=metadata)

    return (
        cr.rarefied_table, faith_pd_vector, cr.observed_otus_vector,
        cr.shannon_vector, cr.evenness_vector, *dms,
        cr.jaccard_distance_matrix, cr.bray_curtis_distance_matrix,
        *pcoas, cr.jaccard_pcoa_results, cr.bray_curtis_pcoa_results,
        *plots, cr.jaccard_emperor, cr.bray_curtis_emperor)
