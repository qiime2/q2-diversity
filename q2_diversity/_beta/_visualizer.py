# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import collections
import urllib.parse
import pkg_resources
import itertools
import functools

import skbio
import skbio.diversity
import biom
import numpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests
import qiime2
import q2templates
import q2_feature_table

from ._method import beta, beta_phylogenetic, phylogenetic_metrics


TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')


def bioenv(output_dir: str, distance_matrix: skbio.DistanceMatrix,
           metadata: qiime2.Metadata) -> None:
    # convert metadata to numeric values where applicable, drop the non-numeric
    # values, and then drop samples that contain NaNs
    df = metadata.to_dataframe()
    df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))

    # filter categorical columns
    pre_filtered_cols = set(df.columns)
    df = df.select_dtypes([numpy.number]).dropna()
    filtered_categorical_cols = pre_filtered_cols - set(df.columns)

    # filter 0 variance numerical columns
    pre_filtered_cols = set(df.columns)
    df = df.loc[:, df.var() != 0]
    filtered_zero_variance_cols = pre_filtered_cols - set(df.columns)

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(df.index, strict=False)
    filtered_dm_length = distance_matrix.shape[0]

    result = skbio.stats.distance.bioenv(distance_matrix, df)
    result = q2templates.df_to_html(result)

    index = os.path.join(TEMPLATES, 'bioenv_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'filtered_categorical_cols': ', '.join(filtered_categorical_cols),
        'filtered_zero_variance_cols': ', '.join(filtered_zero_variance_cols),
        'result': result})


_beta_group_significance_fns = {'permanova': skbio.stats.distance.permanova,
                                'anosim': skbio.stats.distance.anosim}


def _get_distance_boxplot_data(distance_matrix, group_id, groupings):
    x_ticklabels = []
    all_group_distances = []

    # extract the within group distances
    within_group_distances = []
    group = groupings[group_id]
    for i, sid1 in enumerate(group):
        for sid2 in group[:i]:
            within_group_distances.append(distance_matrix[sid1, sid2])
    x_ticklabels.append('%s (n=%d)' %
                        (group_id, len(within_group_distances)))
    all_group_distances.append(within_group_distances)

    # extract between group distances for group to each other group
    for other_group_id, other_group in groupings.items():
        between_group_distances = []
        if group_id == other_group_id:
            continue
        for sid1 in group:
            for sid2 in other_group:
                between_group_distances.append(distance_matrix[sid1, sid2])
        x_ticklabels.append('%s (n=%d)' %
                            (other_group_id, len(between_group_distances)))
        all_group_distances.append(between_group_distances)
    return all_group_distances, x_ticklabels


def _get_pairwise_group_significance_stats(
        distance_matrix, group1_id, group2_id, groupings, metadata,
        beta_group_significance_fn, permutations):
    group1_group2_samples = groupings[group1_id] + groupings[group2_id]
    metadata = metadata[group1_group2_samples]
    distance_matrix = distance_matrix.filter(group1_group2_samples)
    return beta_group_significance_fn(distance_matrix, metadata,
                                      permutations=permutations)


def beta_group_significance(output_dir: str,
                            distance_matrix: skbio.DistanceMatrix,
                            metadata: qiime2.MetadataCategory,
                            method: str='permanova',
                            pairwise: bool=False,
                            permutations: int=999) -> None:
    try:
        beta_group_significance_fn = _beta_group_significance_fns[method]
    except KeyError:
        raise ValueError('Unknown group significance method %s. The available '
                         'options are %s.' %
                         (method,
                          ', '.join(_beta_group_significance_fns)))

    # Cast metadata to numeric (if applicable), which gives better sorting
    # in boxplots. Then filter any samples that are not in the distance matrix,
    # and drop samples with have no data for this metadata
    # category, including those with empty strings as values.
    metadata = pd.to_numeric(metadata.to_series(), errors='ignore')
    metadata = metadata.loc[list(distance_matrix.ids)]
    metadata = metadata.replace(r'', numpy.nan).dropna()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(metadata.index)
    filtered_dm_length = distance_matrix.shape[0]

    # Run the significance test
    result = beta_group_significance_fn(distance_matrix, metadata,
                                        permutations=permutations)

    # Generate distance boxplots
    sns.set_style("white")
    # Identify the groups, then compute the within group distances and the
    # between group distances, and generate one boxplot per group.
    # groups will be an OrderedDict mapping group id to the sample ids in that
    # group. The order is used both on the x-axis, and in the layout of the
    # boxplots in the visualization.
    groupings = collections.OrderedDict(
        [(id, list(series.index))
         for id, series in sorted(metadata.groupby(metadata))])

    for group_id in groupings:
        group_distances, x_ticklabels = \
            _get_distance_boxplot_data(distance_matrix, group_id, groupings)

        ax = sns.boxplot(data=group_distances, flierprops={
            'marker': 'o', 'markeredgecolor': 'black', 'markeredgewidth': 0.5,
            'alpha': 0.5})
        ax.set_xticklabels(x_ticklabels, rotation=90)
        ax.set_xlabel('Group')
        ax.set_ylabel('Distance')
        ax.set_title('Distances to %s' % group_id)
        # change the color of the boxes to white
        for box in ax.artists:
            box.set_facecolor('white')
        sns.despine()
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir, '%s-boxplots.png' %
                                 urllib.parse.quote_plus(str(group_id))))
        fig.savefig(os.path.join(output_dir, '%s-boxplots.pdf' %
                                 urllib.parse.quote_plus(str(group_id))))
        fig.clear()

    result_html = q2templates.df_to_html(result.to_frame())

    if pairwise:
        pairwise_results = []
        for group1_id, group2_id in itertools.combinations(groupings, 2):
            pairwise_result = \
                _get_pairwise_group_significance_stats(
                    distance_matrix=distance_matrix,
                    group1_id=group1_id,
                    group2_id=group2_id,
                    groupings=groupings,
                    metadata=metadata,
                    beta_group_significance_fn=beta_group_significance_fn,
                    permutations=permutations)
            pairwise_results.append([group1_id,
                                     group2_id,
                                     pairwise_result['sample size'],
                                     permutations,
                                     pairwise_result['test statistic'],
                                     pairwise_result['p-value']])
        columns = ['Group 1', 'Group 2', 'Sample size', 'Permutations',
                   result['test statistic name'], 'p-value']
        pairwise_results = pd.DataFrame(pairwise_results, columns=columns)
        pairwise_results.set_index(['Group 1', 'Group 2'], inplace=True)
        pairwise_results['q-value'] = multipletests(
            pairwise_results['p-value'], method='fdr_bh')[1]
        pairwise_results.sort_index(inplace=True)
        pairwise_path = os.path.join(
            output_dir, '%s-pairwise.csv' % method)
        pairwise_results.to_csv(pairwise_path)

        pairwise_results_html = q2templates.df_to_html(pairwise_results)
    else:
        pairwise_results_html = None

    index = os.path.join(
        TEMPLATES, 'beta_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'method': method,
        'groupings': groupings,
        'result': result_html,
        'pairwise_results': pairwise_results_html
    })


def beta_rarefaction(output_dir: str, table: biom.Table, metric: str,
                     sampling_depth: int, iterations: int=10,
                     phylogeny: skbio.TreeNode=None,
                     correlation_method: str='spearman',
                     color_scheme: str='BrBG') -> None:
    if metric in phylogenetic_metrics():
        if phylogeny is None:
            raise ValueError("A phylogenetic metric (%s) was requested, "
                             "but a phylogenetic tree was not provided. "
                             "Phylogeny must be provided when using a "
                             "phylogenetic diversity metric." % metric)
        beta_func = functools.partial(beta_phylogenetic, phylogeny=phylogeny)
    else:
        beta_func = beta

    distance_matrices = _get_multiple_rarefaction(
        beta_func, metric, iterations, table, sampling_depth)

    sm_df = skbio.stats.distance.pwmantel(
        distance_matrices, method=correlation_method, permutations=0,
        strict=True)
    sm = sm_df[['statistic']]  # Drop all other DF columns
    sm = sm.unstack(level=0)  # Reshape for seaborn

    test_statistics = {'spearman': "Spearman's rho", 'pearson': "Pearson's r"}
    ax = sns.heatmap(
        sm, cmap=color_scheme, vmin=-1.0, vmax=1.0, center=0.0, annot=False,
        square=True, xticklabels=False, yticklabels=False,
        cbar_kws={'ticks': [1, 0.5, 0, -0.5, -1],
                  'label': test_statistics[correlation_method]})
    ax.set(xlabel='Iteration', ylabel='Iteration',
           title='Mantel correlation between iterations')
    ax.get_figure().savefig(os.path.join(output_dir, 'heatmap.svg'))

    similarity_mtx_fp = os.path.join(output_dir,
                                     'rarefaction-iteration-correlation.tsv')
    sm_df.to_csv(similarity_mtx_fp, sep='\t')

    index_fp = os.path.join(TEMPLATES, 'beta_rarefaction_assets', 'index.html')
    q2templates.render(index_fp, output_dir)


def _get_multiple_rarefaction(beta_func, metric, iterations, table,
                              sampling_depth):
    distance_matrices = []
    for _ in range(iterations):
        rarefied_table = q2_feature_table.rarefy(table, sampling_depth)
        distance_matrix = beta_func(table=rarefied_table, metric=metric)
        distance_matrices.append(distance_matrix)
    return distance_matrices


def mantel(output_dir: str, distance_matrix1: skbio.DistanceMatrix,
           distance_matrix2: skbio.DistanceMatrix, method: str='spearman',
           permutations: int=999, intersect_ids: bool=False,
           label1: str='Distance Matrix 1',
           label2: str='Distance Matrix 2') -> None:
    test_statistics = {'spearman': 'rho', 'pearson': 'r'}
    alt_hypothesis = 'two-sided'

    # The following code to handle mismatched IDs, and subsequently filter the
    # distance matrices, is not technically necessary because skbio's mantel
    # function will raise an error on mismatches with `strict=True`, and will
    # handle intersection if `strict=False`. However, we need to handle the ID
    # matching explicitly to find *which* IDs are mismatched -- the error
    # message coming from scikit-bio doesn't describe those. We also need to
    # have the mismatched IDs to display as a warning in the viz if
    # `intersect_ids=True`. Finally, the distance matrices are explicitly
    # filtered to matching IDs only because their data are used elsewhere in
    # this function (e.g. extracting scatter plot data).

    # Find the symmetric difference between ID sets.
    ids1 = set(distance_matrix1.ids)
    ids2 = set(distance_matrix2.ids)
    mismatched_ids = ids1 ^ ids2

    if not intersect_ids and mismatched_ids:
        raise ValueError(
            "The following ID(s) are not contained in both distance "
            "matrices. Use `intersect_ids` to discard these mismatches "
            "and apply the Mantel test to only those IDs that are found "
            "in both distance matrices.\n\n%s"
            % ', '.join(sorted(mismatched_ids)))

    if mismatched_ids:
        matched_ids = ids1 & ids2
        # Run in `strict` mode because the matches should all be found in both
        # matrices.
        distance_matrix1 = distance_matrix1.filter(matched_ids, strict=True)
        distance_matrix2 = distance_matrix2.filter(matched_ids, strict=True)

    # Run in `strict` mode because all IDs should be matched at this point.
    r, p, sample_size = skbio.stats.distance.mantel(
            distance_matrix1, distance_matrix2, method=method,
            permutations=permutations, alternative=alt_hypothesis,
            strict=True)

    result = pd.Series([method.title(), sample_size, permutations,
                       alt_hypothesis, r, p],
                       index=['Method', 'Sample size', 'Permutations',
                              'Alternative hypothesis',
                              '%s %s' % (method.title(),
                                         test_statistics[method]),
                              'p-value'],
                       name='Mantel test results')
    table_html = q2templates.df_to_html(result.to_frame())

    # We know the distance matrices have matching ID sets at this point, so we
    # can safely generate all pairs of IDs using one of the matrices' ID sets
    # (it doesn't matter which one).
    scatter_data = []
    for id1, id2 in itertools.combinations(distance_matrix1.ids, 2):
        scatter_data.append((distance_matrix1[id1, id2],
                             distance_matrix2[id1, id2]))

    plt.figure()
    x = 'Pairwise Distance (%s)' % label1
    y = 'Pairwise Distance (%s)' % label2
    scatter_data = pd.DataFrame(scatter_data, columns=[x, y])
    sns.regplot(x=x, y=y, data=scatter_data, fit_reg=False)
    plt.savefig(os.path.join(output_dir, 'mantel-scatter.svg'))

    context = {
        'table': table_html,
        'sample_size': sample_size,
        'mismatched_ids': mismatched_ids
    }
    index = os.path.join(
        TEMPLATES, 'mantel_assets', 'index.html')
    q2templates.render(index, output_dir, context=context)
