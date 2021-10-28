# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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
import tempfile
import subprocess

import skbio
import skbio.diversity
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.sandbox.stats.multicomp import multipletests
import qiime2
import q2templates
from natsort import natsorted
from patsy import ModelDesc


TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')


def bioenv(output_dir: str, distance_matrix: skbio.DistanceMatrix,
           metadata: qiime2.Metadata) -> None:
    # Filter metadata to only include IDs present in the distance matrix.
    # Also ensures every distance matrix ID is present in the metadata.
    metadata = metadata.filter_ids(distance_matrix.ids)

    # drop non-numeric columns and empty columns
    pre_filtered_cols = set(metadata.columns)
    metadata = metadata.filter_columns(column_type='numeric')
    non_numeric_cols = pre_filtered_cols - set(metadata.columns)

    # Drop samples that have any missing values.
    # TODO use Metadata API if more filtering is supported in the future.
    df = metadata.to_dataframe()
    df = df.dropna()
    metadata = qiime2.Metadata(df)

    # filter 0 variance numerical columns and empty columns
    pre_filtered_cols = set(metadata.columns)
    metadata = metadata.filter_columns(drop_zero_variance=True,
                                       drop_all_missing=True)
    zero_variance_cols = pre_filtered_cols - set(metadata.columns)

    df = metadata.to_dataframe()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(df.index)
    filtered_dm_length = distance_matrix.shape[0]

    result = skbio.stats.distance.bioenv(distance_matrix, df)
    result = q2templates.df_to_html(result)

    index = os.path.join(TEMPLATES, 'bioenv_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'non_numeric_cols': ', '.join(sorted(non_numeric_cols)),
        'zero_variance_cols': ', '.join(sorted(zero_variance_cols)),
        'result': result})


_beta_group_significance_fns = {'permanova': skbio.stats.distance.permanova,
                                'anosim': skbio.stats.distance.anosim,
                                'permdisp': skbio.stats.distance.permdisp}


def _get_distance_boxplot_data(distance_matrix, group_id, groupings):
    x_ticklabels = []
    all_group_distances = []

    # extract the within group distances
    within_group_distances = []
    pairs_summary = []
    group = groupings[group_id]
    for i, sid1 in enumerate(group):
        for sid2 in group[:i]:
            dist = distance_matrix[sid1, sid2]
            within_group_distances.append(dist)
            pairs_summary.append((sid1, sid2, group_id, group_id, dist))
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
                dist = distance_matrix[sid1, sid2]
                between_group_distances.append(dist)
                pairs_summary.append(
                    (sid1, sid2, group_id, other_group_id, dist))
        x_ticklabels.append('%s (n=%d)' %
                            (other_group_id, len(between_group_distances)))
        all_group_distances.append(between_group_distances)
    return all_group_distances, x_ticklabels, pairs_summary


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
                            metadata: qiime2.CategoricalMetadataColumn,
                            method: str = 'permanova',
                            pairwise: bool = False,
                            permutations: int = 999) -> None:
    try:
        beta_group_significance_fn = _beta_group_significance_fns[method]
    except KeyError:
        raise ValueError('Unknown group significance method %s. The available '
                         'options are %s.' %
                         (method,
                          ', '.join(_beta_group_significance_fns)))

    # Filter metadata to only include IDs present in the distance matrix.
    # Also ensures every distance matrix ID is present in the metadata.
    metadata = metadata.filter_ids(distance_matrix.ids)
    metadata = metadata.drop_missing_values()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata due to missing values, and keep track of how many samples
    # survived the filtering so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(metadata.ids)
    filtered_dm_length = distance_matrix.shape[0]

    metadata = metadata.to_series()

    # Run the significance test
    result = beta_group_significance_fn(distance_matrix, metadata,
                                        permutations=permutations)

    # Generate distance boxplots
    sns.set_style('white')
    # Identify the groups, then compute the within group distances and the
    # between group distances, and generate one boxplot per group.
    # groups will be an OrderedDict mapping group id to the sample ids in that
    # group. The order is used both on the x-axis, and in the layout of the
    # boxplots in the visualization.
    # TODO: update to use a grouping API and natsort API on
    # CategoricalMetadataColumn, if those become available.
    groupings = collections.OrderedDict(
        [(id, list(series.index))
         for id, series in natsorted(metadata.groupby(metadata))])

    pairs_summary = pd.DataFrame(columns=['SubjectID1', 'SubjectID2', 'Group1',
                                          'Group2', 'Distance'])
    for group_id in groupings:
        group_distances, x_ticklabels, group_pairs_summary = \
            _get_distance_boxplot_data(distance_matrix, group_id, groupings)

        group_pairs_summary = pd.DataFrame(
            group_pairs_summary, columns=['SubjectID1', 'SubjectID2',
                                          'Group1', 'Group2', 'Distance'])

        pairs_summary = pd.concat([pairs_summary, group_pairs_summary])

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
                                 urllib.parse.quote(str(group_id))))
        fig.savefig(os.path.join(output_dir, '%s-boxplots.pdf' %
                                 urllib.parse.quote(str(group_id))))
        fig.clear()

    pairs_summary.to_csv(os.path.join(output_dir, 'raw_data.tsv'), sep='\t')

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

    # repartition groupings for rendering
    group_ids = [
        # We have to DOUBLE encode this, as the file/resource name is a literal
        # URI-encoded string, we do this to prevent issues with the filesystem
        # however, as a result, our links need to escape % so that the browser
        # asks for the right escaped name (instead of the original name, which
        # doesn't exist inside the visualization).
        urllib.parse.quote(urllib.parse.quote(k))
        for k in groupings.keys()
    ]
    row_count, group_count = 3, len(group_ids)  # Start at three plots per row
    while group_count % row_count != 0:
        row_count = row_count - 1

    group_rows = [group_ids[g:g+row_count] for g in range(0, group_count,
                                                          row_count)]

    index = os.path.join(
        TEMPLATES, 'beta_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'method': method,
        'group_rows': group_rows,
        'bootstrap_group_col_size': int(12 / row_count),
        'result': result_html,
        'pairwise_results': pairwise_results_html
    })


def mantel(output_dir: str, dm1: skbio.DistanceMatrix,
           dm2: skbio.DistanceMatrix, method: str = 'spearman',
           permutations: int = 999, intersect_ids: bool = False,
           metadata: qiime2.CategoricalMetadataColumn = None,
           label1: str = 'Distance Matrix 1',
           label2: str = 'Distance Matrix 2') -> None:
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
    ids1 = set(dm1.ids)
    ids2 = set(dm2.ids)
    mismatched_ids = ids1 ^ ids2

    if not intersect_ids and mismatched_ids:
        raise ValueError(
            'The following ID(s) are not contained in both distance matrices. '
            'This sometimes occurs when mismatched files are passed. If this '
            'is not the case, you can use `intersect_ids` to discard these '
            'mismatches and apply the Mantel test to only those IDs that are '
            'found in both distance matrices.\n\n%s'
            % ', '.join(sorted(mismatched_ids)))

    if mismatched_ids:
        matched_ids = ids1 & ids2
        # Run in `strict` mode because the matches should all be found in both
        # matrices.
        dm1 = dm1.filter(matched_ids, strict=True)
        dm2 = dm2.filter(matched_ids, strict=True)

    # Run in `strict` mode because all IDs should be matched at this point.
    r, p, sample_size = skbio.stats.distance.mantel(
            dm1, dm2, method=method, permutations=permutations,
            alternative=alt_hypothesis, strict=True)

    result = pd.Series([method.title(), sample_size, permutations,
                       alt_hypothesis, r, p],
                       index=['Method', 'Sample size', 'Permutations',
                              'Alternative hypothesis',
                              '%s %s' % (method.title(),
                                         test_statistics[method]),
                              'p-value'],
                       name='Mantel test results')
    table_html = q2templates.df_to_html(result.to_frame())

    # if we have metadata, construct a lookup to easily determine the category
    # a sample is associated with
    if metadata is None:
        metadata = pd.Series(['single-group'] * len(dm1.ids),
                             index=list(dm1.ids))
    else:
        metadata = metadata.filter_ids(dm1.ids)
        metadata = metadata.drop_missing_values()
        metadata = metadata.to_series()
    metadata = metadata.to_dict()

    # We know the distance matrices have matching ID sets at this point, so we
    # can safely generate all pairs of IDs using one of the matrices' ID sets
    # (it doesn't matter which one).

    scatter_data = []
    for id1, id2 in itertools.combinations(dm1.ids, 2):
        scatter_data.append((dm1[id1, id2], dm2[id1, id2], metadata[id1]))

    plt.figure()
    x = 'Pairwise Distance (%s)' % label1
    y = 'Pairwise Distance (%s)' % label2
    hue = 'Hue from (%s)' % label1
    scatter_data = pd.DataFrame(scatter_data, columns=[x, y, hue])

    # if we do not have variation (or no metadata was provided), do a simple
    # regplot as prior, otherwise allow for depicting hue
    if len(scatter_data[hue].unique()) == 1:
        sns.regplot(x=x, y=y, data=scatter_data, fit_reg=False)
    else:
        g = sns.lmplot(x=x, y=y, data=scatter_data, hue=hue, fit_reg=False,
                       scatter_kws={"s": 1, 'alpha': 0.65})

        # based off of
        # https://stackoverflow.com/a/48694686
        for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh._sizes = [50]

    plt.savefig(os.path.join(output_dir, 'mantel-scatter.svg'))
    plt.close()

    context = {
        'table': table_html,
        'sample_size': sample_size,
        'mismatched_ids': mismatched_ids
    }
    index = os.path.join(
        TEMPLATES, 'mantel_assets', 'index.html')
    q2templates.render(index, output_dir, context=context)


def adonis(output_dir: str,
           distance_matrix: skbio.DistanceMatrix,
           metadata: qiime2.Metadata,
           formula: str,
           permutations: int = 999,
           n_jobs: int = 1) -> None:
    # Validate sample metadata is superset et cetera
    metadata_ids = set(metadata.ids)
    dm_ids = distance_matrix.ids
    _validate_metadata_is_superset(metadata_ids, set(dm_ids))
    # filter ids. ids must be in same order as dm
    filtered_md = metadata.to_dataframe().reindex(dm_ids)
    filtered_md.index.name = 'sample-id'
    metadata = qiime2.Metadata(filtered_md)

    # Validate formula
    terms = ModelDesc.from_formula(formula)
    for t in terms.rhs_termlist:
        for i in t.factors:
            column = metadata.get_column(i.name())
            if column.has_missing_values():
                raise ValueError('adonis requires metadata columns with no '
                                 'NaN values (missing values in column `%s`.)'
                                 % (column.name, ))

    # Run adonis
    results_fp = os.path.join(output_dir, 'adonis.tsv')
    with tempfile.TemporaryDirectory() as temp_dir_name:
        dm_fp = os.path.join(temp_dir_name, 'dm.tsv')
        distance_matrix.write(dm_fp)
        md_fp = os.path.join(temp_dir_name, 'md.tsv')
        metadata.save(md_fp)
        cmd = ['run_adonis.R', dm_fp, md_fp, formula, str(permutations),
               str(n_jobs), results_fp]
        _run_command(cmd)

    # Visualize results
    results = pd.read_csv(results_fp, sep='\t')
    results = q2templates.df_to_html(results)
    index = os.path.join(TEMPLATES, 'adonis_assets', 'index.html')
    q2templates.render(index, output_dir, context={'results': results})


def _validate_metadata_is_superset(metadata_ids, other_ids):
    missing_ids = other_ids.difference(metadata_ids)
    if len(missing_ids) > 0:
        raise ValueError('Missing samples in metadata: %r' % missing_ids)


# Replace this function with QIIME2 API for wrapping commands/binaries,
# pending https://github.com/qiime2/qiime2/issues/224
def _run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)
