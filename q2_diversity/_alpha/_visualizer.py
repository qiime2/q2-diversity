# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import os
import pkg_resources
import shutil
from urllib.parse import quote
import functools

import scipy
import numpy as np
import pandas as pd
import qiime2
from statsmodels.sandbox.stats.multicomp import multipletests
import q2templates
import biom
import itertools

from . import METRICS
from q2_types.tree import NewickFormat

TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_alpha')


def alpha_group_significance(output_dir: str, alpha_diversity: pd.Series,
                             metadata: qiime2.Metadata) -> None:
    # Filter metadata to only include IDs present in the alpha diversity data.
    # Also ensures every alpha diversity ID is present in the metadata.
    metadata = metadata.filter_ids(alpha_diversity.index)

    # Metadata column filtering could be done in one pass, but this visualizer
    # displays separate warnings for non-categorical columns, and categorical
    # columns that didn't satisfy the requirements of the statistics being
    # computed.
    pre_filtered_cols = set(metadata.columns)
    metadata = metadata.filter_columns(column_type='categorical')
    non_categorical_columns = pre_filtered_cols - set(metadata.columns)

    pre_filtered_cols = set(metadata.columns)
    metadata = metadata.filter_columns(
        drop_all_unique=True, drop_zero_variance=True, drop_all_missing=True)
    filtered_columns = pre_filtered_cols - set(metadata.columns)

    if len(metadata.columns) == 0:
        raise ValueError(
            "Metadata does not contain any columns that satisfy this "
            "visualizer's requirements. There must be at least one metadata "
            "column that contains categorical data, isn't empty, doesn't "
            "consist of unique values, and doesn't consist of exactly one "
            "value.")

    metric_name = alpha_diversity.name

    # save out metadata for download in viz
    alpha_diversity.index.name = 'id'
    alpha = qiime2.Metadata(alpha_diversity.to_frame())
    md = metadata.merge(alpha)
    md.save(os.path.join(output_dir, 'metadata.tsv'))

    filenames = []
    filtered_group_comparisons = []
    for column in metadata.columns:
        metadata_column = metadata.get_column(column)
        metadata_column = metadata_column.drop_missing_values()

        initial_data_length = alpha_diversity.shape[0]
        data = pd.concat([alpha_diversity, metadata_column.to_series()],
                         axis=1, join='inner')
        filtered_data_length = data.shape[0]

        names = []
        groups = []
        for name, group in data.groupby(metadata_column.name):
            names.append('%s (n=%d)' % (name, len(group)))
            groups.append(list(group[metric_name]))

        escaped_column = quote(column)
        escaped_column = escaped_column.replace('/', '%2F')
        filename = 'column-%s.jsonp' % escaped_column
        filenames.append(filename)

        # perform Kruskal-Wallis across all groups
        kw_H_all, kw_p_all = scipy.stats.mstats.kruskalwallis(*groups)

        # perform pairwise Kruskal-Wallis across all pairs of groups and
        # correct for multiple comparisons
        kw_H_pairwise = []
        for i in range(len(names)):
            for j in range(i):
                try:
                    H, p = scipy.stats.mstats.kruskalwallis(groups[i],
                                                            groups[j])
                    kw_H_pairwise.append([names[j], names[i], H, p])
                except ValueError:
                    filtered_group_comparisons.append(
                        ['%s:%s' % (column, names[i]),
                         '%s:%s' % (column, names[j])])
        kw_H_pairwise = pd.DataFrame(
            kw_H_pairwise, columns=['Group 1', 'Group 2', 'H', 'p-value'])
        kw_H_pairwise.set_index(['Group 1', 'Group 2'], inplace=True)
        kw_H_pairwise['q-value'] = multipletests(
            kw_H_pairwise['p-value'], method='fdr_bh')[1]
        kw_H_pairwise.sort_index(inplace=True)
        pairwise_fn = 'kruskal-wallis-pairwise-%s.csv' % escaped_column
        pairwise_path = os.path.join(output_dir, pairwise_fn)
        kw_H_pairwise.to_csv(pairwise_path)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            series = pd.Series(groups, index=names)

            fh.write("load_data('%s'," % column)
            series.to_json(fh, orient='split')
            fh.write(",")
            json.dump({'initial': initial_data_length,
                       'filtered': filtered_data_length}, fh)
            fh.write(",")
            json.dump({'H': kw_H_all, 'p': kw_p_all}, fh)
            fh.write(",'")
            table = q2templates.df_to_html(kw_H_pairwise)
            fh.write(table.replace('\n', '').replace("'", "\\'"))
            fh.write("','%s', '%s');" % (quote(pairwise_fn), metric_name))

    index = os.path.join(
        TEMPLATES, 'alpha_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames],
        'non_categorical_columns': ', '.join(sorted(non_categorical_columns)),
        'filtered_columns': ', '.join(sorted(filtered_columns)),
        'filtered_group_comparisons':
            '; '.join([' vs '.join(e) for e in filtered_group_comparisons])})

    shutil.copytree(
        os.path.join(TEMPLATES, 'alpha_group_significance_assets', 'dist'),
        os.path.join(output_dir, 'dist'))


_alpha_correlation_fns = {'spearman': scipy.stats.spearmanr,
                          'pearson': scipy.stats.pearsonr}


def alpha_correlation(output_dir: str,
                      alpha_diversity: pd.Series,
                      metadata: qiime2.Metadata,
                      method: str = 'spearman',
                      intersect_ids: bool = False) -> None:
    try:
        alpha_correlation_fn = _alpha_correlation_fns[method]
    except KeyError:
        raise ValueError('Unknown alpha correlation method %s. The available '
                         'options are %s.' %
                         (method, ', '.join(_alpha_correlation_fns.keys())))

    # Either filter metadata and alpha_diversity values to match each other
    # (intersect_ids = True) OR
    # Filter metadata to only include IDs present in the alpha diversity data.
    # Also ensures every alpha diversity ID is present in the metadata.
    if intersect_ids:
        ids1 = set(metadata.ids)
        ids2 = set(alpha_diversity.index)
        matched_ids = ids1 & ids2
        metadata = metadata.filter_ids(matched_ids)
        alpha_diversity = alpha_diversity[matched_ids]
    else:
        metadata = metadata.filter_ids(alpha_diversity.index)

    pre_filtered_cols = set(metadata.columns)
    metadata = metadata.filter_columns(column_type='numeric',
                                       drop_all_missing=True)
    filtered_columns = pre_filtered_cols - set(metadata.columns)

    if len(metadata.columns) == 0:
        raise ValueError(
            "Metadata contains only non-numeric or empty columns. This "
            "visualizer requires at least one numeric metadata column to "
            "execute.")

    # save out metadata for download in viz
    alpha_diversity.index.name = 'id'
    alpha = qiime2.Metadata(alpha_diversity.to_frame())
    md = metadata.merge(alpha)
    md.save(os.path.join(output_dir, 'metadata.tsv'))

    filenames = []
    for column in metadata.columns:
        metadata_column = metadata.get_column(column)
        metadata_column = metadata_column.drop_missing_values()

        # create a dataframe containing the data to be correlated, and drop
        # any samples that have no data in either column
        df = pd.concat([metadata_column.to_series(), alpha_diversity], axis=1,
                       join='inner')

        # compute correlation
        correlation_result = alpha_correlation_fn(df[metadata_column.name],
                                                  df[alpha_diversity.name])

        warning = None
        if alpha_diversity.shape[0] != df.shape[0]:
            warning = {'initial': alpha_diversity.shape[0],
                       'method': method.title(),
                       'filtered': df.shape[0]}

        escaped_column = quote(column)
        filename = 'column-%s.jsonp' % escaped_column
        filenames.append(filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            fh.write("load_data('%s'," % column)
            df.to_json(fh, orient='split')
            fh.write(",")
            json.dump(warning, fh)
            fh.write(",")
            json.dump({
                'method': method.title(),
                'testStat': '%1.4f' % correlation_result[0],
                'pVal': '%1.4f' % correlation_result[1],
                'sampleSize': df.shape[0]}, fh)
            fh.write(");")

    index = os.path.join(TEMPLATES, 'alpha_correlation_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'columns': [quote(fn) for fn in filenames],
        'filtered_columns': ', '.join(sorted(filtered_columns))})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_correlation_assets',
                                 'dist'),
                    os.path.join(output_dir, 'dist'))


def _reindex_with_metadata(column, columns, merged):
    reindexed = merged.set_index(column)
    reindexed.sort_index(axis=0, ascending=True, inplace=True)
    grouped = reindexed.groupby(level=[column])
    counts = grouped.count()
    # Removes the column name used to set the index of `merged` above
    col_diff = set(columns) - set([column])
    if col_diff:
        counts.drop(col_diff, axis=1, inplace=True, level=0)
    median_ = grouped.median()
    return median_, counts


def _compute_summary(data, id_label, counts=None):
    perc = [0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]
    describer = functools.partial(pd.DataFrame.describe, percentiles=perc)
    summary_df = data.stack(level=0)
    summary_df = summary_df.apply(describer, axis=1)
    summary_df.drop(['std', 'mean'], axis=1, inplace=True)
    if counts is not None:
        summary_df.drop('count', axis=1, inplace=True)
        stacked_counts = counts.stack(level=0)
        # There will always be at least one iteration, so we grab the first
        stacked_counts = stacked_counts[[1]]
        stacked_counts.rename(columns={1: 'count'}, inplace=True)
        summary_df = summary_df.join(stacked_counts, how='inner')
    else:
        # Reset count (this should always be one if we weren't explicitly
        # passed counts)
        summary_df['count'] = 1
    summary_df = summary_df.reset_index()
    summary_df.rename(columns={'level_0': id_label}, inplace=True)
    return summary_df


def _alpha_rarefaction_jsonp(output_dir, filename, metric, data, column):
    with open(os.path.join(output_dir, filename), 'w') as fh:
        fh.write("load_data('%s', '%s'," % (metric, column))
        data.to_json(fh, orient='split')
        fh.write(");")


def _compute_rarefaction_data(feature_table, min_depth, max_depth, steps,
                              iterations, phylogeny, metrics):
    depth_range = np.linspace(min_depth, max_depth, num=steps, dtype=int)
    iter_range = range(1, iterations + 1)

    rows = feature_table.ids(axis='sample')
    cols = pd.MultiIndex.from_product(
        [list(depth_range), list(iter_range)],
        names=['_alpha_rarefaction_depth_column_', 'iter'])
    data = {k: pd.DataFrame(np.NaN, index=rows, columns=cols)
            for k in metrics}

    with qiime2.sdk.Context() as scope:
        feature_table = scope.ctx.make_artifact(
                'FeatureTable[Frequency]', feature_table)

        if phylogeny:
            phylogeny = scope.ctx.make_artifact('Phylogeny[Rooted]', phylogeny)

        for depth, i in itertools.product(depth_range, iter_range):
            rarefy_method = scope.ctx.get_action('feature_table', 'rarefy')
            rt, = rarefy_method(feature_table, depth)

            for metric in metrics:
                if metric in (METRICS['PHYLO']['IMPL'] |
                              METRICS['PHYLO']['UNIMPL']):
                    alpha_phylo = scope.ctx.get_action('diversity',
                                                       'alpha_phylogenetic')
                    vector, = alpha_phylo(table=rt, metric=metric,
                                          phylogeny=phylogeny)
                else:
                    alpha = scope.ctx.get_action('diversity', 'alpha')
                    vector, = alpha(table=rt, metric=metric)

                vector = vector.view(pd.Series)
                data[metric][(depth, i)] = vector
        return data


def alpha_rarefaction(output_dir: str, table: biom.Table, max_depth: int,
                      phylogeny: NewickFormat = None, metrics: set = None,
                      metadata: qiime2.Metadata = None, min_depth: int = 1,
                      steps: int = 10, iterations: int = 10) -> None:

    if metrics is None:
        metrics = {'observed_features', 'shannon'}
        if phylogeny is not None:
            metrics.add('faith_pd')
    elif not metrics:
        raise ValueError('`metrics` was given an empty set.')
    else:
        phylo_overlap = ((METRICS['PHYLO']['IMPL'] |
                          METRICS['PHYLO']['UNIMPL']) &
                         metrics)
        if phylo_overlap and phylogeny is None:
            raise ValueError('Phylogenetic metric %s was requested but '
                             'phylogeny was not provided.' % phylo_overlap)

    if max_depth <= min_depth:
        raise ValueError('Provided max_depth of %d must be greater than '
                         'provided min_depth of %d.' % (max_depth, min_depth))
    possible_steps = max_depth - min_depth
    if possible_steps < steps:
        raise ValueError('Provided number of steps (%d) is greater than the '
                         'steps possible between min_depth and '
                         'max_depth (%d).' % (steps, possible_steps))
    if table.is_empty():
        raise ValueError('Provided table is empty.')
    max_frequency = max(table.sum(axis='sample'))
    if max_frequency < max_depth:
        raise ValueError('Provided max_depth of %d is greater than '
                         'the maximum sample total frequency of the '
                         'feature_table (%d).' % (max_depth, max_frequency))

    if metadata is None:
        columns, filtered_columns = set(), set()
    else:
        # Filter metadata to only include sample IDs present in the feature
        # table. Also ensures every feature table sample ID is present in the
        # metadata.
        metadata = metadata.filter_ids(table.ids(axis='sample'))

        # Drop metadata columns that aren't categorical, or consist solely of
        # missing values.
        pre_filtered_cols = set(metadata.columns)
        metadata = metadata.filter_columns(column_type='categorical',
                                           drop_all_missing=True)
        filtered_columns = pre_filtered_cols - set(metadata.columns)

        metadata_df = metadata.to_dataframe()
        if metadata_df.empty or len(metadata.columns) == 0:
            raise ValueError("All metadata filtered after dropping columns "
                             "that contained non-categorical data.")
        metadata_df.columns = pd.MultiIndex.from_tuples(
            [(c, '') for c in metadata_df.columns])
        columns = metadata_df.columns.get_level_values(0)

    data = _compute_rarefaction_data(table, min_depth, max_depth,
                                     steps, iterations, phylogeny, metrics)

    filenames = []
    for m, data in data.items():
        metric_name = quote(m)
        filename = '%s.csv' % metric_name

        if metadata is None:
            n_df = _compute_summary(data, 'sample-id')
            jsonp_filename = '%s.jsonp' % metric_name
            _alpha_rarefaction_jsonp(output_dir, jsonp_filename, metric_name,
                                     n_df, '')
            filenames.append(jsonp_filename)
        else:
            merged = data.join(metadata_df, how='left')
            for column in columns:
                column_name = quote(column)
                reindexed_df, counts = _reindex_with_metadata(column,
                                                              columns,
                                                              merged)
                c_df = _compute_summary(reindexed_df, column, counts=counts)
                jsonp_filename = "%s-%s.jsonp" % (metric_name, column_name)
                _alpha_rarefaction_jsonp(output_dir, jsonp_filename,
                                         metric_name, c_df, column)
                filenames.append(jsonp_filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            data.columns = [
                'depth-%d_iter-%d' % (t[0], t[1])
                for t in data.columns.values]
            if metadata is not None:
                data = data.join(metadata.to_dataframe(), how='left')
            data.to_csv(fh, index_label=['sample-id'])

    index = os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'index.html')
    q2templates.render(index, output_dir,
                       context={'metrics': list(metrics),
                                'filenames': [quote(f) for f in filenames],
                                'columns': list(columns),
                                'steps': steps,
                                'filtered_columns': sorted(filtered_columns)})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_rarefaction_assets',
                                 'dist'),
                    os.path.join(output_dir, 'dist'))


alpha_rarefaction_unsupported_metrics = {'osd', 'lladser_ci', 'strong',
                                         'esty_ci', 'kempton_taylor_q',
                                         'chao1_ci'}
