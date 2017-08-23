# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
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
import skbio
import itertools
from q2_feature_table import rarefy

from ._method import (non_phylogenetic_metrics, phylogenetic_metrics,
                      alpha, alpha_phylogenetic)


TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_alpha')


def alpha_group_significance(output_dir: str, alpha_diversity: pd.Series,
                             metadata: qiime2.Metadata) -> None:
    metadata_df = metadata.to_dataframe()
    metadata_df = metadata_df.apply(pd.to_numeric, errors='ignore')
    pre_filtered_cols = set(metadata_df.columns)
    metadata_df = metadata_df.select_dtypes(exclude=[np.number])
    post_filtered_cols = set(metadata_df.columns)
    filtered_numeric_categories = pre_filtered_cols - post_filtered_cols
    filtered_group_comparisons = []

    categories = metadata_df.columns
    metric_name = alpha_diversity.name

    if len(categories) == 0:
        raise ValueError('Only numeric data is present in metadata file.')

    filenames = []
    filtered_categories = []
    for category in categories:
        metadata_category = metadata.get_category(category).to_series()
        metadata_category = metadata_category.loc[alpha_diversity.index]
        metadata_category = metadata_category.replace(r'', np.nan).dropna()

        initial_data_length = alpha_diversity.shape[0]
        data = pd.concat([alpha_diversity, metadata_category], axis=1,
                         join='inner')
        filtered_data_length = data.shape[0]

        names = []
        groups = []
        for name, group in data.groupby(metadata_category.name):
            names.append('%s (n=%d)' % (name, len(group)))
            groups.append(list(group[alpha_diversity.name]))

        if (len(groups) > 1 and len(groups) != len(data.index)):
            escaped_category = quote(category)
            filename = 'category-%s.jsonp' % escaped_category
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
                            ['%s:%s' % (category, names[i]),
                             '%s:%s' % (category, names[j])])
            kw_H_pairwise = pd.DataFrame(
                kw_H_pairwise, columns=['Group 1', 'Group 2', 'H', 'p-value'])
            kw_H_pairwise.set_index(['Group 1', 'Group 2'], inplace=True)
            kw_H_pairwise['q-value'] = multipletests(
                kw_H_pairwise['p-value'], method='fdr_bh')[1]
            kw_H_pairwise.sort_index(inplace=True)
            pairwise_fn = 'kruskal-wallis-pairwise-%s.csv' % escaped_category
            pairwise_path = os.path.join(output_dir, pairwise_fn)
            kw_H_pairwise.to_csv(pairwise_path)

            with open(os.path.join(output_dir, filename), 'w') as fh:
                df = pd.Series(groups, index=names)

                fh.write("load_data('%s'," % category)
                df.to_json(fh, orient='split')
                fh.write(",")
                json.dump({'initial': initial_data_length,
                           'filtered': filtered_data_length}, fh)
                fh.write(",")
                json.dump({'H': kw_H_all, 'p': kw_p_all}, fh)
                fh.write(",'")
                table = kw_H_pairwise.to_html(classes="table table-striped "
                                              "table-hover")
                table = table.replace('border="1"', 'border="0"')
                fh.write(table.replace('\n', '').replace("'", "\\'"))
                fh.write("','%s', '%s');" % (quote(pairwise_fn), metric_name))
        else:
            filtered_categories.append(category)

    index = os.path.join(
        TEMPLATES, 'alpha_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'categories': [quote(fn) for fn in filenames],
        'filtered_numeric_categories': ', '.join(filtered_numeric_categories),
        'filtered_categories': ', '.join(filtered_categories),
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
                      method: str='spearman') -> None:
    try:
        alpha_correlation_fn = _alpha_correlation_fns[method]
    except KeyError:
        raise ValueError('Unknown alpha correlation method %s. The available '
                         'options are %s.' %
                         (method, ', '.join(_alpha_correlation_fns.keys())))
    metadata_df = metadata.to_dataframe()
    metadata_df = metadata_df.apply(pd.to_numeric, errors='ignore')
    pre_filtered_cols = set(metadata_df.columns)
    metadata_df = metadata_df.select_dtypes(include=[np.number])
    post_filtered_cols = set(metadata_df.columns)
    filtered_categories = pre_filtered_cols - post_filtered_cols

    categories = metadata_df.columns

    if len(categories) == 0:
        raise ValueError('Only non-numeric data is present in metadata file.')

    filenames = []
    for category in categories:
        metadata_category = metadata_df[category]
        metadata_category = metadata_category.loc[alpha_diversity.index]
        metadata_category = metadata_category.dropna()

        # create a dataframe containing the data to be correlated, and drop
        # any samples that have no data in either column
        df = pd.concat([metadata_category, alpha_diversity], axis=1,
                       join='inner')

        # compute correlation
        correlation_result = alpha_correlation_fn(df[metadata_category.name],
                                                  df[alpha_diversity.name])

        warning = None
        if alpha_diversity.shape[0] != df.shape[0]:
            warning = {'initial': alpha_diversity.shape[0],
                       'method': method.title(),
                       'filtered': df.shape[0]}

        escaped_category = quote(category)
        filename = 'category-%s.jsonp' % escaped_category
        filenames.append(filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            fh.write("load_data('%s'," % category)
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
        'categories': [quote(fn) for fn in filenames],
        'filtered_categories': ', '.join(filtered_categories)})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_correlation_assets',
                                 'dist'),
                    os.path.join(output_dir, 'dist'))


def _reindex_with_metadata(category, categories, merged):
    merged.set_index(category, inplace=True)
    merged.sort_index(axis=0, ascending=True, inplace=True)
    merged = merged.groupby(level=[category])
    counts = merged.count()
    counts.drop(categories, axis=1, inplace=True, level=0)
    sum_ = merged.sum()
    return sum_, counts


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


def _beta_rarefaction_jsonp(output_dir, filename, metric, data, category):
    with open(os.path.join(output_dir, filename), 'w') as fh:
        fh.write("load_data('%s', '%s'," % (metric, category))
        data.to_json(fh, orient='split')
        fh.write(");")


def _compute_rarefaction_data(feature_table, min_depth, max_depth, steps,
                              iterations, phylogeny, metrics):
    depth_range = np.linspace(min_depth, max_depth, num=steps, dtype=int)
    iter_range = range(1, iterations + 1)

    rows = feature_table.ids(axis='sample')
    cols = pd.MultiIndex.from_product([list(depth_range), list(iter_range)],
                                      names=['depth', 'iter'])
    data = {k: pd.DataFrame(np.NaN, index=rows, columns=cols)
            for k in metrics}

    for d, i in itertools.product(depth_range, iter_range):
        rt = rarefy(feature_table, d)
        for m in metrics:
            if m in phylogenetic_metrics():
                vector = alpha_phylogenetic(table=rt, metric=m,
                                            phylogeny=phylogeny)
            else:
                vector = alpha(table=rt, metric=m)
            data[m][(d, i)] = vector
    return data


def alpha_rarefaction(output_dir: str, table: biom.Table, max_depth: int,
                      phylogeny: skbio.TreeNode=None, metric: str=None,
                      metadata: qiime2.Metadata=None, min_depth: int=1,
                      steps: int=10, iterations: int=10) -> None:
    if steps < 2:
        raise ValueError('Provided steps of %d must greater than 1.' % steps)
    if iterations < 1:
        raise ValueError('Provided iterations of %d must greater than 0.'
                         % steps)
    if metric is None:
        metrics = ['observed_otus', 'shannon']
        if phylogeny is not None:
            metrics.append('faith_pd')
    else:
        if metric == 'faith_pd' and phylogeny is None:
            raise ValueError('Phylogenetic metric %s was requested but '
                             'phylogeny was not provided.' % metric)
        metrics = [metric]

    if max_depth <= min_depth:
        raise ValueError('Provided max_depth of %d must be greater than '
                         'provided min_depth of %d.' % (max_depth, min_depth))
    max_frequency = max(table.sum(axis='sample'))
    if max_frequency < max_depth:
        raise ValueError('Provided max_depth of %d is greater than '
                         'the maximum sample total frequency of the '
                         'feature_table (%d).' % (max_depth, max_frequency))
    filenames, categories = [], []
    data = _compute_rarefaction_data(table, min_depth, max_depth,
                                     steps, iterations, phylogeny, metrics)
    for (m, data) in data.items():
        metric_name = quote(m)
        filename = '%s.csv' % metric_name

        if metadata is None:
            n_df = _compute_summary(data, 'sample-id')
            jsonp_filename = '%s.jsonp' % metric_name
            _beta_rarefaction_jsonp(output_dir, jsonp_filename, metric_name,
                                    n_df, '')
            filenames.append(jsonp_filename)
        else:
            metadata_df = metadata.to_dataframe()
            metadata_df.columns = pd.MultiIndex.from_tuples(
                [(c, '') for c in metadata_df.columns])
            merged = data.join(metadata_df, how='left')
            categories = metadata_df.columns.get_level_values(0)
            for category in categories:
                category_name = quote(category)
                reindexed_df, counts = _reindex_with_metadata(category,
                                                              categories,
                                                              merged)
                c_df = _compute_summary(reindexed_df, category, counts=counts)
                jsonp_filename = "%s-%s.jsonp" % (metric_name, category_name)
                _beta_rarefaction_jsonp(output_dir, jsonp_filename,
                                        metric_name, c_df, category_name)
                filenames.append(jsonp_filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            data.columns = ['depth-%d_iter-%d' % (t[0], t[1])
                            for t in data.columns.values]
            if metadata is not None:
                data = data.join(metadata.to_dataframe(), how='left')
            data.to_csv(fh, index_label=['sample-id'])

    index = os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'index.html')
    q2templates.render(index, output_dir,
                       context={'metrics': list(metrics),
                                'filenames': filenames,
                                'categories': list(categories)})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_rarefaction_assets',
                                 'dist'),
                    os.path.join(output_dir, 'dist'))


alpha_rarefaction_supported_methods = ((non_phylogenetic_metrics()
                                       | phylogenetic_metrics())
                                       - {'osd', 'lladser_ci', 'strong',
                                          'esty_ci', 'kempton_taylor_q',
                                          'chao1_ci'})
