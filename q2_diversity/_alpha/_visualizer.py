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

import scipy
import numpy as np
import pandas as pd
import qiime2
from statsmodels.sandbox.stats.multicomp import multipletests
import q2templates

import biom
import skbio
from itertools import product
from ._method import (non_phylogenetic_metrics, phylogenetic_metrics,
                      alpha, alpha_phylogenetic)
from q2_feature_table import rarefy

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


def _seven_number_summary(g):
    # this should probably be publicly accessible throughout QIIME 2 - it's
    # also currently implemented in q2-demux summarize
    stats = g.describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = stats.index.isin(['std', 'mean'])
    stats = stats[~drop_cols]
    return stats


def _with_metadata_df(category, metadata_df, data, iterations, depth_range):
    rows = []
    # TODO: replace below with metadata API call
    categorical_values = metadata_df[category]
    newData = data.copy()
    newData[category] = categorical_values
    groups = newData.groupby(category)
    for name, group in groups:
        unstacked = group.unstack(level='sample-id')
        unstacked = unstacked.unstack(level='iter')
        unstacked.index = unstacked.index.droplevel(1)
        unstacked[category] = name
        unstacked = unstacked.drop(unstacked.index[len(unstacked) - 1])
        for d in depth_range:
            s = unstacked.loc[d].dropna().iloc[1:iterations]
            rows.append({category: name,
                         'depth': d,
                         **_seven_number_summary(s)})
    return pd.DataFrame(rows)


def _without_metadata_df(data, depth_range):
    rows = []
    for (sample, d) in product(data.index, depth_range):
        s = data.loc[sample,d]
        # TODO: don't hard-code sample-id
        rows.append({'sample-id': sample,
                     'depth': d,
                     **_seven_number_summary(s)})
    return pd.DataFrame(rows)


def write_jsonp(output_dir, filename, metric, data, warnings, category):
    with open(os.path.join(output_dir, filename), 'w') as fh:
        fh.write("load_data('%s', '%s'," % (metric, category))
        data.to_json(fh, orient='split')
        fh.write(",")
        json.dump(list(set(warnings)), fh)
        fh.write(");")


def alpha_rarefaction(output_dir: str,
                      feature_table: biom.Table,
                      max_depth: int,
                      phylogeny: skbio.TreeNode=None,
                      metric: str=None,
                      metadata: qiime2.Metadata=None,
                      min_depth: int=1,
                      steps: int=10,
                      iterations: int=10) -> None:
    if metric is None:
        metrics = ['observed_otus', 'chao1', 'shannon']
        if phylogeny is not None:
            metrics.append('faith_pd')
    else:
        if metric == 'faith_pd' and phylogeny is None:
            raise ValueError('Phylogenetic metric was requested but '
                             'phylogeny was not provided.')
        metrics = [metric]
    
    if max_depth <= min_depth:
        raise ValueError('Provided max_depth of %d must be greater than '
                         'provided min_depth of %d.' % (max_depth, min_depth))
    max_frequency = max(feature_table.sum(axis='sample'))
    if max_frequency < max_depth:
        raise ValueError('Provided max_depth of %d is greater than '
                         'the maximum sample total frequency of the '
                         'feature_table (%d).' % (max_depth, max_frequency))
    warnings = []
    filenames = []
    categories = []    
    
    depth_range = np.linspace(min_depth, max_depth, num=steps, dtype=int)
    iter_range = range(1, iterations + 1)

    rows = feature_table.ids(axis='sample')
    cols = pd.MultiIndex.from_product([list(depth_range), list(iter_range)],
                                      names=['depth', 'iter'])
    data = {k: pd.DataFrame(np.NaN, index=rows, columns=cols)
            for k in metrics}

    for d, i in product(depth_range, iter_range):
        rt = rarefy(feature_table, d)
        for m in metrics:
            if m in phylogenetic_metrics():
                vector = alpha_phylogenetic(table=rt, metric=m,
                                            phylogeny=phylogeny)
            else:
                vector = alpha(table=rt, metric=m)
            data[m][(d, i)] = vector

    for (m, data) in data.items():
        metric_name = quote(m)
        filename = '%s.csv' % metric_name

        if metadata is None:
            jsonp_filename = '%s.jsonp' % metric_name
            n_df = _without_metadata_df(data, depth_range)
            write_jsonp(output_dir, jsonp_filename, metric_name, n_df,
                        warnings, '')
            filenames.append(jsonp_filename)
        else:
            metadata_df = metadata.to_dataframe()
            metadata_df = metadata_df.select_dtypes(exclude=[np.number])
            categories = metadata_df.columns
            for category in categories:
                category_name = quote(category)
                jsonp_filename = "%s-%s.jsonp" % (metric_name, category_name)
                c_df = _with_metadata_df(category_name, metadata_df, data,
                                      iterations, depth_range)
                write_jsonp(output_dir, jsonp_filename, metric_name,
                            c_df, warnings, category_name)
                filenames.append(jsonp_filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            data.to_csv(fh, index_label=['sample-id', 'depth'])

    index = os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'index.html')
    q2templates.render(index, output_dir,
                       context={'metrics': list(metrics),
                                'filenames': filenames,
                                'categories': list(categories)})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'dst'),
                    os.path.join(output_dir, 'dist'))


alpha_rarefaction_supported_methods = (non_phylogenetic_metrics()
                                        & phylogenetic_metrics()
                                        - {'osd', 'lladser_ci', 'strong',
                                           'esty_ci', 'kempton_taylor_q',
                                           'chao1_ci'})
