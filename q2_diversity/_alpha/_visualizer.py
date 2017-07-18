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
from ._method import non_phylogenetic_metrics, alpha
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


def get_stats(group):
    return {'min': group.min(),
            'tf': group.quantile(q=0.25),
            'median': group.median(),
            'sf': group.quantile(q=0.75),
            'max': group.max()}


def alpha_rarefaction(output_dir: str,
                      feature_table: biom.Table,
                      phylogeny: skbio.TreeNode=None,
                      metrics: set=None,
                      metadata: qiime2.Metadata=None,
                      min_depth: int=1,
                      max_depth: int=100,
                      steps: int=10,
                      iterations: int=10) -> None:

    warnings = []

    for m in metrics:
        if m not in non_phylogenetic_metrics():
            warnings.append("Warning: requested metric %s "
                            "not a known metric." % m)

    max_depth = int(min(max_depth, feature_table.nnz))
    min_depth = int(min_depth)
    step_size = int(max((max_depth - min_depth) / steps, 1))

    depth_range = range(min_depth, max_depth, step_size)
    iter_range = range(1, iterations)

    rows = feature_table.ids()
    cols = pd.MultiIndex.from_product([list(depth_range), list(iter_range)],
                                      names=['depth', 'iter'])

    data = {k: pd.DataFrame(np.NaN, rows, cols) for k in metrics}

    for d, i in product(depth_range, iter_range):
        rt = rarefy(feature_table, d)
        for m in metrics:
            try:
                vector = alpha(table=rt, metric=m)
                data[m][(d, i)] = vector
            except Exception as e:
                warnings.append(str(e))
                pass

    filenames = []
    for (k, v) in data.items():
        metric_name = quote(k)
        filename = 'metric-%s.csv' % metric_name
        with open(os.path.join(output_dir, filename), 'w') as fh:
            # I think move some collation stats into here probably
            v = v.stack('depth')
            v.to_csv(fh, index_label=['sample-id', 'depth'])

        jsonp_filename = '%s.jsonp' % metric_name
        if metadata is None:
            filenames.append(jsonp_filename)

            # TODO: calculate five figure summary <-----

            with open(os.path.join(output_dir, jsonp_filename), 'w') as fh:
                fh.write("load_data('%s'," % metric_name)
                v.to_json(fh, orient='split')
                fh.write(",")
                json.dump(warnings, fh)
                fh.write(",")
                fh.write(");")
        else:
            metadata_df = metadata.to_dataframe()
            categories = metadata_df.columns
            for category in categories:
                metadata_category = metadata_df[category]
                metadata_category = metadata_category.loc[v.index.levels[0]]
                metadata_category = metadata_category.dropna()
                v[category] = [metadata_category[row.name[0]]
                               for _, row in v.iterrows()]
                print(v.iloc[:, 0:iterations-1])
                # vc = v.groupby([category, 'depth'])[[:, 1:iterations]]
                # vc = vc.apply(get_stats).unstack()

                # TODO: make not broken, and stick in jsonp <-----

    index = os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'index.html')
    q2templates.render(index, output_dir,
                       context={'metrics': metrics, 'filenames': filenames})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'dst'),
                    os.path.join(output_dir, 'dist'))
