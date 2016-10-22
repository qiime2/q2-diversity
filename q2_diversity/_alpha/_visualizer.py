# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
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
import qiime
from statsmodels.sandbox.stats.multicomp import multipletests
import q2templates


TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_alpha')


def alpha_group_significance(output_dir: str, alpha_diversity: pd.Series,
                             metadata: qiime.Metadata) -> None:
    metadata_df = metadata.to_dataframe()
    metadata_df = metadata_df.apply(pd.to_numeric, errors='ignore')
    pre_filtered_cols = set(metadata_df.columns)
    metadata_df = metadata_df.select_dtypes(exclude=[np.number])
    post_filtered_cols = set(metadata_df.columns)
    filtered_numeric_categories = pre_filtered_cols - post_filtered_cols

    categories = metadata_df.columns
    metric_name = alpha_diversity.name

    if len(categories) == 0:
        raise ValueError('Only numeric data is present in metadata file.')

    filenames = []
    filtered_categories = []
    for category in categories:
        metadata_category = metadata.get_category(category).to_series()
        metadata_category = metadata_category[alpha_diversity.index]
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
                    H, p = scipy.stats.mstats.kruskalwallis(groups[i],
                                                            groups[j])
                    kw_H_pairwise.append([names[j], names[i], H, p])
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
                fh.write(table.replace('\n', ''))
                fh.write("','%s', '%s');" % (quote(pairwise_fn), metric_name))
        else:
            filtered_categories.append(category)

    index = os.path.join(
        TEMPLATES, 'alpha_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'categories': [quote(fn) for fn in filenames],
        'filtered_numeric_categories': ', '.join(filtered_numeric_categories),
        'filtered_categories': ', '.join(filtered_categories)})

    shutil.copytree(os.path.join(TEMPLATES, 'dst'),
                    os.path.join(output_dir, 'dist'))


_alpha_correlation_fns = {'spearman': scipy.stats.spearmanr,
                          'pearson': scipy.stats.pearsonr}


def alpha_correlation(output_dir: str,
                      alpha_diversity: pd.Series,
                      metadata: qiime.Metadata,
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
        metadata_category = metadata_category[alpha_diversity.index]
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

    shutil.copytree(os.path.join(TEMPLATES, 'dst'),
                    os.path.join(output_dir, 'dist'))
