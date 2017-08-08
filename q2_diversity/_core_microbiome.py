# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources

import biom
import qiime2
import q2templates
import numpy as np
import pandas as pd
import seaborn as sns


TEMPLATES = pkg_resources.resource_filename('q2_diversity',
                                            'core_microbiome_assets')


def core_microbiome(output_dir, table: biom.Table,
                    min_fraction: float=0.5,
                    max_fraction: float=1.0,
                    steps: int=11,
                    where: str=None,
                    metadata: qiime2.Metadata=None) -> None:

    if max_fraction < min_fraction:
        raise ValueError('min_fraction (%f) parameter must be less than '
                         'max_fraction (%f) parameter.' %
                         (min_fraction, max_fraction))

    index_fp = os.path.join(TEMPLATES, 'index.html')
    context = {}

    if where is not None:
        if metadata is None:
            raise ValueError("Metadata must be provided if 'where' is "
                             "specified.")
        else:
            context['filtered'] = True
            context['samples_pre_filter'] = table.shape[1]
            context['features_pre_filter'] = table.shape[0]
            context['filter'] = where

            sample_ids = metadata.ids(where=where)
            if len(sample_ids) == 0:
                raise ValueError('Where clause (%s) matched no samples: '
                                 % where)
            ids_to_keep = sample_ids & set(table.ids(axis='sample'))
            table.filter(ids_to_keep=ids_to_keep, axis='sample')
            # adapted from biom-format 2.1.6 Table.remove_empty()
            # this should be replaced with table.remove_empty() when we
            # depend on biom-format 2.1.6.
            table.filter(lambda v, i, md: (v > 0).sum(), axis='observation')
            context['samples_post_filter'] = table.shape[1]
            context['features_post_filter'] = table.shape[0]
    else:
        context['filtered'] = False
        context['samples'] = table.shape[1]
        context['features'] = table.shape[0]

    fractions = np.linspace(min_fraction, max_fraction, steps)
    data = []
    file_links = []
    for fraction in fractions:
        core_features = _get_core_features(table, fraction)
        core_feature_count = len(core_features)
        data.append([fraction, core_feature_count])

        if core_feature_count > 0:
            core_feature_fn = 'core-features-%1.3f.tsv' % fraction
            core_feature_fp = os.path.join(output_dir, core_feature_fn)

            file_links.append(
                "<a href='./%s'>TSV</a>" % core_feature_fn)

            core_features.to_csv(core_feature_fp, sep='\t',
                                 index_label='Feature ID')
        else:
            file_links.append('No core features')

    df = pd.DataFrame(data, columns=['Fraction of samples', 'Feature count'])
    df['Fraction of features'] = df['Feature count'] / table.shape[0]
    df['Feature list'] = file_links

    ax = sns.regplot(x='Fraction of samples', y='Feature count',
                     data=df, fit_reg=False)
    ax.set_xbound(min(fractions), max(fractions))
    ax.set_ybound(0, max(df['Feature count']) + 1)
    ax.get_figure().savefig(os.path.join(output_dir,
                                         'core-feature-counts.svg'))

    table_html = df.to_html(
        index=False, escape=False, classes=("table table-striped table-hover"))
    table_html = table_html.replace('border="1"', 'border="0"')
    context['table_html'] = table_html

    q2templates.render(index_fp, output_dir, context=context)


def _get_filter_to_core_f(table, fraction):
    if not (0.0 <= fraction <= 1.0):
        raise ValueError("Invalid fraction passed to core filter: %1.2f is "
                         "outside of range [0,1]." % fraction)

    # determine the number of samples that must have a non-zero
    # value for an OTU to be considered part of the core
    min_count = fraction * len(table.ids(axis='sample'))

    def f(values, obs_ids, obs_md):
        return np.count_nonzero(values) >= min_count
    return f


def _seven_number_summary(a):
    # this should probably be publicly accessible throughout QIIME 2 - it's
    # also currently implemented in q2-demux summarize
    stats = pd.Series(a).describe(
        percentiles=[0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98])
    drop_cols = stats.index.isin(['std', 'mean', 'min', 'max', 'count'])
    stats = stats[~drop_cols]
    return stats


def _get_core_features(table, fraction):
    filter_f = _get_filter_to_core_f(table, fraction)
    feature_filtered_table = table.filter(
        filter_f, axis='observation', inplace=False)
    index = []
    data = []
    for values, id_, _ in feature_filtered_table.iter(axis='observation'):
        index.append(id_)
        data.append(_seven_number_summary(values))
    if len(data) > 0:
        return pd.DataFrame(data, index=index).sort_values(
                by='50%', ascending=False)
    else:
        return pd.DataFrame()
