# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import collections
import urllib.parse

import numpy as np
import qiime
import biom
import skbio
import skbio.diversity
import numpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'unweighted_unifrac', 'weighted_unifrac'}


def non_phylogenetic_metrics():
    return {'cityblock', 'euclidean', 'seuclidean', 'sqeuclidean', 'cosine',
            'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
            'braycurtis', 'mahalanobis', 'yule', 'matching', 'dice',
            'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
            'sokalsneath', 'wminkowski'}


def beta_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                      metric: str)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        otu_ids=feature_ids,
        tree=phylogeny
    )


def beta(table: biom.Table, metric: str)-> skbio.DistanceMatrix:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids
    )


def bioenv(output_dir: str, distance_matrix: skbio.DistanceMatrix,
           metadata: qiime.Metadata) -> None:
    # convert metadata to numeric values where applicable, drop the non-numeric
    # values, and then drop samples that contain NaNs
    df = metadata.to_dataframe()
    df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    df = df.select_dtypes([numpy.number]).dropna()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(df.index, strict=False)
    filtered_dm_length = distance_matrix.shape[0]

    result = skbio.stats.distance.bioenv(distance_matrix, df)
    index_fp = os.path.join(output_dir, 'index.html')
    with open(index_fp, 'w') as fh:
        fh.write('<html><body>')
        if initial_dm_length != filtered_dm_length:
            fh.write("<b>Warning</b>: Some samples were filtered from the "
                     "input distance matrix because they were missing "
                     "metadata values.<br><b>The input distance matrix "
                     "contained %d samples but bioenv was computed on "
                     "only %d samples.</b><p>"
                     % (initial_dm_length, filtered_dm_length))
        fh.write(result.to_html())
        fh.write('</body></html>')


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


def beta_group_significance(output_dir: str,
                            distance_matrix: skbio.DistanceMatrix,
                            metadata: qiime.MetadataCategory,
                            method: str='permanova',
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
    metadata = metadata.replace(r'', np.nan).dropna()

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

        ax = sns.boxplot(data=group_distances)
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

    index_fp = os.path.join(output_dir, 'index.html')
    with open(index_fp, 'w') as fh:
        fh.write('<html><body>')
        if initial_dm_length != filtered_dm_length:
            fh.write("<b>Warning</b>: Some samples were filtered from the "
                     "input distance matrix because they were missing "
                     "metadata values.<br><b>The input distance matrix "
                     "contained %d samples but %s was computed on "
                     "only %d samples.</b><p>"
                     % (initial_dm_length, method, filtered_dm_length))
        fh.write(result.to_frame().to_html())
        for group_id in groupings:
            fh.write('<p>\n')
            fh.write('<a href="%s-boxplots.pdf">\n' %
                     urllib.parse.quote_plus(str(group_id)))
            fh.write(' <img src="%s-boxplots.png">' %
                     urllib.parse.quote_plus(str(group_id)))
            fh.write(' <p>Download as PDF</p>\n')
            fh.write('</a>\n\n')
            fh.write('</body></html>')
