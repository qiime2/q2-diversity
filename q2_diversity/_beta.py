# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path

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
    distance_matrix = distance_matrix.filter(df.index)
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
                         (method, ', '.join(_beta_group_significance_fns)))

    # Cast metadata to numeric (if applicable), which gives better sorting
    # in boxplots. Then drop samples with have no data for this metadata
    # category, including those with empty strings as values.
    metadata = pd.to_numeric(metadata.to_series(), errors='ignore')
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
    # Identify the groups (it would be nice if skbio had a public API for
    # this, so we'd be sure to have the same groups as as used by the test).
    grouping = metadata.groupby(metadata)
    group_distances = []
    group_ids = []
    for group_id, group in grouping:
        group_ids.append(group_id)
        d = []
        for i, sid1 in enumerate(group.index):
            for sid2 in group.index[:i]:
                d.append(distance_matrix[sid1, sid2])
        group_distances.append(d)
    group_distances = list(zip(group_ids, group_distances))
    group_distances.sort()

    # Plot the within group distances by group
    ax = sns.boxplot(data=[e[1] for e in group_distances])
    ax.set_xticklabels(['%s (n=%d)' % (e[0], len(e[1]))
                        for e in group_distances], rotation=90)
    ax.set_xlabel(metadata.name)
    ax.set_ylabel('Within group distances')
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(os.path.join(output_dir, 'boxplots.png'))
    fig.savefig(os.path.join(output_dir, 'boxplots.pdf'))

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
        fh.write('<p>\n')
        fh.write('<a href="boxplots.pdf">\n')
        fh.write(' <img src="boxplots.png">')
        fh.write(' <p>Download as PDF</p>\n')
        fh.write('</a>\n\n')
        fh.write('</body></html>')
