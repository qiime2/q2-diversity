# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import os.path
import functools

import qiime2
import biom
import skbio
import seaborn as sns
import scipy
from emperor import Emperor

import q2templates

from . import METRICS
from .._ordination import pcoa

TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')


def beta_rarefaction(output_dir: str, table: biom.Table, metric: str,
                     clustering_method: str, metadata: qiime2.Metadata,
                     sampling_depth: int, iterations: int = 10,
                     phylogeny: skbio.TreeNode = None,
                     correlation_method: str = 'spearman',
                     color_scheme: str = 'BrBG') -> None:
    with qiime2.sdk.Context() as scope:
        if table.is_empty():
            raise ValueError("Input feature table is empty.")

        # Filter metadata to only include sample IDs present in the feature
        # table. Also ensures every feature table sample ID is present in the
        # metadata.
        metadata = metadata.filter_ids(table.ids(axis='sample'))

        table = qiime2.Artifact.import_data('FeatureTable[Frequency]', table)

        if metric in METRICS['PHYLO']['IMPL'] | METRICS['PHYLO']['UNIMPL']:
            if phylogeny is None:
                raise ValueError("A phylogenetic metric (%s) was requested, "
                                 "but a phylogenetic tree was not provided. "
                                 "Phylogeny must be provided when using a "
                                 "phylogenetic diversity metric." % metric)

            phylogeny = qiime2.Artifact.import_data('Phylogeny[Rooted]',
                                                    phylogeny)
            api_method = scope.ctx.get_action('diversity', 'beta_phylogenetic')
            beta_func = functools.partial(api_method, phylogeny=phylogeny)
        else:
            beta_func = scope.ctx.get_action('diversity', 'beta')

        rare_func = scope.ctx.get_action('feature-table', 'rarefy')

        distance_matrices = _get_multiple_rarefaction(
            beta_func, rare_func, metric, iterations, table, sampling_depth)

    primary = distance_matrices[0]
    support = distance_matrices[1:]

    heatmap_fig, similarity_df = _make_heatmap(
        distance_matrices, metric, correlation_method, color_scheme)
    heatmap_fig.savefig(os.path.join(output_dir, 'heatmap.svg'))
    similarity_df.to_csv(
        os.path.join(output_dir, 'rarefaction-iteration-correlation.tsv'),
        sep='\t')

    tree = _cluster_samples(primary, support, clustering_method)
    tree.write(os.path.join(output_dir,
                            'sample-clustering-%s.tre' % clustering_method))

    emperor = _jackknifed_emperor(primary, support, metadata)
    emperor_dir = os.path.join(output_dir, 'emperor')
    emperor.copy_support_files(emperor_dir)
    with open(os.path.join(emperor_dir, 'index.html'), 'w') as fh:
        fh.write(emperor.make_emperor(standalone=True))

    templates = list(map(
        lambda page: os.path.join(TEMPLATES, 'beta_rarefaction_assets', page),
        ['index.html', 'heatmap.html', 'tree.html', 'emperor.html']))

    context = {
        'metric': metric,
        'clustering_method': clustering_method,
        'tabs': [{'url': 'emperor.html',
                  'title': 'PCoA'},
                 {'url': 'heatmap.html',
                  'title': 'Heatmap'},
                 {'url': 'tree.html',
                  'title': 'Clustering'}]
    }

    q2templates.render(templates, output_dir, context=context)


def _get_multiple_rarefaction(beta_func, rare_func, metric, iterations, table,
                              sampling_depth):
    distance_matrices = []
    for _ in range(iterations):
        rarefied_table, = rare_func(table=table, sampling_depth=sampling_depth)
        distance_matrix, = beta_func(table=rarefied_table, metric=metric)
        distance_matrices.append(distance_matrix.view(skbio.DistanceMatrix))
    return distance_matrices


def _make_heatmap(distance_matrices, metric, correlation_method, color_scheme):
    test_statistics = {'spearman': "Spearman's rho", 'pearson': "Pearson's r"}

    sm_df = skbio.stats.distance.pwmantel(
        distance_matrices, method=correlation_method, permutations=0,
        strict=True)
    sm = sm_df[['statistic']]  # Drop all other DF columns
    sm = sm.unstack(level=0)  # Reshape for seaborn

    ax = sns.heatmap(
        sm, cmap=color_scheme, vmin=-1.0, vmax=1.0, center=0.0, annot=False,
        square=True, xticklabels=False, yticklabels=False,
        cbar_kws={'ticks': [1, 0.5, 0, -0.5, -1],
                  'label': test_statistics[correlation_method]})
    ax.set(xlabel='Iteration', ylabel='Iteration',
           title='%s - Mantel correlation between iterations' % metric)

    return ax.get_figure(), sm_df


def _cluster_samples(primary, support, clustering_method):
    cluster = {'nj': _nj, 'upgma': _upgma}[clustering_method]

    primary = cluster(primary)
    primary_internal_nodes = list(primary.non_tips())
    support_total = len(support)

    for n in primary_internal_nodes:
        n.support_count = 0

    for dm in support:
        _add_support_count(primary_internal_nodes, cluster(dm))

    for n in primary_internal_nodes:
        n.name = str(n.support_count / support_total)
        del n.support_count

    return primary


def _upgma(dm):
    upper_triangle = dm.condensed_form()
    linkage = scipy.cluster.hierarchy.average(upper_triangle)
    tree = skbio.TreeNode.from_linkage_matrix(linkage, dm.ids)
    tree.name = "root"  # root_at_midpoint for _nj labels the root
    return tree


def _nj(dm):
    # Negative branch lengths are strange, BUT we are clustering, not modeling
    # evolution, so it's not necessarily a problem
    nj = skbio.tree.nj(dm, disallow_negative_branch_length=False)
    return nj.root_at_midpoint()


def _add_support_count(nodes, support):
    # This isn't a fast or really good way to compute this, but it is obvious.
    for n in nodes:
        n_tips = {t.name for t in n.tips()}
        corresponding_node = support.lca(n_tips)
        if {t.name for t in corresponding_node.tips()} == n_tips:
            # This node has the same tips as the lca of the support tree with
            # the same tips, so there aren't any missing or extra nodes in the
            # support's subtree. (Though the subtree's topology may differ.)
            n.support_count += 1


def _jackknifed_emperor(primary_matrix, support_matrices, metadata):
    primary_pcoa = pcoa(primary_matrix)
    jackknifed_pcoa = list(map(pcoa, support_matrices))
    df = metadata.to_dataframe()
    return Emperor(primary_pcoa, df, jackknifed=jackknifed_pcoa, remote='.')
