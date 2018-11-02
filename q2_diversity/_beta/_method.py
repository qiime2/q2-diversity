# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import h5py
import biom
import skbio
import skbio.diversity
import skbio.tree
import sklearn.metrics
import unifrac
import psutil

from q2_types.feature_table import BIOMV210Format
from q2_types.tree import NewickFormat

from functools import partial

from skbio.stats.composition import clr
from scipy.spatial.distance import euclidean


def phylogenetic_metrics_dict():
    return {'unweighted_unifrac': unifrac.unweighted,
            'weighted_unifrac': unifrac.weighted_unnormalized,
            'weighted_normalized_unifrac': unifrac.weighted_normalized,
            'generalized_unifrac': unifrac.generalized}


def phylogenetic_metrics():
    return set(phylogenetic_metrics_dict())


def non_phylogenetic_metrics():
    return {'cityblock', 'euclidean', 'seuclidean', 'sqeuclidean', 'cosine',
            'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
            'braycurtis', 'mahalanobis', 'yule', 'matching', 'dice',
            'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
            'sokalsneath', 'wminkowski', 'aitchison'}


def all_metrics():
    return phylogenetic_metrics() | non_phylogenetic_metrics()


def beta_phylogenetic(table: BIOMV210Format, phylogeny: NewickFormat,
                      metric: str, n_jobs: int = 1,
                      variance_adjusted: bool = False,
                      alpha: float = None,
                      bypass_tips: bool = False) -> skbio.DistanceMatrix:

    metrics = phylogenetic_metrics_dict()
    generalized_unifrac = 'generalized_unifrac'

    if metric not in metrics:
        raise ValueError("Unknown metric: %s" % metric)

    if alpha is not None and metric != generalized_unifrac:
        raise ValueError('The alpha parameter is only allowed when the choice'
                         ' of metric is generalized_unifrac')

    # this behaviour is undefined, so let's avoid a seg fault
    cpus = psutil.cpu_count(logical=False)
    if n_jobs > cpus:
        raise ValueError('The value of n_jobs cannot exceed the number of '
                         'processors (%d) available in this system.' % cpus)

    if metric == generalized_unifrac:
        alpha = 1.0 if alpha is None else alpha
        f = partial(metrics[metric], alpha=alpha)
    else:
        f = metrics[metric]

    # avoid a full parse of the table just to check shape and IDs
    fp_path = str(table)
    h5 = h5py.File(fp_path)
    n_features, n_samples = h5.attrs['shape']
    if n_features == 0 or n_samples == 0:
        raise ValueError('The table appears to be empty.')

    dataset = h5['observation/ids']
    if isinstance(dataset[0], bytes):
        obs_ids = {i.decode('ascii') for i in dataset}
    else:
        obs_ids = {i for i in dataset}

    tmp_tree = skbio.TreeNode.read(str(phylogeny), convert_underscores=False)
    if not obs_ids.issubset({n.name for n in tmp_tree.tips()}):
        raise ValueError("Table does not appear to be completed represented "
                         "by the phylogeny.")

    # unifrac processes tables and trees should be filenames
    return f(str(table), str(phylogeny), threads=n_jobs,
             variance_adjusted=variance_adjusted, bypass_tips=bypass_tips)


def beta(table: biom.Table, metric: str,
         pseudocount: int = 1, n_jobs: int = 1)-> skbio.DistanceMatrix:

    if not (metric in non_phylogenetic_metrics()):
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().T

    def aitchison(x, y, **kwds):
        return euclidean(clr(x), clr(y))

    if metric == 'aitchison':
        counts += pseudocount
        metric = aitchison

    if table.is_empty():
        raise ValueError("The provided table object is empty")

    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        validate=True,
        pairwise_func=sklearn.metrics.pairwise_distances,
        n_jobs=n_jobs
    )
