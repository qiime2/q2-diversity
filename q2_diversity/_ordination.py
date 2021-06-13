# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from sklearn.manifold import TSNE
import skbio.stats.ordination
import pandas as pd
import numpy as np
import umap as up


def pcoa(distance_matrix: skbio.DistanceMatrix,
         number_of_dimensions: int = None) -> skbio.OrdinationResults:
    if number_of_dimensions is None:
        # calculate full decomposition using eigh
        return skbio.stats.ordination.pcoa(distance_matrix, method='eigh',
                                           inplace=False)
    else:
        # calculate the decomposition only for the `number_of_dimensions`
        # using fast heuristic eigendecomposition (fsvd)
        return skbio.stats.ordination.pcoa(
            distance_matrix, method='fsvd',
            number_of_dimensions=number_of_dimensions,
            inplace=True)


def pcoa_biplot(pcoa: skbio.OrdinationResults,
                features: pd.DataFrame) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa_biplot(pcoa, features)


def tsne(distance_matrix: skbio.DistanceMatrix,
         number_of_dimensions: int = 2,
         perplexity: float = 25.0,
         n_iter: int = 1000,
         learning_rate: float = 200.0,
         early_exaggeration: float = 12.0) -> skbio.OrdinationResults:

    data = distance_matrix.data
    ids = distance_matrix.ids

    tsne = TSNE(number_of_dimensions, perplexity=perplexity,
                learning_rate=learning_rate,
                n_iter=n_iter,
                early_exaggeration=early_exaggeration).fit_transform(data)

    if number_of_dimensions == 2:
        number_of_dimensions = 3
        add_zeros = np.zeros((tsne.shape[0], 1), dtype=np.int64)
        tsne = np.append(tsne, add_zeros, axis=1)

    axis_labels = ["TSNE%d" % i for i in range(1, number_of_dimensions + 1)]
    eigenvalues = [0 for i in axis_labels]

    return skbio.OrdinationResults(
        short_method_name="T-SNE",
        long_method_name="t-distributed stochastic neighbor embedding",
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        proportion_explained=pd.Series(None, index=axis_labels),
        samples=pd.DataFrame(tsne, index=ids, columns=axis_labels),
    )


def umap(distance_matrix: skbio.DistanceMatrix,
         number_of_dimensions: int = 2,
         n_neighbors: int = 15,
         min_dist: float = 0.4) -> skbio.OrdinationResults:

    data = distance_matrix.data
    ids = distance_matrix.ids

    umap = up.UMAP(n_components=number_of_dimensions,
                   n_neighbors=n_neighbors,
                   min_dist=min_dist).fit_transform(data)

    if number_of_dimensions == 2:
        number_of_dimensions = 3
        add_zeros = np.zeros((umap.shape[0], 1), dtype=np.int64)
        umap = np.append(umap, add_zeros, axis=1)

    axis_labels = ["UMAP%d" % i for i in range(1, number_of_dimensions + 1)]
    eigenvalues = [0 for i in axis_labels]

    return skbio.OrdinationResults(
        short_method_name="UMAP",
        long_method_name="Uniform Manifold Approximation and Projection",
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        proportion_explained=pd.Series(None, index=axis_labels),
        samples=pd.DataFrame(umap, index=ids, columns=axis_labels),
        )
