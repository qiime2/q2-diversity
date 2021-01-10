# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from sklearn.manifold import TSNE
import skbio.stats.ordination
import pandas as pd
import numpy as np
import umap


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
         number_of_dimensions: int = 3) -> skbio.OrdinationResults:

    data = distance_matrix.data
    ids = distance_matrix.ids

    if number_of_dimensions == 2:
        number_of_dimensions = 3
        tsne = TSNE(2).fit_transform(data)
        add_zeros = np.zeros((tsne.shape[0], 1), dtype=np.int64)
        tsneData = np.append(tsne, add_zeros, axis=1)

    else:
        tsneData = TSNE().fit_transform(data)

    axis_labels = ["TSNE%d" % i for i in range(1, number_of_dimensions + 1)]
    eigenvalues = [0 for i in axis_labels]

    return skbio.OrdinationResults(
        short_method_name="T-SNE",
        long_method_name="t-distributed stochastic neighbor embedding",
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        samples=pd.DataFrame(tsneData, index=ids, columns=axis_labels),
    )


def uMAP(distance_matrix: skbio.DistanceMatrix,
         number_of_dimensions: int = 3) -> skbio.OrdinationResults:

    data = distance_matrix.data
    ids = distance_matrix.ids

    if number_of_dimensions == 2:
        number_of_dimensions = 3
        reducer = umap.UMAP(n_components=2).fit_transform(data)
        add_zeros = np.zeros((reducer.shape[0], 1), dtype=np.int64)
        umap_data = np.append(reducer, add_zeros, axis=1)
    else:
        reducer = umap.UMAP(n_components=number_of_dimensions)
        umap_data = reducer.fit_transform(data)

    axis_labels = ["UMAP%d" % i for i in range(1, number_of_dimensions + 1)]
    eigenvalues = [0 for i in axis_labels]

    return skbio.OrdinationResults(
        short_method_name="UMAP",
        long_method_name="Uniform Manifold Approximation and Projection",
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        samples=pd.DataFrame(umap_data, index=ids, columns=axis_labels),
    )
