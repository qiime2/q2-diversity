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
import umap


def pcoa(
    distance_matrix: skbio.DistanceMatrix, number_of_dimensions: int = None
) -> skbio.OrdinationResults:
    if number_of_dimensions is None:
        # calculate full decomposition using eigh
        return skbio.stats.ordination.pcoa(
            distance_matrix, method="eigh", inplace=False
        )
    else:
        # calculate the decomposition only for the `number_of_dimensions`
        # using fast heuristic eigendecomposition (fsvd)
        return skbio.stats.ordination.pcoa(
            distance_matrix,
            method="fsvd",
            number_of_dimensions=number_of_dimensions,
            inplace=True,
        )


def pcoa_biplot(
    pcoa: skbio.OrdinationResults, features: pd.DataFrame
) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa_biplot(pcoa, features)


def tsne(
    distance_matrix: skbio.DistanceMatrix, **kwargs
) -> skbio.OrdinationResults:

    p = {
        "n_components": 3,
        "perplexity": 30.0,
        "early_exaggeration": 12.0,
        "learning_rate": 200.0,
        "n_iter": 1000,
        "n_iter_without_progress": 300,
        "min_grad_norm": 1e-07,
        "metric": "euclidean",
        "init": "random",
        "verbose": 0,
        "random_state": None,
        "method": "barnes_hut",
        "angle": 0.5,
        "n_jobs": None,
    }

    for key, value in kwargs.items():
        if key not in p:
            continue
            p[key] = value

    data = distance_matrix.data
    ids = distance_matrix.ids
    long_method_name = "t-distributed stochastic neighbor embedding"

    tsneData = TSNE(
        p["n_components"],
        p["perplexity"],
        p["early_exaggeration"],
        p["learning_rate"],
        p["n_iter"],
        p["n_iter_without_progress"],
        p["min_grad_norm"],
        p["metric"],
        p["init"],
        p["verbose"],
        p["random_state"],
        p["method"],
        p["angle"],
        p["n_jobs"],
    ).fit_transform(data)

    axis_labels = ["TSNE%d" % i for i in range(1, p["n_components"] + 1)]
    eigenvalues = [0 for i in axis_labels]

    return skbio.OrdinationResults(
        short_method_name="T-SNE",
        long_method_name=long_method_name,
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        samples=pd.DataFrame(tsneData, index=ids, columns=axis_labels),
    )


def uMAP(
    distance_matrix: skbio.DistanceMatrix, n_components=3
) -> skbio.OrdinationResults:

    reducer = umap.UMAP(n_components=n_components)

    data = distance_matrix.data
    ids = distance_matrix.ids
    umap_data = reducer.fit_transform(data)
    axis_labels = ["UMAP%d" % i for i in range(1, n_components + 1)]
    eigenvalues = [0 for i in axis_labels]
    long_method_name = "Uniform Manifold Approximation and Projection"

    return skbio.OrdinationResults(
        short_method_name="UMAP",
        long_method_name=long_method_name,
        eigvals=pd.Series(eigenvalues, index=axis_labels),
        samples=pd.DataFrame(umap_data, index=ids, columns=axis_labels),
    )
