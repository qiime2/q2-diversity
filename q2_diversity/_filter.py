# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import qiime2

import pandas as pd


def filter_distance_matrix(distance_matrix: skbio.DistanceMatrix,
                           metadata: qiime2.Metadata,
                           where: str = None,
                           exclude_ids: bool = False) -> skbio.DistanceMatrix:
    ids_to_keep = metadata.get_ids(where=where)
    if exclude_ids:
        ids_to_keep = set(distance_matrix.ids) - set(ids_to_keep)
    # NOTE: there is no guaranteed ordering to output distance matrix because
    # `ids_to_keep` is a set, and `DistanceMatrix.filter` uses its iteration
    # order.
    try:
        return distance_matrix.filter(ids_to_keep, strict=False)
    except skbio.stats.distance.DissimilarityMatrixError:
        raise ValueError(
            "All samples were filtered out of the distance matrix.")


def filter_alpha_diversity(alpha_diversity: pd.Series,
                           metadata: qiime2.Metadata,
                           where: str = None,
                           exclude_ids: bool = False) -> pd.Series:
    """
    Filters SampleData[AlphaDiversity] using `metadata`.

    Parameters
    ----------
    alpha_diversity : pd.Series
        The alpha diversity metrics, indexed by sample.
    metadata : qiime2.Metadata
        The metadata object to be used for filtering.
    where : str, optional
        A SQLite WHERE clause specifying which samples to select from the
        metadata.
    exclude_ids : bool, optional
        Whether to keep (default) or exclude the selected IDs in the metadata.

    Returns
    -------
    pd.Series
        The filtered alpha diversity values.

    """
    ids_to_keep = metadata.get_ids(where=where)
    if exclude_ids:
        ids_to_keep = set(alpha_diversity.index) - set(ids_to_keep)
    filtered_metric = alpha_diversity[alpha_diversity.index.isin(ids_to_keep)]
    if filtered_metric.empty:
        raise ValueError(
            "All samples were filtered out of the alpha diversity artifact.")
    return filtered_metric
