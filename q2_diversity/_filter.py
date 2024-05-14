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


# This function filters the SampleData[AlphaDiversity] table by the metadata
# only sample present in the metadata will remain
def filter_alpha_diversity_artifact(alpha_diversity: pd.Series,
                                    metadata: qiime2.Metadata,
                                    exclude_ids: bool = False) -> pd.Series:
    ids_to_keep = metadata.get_ids()
    if exclude_ids:
        ids_to_keep = set(alpha_diversity.index) - set(ids_to_keep)
    filtered_table = alpha_diversity[alpha_diversity.index.isin(ids_to_keep)]
    if filtered_table.empty:
        raise ValueError(
            "All samples were filtered out of the alpha diversity artifact.")
    return filtered_table
