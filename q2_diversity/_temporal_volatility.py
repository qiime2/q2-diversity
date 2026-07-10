# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2

import numpy as np
import pandas as pd
from skbio import OrdinationResults


def pcoa_centroid_temporal_volatility(
    pcoa: OrdinationResults,
    metadata: qiime2.Metadata,
    subject_column: str,
    group_column: str = "",
    dimensions: int = 3,
    min_group_size: int = 2,
    metric: str = "median",
) -> pd.DataFrame:
    if dimensions < 1:
        raise ValueError("dimensions must be a positive integer.")
    if min_group_size < 2:
        raise ValueError(
            "min_group_size must be at least 2: a centroid distance and "
            "its standard error cannot be computed from a single sample."
        )
    if pcoa.samples.shape[1] < dimensions:
        raise ValueError(
            f"Cannot compute centroids in {dimensions} dimensions; the "
            f"PCoA only has {pcoa.samples.shape[1]} axes."
        )

    md_df = metadata.to_dataframe()
    for column in (subject_column, group_column):
        if column and column not in md_df.columns:
            raise ValueError(f'"{column}" was not found in the metadata.')

    axes = list(pcoa.samples.columns[:dimensions])
    shared_ids = pcoa.samples.index.intersection(md_df.index)
    data = pcoa.samples.loc[shared_ids, axes].copy()
    data["subject"] = md_df.loc[shared_ids, subject_column]
    data = data[~data["subject"].isnull()]

    if group_column:
        data["group"] = md_df.loc[data.index, group_column]
        data = data[~data["group"].isnull()]
        group_cols = ["subject", "group"]
    else:
        group_cols = ["subject"]

    if data.empty:
        raise ValueError(
            "No samples are shared between the PCoA and the provided "
            "metadata columns."
        )

    centroids = data.groupby(group_cols)[axes].transform("mean")
    data["distance_to_centroid"] = np.sqrt(((data[axes] - centroids) ** 2).sum(axis=1))

    grouped = data.groupby(group_cols)["distance_to_centroid"]
    n_samples = grouped.size().rename("n_samples")
    measure = grouped.agg(metric).rename("measure")
    # standard error of the per-group distances: SD / sqrt(n)
    error = (grouped.std() / np.sqrt(n_samples)).rename("error")

    volatility = pd.concat([measure, error, n_samples], axis=1)
    volatility = volatility[volatility["n_samples"] >= min_group_size]
    if volatility.empty:
        group_desc = "subject/group" if group_column else "subject"
        raise ValueError(
            f"No {group_desc} group had at least min_group_size="
            f"{min_group_size} samples; no volatility scores could be "
            "computed."
        )

    volatility = volatility.reset_index()
    if not group_column:
        volatility["group"] = "all"

    volatility["id"] = (
        volatility["subject"].astype(str) + ":" + volatility["group"].astype(str)
    )
    volatility = volatility.set_index("id", drop=False)

    return volatility[["id", "measure", "error", "group", "subject", "n_samples"]]
