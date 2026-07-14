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


def beta_dispersion(
    pcoa: OrdinationResults,
    group: qiime2.CategoricalMetadataColumn,
    dimensions: int = 3,
    min_group_size: int = 2,
    metric: str = "median",
) -> qiime2.Metadata:
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
    if metric not in ("mean", "sum", "median"):
        raise ValueError(
            f'metric must be one of "mean", "sum", or "median"; received '
            f'"{metric}".'
        )

    group_series = group.drop_missing_values().to_series()

    axes = list(pcoa.samples.columns[:dimensions])
    shared_ids = pcoa.samples.index.intersection(group_series.index)
    data = pcoa.samples.loc[shared_ids, axes].copy()
    data["group"] = group_series.loc[shared_ids]

    if data.empty:
        raise ValueError(
            "No samples are shared between the PCoA and the provided "
            "group metadata column."
        )

    centroids = data.groupby("group")[axes].transform("mean")
    data["distance_to_centroid"] = np.sqrt(((data[axes] - centroids) ** 2).sum(axis=1))

    grouped = data.groupby("group")["distance_to_centroid"]
    n_samples = grouped.size().rename("n_samples")
    measure = grouped.agg(metric).rename("measure")
    # standard error of the per-group distances: SD / sqrt(n)
    error = (grouped.std() / np.sqrt(n_samples)).rename("error")

    dispersion = pd.concat([measure, error, n_samples], axis=1)
    dispersion = dispersion[dispersion["n_samples"] >= min_group_size]
    if dispersion.empty:
        raise ValueError(
            f"No group had at least min_group_size={min_group_size} "
            "samples; no dispersion scores could be computed."
        )

    dispersion = dispersion.reset_index()
    dispersion["id"] = dispersion["group"].astype(str)
    dispersion = dispersion.set_index("id")

    return qiime2.Metadata(dispersion[["measure", "error", "group", "n_samples"]])
