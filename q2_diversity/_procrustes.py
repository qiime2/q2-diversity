# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio import OrdinationResults
from scipy.spatial import procrustes


def procrustes_analysis(reference: OrdinationResults, other: OrdinationResults,
                        dimensions: int = 5) -> (OrdinationResults,
                                                 OrdinationResults):

    if reference.samples.shape != other.samples.shape:
        raise ValueError('The matrices cannot be fitted unless they have the '
                         'same dimensions')

    if reference.samples.shape[1] < dimensions:
        raise ValueError('Cannot fit fewer dimensions than available')

    # fail if there are any elements in the symmetric difference
    if not (reference.samples.index ^ other.samples.index).empty:
        raise ValueError('The ordinations represent two different sets of '
                         'samples')

    # make the matrices be comparable
    other.samples = other.samples.reindex(index=reference.samples.index)

    mtx1, mtx2, _ = procrustes(reference.samples.values[:, :dimensions],
                               other.samples.values[:, :dimensions])

    axes = reference.samples.columns[:dimensions]
    samples1 = pd.DataFrame(data=mtx1,
                            index=reference.samples.index.copy(),
                            columns=axes.copy())
    samples2 = pd.DataFrame(data=mtx2,
                            index=reference.samples.index.copy(),
                            columns=axes.copy())

    out1 = OrdinationResults(
            short_method_name=reference.short_method_name,
            long_method_name=reference.long_method_name,
            eigvals=reference.eigvals[:dimensions].copy(),
            samples=samples1,
            features=reference.features,
            biplot_scores=reference.biplot_scores,
            sample_constraints=reference.sample_constraints,
            proportion_explained=reference.proportion_explained[:dimensions]
            .copy())
    out2 = OrdinationResults(
            short_method_name=other.short_method_name,
            long_method_name=other.long_method_name,
            eigvals=other.eigvals[:dimensions].copy(),
            samples=samples2,
            features=other.features,
            biplot_scores=other.biplot_scores,
            sample_constraints=other.sample_constraints,
            proportion_explained=other.proportion_explained[:dimensions]
            .copy())
    return out1, out2
