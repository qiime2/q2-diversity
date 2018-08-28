# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio.stats.ordination
import pandas as pd


def pcoa(distance_matrix: skbio.DistanceMatrix,
         number_of_dimensions: int = 0) -> skbio.OrdinationResults:
    if number_of_dimensions == 0:
        # calculate full decomposition using eigh
        method = 'eigh'
    else:
        # calculate the decomposition only for the `number_of_dimensions`
        # using fast heuristic eigendecomposition (fsvd)
        method = 'fsvd'

    return skbio.stats.ordination.pcoa(
        distance_matrix, method=method,
        number_of_dimensions=number_of_dimensions,
        inplace=False)


def pcoa_biplot(pcoa: skbio.OrdinationResults,
                features: pd.DataFrame) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa_biplot(pcoa, features)
