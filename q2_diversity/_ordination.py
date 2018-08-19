# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import typing

import skbio.stats.ordination
import pandas as pd


def pcoa(distance_matrix: skbio.DistanceMatrix, method: str = 'eigh',
         number_of_dimensions: typing.Optional[int] = None,
         inplace: bool = False) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa(distance_matrix, method=method,
                                       number_of_dimensions=
                                       number_of_dimensions,
                                       inplace=inplace)


def pcoa_biplot(pcoa: skbio.OrdinationResults,
                features: pd.DataFrame) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa_biplot(pcoa, features)
