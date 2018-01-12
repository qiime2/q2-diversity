# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio.stats.ordination


def pcoa(distance_matrix: skbio.DistanceMatrix) -> skbio.OrdinationResults:
    return skbio.stats.ordination.pcoa(distance_matrix)
