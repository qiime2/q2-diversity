# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio

from q2_diversity import pcoa


class PCoATests(unittest.TestCase):

    def test_pcoa(self):
        dm = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                   [0.3333333, 0.0000000, 0.4285714],
                                   [0.6666667, 0.4285714, 0.0000000]],
                                  ids=['S1', 'S2', 'S3'])
        actual = pcoa(dm)
        expected = skbio.stats.ordination.pcoa(dm)
        skbio.util.assert_ordination_results_equal(
            actual, expected, ignore_directionality=True)
