# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import io

import skbio
import numpy as np
import numpy.testing as npt
from biom.table import Table

from q2_diversity import beta_diversity


class BetaDiversityTests(TestCase):

    def test_non_phylogenetic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        actual = beta_diversity('braycurtis', t)
        # expected computed with scipy.spatial.distance.braycurtis
        expected = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                         [0.3333333, 0.0000000, 0.4285714],
                                         [0.6666667, 0.4285714, 0.0000000]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_phylogenetic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        actual = beta_diversity('unweighted_unifrac', t, phylogeny=tree)
        # expected computed with skbio.diversity.beta_diversity
        expected = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                         [0.25, 0.00, 0.00],
                                         [0.25, 0.00, 0.00]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_non_phylogenetic_invalid_input(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.assertRaises(TypeError, beta_diversity, 'braycurtis', t,
                          phylogeny=tree)

    def test_phylogenetic_invalid_input(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        self.assertRaises(TypeError, beta_diversity, 'unweighted_unifrac', t)

if __name__ == "__main__":
    main()
