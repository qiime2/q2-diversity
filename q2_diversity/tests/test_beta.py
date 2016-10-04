# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io
import os
import tempfile

import skbio
import numpy as np
import numpy.testing as npt
from biom.table import Table
import pandas as pd
import qiime

from q2_diversity import beta, beta_phylogenetic, bioenv


class BetaDiversityTests(unittest.TestCase):

    def test_beta(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        actual = beta(table=t, metric='braycurtis')
        # expected computed with scipy.spatial.distance.braycurtis
        expected = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                         [0.3333333, 0.0000000, 0.4285714],
                                         [0.6666667, 0.4285714, 0.0000000]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_phylo_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            beta(table=t, metric='unweighted_unifrac')

    def test_beta_unknown_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            beta(table=t, metric='not-a-metric')

    def test_beta_phylogenetic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        actual = beta_phylogenetic(
            table=t, phylogeny=tree, metric='unweighted_unifrac')
        # expected computed with skbio.diversity.beta_diversity
        expected = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                         [0.25, 0.00, 0.00],
                                         [0.25, 0.00, 0.00]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_phylogenetic_non_phylo_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            beta_phylogenetic(table=t, phylogeny=tree, metric='braycurtis')

    def test_beta_phylogenetic_unknown_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            beta_phylogenetic(table=t, phylogeny=tree, metric='not-a-metric')


class BioenvTests(unittest.TestCase):

    def test_bioenv(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.Metadata(
            pd.DataFrame([['1.0', 'a'], ['2.0', 'b'], ['3.0', 'c']],
                         index=['sample1', 'sample2', 'sample3'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertFalse('metadata2' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())

    def test_bioenv_exclude_missing_data(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.Metadata(
            pd.DataFrame([['1.0', '2.0'], ['2.0', ''], ['3.0', '42.0']],
                         index=['sample1', 'sample2', 'sample3'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertTrue('metadata2' in open(index_fp).read())
            self.assertTrue('Warning' in open(index_fp).read())
            self.assertTrue('contained 3 samples' in open(index_fp).read())
            self.assertTrue('only 2 sample' in open(index_fp).read())


if __name__ == "__main__":
    unittest.main()
