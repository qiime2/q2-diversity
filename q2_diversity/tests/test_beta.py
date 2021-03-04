# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io
import os
import tempfile
import glob
import collections

import skbio
import numpy as np
import numpy.testing as npt
from biom.table import Table
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase


from qiime2 import Artifact
from q2_diversity import (bioenv, beta_group_significance, mantel)
from q2_diversity._beta._visualizer import _get_distance_boxplot_data


class BetaDiversityTests(TestPluginBase):
    # Note that some of these tests replicate the cases in biocore/unifrac
    package = 'q2_diversity.tests'

    def setUp(self):
        super().setUp()
        self.beta = self.plugin.pipelines['beta']
        self.beta_phylogenetic = self.plugin.pipelines['beta_phylogenetic']

        two_feature_table = self.get_data_path('two_feature_table.biom')
        self.two_feature_table = Artifact.import_data(
                'FeatureTable[Frequency]',
                two_feature_table)

        three_feature_tree = self.get_data_path('three_feature.tree')
        self.three_feature_tree = Artifact.import_data('Phylogeny[Rooted]',
                                                       three_feature_tree)

        crawford_table = self.get_data_path('crawford.biom')
        self.crawford_table = Artifact.import_data('FeatureTable[Frequency]',
                                                   crawford_table)
        crawford_tree = self.get_data_path('crawford.nwk')
        self.crawford_tree = Artifact.import_data('Phylogeny[Rooted]',
                                                  crawford_tree)

        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        self.t = Artifact.import_data('FeatureTable[Frequency]', t)
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.tree = Artifact.import_data('Phylogeny[Rooted]', tree)

    def test_beta(self):
        actual = self.beta(table=self.t, metric='braycurtis')
        actual = actual[0].view(skbio.DistanceMatrix)
        # expected computed with scipy.spatial.distance.braycurtis
        expected = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                         [0.3333333, 0.0000000, 0.4285714],
                                         [0.6666667, 0.4285714, 0.0000000]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_parallel_beta(self):
        parallel = self.beta(table=self.t, metric='braycurtis', n_jobs='auto')
        parallel = parallel[0].view(skbio.DistanceMatrix)
        single_thread = self.beta(table=self.t, metric='braycurtis', n_jobs=1)
        single_thread = single_thread[0].view(skbio.DistanceMatrix)
        # expected computed with scipy.spatial.distance.braycurtis
        expected = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                         [0.3333333, 0.0000000, 0.4285714],
                                         [0.6666667, 0.4285714, 0.0000000]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(parallel.ids, expected.ids)
        self.assertEqual(single_thread.ids, expected.ids)
        for id1 in parallel.ids:
            for id2 in parallel.ids:
                npt.assert_almost_equal(parallel[id1, id2], expected[id1, id2])
        for id1 in single_thread.ids:
            for id2 in single_thread.ids:
                npt.assert_almost_equal(single_thread[id1, id2],
                                        expected[id1, id2])

    def test_beta_phylo_metric(self):
        with self.assertRaisesRegex(TypeError,
                                    'received \'unweighted_unifrac\''):
            self.beta(table=self.t, metric='unweighted_unifrac')

    def test_beta_unknown_metric(self):
        with self.assertRaisesRegex(TypeError,
                                    'received \'not-a-metric\''):
            self.beta(table=self.t, metric='not-a-metric')

    def test_beta_empty_table(self):
        t = Table(np.array([]), [], [])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        with self.assertRaisesRegex(ValueError, 'empty'):
            self.beta(table=t, metric='braycurtis')

    def test_beta_phylogenetic(self):
        t = self.two_feature_table
        tree = self.three_feature_tree
        actual = self.beta_phylogenetic(
            table=t, phylogeny=tree, metric='unweighted_unifrac')

        self.assertEqual(len(actual), 1)

        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')

        # expected computed with skbio.diversity.beta_diversity
        expected = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                         [0.25, 0.00, 0.00],
                                         [0.25, 0.00, 0.00]],
                                        ids=['S1', 'S2', 'S3'])

        actual = actual[0].view(skbio.DistanceMatrix)
        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_phylogenetic_non_phylo_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)
        with self.assertRaisesRegex(TypeError, 'received \'braycurtis'):
            self.beta_phylogenetic(table=t, phylogeny=tree,
                                   metric='braycurtis')

    def test_beta_phylogenetic_unknown_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)
        with self.assertRaisesRegex(TypeError, 'received \'not-a-metric\''):
            self.beta_phylogenetic(table=t, phylogeny=tree,
                                   metric='not-a-metric')

    def test_beta_phylogenetic_empty_table(self):
        t = self.get_data_path('empty.biom')
        t = Artifact.import_data('FeatureTable[Frequency]', t)
        tree = self.get_data_path('three_feature.tree')
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        with self.assertRaisesRegex(ValueError, 'empty'):
            self.beta_phylogenetic(table=t, phylogeny=tree,
                                   metric='unweighted_unifrac')

    def test_beta_unweighted(self):
        actual = self.beta_phylogenetic(table=self.crawford_table,
                                        phylogeny=self.crawford_tree,
                                        metric='unweighted_unifrac')

        # computed with beta-phylogenetic
        data = np.array([0.71836067, 0.71317361, 0.69746044, 0.62587207,
                         0.72826674, 0.72065895, 0.72640581, 0.73606053,
                         0.70302967, 0.73407301, 0.6548042, 0.71547381,
                         0.78397813, 0.72318399, 0.76138933, 0.61041275,
                         0.62331299, 0.71848305, 0.70416337, 0.75258475,
                         0.79249029, 0.64392779, 0.70052733, 0.69832716,
                         0.77818938, 0.72959894, 0.75782689, 0.71005144,
                         0.75065046, 0.78944369, 0.63593642, 0.71283615,
                         0.58314638, 0.69200762, 0.68972056, 0.71514083])
        ids = ('10084.PC.481', '10084.PC.593', '10084.PC.356', '10084.PC.355',
               '10084.PC.354', '10084.PC.636', '10084.PC.635', '10084.PC.607',
               '10084.PC.634')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_unweighted_parallel(self):
        bt_fp = self.get_data_path('crawford.biom')
        bt = Artifact.import_data('FeatureTable[Frequency]', bt_fp)
        tree_fp = self.get_data_path('crawford.nwk')
        tree = Artifact.import_data('Phylogeny[Rooted]', tree_fp)

        actual = self.beta_phylogenetic(table=bt,
                                        phylogeny=tree,
                                        metric='unweighted_unifrac',
                                        threads=2)

        # computed with beta-phylogenetic
        data = np.array([0.71836067, 0.71317361, 0.69746044, 0.62587207,
                         0.72826674, 0.72065895, 0.72640581, 0.73606053,
                         0.70302967, 0.73407301, 0.6548042, 0.71547381,
                         0.78397813, 0.72318399, 0.76138933, 0.61041275,
                         0.62331299, 0.71848305, 0.70416337, 0.75258475,
                         0.79249029, 0.64392779, 0.70052733, 0.69832716,
                         0.77818938, 0.72959894, 0.75782689, 0.71005144,
                         0.75065046, 0.78944369, 0.63593642, 0.71283615,
                         0.58314638, 0.69200762, 0.68972056, 0.71514083])
        ids = ('10084.PC.481', '10084.PC.593', '10084.PC.356', '10084.PC.355',
               '10084.PC.354', '10084.PC.636', '10084.PC.635', '10084.PC.607',
               '10084.PC.634')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_weighted(self):
        actual = self.beta_phylogenetic(table=self.crawford_table,
                                        phylogeny=self.crawford_tree,
                                        metric='weighted_unifrac')

        # computed with beta-phylogenetic (weighted_unifrac)
        data = np.array([0.44656238, 0.23771096, 0.30489123, 0.23446002,
                         0.65723575, 0.44911772, 0.381904, 0.69144829,
                         0.39611776, 0.36568012, 0.53377975, 0.48908025,
                         0.35155196, 0.28318669, 0.57376916, 0.23395746,
                         0.24658122, 0.60271637, 0.39802552, 0.36567394,
                         0.68062701, 0.36862049, 0.48350632, 0.33024631,
                         0.33266697, 0.53464744, 0.74605075, 0.53951035,
                         0.49680733, 0.79178838, 0.37109012, 0.52629343,
                         0.22118218, 0.32400805, 0.43189708, 0.59705893])
        ids = ('10084.PC.481', '10084.PC.593', '10084.PC.356', '10084.PC.355',
               '10084.PC.354', '10084.PC.636', '10084.PC.635', '10084.PC.607',
               '10084.PC.634')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_variance_adjusted_normalized(self):
        bt_fp = self.get_data_path('vaw.biom')
        bt = Artifact.import_data('FeatureTable[Frequency]', bt_fp)
        tree_fp = self.get_data_path('vaw.nwk')
        tree = Artifact.import_data('Phylogeny[Rooted]', tree_fp)

        actual = self.beta_phylogenetic(table=bt,
                                        phylogeny=tree,
                                        metric='weighted_normalized_unifrac',
                                        variance_adjusted=True)

        data = np.array([[0.0000000, 0.4086040, 0.6240185, 0.4639481,
                          0.2857143, 0.2766318],
                         [0.4086040, 0.0000000, 0.3798594, 0.6884992,
                          0.6807616, 0.4735781],
                         [0.6240185, 0.3798594, 0.0000000, 0.7713254,
                          0.8812897, 0.5047114],
                         [0.4639481, 0.6884992, 0.7713254, 0.0000000,
                          0.6666667, 0.2709298],
                         [0.2857143, 0.6807616, 0.8812897, 0.6666667,
                          0.0000000, 0.4735991],
                         [0.2766318, 0.4735781, 0.5047114, 0.2709298,
                          0.4735991, 0.0000000]])
        ids = ('Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5',
               'Sample6')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_generalized_unifrac(self):
        bt_fp = self.get_data_path('vaw.biom')
        bt = Artifact.import_data('FeatureTable[Frequency]', bt_fp)
        tree_fp = self.get_data_path('vaw.nwk')
        tree = Artifact.import_data('Phylogeny[Rooted]', tree_fp)

        actual = self.beta_phylogenetic(table=bt,
                                        phylogeny=tree,
                                        metric='generalized_unifrac',
                                        alpha=0.5)

        data = np.array([[0.0000000, 0.4040518, 0.6285560, 0.5869439,
                          0.4082483, 0.2995673],
                         [0.4040518, 0.0000000, 0.4160597, 0.7071068,
                          0.7302479, 0.4860856],
                         [0.6285560, 0.4160597, 0.0000000, 0.8005220,
                          0.9073159, 0.5218198],
                         [0.5869439, 0.7071068, 0.8005220, 0.0000000,
                          0.4117216, 0.3485667],
                         [0.4082483, 0.7302479, 0.9073159, 0.4117216,
                          0.0000000, 0.6188282],
                         [0.2995673, 0.4860856, 0.5218198, 0.3485667,
                          0.6188282, 0.0000000]])
        ids = ('Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5',
               'Sample6')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_generalized_unifrac_no_alpha(self):
        actual = self.beta_phylogenetic(table=self.crawford_table,
                                        phylogeny=self.crawford_tree,
                                        metric='generalized_unifrac',
                                        alpha=None)

        # alpha=1 should be equal to weighted normalized UniFrac
        data = np.array([0.2821874, 0.16148405, 0.20186143, 0.1634832,
                         0.40351108, 0.29135056, 0.24790944, 0.41967404,
                         0.24642185, 0.22218489, 0.34007547, 0.27722011,
                         0.20963881, 0.16897221, 0.3217958, 0.15237816,
                         0.16899207, 0.36445044, 0.25408941, 0.23358681,
                         0.4069374, 0.24615927, 0.28573888, 0.20578184,
                         0.20742006, 0.31249151, 0.46169893, 0.35294595,
                         0.32522355, 0.48437103, 0.21534558, 0.30558908,
                         0.12091004, 0.19817777, 0.24792853, 0.34293674])
        ids = ('10084.PC.481', '10084.PC.593', '10084.PC.356', '10084.PC.355',
               '10084.PC.354', '10084.PC.636', '10084.PC.635', '10084.PC.607',
               '10084.PC.634')
        expected = skbio.DistanceMatrix(data, ids=ids)

        self.assertEqual(len(actual), 1)
        self.assertEqual(repr(actual.distance_matrix.type), 'DistanceMatrix')
        actual = actual[0].view(skbio.DistanceMatrix)

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_not_generalized_passed_alpha(self):
        with self.assertRaisesRegex(ValueError,
                                    "alpha.*only allowed.*when.*generalized"):
            self.beta_phylogenetic(table=self.crawford_table,
                                   phylogeny=self.crawford_tree,
                                   metric='unweighted_unifrac',
                                   alpha=0.5)

    def test_beta_phylogenetic_too_many_jobs(self):
        with self.assertRaises(ValueError):
            # cannot guarantee that this will always be true, but it would be
            # odd to see a machine with these many CPUs
            self.beta_phylogenetic(table=self.crawford_table,
                                   phylogeny=self.crawford_tree,
                                   metric='unweighted_unifrac', threads=11117)


class BioenvTests(TestPluginBase):
    package = 'q2_diversity.tests'

    def test_bioenv(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                [[1.0, 'a'], [2.0, 'b'], [3.0, 'c']],
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id'),
                columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())

            self.assertTrue('not numeric:' in open(index_fp).read())
            self.assertTrue('<strong>metadata2' in open(index_fp).read())

            self.assertFalse('Warning' in open(index_fp).read())

    def test_bioenv_exclude_missing_data(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                [[1.0, 2.0], [2.0, np.nan], [3.0, 42.0]],
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id'),
                columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertTrue('metadata2' in open(index_fp).read())
            self.assertTrue('Warning' in open(index_fp).read())
            self.assertTrue('contained 3 samples' in open(index_fp).read())
            self.assertTrue('2 samples' in open(index_fp).read())

    def test_bioenv_extra_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                [[1.0, 'a'], [2.0, 'b'], [3.0, 'c'], [4.0, 'd']],
                index=pd.Index(['sample1', 'sample2', 'sample3', 'sample4'],
                               name='id'),
                columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())

            self.assertTrue('not numeric:' in open(index_fp).read())
            self.assertTrue('<strong>metadata2' in open(index_fp).read())

            self.assertFalse('Warning' in open(index_fp).read())

    def test_bioenv_zero_variance_column(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                [[1.0, 2.0], [2.0, 2.0], [3.0, 2.0]],
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id'),
                columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('metadata1' in open(index_fp).read())

            self.assertTrue('no variance' in open(index_fp).read())
            self.assertTrue('<strong>metadata2' in open(index_fp).read())

            self.assertFalse('Warning' in open(index_fp).read())


class BetaGroupSignificanceTests(unittest.TestCase):

    def test_permanova(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('PERMANOVA results' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())
            self.assertFalse('Pairwise permanova' in open(index_fp).read())

    def test_anosim(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, method='anosim',
                                    permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('ANOSIM results' in open(index_fp).read())
            self.assertTrue('<td>42</td>' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())
            self.assertFalse('Pairwise anosim' in open(index_fp).read())

    def test_permanova_pairwise(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, pairwise=True)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('PERMANOVA results' in open(index_fp).read())
            self.assertTrue('Pairwise permanova' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())

    def test_anosim_pairwise(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, method='anosim',
                                    permutations=42, pairwise=True)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('ANOSIM results' in open(index_fp).read())
            self.assertTrue('<td>42</td>' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())
            self.assertTrue('Pairwise anosim' in open(index_fp).read())

    def test_alt_permutations(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('<td>42</td>' in open(index_fp).read())

    def test_invalid_method(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

        with self.assertRaises(ValueError):
            with tempfile.TemporaryDirectory() as output_dir:
                beta_group_significance(output_dir, dm, md, method='bad!')

    def test_filtered_samples_numeric_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.66, 0.66, 0.66, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3',
                                       'sample4'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['1.0', '2.0', '2.0', np.nan], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3',
                                      'sample4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('Warning' in open(index_fp).read())

    def test_filtered_samples_str_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.66, 0.66, 0.66, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3',
                                       'sample4'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b', np.nan], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3',
                                      'sample4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('Warning' in open(index_fp).read())

    def test_extra_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime2.CategoricalMetadataColumn(
            pd.Series(['a', 'b', 'b', 'c'], name='a or b',
                      index=pd.Index(['sample1', 'sample2', 'sample3',
                                      'sample4'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('<td>2</td>' in open(index_fp).read())

    def test_get_distance_boxplot_data_two_groups(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g2', ['s3', 's4', 's5'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.14, 0.15, 0.22, 0.23, 0.24]]
        exp_labels = ['g1 (n=1)', 'g2 (n=6)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

    def test_get_distance_boxplot_data_within_always_first(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g2', ['s3', 's4', 's5']), ('g1', ['s1', 's2'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.14, 0.15, 0.22, 0.23, 0.24]]
        exp_labels = ['g1 (n=1)', 'g2 (n=6)']
        exp_summary = [('s2', 's1', 'g1', 'g1', 0.12),
                       ('s1', 's3', 'g1', 'g2', 0.13),
                       ('s1', 's4', 'g1', 'g2', 0.14000000000000001),
                       ('s1', 's5', 'g1', 'g2', 0.14999999999999999),
                       ('s2', 's3', 'g1', 'g2', 0.22),
                       ('s2', 's4', 'g1', 'g2', 0.23000000000000001),
                       ('s2', 's5', 'g1', 'g2', 0.23999999999999999)]
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)
        self.assertEqual(obs[2], exp_summary)

    def test_get_distance_boxplot_data_three_groups(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g2', ['s3', 's5']), ('g3', ['s4'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.15, 0.22, 0.24], [0.14, 0.23]]
        exp_labels = ['g1 (n=1)', 'g2 (n=4)', 'g3 (n=2)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

    def test_get_distance_boxplot_data_between_order_retained(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g3', ['s4']), ('g2', ['s3', 's5'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.14, 0.23], [0.13, 0.15, 0.22, 0.24]]
        exp_labels = ['g1 (n=1)', 'g3 (n=2)', 'g2 (n=4)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)


class TestMantel(unittest.TestCase):
    def setUp(self):
        self.dm1 = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                         [0.25, 0.00, 0.00],
                                         [0.25, 0.00, 0.00]],
                                        ids=['sample1', 'sample2', 'sample3'])

        # Positive correlation with `dm1`
        self.dm2 = skbio.DistanceMatrix([[0.00, 1.00, 2.00],
                                         [1.00, 0.00, 1.00],
                                         [2.00, 1.00, 0.00]],
                                        ids=['sample1', 'sample2', 'sample3'])

        # Perfect negative correlation with `dm1`
        self.dm3 = skbio.DistanceMatrix([[0.00, 0.00, 0.00],
                                         [0.00, 0.00, 0.25],
                                         [0.00, 0.25, 0.00]],
                                        ids=['sample1', 'sample2', 'sample3'])

        self.dm2_reordered = skbio.DistanceMatrix(
            [[0.00, 2.00, 1.00],
             [2.00, 0.00, 1.00],
             [1.00, 1.00, 0.00]],
            ids=['sample3', 'sample1', 'sample2'])

        self.mismatched_dm = skbio.DistanceMatrix(
            [[0.0, 0.0, 0.0, 0.0, 0.0],
             [0.0, 0.0, 1.0, 2.0, 3.0],
             [0.0, 1.0, 0.0, 1.0, 2.0],
             [0.0, 2.0, 1.0, 0.0, 1.0],
             [0.0, 3.0, 2.0, 1.0, 0.0]],
            ids=['foo', 'sample1', 'sample2', 'sample3', 'x'])

        self.output_dir_obj = tempfile.TemporaryDirectory(
                prefix='q2-diversity-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertBasicVizValidity(self, viz_dir, sample_size, method='spearman',
                               permutations=999, label1='Distance Matrix 1',
                               label2='Distance Matrix 2',
                               mismatched_ids=None, exp_test_stat=None,
                               exp_p_value=None):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))

        with open(index_fp, 'r') as fh:
            index_contents = fh.read()

        self.assertIn('Mantel test results', index_contents)
        self.assertIn('<td>%d</td>' % sample_size, index_contents)
        self.assertIn('<td>%d</td>' % permutations, index_contents)

        method_labels = {'spearman': "Spearman rho", 'pearson': "Pearson r"}
        self.assertIn(method_labels[method], index_contents)

        if mismatched_ids is None:
            self.assertNotIn('Warning:', index_contents)
        else:
            self.assertIn('Warning', index_contents)
            self.assertIn('%d ID(s)' % len(mismatched_ids), index_contents)
            self.assertIn('remaining <strong>%d IDs</strong>' % sample_size,
                          index_contents)
            self.assertIn(', '.join(sorted(mismatched_ids)), index_contents)

        if exp_test_stat is not None:
            self.assertIn('<td>%r</td>' % exp_test_stat, index_contents)

        if exp_p_value is not None:
            self.assertIn('<td>%s</td>' % exp_p_value, index_contents)

        svg_fp = os.path.join(viz_dir, 'mantel-scatter.svg')
        self.assertTrue(os.path.exists(svg_fp))

        with open(svg_fp, 'r') as fh:
            svg_contents = fh.read()

        self.assertIn('Pairwise Distance (%s)' % label1, svg_contents)
        self.assertIn('Pairwise Distance (%s)' % label2, svg_contents)

    def test_defaults_positive_correlation(self):
        mantel(self.output_dir, self.dm1, self.dm2)

        self.assertBasicVizValidity(self.output_dir, 3, exp_test_stat=0.5,
                                    exp_p_value=1)

    def test_defaults_negative_correlation(self):
        mantel(self.output_dir, self.dm1, self.dm3)

        # p-value will be stochastic with this dataset so not asserting its
        # value.
        self.assertBasicVizValidity(self.output_dir, 3, exp_test_stat=-1)

    def test_defaults_reverse_comparison(self):
        # Comparing X to Y should be the same as comparing Y to X.
        mantel(self.output_dir, self.dm2, self.dm1)

        self.assertBasicVizValidity(self.output_dir, 3, exp_test_stat=0.5,
                                    exp_p_value=1)

    def test_defaults_reordered(self):
        # Order of IDs in distance matrices shouldn't change the results.
        mantel(self.output_dir, self.dm1, self.dm2_reordered)

        self.assertBasicVizValidity(self.output_dir, 3, exp_test_stat=0.5,
                                    exp_p_value=1)

    def test_pearson(self):
        mantel(self.output_dir, self.dm1, self.dm2, method='pearson')

        self.assertBasicVizValidity(self.output_dir, 3, method='pearson',
                                    exp_test_stat=0.5, exp_p_value=1)

    def test_alt_permutations(self):
        mantel(self.output_dir, self.dm1, self.dm2, permutations=42)

        self.assertBasicVizValidity(self.output_dir, 3, permutations=42,
                                    exp_test_stat=0.5, exp_p_value=1)

    def test_zero_permutations(self):
        mantel(self.output_dir, self.dm1, self.dm2, permutations=0)

        self.assertBasicVizValidity(self.output_dir, 3, permutations=0,
                                    exp_test_stat=0.5, exp_p_value='NaN')

    def test_alt_labels(self):
        mantel(self.output_dir, self.dm1, self.dm2, label1='Peanut',
               label2='Milo')

        self.assertBasicVizValidity(self.output_dir, 3, label1='Peanut',
                                    label2='Milo', exp_test_stat=0.5,
                                    exp_p_value=1)

    def test_error_on_sample_mismatch(self):
        with self.assertRaisesRegex(ValueError,
                                    'intersect_ids.*mismatches.*\n\nfoo, x'):
            mantel(self.output_dir, self.dm1, self.mismatched_dm)

    def test_warn_on_sample_mismatch(self):
        mantel(self.output_dir, self.dm1, self.mismatched_dm,
               intersect_ids=True)

        self.assertBasicVizValidity(self.output_dir, 3,
                                    mismatched_ids={'foo', 'x'},
                                    exp_test_stat=0.5, exp_p_value=1)

    def test_warn_on_sample_mismatch_reverse_comparison(self):
        # Comparing X to Y should be the same as comparing Y to X.
        mantel(self.output_dir, self.mismatched_dm, self.dm1,
               intersect_ids=True)

        self.assertBasicVizValidity(self.output_dir, 3,
                                    mismatched_ids={'foo', 'x'},
                                    exp_test_stat=0.5, exp_p_value=1)

    def test_same_matrices(self):
        mantel(self.output_dir, self.dm1, self.dm1)

        self.assertBasicVizValidity(self.output_dir, 3, exp_test_stat=1,
                                    exp_p_value=1)

    def test_same_filtered_matrices(self):
        # These two matrices filter down to the same data.
        mantel(self.output_dir, self.dm2, self.mismatched_dm,
               intersect_ids=True)

        self.assertBasicVizValidity(self.output_dir, 3,
                                    mismatched_ids={'foo', 'x'},
                                    exp_test_stat=1, exp_p_value=1)


if __name__ == '__main__':
    unittest.main()
