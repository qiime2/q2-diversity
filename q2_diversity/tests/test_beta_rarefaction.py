# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import functools
import tempfile
import os

import skbio
import qiime2
from biom import Table
import numpy as np
import numpy.testing as npt
import pandas as pd
import scipy

from qiime2.plugin.testing import TestPluginBase
from q2_diversity import beta_rarefaction
from q2_diversity._beta._beta_rarefaction import (
    _get_multiple_rarefaction, _upgma, _cluster_samples, _add_support_count,
    _jackknifed_emperor)


class SharedSetup:
    def setUp(self):
        self.table = Table(np.array([[0, 1, 3],
                                     [1, 1, 2],
                                     [2, 1, 0]]),
                           ['O1', 'O2', 'O3'], ['S1', 'S2', 'S3'])
        self.tree = skbio.TreeNode.read([
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'])

        self.md = qiime2.Metadata(
            pd.DataFrame({'a': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        self.output_dir_obj = tempfile.TemporaryDirectory(
                prefix='q2-diversity-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()


class BetaRarefactionTests(SharedSetup, TestPluginBase):
    package = 'q2_diversity.tests'

    def check_heatmap(self, viz_dir, iterations, correlation_method):
        heatmap_fp = os.path.join(viz_dir, 'heatmap.html')
        self.assertTrue(os.path.exists(heatmap_fp))
        with open(heatmap_fp, 'r') as fh:
            heatmap_contents = fh.read()
        self.assertIn('Heatmap -', heatmap_contents)

        svg_fp = os.path.join(viz_dir, 'heatmap.svg')
        self.assertTrue(os.path.exists(svg_fp))
        with open(svg_fp, 'r') as fh:
            svg_contents = fh.read()

        self.assertIn('Mantel correlation', svg_contents)
        self.assertIn('Iteration', svg_contents)
        test_statistics = {'spearman': "Spearman's rho",
                           'pearson': "Pearson's r"}
        self.assertIn(test_statistics[correlation_method], svg_contents)

        tsv_fp = os.path.join(viz_dir, 'rarefaction-iteration-correlation.tsv')
        self.assertTrue(os.path.exists(tsv_fp))
        with open(tsv_fp, 'r') as fh:
            tsv_contents = fh.read()
        self.assertIn(correlation_method, tsv_contents)

        # TSV has a header line and trailing newline, substract 2 to get the
        # number of data lines in the file. Each data line represents a
        # pairwise comparison; assert the number of comparisons is equal to
        # nCr (n Choose r), where n=`iterations` and r=2.
        self.assertEqual(len(tsv_contents.split('\n')) - 2,
                         scipy.special.comb(iterations, 2))

    def check_clustering(self, viz_dir, clustering_method):
        cluster_fp = os.path.join(viz_dir, 'tree.html')
        self.assertTrue(os.path.exists(cluster_fp))

        newick_fp = os.path.join(
            viz_dir, 'sample-clustering-%s.tre' % clustering_method)
        self.assertTrue(os.path.exists(newick_fp))

        tree = skbio.TreeNode.read(newick_fp)
        self.assertEqual({t.name for t in tree.tips()}, {'S1', 'S2', 'S3'})

    def check_emperor(self, viz_dir):
        emperor_tab_fp = os.path.join(viz_dir, 'emperor.html')
        self.assertTrue(os.path.exists(emperor_tab_fp))

        emperor_dir = os.path.join(viz_dir, 'emperor')
        self.assertTrue(os.path.isdir(emperor_dir))

    def assertBetaRarefactionValidity(self, viz_dir, iterations,
                                      correlation_method, clustering_method):
        self.check_heatmap(viz_dir, iterations, correlation_method)
        self.check_clustering(viz_dir, clustering_method)
        self.check_emperor(viz_dir)

    def test_beta_rarefaction_with_phylogeny(self):
        beta_rarefaction(self.output_dir, self.table,
                         'weighted_unifrac',
                         'upgma', self.md, 2, phylogeny=self.tree)

        self.assertBetaRarefactionValidity(
            self.output_dir, 10, 'spearman', 'upgma')

    def test_beta_rarefaction_without_phylogeny(self):
        beta_rarefaction(self.output_dir, self.table, 'braycurtis', 'upgma',
                         self.md, 2)

        self.assertBetaRarefactionValidity(
            self.output_dir, 10, 'spearman', 'upgma')

    def test_beta_rarefaction_minimum_iterations(self):
        beta_rarefaction(self.output_dir, self.table, 'braycurtis', 'upgma',
                         self.md, 2, iterations=2)

        self.assertBetaRarefactionValidity(
            self.output_dir, 2, 'spearman', 'upgma')

    def test_beta_rarefaction_pearson_correlation(self):
        beta_rarefaction(self.output_dir, self.table, 'jaccard', 'upgma',
                         self.md, 2, iterations=7, correlation_method='pearson'
                         )

        self.assertBetaRarefactionValidity(
            self.output_dir, 7, 'pearson', 'upgma')

    def test_beta_rarefaction_non_default_color_scheme(self):
        beta_rarefaction(self.output_dir, self.table, 'euclidean', 'upgma',
                         self.md, 3, iterations=5, color_scheme='PiYG')

        self.assertBetaRarefactionValidity(
            self.output_dir, 5, 'spearman', 'upgma')

    def test_beta_rarefaction_neighbor_joining(self):
        beta_rarefaction(self.output_dir, self.table, 'euclidean', 'nj',
                         self.md, 3, iterations=5, color_scheme='PiYG')

        self.assertBetaRarefactionValidity(
            self.output_dir, 5, 'spearman', 'nj')

    def test_beta_rarefaction_empty_table(self):
        table = Table(np.array([[]]), [], [])
        with self.assertRaisesRegex(ValueError, 'feature table is empty'):
            beta_rarefaction(self.output_dir, table, 'braycurtis', 'upgma',
                             self.md, 1)

    def test_beta_rarefaction_all_samples_dropped(self):
        with self.assertRaisesRegex(ValueError,
                                    'shallow enough sampling depth'):
            beta_rarefaction(self.output_dir, self.table, 'braycurtis',
                             'upgma', self.md, 100)

    def test_beta_rarefaction_too_many_samples_dropped(self):
        # mantel needs 3x3 or larger distance matrix
        table = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                      ['O1', 'O2'], ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, '3x3 in size'):
            beta_rarefaction(self.output_dir, table, 'braycurtis', 'upgma',
                             self.md, 2)

    def test_beta_rarefaction_missing_phylogeny(self):
        with self.assertRaisesRegex(ValueError, 'Phylogeny must be provided'):
            beta_rarefaction(self.output_dir, self.table,
                             'weighted_unifrac',
                             'upgma', self.md, 2)


class GetMultipleRarefactionTests(SharedSetup, TestPluginBase):
    package = 'q2_diversity.tests'

    def test_with_phylogeny(self):
        with qiime2.sdk.Context() as scope:
            table = qiime2.Artifact.import_data('FeatureTable[Frequency]',
                                                self.table)
            tree = qiime2.Artifact.import_data('Phylogeny[Rooted]',
                                               self.tree)
            api_method = scope.ctx.get_action('diversity', 'beta_phylogenetic')
            beta_func = functools.partial(api_method, phylogeny=tree)
            rare_func = scope.ctx.get_action('feature-table', 'rarefy')

            for iterations in range(1, 4):
                obs_dms = _get_multiple_rarefaction(
                        beta_func, rare_func, 'weighted_unifrac',
                        iterations, table, 2)

                self.assertEqual(len(obs_dms), iterations)
                for obs in obs_dms:
                    self.assertEqual(obs.shape, (3, 3))
                    self.assertEqual(set(obs.ids), set(['S1', 'S2', 'S3']))

    def test_without_phylogeny(self):
        with qiime2.sdk.Context() as scope:
            table = qiime2.Artifact.import_data('FeatureTable[Frequency]',
                                                self.table)
            beta_func = scope.ctx.get_action('diversity', 'beta')
            rare_func = scope.ctx.get_action('feature-table', 'rarefy')
            for iterations in range(1, 4):
                obs_dms = _get_multiple_rarefaction(beta_func, rare_func,
                                                    'braycurtis', iterations,
                                                    table, 2)

                self.assertEqual(len(obs_dms), iterations)
                for obs in obs_dms:
                    self.assertEqual(obs.shape, (3, 3))
                    self.assertEqual(set(obs.ids), set(['S1', 'S2', 'S3']))


class UPGMATests(unittest.TestCase):
    # The translation between skbio and scipy is a little spooky, so these
    # tests just confirm that the ids don't get jumbled along the way
    def test_ids_retained(self):
        # This makes a very simple (and comb-like) UPGMA tree
        dm = skbio.DistanceMatrix(
            [[0, 1, 3, 5],
             [1, 0, 7, 9],
             [3, 7, 0, 11],
             [5, 9, 11, 0]],
            ids=['a', 'b', 'c', 'd'])

        tree = _upgma(dm)

        # Nodes exist
        a = tree.find('a')
        b = tree.find('b')
        c = tree.find('c')
        d = tree.find('d')

        # Check topology quickly. If the IDs were flipped or wrong, these
        # checks would fail, as we're starting from the tips and working to the
        # root.
        self.assertIs(a.parent, b.parent)
        self.assertIs(b.parent.parent, c.parent)
        self.assertIs(c.parent.parent, d.parent)
        self.assertIs(d.parent.parent, None)


class ClusterSamplesTests(unittest.TestCase):
    def setUp(self):
        self.dm = skbio.DistanceMatrix([[0, 1, 2.1], [1, 0, 3], [2.1, 3, 0]],
                                       ids=['S1', 'S2', 'S3'])

        # Since support is traditionally held as the name, we'll use only two
        # trees since 1/2 has an exact floating point representation and will
        # look like `"0.5"` on any machine.
        self.support = [
            skbio.DistanceMatrix([[0, 1.1, 2], [1.1, 0, 3], [2, 3, 0]],
                                 ids=['S1', 'S2', 'S3']),
            skbio.DistanceMatrix([[0, 2, 3.1], [2, 0, 1], [3.1, 1, 0]],
                                 ids=['S1', 'S2', 'S3'])
        ]

    def test_nj_support(self):
        result = _cluster_samples(self.dm, self.support, 'nj')

        s1 = result.find('S1')
        s2 = result.find('S2')
        s3 = result.find('S3')

        npt.assert_almost_equal(s1.length, 0.05)
        npt.assert_almost_equal(s2.length, 0.95)
        self.assertIs(s1.parent, s2.parent)
        s1_s2 = s1.parent

        self.assertEqual(s1_s2.name, '0.5')  # half of the support agrees
        npt.assert_almost_equal(s1_s2.length, 1)

        npt.assert_almost_equal(s3.length, 1.05)
        self.assertIs(s3.parent, s1_s2.parent)
        self.assertIs(s3.parent.name, 'root')  # root support is pointless

    def test_upgma_support(self):
        result = _cluster_samples(self.dm, self.support, 'upgma')

        s1 = result.find('S1')
        s2 = result.find('S2')
        s3 = result.find('S3')

        npt.assert_almost_equal(s1.length, 0.5)
        npt.assert_almost_equal(s2.length, 0.5)  # unimetric tree!
        self.assertIs(s1.parent, s2.parent)
        s1_s2 = s1.parent

        self.assertEqual(s1_s2.name, '0.5')  # half of the support agrees
        npt.assert_almost_equal(s1_s2.length, 0.775)

        npt.assert_almost_equal(s3.length, 1.275)
        self.assertIs(s3.parent, s1_s2.parent)
        self.assertIs(s3.parent.name, 'root')  # root support is pointless


class AddSupportCountTests(unittest.TestCase):
    def test_same_topology(self):
        p = skbio.TreeNode.read(['((a,b),c);'])
        s = skbio.TreeNode.read(['((a,b),c);'])

        internal = list(p.non_tips())
        for n in internal:
            n.support_count = 0

        _add_support_count(internal, s)

        for n in internal:
            self.assertEqual(n.support_count, 1)

    def test_differing_topology(self):
        p = skbio.TreeNode.read(['(((a,b),c),d);'])
        s = skbio.TreeNode.read(['(((a,c),b),d);'])

        internal = list(p.non_tips())
        for n in internal:
            n.support_count = 0

        _add_support_count(internal, s)

        a_b = p.find('a').parent
        a_b_c = a_b.parent

        self.assertEqual(a_b.support_count, 0)
        self.assertEqual(a_b_c.support_count, 1)

    def test_extra_node(self):
        p = skbio.TreeNode.read(['(((a,b),c),d);'])
        s = skbio.TreeNode.read(['((a,b),(c,d));'])
        # The first node that has a, b, and c, also has d as a descendant

        internal = list(p.non_tips())
        for n in internal:
            n.support_count = 0

        _add_support_count(internal, s)

        a_b = p.find('a').parent
        a_b_c = a_b.parent

        self.assertEqual(a_b.support_count, 1)
        self.assertEqual(a_b_c.support_count, 0)

    def test_multiple_calls_with_outgroup(self):
        p = skbio.TreeNode.read(['((((a,b),c),d),e);'])
        s1 = skbio.TreeNode.read(['(((c,b),(a,d)),e);'])
        s2 = skbio.TreeNode.read(['((((a,c),b),d),e);'])
        # e is the outgroup here

        internal = list(p.non_tips())
        for n in internal:
            n.support_count = 0

        _add_support_count(internal, s1)
        _add_support_count(internal, s2)

        a_b = p.find('a').parent
        a_b_c = a_b.parent
        a_b_c_d = a_b_c.parent

        self.assertEqual(a_b.support_count, 0)
        self.assertEqual(a_b_c.support_count, 1)
        self.assertEqual(a_b_c_d.support_count, 2)


class JackknifedEmperorTests(SharedSetup, unittest.TestCase):
    def test_simple(self):
        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ids=['S1', 'S2', 'S3'])

        j1 = skbio.DistanceMatrix([[0, 1.1, 2], [1.1, 0, 3], [2, 3, 0]],
                                  ids=['S1', 'S2', 'S3'])
        j2 = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3.1], [2, 3.1, 0]],
                                  ids=['S1', 'S2', 'S3'])
        j3 = skbio.DistanceMatrix([[0, 1.1, 1.9], [1.1, 0, 3], [1.9, 3, 0]],
                                  ids=['S1', 'S2', 'S3'])

        e = _jackknifed_emperor(dm, [j1, j2, j3], self.md)

        self.assertEqual(len(e.jackknifed), 3)


if __name__ == "__main__":
    unittest.main()
