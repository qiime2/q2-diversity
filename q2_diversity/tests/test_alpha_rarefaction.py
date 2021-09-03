# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import tempfile
import unittest

import biom
import numpy as np
import pandas.testing as pdt
import qiime2
import skbio
import pandas as pd

from qiime2.plugin.util import transform
from q2_types.tree import NewickFormat
from q2_diversity import alpha_rarefaction
from q2_diversity._alpha._visualizer import (
    _compute_rarefaction_data, _compute_summary, _reindex_with_metadata,
    _alpha_rarefaction_jsonp)


class AlphaRarefactionTests(unittest.TestCase):

    @staticmethod
    def _to_newick(tree: skbio.TreeNode):
        return transform(tree, from_type=skbio.TreeNode,
                         to_type=NewickFormat)

    def test_alpha_rarefaction_without_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)

    def test_alpha_rarefaction_with_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)

    def test_alpha_rarefaction_with_superset_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut', 'summer']},
                         index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)
            metric_fp = os.path.join(output_dir, 'shannon-pet.jsonp')
            with open(metric_fp) as metric_fh:
                self.assertTrue('summer' not in metric_fh.read())

    def test_alpha_rarefaction_with_filtered_metadata_columns(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        # Empty column and numeric column should both be filtered.
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut', 'summer'],
                          'foo': [np.nan, np.nan, np.nan, 'bar'],
                          'bar': [42, 4.2, 99.9, 100.0]},
                         index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp, 'r') as fh:
                contents = fh.read()
            self.assertTrue('observed_features' in contents)
            self.assertTrue('shannon' in contents)
            self.assertTrue('didn\'t contain categorical data' in contents)
            self.assertTrue('consisted only of missing values:' in contents)
            self.assertTrue('<strong>bar, foo' in contents)

            metric_fp = os.path.join(output_dir, 'shannon-pet.jsonp')
            with open(metric_fp) as metric_fh:
                self.assertTrue('summer' not in metric_fh.read())
            self.assertFalse(
                os.path.exists(os.path.join(output_dir,
                               'shannon_entropy-foo.jsonp')))
            self.assertFalse(
                os.path.exists(os.path.join(output_dir,
                               'shannon_entropy-bar.jsonp')))

    def test_alpha_rarefaction_with_depth_column_in_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'depth': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)

    def test_alpha_rarefaction_with_phylogeny(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        p = self._to_newick(skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, phylogeny=p)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)
            self.assertTrue('faith_pd' in index_content)

    def test_alpha_rarefaction_with_phylogeny_and_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        p = self._to_newick(skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;')))
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, phylogeny=p,
                              metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)
            self.assertTrue('faith_pd' in index_content)

    def test_invalid_invocations(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        empty_table = biom.Table(np.array([]), [], [])

        bad_metadata = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'summer']},
                         index=pd.Index(['S1', 'S2', 'S4'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(ValueError, 'must be greater'):
                alpha_rarefaction(output_dir, t, min_depth=200, max_depth=1,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'phylogeny was not'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metrics=set(['faith_pd']))

            with self.assertRaisesRegex(TypeError, 'pole.*incompatible'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metrics=set(['pole-position']))

            with self.assertRaisesRegex(ValueError, 'max_depth'):
                alpha_rarefaction(output_dir, t, max_depth=1000)

            with self.assertRaisesRegex(ValueError, 'steps'):
                alpha_rarefaction(output_dir, t, max_depth=2)

            with self.assertRaisesRegex(ValueError, 'empty'):
                alpha_rarefaction(output_dir, empty_table, max_depth=200)

            with self.assertRaisesRegex(ValueError, 'not present.*S3'):
                alpha_rarefaction(output_dir, t, metadata=bad_metadata,
                                  max_depth=200)

            with self.assertRaisesRegex(ValueError, 'empty set'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metrics=set())

    def test_alpha_rarefaction_with_metric_set(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        metrics = set(['observed_features', 'shannon', 'pielou_e'])
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, metrics=metrics,
                              max_depth=200)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('observed_features' in index_content)
            self.assertTrue('shannon' in index_content)
            self.assertTrue('pielou_e' in index_content)

    def test_alpha_rarefaction_with_metadata_column_with_spaces(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet name': ['russ', 'milo', 'peanut']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            with open(index_fp) as index_fh:
                self.assertTrue('pet%2520name' in index_fh.read())
            jsonp_fp = os.path.join(output_dir,
                                    'shannon-pet%20name.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))
            with open(jsonp_fp) as jsonp_fh:
                self.assertTrue('pet name' in jsonp_fh.read())

    def test_alpha_rarefaction_with_fully_filtered_metadata_columns(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        # All columns should both be filtered and raise an error
        md = qiime2.Metadata(
            pd.DataFrame({'foo': [np.nan, np.nan, np.nan, 'bar'],
                          'bar': [42, 4.2, 99.9, 100.0]},
                         index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(ValueError, "non-categorical"):
                alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)


class ComputeRarefactionDataTests(unittest.TestCase):
    def setUp(self):
        np.random.seed(0)

    @staticmethod
    def _to_newick(tree: skbio.TreeNode):
        return transform(tree, from_type=skbio.TreeNode,
                         to_type=NewickFormat)

    def test_observed_features(self):
        t = biom.Table(np.array([[150, 100, 100], [50, 100, 100]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        obs = _compute_rarefaction_data(feature_table=t,
                                        min_depth=1,
                                        max_depth=200,
                                        steps=2,
                                        iterations=1,
                                        phylogeny=None,
                                        metrics=['observed_features'])

        exp_ind = pd.MultiIndex.from_product(
            [[1, 200], [1]],
            names=['_alpha_rarefaction_depth_column_', 'iter'])
        exp = pd.DataFrame(data=[[1, 2], [1, 2], [1, 2]],
                           columns=exp_ind,
                           index=['S1', 'S2', 'S3'])
        pdt.assert_frame_equal(obs['observed_features'], exp)

    def test_faith_pd(self):
        t = biom.Table(np.array([[150, 100, 100], [50, 100, 100]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        p = self._to_newick(skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;')))

        obs = _compute_rarefaction_data(feature_table=t,
                                        min_depth=1,
                                        max_depth=200,
                                        steps=2,
                                        iterations=1,
                                        phylogeny=p,
                                        metrics=['faith_pd'])

        self.assertTrue('faith_pd' in obs)

    def test_multiple_metrics(self):
        t = biom.Table(np.array([[150, 100, 100], [50, 100, 100]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        obs = _compute_rarefaction_data(feature_table=t,
                                        min_depth=1,
                                        max_depth=200,
                                        steps=2,
                                        iterations=1,
                                        phylogeny=None,
                                        metrics=['observed_features',
                                                 'shannon'])

        exp_ind = pd.MultiIndex.from_product(
            [[1, 200], [1]],
            names=['_alpha_rarefaction_depth_column_', 'iter'])
        exp = pd.DataFrame(data=[[1, 2], [1, 2], [1, 2]],
                           columns=exp_ind,
                           index=['S1', 'S2', 'S3'])
        pdt.assert_frame_equal(obs['observed_features'], exp)

        exp = pd.DataFrame(data=[[0., 0.811278124459], [0., 1.], [0., 1.]],
                           columns=exp_ind,
                           index=['S1', 'S2', 'S3'])
        pdt.assert_frame_equal(obs['shannon'], exp)


class ComputeSummaryTests(unittest.TestCase):
    def test_one_iteration_no_metadata(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2], [1, 2], [1, 2]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        # No counts provided because no metadata
        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   1, 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S1', 200, 1, 2., 2., 2., 2., 2., 2., 2., 2., 2.],
             ['S2', 1,   1, 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S2', 200, 1, 2., 2., 2., 2., 2., 2., 2., 2., 2.],
             ['S3', 1,   1, 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S3', 200, 1, 2., 2., 2., 2., 2., 2., 2., 2., 2.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_two_iterations_no_metadata(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        # No counts provided because no metadata
        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   1, 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S1', 200, 1, 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
             ['S2', 1,   1, 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S2', 200, 1, 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
             ['S3', 1,   1, 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S3', 200, 1, 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_three_iterations_no_metadata(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2, 3]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6],
                                  [1, 2, 3, 4, 5, 6]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        # No counts provided because no metadata
        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   1, 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S1', 200, 1, 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.],
             ['S2', 1,   1, 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S2', 200, 1, 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.],
             ['S3', 1,   1, 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S3', 200, 1, 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_two_iterations_with_metadata_were_values_are_unique(self):
        # This should be identical to test_without_metadata_df_two_iterations,
        # with just the `sample-id` replaced with `pet`.
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                            columns=columns, index=['russ', 'milo', 'pea'])

        counts = pd.DataFrame(data=[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
                              columns=columns, index=['russ', 'milo', 'pea'])

        obs = _compute_summary(data, 'pet', counts=counts)

        d = [
            ['russ', 1,   1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2., 1],
            ['russ', 200, 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4., 1],
            ['milo', 1,   1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2., 1],
            ['milo', 200, 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4., 1],
            ['pea', 1,    1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2., 1],
            ['pea', 200,  3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4., 1],
        ]
        exp = pd.DataFrame(data=d, columns=['pet', 'depth', 'min', '2%', '9%',
                                            '25%', '50%', '75%', '91%', '98%',
                                            'max', 'count'])
        pdt.assert_frame_equal(exp, obs)

    def test_two_iterations_with_metadata_were_values_are_identical(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[3, 6, 9, 9]], columns=columns,
                            index=['milo'])

        counts = pd.DataFrame(data=[[3, 3, 3, 3]], columns=columns,
                              index=['milo'])

        obs = _compute_summary(data, 'pet', counts=counts)

        d = [
            ['milo', 1,   3., 3.06, 3.27, 3.75, 4.5,  5.25, 5.73, 5.94, 6., 3],
            ['milo', 200, 9.,   9.,   9.,   9.,  9.,    9.,   9.,   9., 9., 3],
        ]
        exp = pd.DataFrame(data=d, columns=['pet', 'depth', 'min', '2%', '9%',
                                            '25%', '50%', '75%', '91%', '98%',
                                            'max', 'count'])
        pdt.assert_frame_equal(exp, obs)


class ReindexWithMetadataTests(unittest.TestCase):
    def test_unique_metadata_groups(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'milo'],
                                  [9, 10, 11, 12, 'peanut']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        median, counts = _reindex_with_metadata('pet', ['pet'], data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'peanut', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [9, 10, 11, 12], [1, 2, 3, 4]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, median)

        exp = pd.DataFrame(data=[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, counts)

    def test_some_duplicates_in_column(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'milo'],
                                  [9, 10, 11, 12, 'russ']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        median, counts = _reindex_with_metadata('pet', ['pet'], data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [5, 6, 7, 8]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, median)

        exp = pd.DataFrame(data=[[1, 1, 1, 1], [2, 2, 2, 2]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, counts)

    def test_all_identical(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'russ'],
                                  [9, 10, 11, 12, 'russ']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        median, counts = _reindex_with_metadata('pet', ['pet'], data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, median)

        exp = pd.DataFrame(data=[[3, 3, 3, 3]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, counts)

    def test_multiple_columns(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', ''),
                                             ('toy', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ', 'stick'],
                                  [5, 6, 7, 8, 'milo', 'yeti'],
                                  [9, 10, 11, 12, 'peanut', 'stick']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        median, counts = _reindex_with_metadata('pet', ['pet', 'toy'], data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet', 'toy'], [1, 2, '']],
                                codes=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'peanut', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [9, 10, 11, 12], [1, 2, 3, 4]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, median)

        exp = pd.DataFrame(data=[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, counts)

        median, counts = _reindex_with_metadata('toy', ['pet', 'toy'], data)

        exp_ind = pd.Index(['stick', 'yeti'], name='toy')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [5, 6, 7, 8]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, median)

        exp = pd.DataFrame(data=[[2, 2, 2, 2], [1, 1, 1, 1]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, counts)


class AlphaRarefactionJSONPTests(unittest.TestCase):
    def test_simple(self):
        d = [[1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S1'],
             [1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S2'],
             [1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S3']]

        data = pd.DataFrame(data=d, columns=['2%', '25%', '50%', '75%', '9%',
                                             '91%', '98%', 'count', 'depth',
                                             'max', 'min', 'sample-id'])

        with tempfile.TemporaryDirectory() as output_dir:
            _alpha_rarefaction_jsonp(output_dir, 'peanut.jsonp', 'shannon',
                                     data, '')

            jsonp_fp = os.path.join(output_dir, 'peanut.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))
            with open(jsonp_fp) as jsonp_fh:
                jsonp_content = jsonp_fh.read()
            self.assertTrue('load_data' in jsonp_content)
            self.assertTrue('columns' in jsonp_content)
            self.assertTrue('index' in jsonp_content)
            self.assertTrue('data' in jsonp_content)
            self.assertTrue('sample-id' in jsonp_content)
            self.assertTrue('shannon' in jsonp_content)


if __name__ == '__main__':
    unittest.main()
