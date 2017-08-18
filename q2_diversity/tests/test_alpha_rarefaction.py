# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
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
from q2_diversity import alpha_rarefaction
from q2_diversity._alpha._visualizer import (
    _compute_rarefaction_data, _compute_summary, _reindex_with_metadata,
    _beta_rarefaction_jsonp)


class AlphaRarefactionTests(unittest.TestCase):

    def test_alpha_rarefaction_without_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('observed_otus' in open(index_fp).read())
            self.assertTrue('shannon' in open(index_fp).read())

    def test_alpha_rarefaction_with_phylogeny(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        p = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, phylogeny=p)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('observed_otus' in open(index_fp).read())
            self.assertTrue('shannon' in open(index_fp).read())
            self.assertTrue('faith_pd' in open(index_fp).read())

    def test_alpha_rarefaction_with_metadata(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut']},
                         index=['S1', 'S2', 'S3']))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_rarefaction(output_dir, t, max_depth=200, metadata=md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('observed_otus' in open(index_fp).read())
            self.assertTrue('shannon' in open(index_fp).read())

    def test_compute_rarefaction_data(self):
        t = biom.Table(np.array([[150, 100, 100], [50, 100, 100]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        obs = _compute_rarefaction_data(feature_table=t,
                                        min_depth=1,
                                        max_depth=200,
                                        steps=2,
                                        iterations=1,
                                        phylogeny=None,
                                        metrics=['observed_otus'])

        exp_ind = pd.MultiIndex.from_product(
            [[1, 200], [1]],
            names=['depth', 'iter'])
        exp = pd.DataFrame(data=[[1, 2], [1, 2], [1, 2]],
                           columns=exp_ind,
                           index=['S1', 'S2', 'S3'])
        pdt.assert_frame_equal(obs['observed_otus'], exp)

    def test_compute_summary_one_iteration(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2], [1, 2], [1, 2]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S1', 200, 1., 2., 2., 2., 2., 2., 2., 2., 2., 2.],
             ['S2', 1,   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S2', 200, 1., 2., 2., 2., 2., 2., 2., 2., 2., 2.],
             ['S3', 1,   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.],
             ['S3', 200, 1., 2., 2., 2., 2., 2., 2., 2., 2., 2.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_compute_summary_two_iterations(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S1', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
             ['S2', 1,   2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S2', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
             ['S3', 1,   2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
             ['S3', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_compute_summary_three_iterations(self):
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2, 3]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6],
                                  [1, 2, 3, 4, 5, 6]],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _compute_summary(data, 'sample-id')

        d = [['S1', 1,   3., 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S1', 200, 3., 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.],
             ['S2', 1,   3., 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S2', 200, 3., 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.],
             ['S3', 1,   3., 1., 1.04, 1.18, 1.5, 2., 2.5, 2.82, 2.96, 3.],
             ['S3', 200, 3., 4., 4.04, 4.18, 4.5, 5., 5.5, 5.82, 5.96, 6.]]
        exp = pd.DataFrame(data=d, columns=['sample-id', 'depth', 'count',
                                            'min', '2%', '9%', '25%', '50%',
                                            '75%', '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_with_metadata_two_iterations_unique_metadata_groups(self):
        # This should be identical to test_without_metadata_df_two_iterations,
        # with just the `sample-id` replaced with `pet`.
        columns = pd.MultiIndex.from_product([[1, 200], [1, 2]],
                                             names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]],
                            columns=columns, index=['russ', 'milo', 'pea'])

        obs = _compute_summary(data, 'pet')

        d = [
            ['russ', 1, 2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
            ['russ', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
            ['milo', 1, 2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
            ['milo', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
            ['pea', 1, 2., 1., 1.02, 1.09, 1.25, 1.5, 1.75, 1.91, 1.98, 2.],
            ['pea', 200, 2., 3., 3.02, 3.09, 3.25, 3.5, 3.75, 3.91, 3.98, 4.],
        ]
        exp = pd.DataFrame(data=d, columns=['pet', 'depth', 'count', 'min',
                                            '2%', '9%', '25%', '50%', '75%',
                                            '91%', '98%', 'max'])
        pdt.assert_frame_equal(exp, obs)

    def test_reindex_with_metadata_unique_metadata_groups(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'milo'],
                                  [9, 10, 11, 12, 'peanut']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _reindex_with_metadata('pet', data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'peanut', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [9, 10, 11, 12], [1, 2, 3, 4]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, obs)

    def test_reindex_with_metadata_some_dupes(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'milo'],
                                  [9, 10, 11, 12, 'russ']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _reindex_with_metadata('pet', data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [10, 12, 14, 16]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, obs)

    def test_reindex_with_metadata_all_dupes(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ'], [5, 6, 7, 8, 'russ'],
                                  [9, 10, 11, 12, 'russ']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _reindex_with_metadata('pet', data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet'], [1, 2, '']],
                                labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['russ'], name='pet')
        exp = pd.DataFrame(data=[[15, 18, 21, 24]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, obs)

    def test_reindex_with_metadata_multiple_categories(self):
        columns = pd.MultiIndex.from_tuples([(1, 1), (1, 2), (200, 1),
                                             (200, 2), ('pet', ''),
                                             ('toy', '')],
                                            names=['depth', 'iter'])
        data = pd.DataFrame(data=[[1, 2, 3, 4, 'russ', 'stick'],
                                  [5, 6, 7, 8, 'milo', 'yeti'],
                                  [9, 10, 11, 12, 'peanut', 'stick']],
                            columns=columns, index=['S1', 'S2', 'S3'])

        obs = _reindex_with_metadata('pet', data)

        exp_col = pd.MultiIndex(levels=[[1, 200, 'pet', 'toy'], [1, 2, '']],
                                labels=[[0, 0, 1, 1], [0, 1, 0, 1]],
                                names=['depth', 'iter'])
        exp_ind = pd.Index(['milo', 'peanut', 'russ'], name='pet')
        exp = pd.DataFrame(data=[[5, 6, 7, 8], [9, 10, 11, 12], [1, 2, 3, 4]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, obs)

        obs = _reindex_with_metadata('toy', data)

        exp_ind = pd.Index(['stick', 'yeti'], name='toy')
        exp = pd.DataFrame(data=[[10, 12, 14, 16], [5, 6, 7, 8]],
                           columns=exp_col, index=exp_ind)

        pdt.assert_frame_equal(exp, obs)

    def test_beta_rarefaction_jsonp(self):
        d = [[1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S1'],
             [1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S2'],
             [1.04, 1.5, 2., 2.5, 1.18, 2.82, 2.96, 3., 1, 3., 1., 'S3']]

        data = pd.DataFrame(data=d, columns=['2%', '25%', '50%', '75%', '9%',
                                             '91%', '98%', 'count', 'depth',
                                             'max', 'min', 'sample-id'])

        with tempfile.TemporaryDirectory() as output_dir:
            _beta_rarefaction_jsonp(output_dir, 'peanut.jsonp', 'braycurtis',
                                    data, '')

            jsonp_fp = os.path.join(output_dir, 'peanut.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))
            jsonp_content = open(jsonp_fp).read()
            self.assertTrue('load_data' in jsonp_content)
            self.assertTrue('columns' in jsonp_content)
            self.assertTrue('index' in jsonp_content)
            self.assertTrue('data' in jsonp_content)
            self.assertTrue('sample-id' in jsonp_content)

    def test_alpha_rarefaction_invalid(self):
        t = biom.Table(np.array([[100, 111, 113], [111, 111, 112]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame({'pet': ['russ', 'milo', 'peanut']},
                         index=['S1', 'S2', 'S3']))

        with tempfile.TemporaryDirectory() as output_dir:

            with self.assertRaisesRegex(ValueError, 'must be greater'):
                alpha_rarefaction(output_dir, t, min_depth=200, max_depth=1,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'Provided steps'):
                alpha_rarefaction(output_dir, t, max_depth=200, steps=1,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'Provided iterations'):
                alpha_rarefaction(output_dir, t, max_depth=200, iterations=0,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'phylogeny was not'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metric='faith_pd')

            with self.assertRaisesRegex(ValueError, 'Unknown metric: pole'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metric='pole-position')

            with self.assertRaisesRegex(ValueError, 'max_depth'):
                alpha_rarefaction(output_dir, t, max_depth=1000)
