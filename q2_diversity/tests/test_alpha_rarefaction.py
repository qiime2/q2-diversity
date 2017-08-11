# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

import biom
import numpy as np
import numpy.testing as npt
import pandas.testing as pdt
import qiime2
import pandas as pd
from q2_diversity import alpha_rarefaction
from q2_diversity._alpha._visualizer import _compute_rarefaction_data

class AlphaRarefactionTests(unittest.TestCase):

    def test_alpha_rarefaction(self):
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
        pdt.assert_frame_equal(obs[0]['observed_otus'], exp)
        npt.assert_array_equal(obs[1], np.array([1, 200]))
        npt.assert_array_equal(obs[2], np.array([1]))

    def test_compute_rarefaction_data_123(self):
        raise NotImplementedError('Many more tests are needed - this is just '
                                  'a stub so far! ')

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

            with self.assertRaisesRegex(ValueError, 'too few steps'):
                alpha_rarefaction(output_dir, t, max_depth=200, steps=1,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'too few iterations'):
                alpha_rarefaction(output_dir, t, max_depth=200, iterations=0,
                                  metadata=md)

            with self.assertRaisesRegex(ValueError, 'no phylogeny'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metrics=['faith_pd'])

            with self.assertRaisesRegex(ValueError, 'unknown metric'):
                alpha_rarefaction(output_dir, t, max_depth=200,
                                  metadata=md, metrics=['pole-position'])
