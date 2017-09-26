# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import io
import unittest

import biom
import skbio
import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from q2_diversity import core_metrics_phylogenetic, core_metrics


class CoreMetricsTests(unittest.TestCase):
    def test_core_metrics_phylogenetic(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        results = core_metrics_phylogenetic(table, tree, 13)

        self.assertEqual(len(results), 13)

        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2},
                             name='observed_otus')
        pdt.assert_series_equal(results[2], expected)

    def test_core_metrics_phylogenetic_rarefy_drops_sample(self):
        table = biom.Table(np.array([[0, 11, 11], [12, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        results = core_metrics_phylogenetic(table, tree, 13)

        self.assertEqual(len(results), 13)

        expected = pd.Series({'S2': 2, 'S3': 2},
                             name='observed_otus')
        pdt.assert_series_equal(results[2], expected)

    def test_core_metrics(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        results = core_metrics(table, 13)

        self.assertEqual(len(results), 8)

        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2}, name='observed_otus')
        pdt.assert_series_equal(results[1], expected)
