# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import tempfile
import os.path

import biom
import qiime2
import numpy as np
import pandas as pd
import pandas.testing as pdt

from q2_diversity import core_microbiome
from q2_diversity._core_microbiome import (
    _get_core_features, _seven_number_summary)


class CoreMicrobiomeTests(unittest.TestCase):

    def test_core_microbiome(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            core_microbiome(output_dir, table)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('Sample count: 3' in open(index_fp).read())
            self.assertTrue('Feature count: 2' in open(index_fp).read())
            self.assertFalse('Applied filter' in open(index_fp).read())

            svg_fp = os.path.join(output_dir, 'core-feature-counts.svg')
            self.assertTrue(os.path.exists(svg_fp))

            core_50_fp = os.path.join(output_dir, 'core-features-0.500.tsv')
            self.assertTrue(os.path.exists(core_50_fp))
            self.assertTrue('O1' in open(core_50_fp).read())
            self.assertTrue('O2' in open(core_50_fp).read())

            core_100_fp = os.path.join(output_dir, 'core-features-1.000.tsv')
            self.assertTrue(os.path.exists(core_100_fp))
            self.assertFalse('O1' in open(core_100_fp).read())
            self.assertTrue('O2' in open(core_100_fp).read())

    def test_core_microbiome_fraction_range(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            core_microbiome(output_dir, table, min_fraction=0.55,
                            max_fraction=0.95, steps=9)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('Sample count: 3' in open(index_fp).read())
            self.assertTrue('Feature count: 2' in open(index_fp).read())
            self.assertFalse('Applied filter' in open(index_fp).read())

            svg_fp = os.path.join(output_dir, 'core-feature-counts.svg')
            self.assertTrue(os.path.exists(svg_fp))

            core_50_fp = os.path.join(output_dir, 'core-features-0.500.tsv')
            self.assertFalse(os.path.exists(core_50_fp))

            core_55_fp = os.path.join(output_dir, 'core-features-0.550.tsv')
            self.assertTrue(os.path.exists(core_55_fp))
            self.assertTrue('O1' in open(core_55_fp).read())
            self.assertTrue('O2' in open(core_55_fp).read())

            core_95_fp = os.path.join(output_dir, 'core-features-0.950.tsv')
            self.assertTrue(os.path.exists(core_95_fp))
            self.assertFalse('O1' in open(core_95_fp).read())
            self.assertTrue('O2' in open(core_95_fp).read())

            core_100_fp = os.path.join(output_dir, 'core-features-1.000.tsv')
            self.assertFalse(os.path.exists(core_100_fp))

    def test_core_microbiome_no_core(self):
        table = biom.Table(np.array([[0, 11, 11], [11, 11, 0]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            core_microbiome(output_dir, table)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('Sample count: 3' in open(index_fp).read())
            self.assertTrue('Feature count: 2' in open(index_fp).read())
            self.assertFalse('Applied filter' in open(index_fp).read())

            svg_fp = os.path.join(output_dir, 'core-feature-counts.svg')
            self.assertTrue(os.path.exists(svg_fp))

            core_50_fp = os.path.join(output_dir, 'core-features-0.500.tsv')
            self.assertTrue(os.path.exists(core_50_fp))
            self.assertTrue('O1' in open(core_50_fp).read())
            self.assertTrue('O2' in open(core_50_fp).read())

            # No core features exist at fraction=1.0
            self.assertTrue('No core features' in open(index_fp).read())
            core_100_fp = os.path.join(output_dir, 'core-features-1.000.tsv')
            self.assertFalse(os.path.exists(core_100_fp))

    def test_core_microbiome_where(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame([['a', 'b'], ['b', 'b'], ['b', 'b']],
                         index=['S1', 'S2', 'S3'],
                         columns=['metadata1', 'metadata2']))

        with tempfile.TemporaryDirectory() as output_dir:
            core_microbiome(output_dir, table, metadata=md,
                            where="metadata1='b'")
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('Sample count pre-filtering: 3'
                            in open(index_fp).read())
            self.assertTrue('Feature count pre-filtering: 2'
                            in open(index_fp).read())
            self.assertTrue("Applied filter: metadata1='b'"
                            in open(index_fp).read())
            self.assertTrue('Sample count post-filtering: 2'
                            in open(index_fp).read())
            self.assertTrue('Feature count post-filtering: 2'
                            in open(index_fp).read())

            svg_fp = os.path.join(output_dir, 'core-feature-counts.svg')
            self.assertTrue(os.path.exists(svg_fp))

            core_50_fp = os.path.join(output_dir, 'core-features-0.500.tsv')
            self.assertTrue(os.path.exists(core_50_fp))
            self.assertTrue('O1' in open(core_50_fp).read())
            self.assertTrue('O2' in open(core_50_fp).read())

            core_100_fp = os.path.join(output_dir, 'core-features-1.000.tsv')
            self.assertTrue(os.path.exists(core_100_fp))
            # O1 is now in 100% of the samples
            self.assertTrue('O1' in open(core_100_fp).read())
            self.assertTrue('O2' in open(core_100_fp).read())

    def test_core_microbiome_invalid(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        md = qiime2.Metadata(
            pd.DataFrame([['a', 'b'], ['b', 'b'], ['b', 'b']],
                         index=['S1', 'S2', 'S3'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:

            with self.assertRaisesRegex(ValueError,
                                        expected_regex='fraction'):
                core_microbiome(output_dir, table, min_fraction=0.75,
                                max_fraction=0.5)

            with self.assertRaisesRegex(ValueError,
                                        expected_regex='Metadata'):
                core_microbiome(output_dir, table, where="name='peanut'")

            with self.assertRaisesRegex(
                    ValueError, expected_regex="matched no samples"):
                core_microbiome(output_dir, table, where="metadata2='peanut'",
                                metadata=md)

            with self.assertRaisesRegex(ValueError,
                                        expected_regex='Invalid fraction'):
                core_microbiome(output_dir, table, min_fraction=-0.01)

            with self.assertRaisesRegex(ValueError,
                                        expected_regex='Invalid fraction'):
                core_microbiome(output_dir, table, max_fraction=1.01)


class CoreMicrobiomePrivateFunctionTests(unittest.TestCase):
    # Most of the work happens in private functions for this visualizer,
    # so it's important to have some tests of the private functionality

    def test_seven_number_summary(self):
        # validated against calling pd.Series.describe directly
        exp = pd.Series([11., 11., 11., 11., 12., 12.64, 12.92],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([13, 11, 11]), exp)

        exp = pd.Series([1., 1., 1., 1., 1., 1., 1.],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([1]), exp)

        exp = pd.Series([np.nan, np.nan, np.nan, np.nan, np.nan,
                         np.nan, np.nan],
                        index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        pdt.assert_series_equal(_seven_number_summary([]), exp)

    def test_get_core_features(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        o2_seven_num = pd.Series(
            [11., 11., 11., 11., 12., 12.64, 12.92],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        exp = pd.DataFrame([o2_seven_num], index=['O2'])

        obs = _get_core_features(table, fraction=0.67)
        pdt.assert_frame_equal(obs, exp)

        obs = _get_core_features(table, fraction=1.0)
        pdt.assert_frame_equal(obs, exp)

    def test_get_core_features_all(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        o1_seven_num = pd.Series(
            [0.44, 1.98, 5.5, 11., 11., 11., 11.],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        o2_seven_num = pd.Series(
            [11., 11., 11., 11., 12., 12.64, 12.92],
            index=['2%', '9%', '25%', '50%', '75%', '91%', '98%'])
        exp = pd.DataFrame([o1_seven_num, o2_seven_num], index=['O1', 'O2'])

        # fraction = 2/3 - 0.006
        obs = _get_core_features(table, fraction=0.660)
        pdt.assert_frame_equal(obs, exp)

        obs = _get_core_features(table, fraction=0.0)
        pdt.assert_frame_equal(obs, exp)
