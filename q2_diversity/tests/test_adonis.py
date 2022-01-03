# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import os
import tempfile

import skbio
import pandas as pd
import qiime2
from qiime2.util import redirected_stdio
import numpy as np
from q2_diversity import adonis
import pandas.testing as pdt

from qiime2.plugin.testing import TestPluginBase


class AdonisTests(TestPluginBase):
    package = 'q2_diversity'

    def setUp(self):
        super().setUp()

        self.dm = skbio.DistanceMatrix(
            [[0, 0.5, 1], [0.5, 0, 0.75], [1, 0.75, 0]],
            ids=['sample1', 'sample2', 'sample3'])

    def test_execute_and_validate_output(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b']], columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        exp = pd.DataFrame(
            [[1.0, 0.322916667, 0.322916667, 0.0, 0.534482759, 1.0],
             [1.0, 0.281250000, 0.281250000, 0.0, 0.465517241, 1.0],
             [0.0, -1.403048e-18, -np.Infinity, np.nan, -2.322286e-18, np.nan],
             [2.0, 0.604166667, np.nan, np.nan, 1.0, np.nan]],
            columns=['Df', 'SumsOfSqs', 'MeanSqs', 'F.Model', 'R2', 'Pr(>F)'],
            index=['letter', 'number', 'Residuals', 'Total'])

        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

            with open(os.path.join(temp_dir_name, 'adonis.tsv'), 'r') as fh:
                res = pd.read_csv(fh, sep='\t')
                pdt.assert_frame_equal(
                    res, exp, check_dtype=False, check_frame_type=False)

    def test_adonis_handles_single_quotes_in_metadata(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a\'s'], [1, 'b\'s'], [2, 'b\'s'], [2, 'a\'s']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3', 'F'], name='id')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_metadata_is_superset(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3', 'F'], name='id')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_metadata_is_subset(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b']], columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'peanuts'], name='id')))
        with self.assertRaisesRegex(ValueError, "Missing samples"):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_invalid_formula(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b']], columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with self.assertRaisesRegex(ValueError, "not a column"):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                adonis(temp_dir_name, self.dm, md, 'letter+fakecolumn')

    def test_metadata_index_rename(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3', 'F'],
                           name='#SampleID')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_nans_in_formula_column(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [np.nan, 'b']], columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with redirected_stdio(stderr=os.devnull):
            with self.assertRaisesRegex(ValueError, "no NaN values.*`number`"):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_nans_in_unused_column(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [np.nan, 'b']], columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with redirected_stdio(stderr=os.devnull):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                adonis(temp_dir_name, self.dm, md, 'letter+letter')


if __name__ == '__main__':
    unittest.main()
