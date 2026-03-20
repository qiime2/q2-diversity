# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
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
            [[0, 0.5, 1, 0.25, 1],
             [0.5, 0, 0.75, 0.25, 0.75],
             [1, 0.75, 0, 0.25, 0],
             [0.25, 0.25, 0.25, 0, 0.25],
             [1, 0.75, 0, 0.25, 0]],
            ids=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])

    # due to the adonis -> adonis2 migration in preparation for 2026.4 release,
    # this test data has been updated to provide a distance matrix with
    # sufficient residuals so as to not overfit the model when running adonis2
    # this confirms that we still get the same terms in our output
    # the expected DF has been validated by using the updated test data
    # and running the prior adonis implementation, providing the same results
    def test_execute_and_validate_output(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))

        exp = pd.DataFrame(
            [[2, 0.41250, 0.568965517241379, -6.6, 0.933333],
             [1, 0.34375, 0.474137931034483, -11.0, 1.0],
             [1, -0.03125, -0.043103448275862, np.nan, np.nan],
             [4, 0.72500, 1.0, np.nan, np.nan]],
            columns=['Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'],
            index=['letter', 'number', 'Residual', 'Total'])

        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

            with open(os.path.join(temp_dir_name, 'adonis.tsv'), 'r') as fh:
                res = pd.read_csv(fh, sep='\t')
                pdt.assert_frame_equal(
                    res, exp, check_dtype=False, check_frame_type=False)

    def test_adonis_handles_single_quotes_in_metadata(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a\'s'], [1, 'b\'s'], [2, 'b\'s'], [2, 'a\'s'], [3, 'c\'s']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_metadata_is_superset(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a'], [3, 'c'], [4, 'd']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5', 'sample6'], name='id')))
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
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with self.assertRaisesRegex(ValueError, "not a column"):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                adonis(temp_dir_name, self.dm, md, 'letter+fakecolumn')

    def test_metadata_index_rename(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='#SampleID')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_nans_in_formula_column(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [np.nan, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with redirected_stdio(stderr=os.devnull):
            with self.assertRaisesRegex(ValueError, "no NaN values.*`number`"):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    adonis(temp_dir_name, self.dm, md, 'letter+number')

    def test_nans_in_unused_column(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [np.nan, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with redirected_stdio(stderr=os.devnull):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                adonis(temp_dir_name, self.dm, md, 'letter+letter')


if __name__ == '__main__':
    unittest.main()
