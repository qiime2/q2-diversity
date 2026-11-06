# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from unittest.mock import patch
import os
import tempfile

import skbio
import pandas as pd
import qiime2
from qiime2.util import redirected_stdio
import numpy as np
from q2_diversity import adonis, blocked_adonis
from q2_diversity._beta._adonis import _permute_target_column
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
            [[1, 'a\'s'], [1, 'b\'s'], [2, 'b\'s'],
             [2, 'a\'s'], [3, 'c\'s']],
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
            with self.assertRaisesRegex(ValueError, "no NaN values.*number"):
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

    # addresses https://github.com/qiime2/q2-diversity/pull/394
    def test_njobs_handled_as_integer(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a\'s'], [1, 'b\'s'], [2, 'b\'s'],
             [2, 'a\'s'], [3, 'c\'s']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, distance_matrix=self.dm, metadata=md,
                   formula='letter+number', n_jobs='8')

    def _unit_permutation_test_data(self):
        ids = ['sample1', 'sample2', 'sample3',
               'sample4', 'sample5', 'sample6']

        dm = skbio.DistanceMatrix(
            [[0, 0.2, 0.3, 0.7, 0.8, 0.9],
             [0.2, 0, 0.25, 0.75, 0.85, 0.95],
             [0.3, 0.25, 0, 0.65, 0.7, 0.8],
             [0.7, 0.75, 0.65, 0, 0.2, 0.3],
             [0.8, 0.85, 0.7, 0.2, 0, 0.25],
             [0.9, 0.95, 0.8, 0.3, 0.25, 0]],
            ids=ids)

        md = qiime2.Metadata(pd.DataFrame(
            [['a', 'unit1'],
             ['a', 'unit1'],
             ['a', 'unit2'],
             ['b', 'unit3'],
             ['b', 'unit3'],
             ['b', 'unit4']],
            columns=['treatment', 'subject_id'],
            index=pd.Index(ids, name='id')))

        return dm, md

    def test_permute_target_column_preserves_shared_unit_labels(self):
        _, md = self._unit_permutation_test_data()
        expected = md.to_dataframe()
        expected.loc[:, 'treatment'] = ['a', 'a', 'b', 'a', 'a', 'b']

        with patch.object(pd.Series, 'sample', autospec=True) as sample:
            sample.side_effect = lambda labels, **kwargs: labels.iloc[
                [0, 2, 1, 3]]
            shuffled = _permute_target_column(md, 'treatment', 'subject_id')

        pdt.assert_frame_equal(shuffled.to_dataframe(), expected)

    def test_blocked_adonis_with_permutation_group_column(self):
        dm, md = self._unit_permutation_test_data()

        distance_matrix = qiime2.Artifact.import_data('DistanceMatrix', dm)
        visualization, = self.plugin.visualizers['blocked_adonis'](
            distance_matrix=distance_matrix, metadata=md, formula='treatment',
            permutations=19, permutation_group_column='subject_id')

        with tempfile.TemporaryDirectory() as temp_dir_name:
            visualization.export_data(temp_dir_name)

            res = pd.read_csv(
                os.path.join(temp_dir_name, 'adonis.tsv'),
                sep='\t',
                index_col=0)

            self.assertIn('treatment', res.index)
            self.assertIn('Pr(>F)', res.columns)
            self.assertGreaterEqual(res.loc['treatment', 'Pr(>F)'], 0)
            self.assertLessEqual(res.loc['treatment', 'Pr(>F)'], 1)

    def test_blocked_adonis_counts_exact_and_numerical_ties(self):
        dm, md = self._unit_permutation_test_data()
        f_statistics = [4.0, 1.0, 4.0, 4.0 - 1e-14, 6.0]
        results = [
            pd.DataFrame({'F': [f_stat], 'Pr(>F)': [np.nan]},
                         index=['treatment'])
            for f_stat in f_statistics
        ]

        with patch('q2_diversity._beta._adonis._run_adonis',
                   side_effect=results):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                blocked_adonis(temp_dir_name, dm, md, 'treatment',
                               permutations=4,
                               permutation_group_column='subject_id')
                res = pd.read_csv(
                    os.path.join(temp_dir_name, 'adonis.tsv'),
                    sep='\t', index_col=0)

        self.assertEqual(res.loc['treatment', 'Pr(>F)'], 0.8)

    def test_unit_permutations_preserve_observed_statistics(self):
        dm, md = self._unit_permutation_test_data()

        with tempfile.TemporaryDirectory() as old_dir:
            adonis(old_dir, dm, md, 'treatment', permutations=19)
            old = pd.read_csv(
                os.path.join(old_dir, 'adonis.tsv'),
                sep='\t',
                index_col=0)

        with tempfile.TemporaryDirectory() as unit_dir:
            blocked_adonis(unit_dir, dm, md, 'treatment',
                           permutations=19,
                           permutation_group_column='subject_id')
            unit = pd.read_csv(
                os.path.join(unit_dir, 'adonis.tsv'),
                sep='\t',
                index_col=0)

        observed_cols = ['Df', 'SumOfSqs', 'R2', 'F']

        pdt.assert_series_equal(
            old.loc['treatment', observed_cols],
            unit.loc['treatment', observed_cols],
            check_dtype=False,
            check_names=False,
            check_exact=False)

    def test_permutation_group_column_must_exist(self):
        dm, md = self._unit_permutation_test_data()

        with self.assertRaisesRegex(
                ValueError, 'fake_subject_id.*not a column'):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                blocked_adonis(temp_dir_name, dm, md, 'treatment',
                               permutations=19,
                               permutation_group_column='fake_subject_id')

    def test_unit_permutations_require_single_formula_column(self):
        dm, md = self._unit_permutation_test_data()

        with self.assertRaisesRegex(ValueError, 'exactly one metadata column'):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                blocked_adonis(temp_dir_name, dm, md, 'treatment+subject_id',
                               permutations=19,
                               permutation_group_column='subject_id')

    def test_unit_permutations_require_constant_label_within_unit(self):
        ids = ['sample1', 'sample2', 'sample3',
               'sample4', 'sample5', 'sample6']

        dm = skbio.DistanceMatrix(
            [[0, 0.2, 0.3, 0.7, 0.8, 0.9],
             [0.2, 0, 0.25, 0.75, 0.85, 0.95],
             [0.3, 0.25, 0, 0.65, 0.7, 0.8],
             [0.7, 0.75, 0.65, 0, 0.2, 0.3],
             [0.8, 0.85, 0.7, 0.2, 0, 0.25],
             [0.9, 0.95, 0.8, 0.3, 0.25, 0]],
            ids=ids)

        md = qiime2.Metadata(pd.DataFrame(
            [['a', 'unit1'],
             ['b', 'unit1'],
             ['a', 'unit2'],
             ['b', 'unit3'],
             ['b', 'unit3'],
             ['b', 'unit4']],
            columns=['treatment', 'subject_id'],
            index=pd.Index(ids, name='id')))

        with self.assertRaisesRegex(ValueError, 'constant within'):
            with tempfile.TemporaryDirectory() as temp_dir_name:
                blocked_adonis(temp_dir_name, dm, md, 'treatment',
                               permutations=19,
                               permutation_group_column='subject_id')


if __name__ == '__main__':
    unittest.main()
