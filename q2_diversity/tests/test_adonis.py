# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from unittest.mock import ANY, patch
import os
import tempfile

import skbio
import pandas as pd
import qiime2
from qiime2.util import redirected_stdio
import numpy as np
from q2_diversity import adonis
import pandas.testing as pdt

from qiime2.plugin import Choices, Str
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

    # addresses https://github.com/qiime2/q2-diversity/pull/394
    def test_njobs_handled_as_integer(self):
        md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a\'s'], [1, 'b\'s'], [2, 'b\'s'], [2, 'a\'s'], [3, 'c\'s']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, distance_matrix=self.dm, metadata=md,
                   formula='letter+number', n_jobs='8')


class AdonisSignificanceByTests(TestPluginBase):
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
        self.md = qiime2.Metadata(pd.DataFrame(
            [[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a'], [3, 'c']],
            columns=['number', 'letter'],
            index=pd.Index(['sample1', 'sample2', 'sample3',
                            'sample4', 'sample5'], name='id')))

        # These five samples allow fewer permutations than are requested, so
        # vegan enumerates all of them and the p-values are deterministic.
        self.exp_terms = pd.DataFrame(
            [[2, 0.41250, 0.568965517241379, -6.6, 0.933333],
             [1, 0.34375, 0.474137931034483, -11.0, 1.0],
             [1, -0.03125, -0.043103448275862, np.nan, np.nan],
             [4, 0.72500, 1.0, np.nan, np.nan]],
            columns=['Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'],
            index=['letter', 'number', 'Residual', 'Total'])

    def _run_adonis(self, formula, **kwargs):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, self.md, formula, **kwargs)
            return pd.read_csv(os.path.join(temp_dir_name, 'adonis.tsv'),
                               sep='\t')

    def _assert_results_equal(self, obs, exp):
        pdt.assert_frame_equal(
            obs, exp, check_dtype=False, check_frame_type=False)

    def test_defaults_to_terms(self):
        self._assert_results_equal(
            self._run_adonis('letter+number'), self.exp_terms)

    def test_terms(self):
        obs = self._run_adonis('letter+number', significance_by='terms')

        self._assert_results_equal(obs, self.exp_terms)

    def test_terms_depends_on_formula_order(self):
        exp = pd.DataFrame(
            [[1, 0.56875, 0.784482758620689, -18.2, 1.0],
             [2, 0.18750, 0.258620689655172, -3.0, 0.9],
             [1, -0.03125, -0.043103448275862, np.nan, np.nan],
             [4, 0.72500, 1.0, np.nan, np.nan]],
            columns=['Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'],
            index=['number', 'letter', 'Residual', 'Total'])

        obs = self._run_adonis('number+letter', significance_by='terms')

        self._assert_results_equal(obs, exp)

    def test_margin(self):
        exp = pd.DataFrame(
            [[2, 0.18750, 0.258620689655173, -3.0, 0.9],
             [1, 0.34375, 0.474137931034482, -11.0, 1.0],
             [1, -0.03125, -0.043103448275862, np.nan, np.nan],
             [4, 0.72500, 1.0, np.nan, np.nan]],
            columns=['Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'],
            index=['letter', 'number', 'Residual', 'Total'])

        obs = self._run_adonis('letter+number', significance_by='margin')

        self._assert_results_equal(obs, exp)

    def test_margin_does_not_depend_on_formula_order(self):
        letter_first = self._run_adonis(
            'letter+number', significance_by='margin')
        number_first = self._run_adonis(
            'number+letter', significance_by='margin')

        self._assert_results_equal(
            letter_first, number_first.loc[letter_first.index])

    def test_null(self):
        exp = pd.DataFrame(
            [[3, 0.75625, 1.043103448275860, -8.066667, 1.0],
             [1, -0.03125, -0.043103448275862, np.nan, np.nan],
             [4, 0.72500, 1.0, np.nan, np.nan]],
            columns=['Df', 'SumOfSqs', 'R2', 'F', 'Pr(>F)'],
            index=['Model', 'Residual', 'Total'])

        obs = self._run_adonis('letter+number', significance_by='null')

        self._assert_results_equal(obs, exp)

    @patch('q2_diversity._beta._visualizer._run_command')
    def test_invalid_value(self, run_command):
        action = self.plugin.visualizers['adonis']
        dm = qiime2.Artifact.import_data('DistanceMatrix', self.dm)

        for value in ('bogus', 'Terms', 'NULL', ''):
            with self.subTest(value=value):
                with self.assertRaisesRegex(TypeError, 'significance_by'):
                    action(dm, self.md, 'letter+number',
                           significance_by=value)

        run_command.assert_not_called()

    @patch('q2_diversity._beta._visualizer._run_command')
    def test_command_argument_order(self, run_command):
        # run_adonis.R reads its arguments by position
        def write_results(cmd):
            self.exp_terms.to_csv(cmd[-1], sep='\t')

        run_command.side_effect = write_results

        with tempfile.TemporaryDirectory() as temp_dir_name:
            adonis(temp_dir_name, self.dm, self.md, 'letter+number',
                   significance_by='margin', permutations=99, n_jobs=2)

            run_command.assert_called_once_with(
                ['run_adonis.R', ANY, ANY, 'letter+number', 'margin', '99',
                 '2', os.path.join(temp_dir_name, 'adonis.tsv')])

    def test_plugin_registration(self):
        parameters = self.plugin.visualizers['adonis'].signature.parameters

        self.assertIn('significance_by', parameters)
        self.assertEqual(parameters['significance_by'].default, 'terms')
        self.assertEqual(parameters['significance_by'].qiime_type,
                         Str % Choices('terms', 'margin', 'onedf', 'null'))


if __name__ == '__main__':
    unittest.main()
