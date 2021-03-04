# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import numpy as np
import pandas as pd

from q2_diversity import procrustes_analysis


class PCoATests(unittest.TestCase):

    def setUp(self):
        axes = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
        eigvals = pd.Series(np.array([1.5, 0.75, 0.3, 0.15, 0.15, 0.15]),
                            index=axes)
        samples = np.array([[0, 3, 4, 4, 0, 0],
                            [1, 2, 1, 4, 3, 3],
                            [2, 3, 1, 0, 0, 1],
                            [0, 3, 2, 4, 3, 0]])

        proportion_explained = pd.Series([0.50, 0.25, 0.10, 0.05, 0.05, 0.05],
                                         index=axes)
        samples_df = pd.DataFrame(samples,
                                  index=['A', 'B', 'C', 'D'],
                                  columns=axes)
        self.reference = skbio.OrdinationResults(
                'PCoA',
                'Principal Coordinate Analysis',
                eigvals,
                samples_df,
                proportion_explained=proportion_explained)

        samples = np.array([[0.7, 3.7, 4.7, 4.7, 0.7, 0.7],
                            [1.7, 2.7, 1.7, 4.7, 3.7, 3.7],
                            [2.7, 3.7, 1.7, 0.7, 0.7, 1.7],
                            [30, 3.7, 2.7, 4.7, 3.7, 0.7]])
        samples_df = pd.DataFrame(samples,
                                  index=['A', 'B', 'C', 'D'],
                                  columns=axes)
        self.other = skbio.OrdinationResults(
                'PCoA',
                'Principal Coordinate Analysis',
                eigvals.copy(),
                samples_df.copy(),
                proportion_explained=proportion_explained.copy())

        S = [[-0.1358036, 0.0452679, 0.3621430, 0.1810715, -0.2716072],
             [0.0452679, -0.1358036, -0.1810715, 0.1810715, 0.2716072],
             [0.2263394, 0.0452679, -0.1810715, -0.5432145, -0.2716072],
             [-0.1358036, 0.0452679, 0.0000000, 0.1810715, 0.2716072]]
        samples_df = pd.DataFrame(np.array(S),
                                  index=['A', 'B', 'C', 'D'],
                                  columns=axes[:5])
        self.expected_ref = skbio.OrdinationResults(
                'PCoA',
                'Principal Coordinate Analysis',
                eigvals[:5].copy(),
                samples_df.copy(),
                proportion_explained=proportion_explained[:5].copy())
        S = [[0.0482731, -0.0324317, 0.0494312, -0.0316828, -0.1584374],
             [0.0803620, -0.0718115, -0.0112234, -0.0171011, -0.1101209],
             [0.0527554, -0.0042753, -0.0126739, -0.0969602, -0.0964822],
             [-0.1813905, 0.1085184, -0.0255339, 0.1457440, 0.3650405]]
        samples_df = pd.DataFrame(np.array(S),
                                  index=['A', 'B', 'C', 'D'],
                                  columns=axes[:5])
        self.expected_other = skbio.OrdinationResults(
                'PCoA',
                'Principal Coordinate Analysis',
                eigvals[:5].copy(),
                samples_df.copy(),
                proportion_explained=proportion_explained[:5].copy())

    def test_procrustes(self):
        ref, other = procrustes_analysis(self.reference, self.other)

        skbio.util.assert_ordination_results_equal(ref, self.expected_ref)
        skbio.util.assert_ordination_results_equal(other, self.expected_other)

    def test_procrustes_bad_dimensions(self):

        self.other.samples = self.other.samples.iloc[:, :4]
        self.other.eigvals = self.other.eigvals[:4]
        self.other.proportion_explained = self.other.proportion_explained[:4]

        with self.assertRaisesRegex(ValueError, 'The matrices cannot be '):
            procrustes_analysis(self.reference, self.other)

    def test_procrustes_over_dimensions(self):
        with self.assertRaisesRegex(ValueError, 'Cannot fit fewer dimensions '
                                    'than available'):
            procrustes_analysis(self.reference, self.other, 11)

    def test_procrustes_id_mismatch(self):
        msg = 'The ordinations represent two different sets of samples'
        self.other.samples.index = pd.Index([':L', ':D', ':)', ':('])
        with self.assertRaisesRegex(ValueError, msg):
            procrustes_analysis(self.reference, self.other)

        self.other.samples.index = pd.Index([':L', 'B', 'C', 'D'])
        with self.assertRaisesRegex(ValueError, msg):
            procrustes_analysis(self.reference, self.other)

        self.other.samples.index = pd.Index(['a', 'b', 'c', 'd'])
        with self.assertRaisesRegex(ValueError, msg):
            procrustes_analysis(self.reference, self.other)
