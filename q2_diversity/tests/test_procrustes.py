# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import numpy as np
import numpy.testing as npt
import pandas as pd

from q2_diversity import procrustes_analysis
from q2_diversity._procrustes import  (_deconstructed_procrustes,
                                       _partial_procrustes)


class PartialProcrustesTests(unittest.TestCase):
    def _complete_procrustes(self, mtx1, mtx2, mtx1_translate, mtx2_translate,
                             norm1, norm2, R, s):
        mtx1 = mtx1.copy()
        mtx2 = mtx2.copy()

        mtx1_translate = np.mean(mtx1, 0)
        mtx2_translate = np.mean(mtx2, 0)
        mtx1 -= mtx1_translate
        mtx2 -= mtx2_translate
        mtx1 /= norm1
        mtx2 /= norm2

        mtx2 = np.dot(mtx2, R.T) * s

        disparity = np.sum(np.square(mtx1 - mtx2))

        return mtx1, mtx2, disparity

    def test_deconstructed_procrustes(self):
        # test adapted from scipy
        # https://github.com/scipy/scipy/blob/d541c752246a9e196034957d3e044950eec75907/scipy/spatial/tests/test__procrustes.py#L37-L54
        # an L
        data1 = np.array([[1, 3], [1, 2], [1, 1], [2, 1]], 'd')

        # a larger, shifted, mirrored L
        data2 = np.array([[4, -2], [4, -4], [4, -6], [2, -6]], 'd')

        # data4, data5 are standardized (trace(A*A') = 1).
        # procrustes should return an identical copy if they are used
        # as the first matrix argument.
        shiftangle = np.pi / 8
        data4 = np.array([[1, 0], [0, 1], [-1, 0],
                          [0, -1]], 'd') / np.sqrt(4)
        data5 = np.array([[np.cos(shiftangle), np.sin(shiftangle)],
                          [np.cos(np.pi / 2 - shiftangle),
                           np.sin(np.pi / 2 - shiftangle)],
                          [-np.cos(shiftangle),
                           -np.sin(shiftangle)],
                          [-np.cos(np.pi / 2 - shiftangle),
                           -np.sin(np.pi / 2 - shiftangle)]],
                          'd') / np.sqrt(4)
        result = _deconstructed_procrustes(data1, data2)
        a, b, disparity = self._complete_procrustes(data1, data2, *result)
        npt.assert_allclose(b, a)
        npt.assert_almost_equal(disparity, 0.)

        result = _deconstructed_procrustes(data4, data5)
        m4, m5, disp45 = self._complete_procrustes(data4, data5, *result)
        npt.assert_equal(m4, data4)

    def test_partial_procrustes(self):
        # an L with an extra point on vertical
        mtx1 = np.array([[1, 4], [1, 3], [1, 2], [1, 1], [2, 1]], 'd')
        data1 = pd.DataFrame(mtx1, index=['foo', 'a1', 'b1', 'c1', 'd1'],
                             columns=[0, 1])
        data1paired = ['a1', 'b1', 'c1', 'd1']

        # a larger L, shifted, mirrored, with an extra point on the horizontal
        mtx2 = np.array([[4, -2], [4, -4], [4, -6], [2, -6], [0, -6]], 'd')
        data2 = pd.DataFrame(mtx2, index=['a2', 'b2', 'c2', 'd2', 'bar'],
                             columns=[0, 1])
        data2paired = ['a2', 'b2', 'c2', 'd2']

        df1_obs, df2_obs = _partial_procrustes(data1, data2, data1paired,
                                               data2paired)

        # we should be anchored
        npt.assert_allclose(df1_obs.loc[data1paired].values,
                            df2_obs.loc[data2paired].values)


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

        noise = [
            [0.04988341, -0.03234447, 0.03177641, -0.03507789, -0.13564394],
            [0.09117347, -0.08318546, -0.02249053, -0.01597601, -0.10901541],
            [0.05077765, -0.003994, -0.00984688, -0.09356729, -0.09648388],
            [-0.19183453, 0.11952393, 0.000561, 0.14462118, 0.34114323]]
        samples_df = pd.DataFrame(np.array(noise),
                                  index=['A', 'B', 'C', 'D'],
                                  columns=axes[:5])
        self.expected_noise = skbio.OrdinationResults(
                'PCoA',
                'Principal Coordinate Analysis',
                eigvals[:5].copy(),
                samples_df.copy(),
                proportion_explained=proportion_explained[:5].copy())

        self.expected_m2 = 0.72240956
        self.expected_p = 0.5

    def test_procrustes(self):
        ref, other, m2_results = procrustes_analysis(self.reference,
                                                     self.other)
        true_m2 = m2_results['true M^2 value'][0]
        true_p_value = m2_results['p-value for true M^2 value'][0]

        skbio.util.assert_ordination_results_equal(ref, self.expected_ref)
        skbio.util.assert_ordination_results_equal(other, self.expected_other)

        self.assertAlmostEqual(true_m2, self.expected_m2)
        self.assertNotAlmostEqual(true_p_value, self.expected_p)

    def test_non_zero_p(self):
        # generated with np.random.seed(3); np.random.randn(4, 6)
        noise = np.array(
            [[1.78862847, 0.43650985, 0.09649747, -1.8634927, -0.2773882,
              -0.35475898],
             [-0.08274148, -0.62700068, -0.04381817, -0.47721803, -1.31386475,
              0.88462238],
             [0.88131804, 1.70957306, 0.05003364, -0.40467741, -0.54535995,
              -1.54647732],
             [0.98236743, -1.10106763, -1.18504653, -0.2056499, 1.48614836,
              0.23671627]])
        self.other.samples += noise

        ref, other, m2_results = procrustes_analysis(self.reference,
                                                     self.other)

        true_m2 = m2_results['true M^2 value'][0]
        true_p_value = m2_results['p-value for true M^2 value'][0]

        skbio.util.assert_ordination_results_equal(ref, self.expected_ref)
        skbio.util.assert_ordination_results_equal(other, self.expected_noise)

        # the p value shouldn't be zero even in the presence of noise
        self.assertAlmostEqual(true_m2, 0.7388121)
        self.assertNotAlmostEqual(true_p_value, 0.001)

    def test_zero_permutations_nan_pvalue(self):
        ref, other, m2_results = procrustes_analysis(self.reference,
                                                     self.other,
                                                     permutations='disable')
        true_m2 = m2_results['true M^2 value'][0]
        true_p_value = m2_results['p-value for true M^2 value'][0]

        skbio.util.assert_ordination_results_equal(ref, self.expected_ref)
        skbio.util.assert_ordination_results_equal(other, self.expected_other)

        self.assertAlmostEqual(true_m2, self.expected_m2)
        self.assertTrue(np.isnan(true_p_value))

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
