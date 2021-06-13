# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
from skbio.util import assert_ordination_results_equal
import pandas as pd
import numpy as np

from q2_diversity import pcoa, pcoa_biplot, tsne, umap


class PCoATests(unittest.TestCase):

    def setUp(self):
        self.dm = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                        [0.3333333, 0.0000000, 0.4285714],
                                        [0.6666667, 0.4285714, 0.0000000]],
                                       ids=['S1', 'S2', 'S3'])
        self.ordination = skbio.stats.ordination.pcoa(self.dm)

    def test_pcoa(self):
        observed = pcoa(self.dm)
        skbio.util.assert_ordination_results_equal(
            observed, self.ordination, ignore_directionality=True)

    def test_pcoa_fsvd(self):
        # Run fsvd, computing all dimensions.
        fsvd_result = pcoa(self.dm,
                           number_of_dimensions=self.dm.data.shape[0])

        # Run eigh, which computes all dimensions by default.
        eigh_result = pcoa(self.dm)

        assert_ordination_results_equal(fsvd_result, eigh_result,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

    def test_pcoa_biplot(self):
        features = pd.DataFrame([[1, 0], [3, 0.1], [8, -0.4]],
                                index=['S1', 'S2', 'S3'],
                                columns=['positive', 'neither'])

        expected = skbio.stats.ordination.pcoa_biplot(self.ordination,
                                                      features)
        observed = pcoa_biplot(self.ordination, features)

        skbio.util.assert_ordination_results_equal(observed, expected,
                                                   ignore_directionality=True)


class tsneTests(unittest.TestCase):

    def setUp(self):
        self.dm = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                        [0.3333333, 0.0000000, 0.4285714],
                                        [0.6666667, 0.4285714, 0.0000000]],
                                       ids=['S1', 'S2', 'S3'])

    def test_tsne_default(self):
        observed = tsne(self.dm).samples.shape
        self.assertEqual(observed, (3, 3))

    def test_tsne_custom_number_of_dims(self):
        observed_second = tsne(self.dm, 3).samples.shape
        self.assertEqual(observed_second, (3, 3))


class umapTests(unittest.TestCase):

    def setUp(self):
        data = np.ones((5, 5))
        np.fill_diagonal(data, 0)
        self.dm = skbio.DistanceMatrix(data,
                                       ids=['S1', 'S2', 'S3', 'S4', 'S5'])

    def test_umap_default(self):
        observed = umap(self.dm).samples.shape
        self.assertEqual(observed, (5, 3))

    def test_umap_custom(self):
        observed_second = umap(self.dm, 3).samples.shape
        self.assertEqual(observed_second, (5, 3))
