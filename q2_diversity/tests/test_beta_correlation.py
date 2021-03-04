# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import skbio
import pandas as pd
import qiime2
from qiime2 import Artifact

from qiime2.plugin.testing import TestPluginBase


class BetaCorrelationTests(TestPluginBase):
    package = 'q2_diversity'

    def setUp(self):
        super().setUp()
        self.beta_correlation = self.plugin.pipelines['beta_correlation']
        dm = skbio.DistanceMatrix([[0, 1, 2],
                                   [1, 0, 1],
                                   [2, 1, 0]],
                                  ids=['sample1', 'sample2', 'sample3'])
        self.dm = Artifact.import_data('DistanceMatrix', dm)

        self.md = qiime2.NumericMetadataColumn(
            pd.Series([1, 2, 3], name='number',
                      index=pd.Index(['sample1', 'sample2', 'sample3'],
                                     name='id')))

    def test_execution(self):
        # does it run?
        self.beta_correlation(self.md, self.dm)

    def test_outputs(self):
        result = self.beta_correlation(self.md, self.dm)
        # correct number of outputs?
        self.assertEqual(2, len(result))
        # correct types?
        self.assertEqual('DistanceMatrix',
                         str(result.metadata_distance_matrix.type))
        self.assertEqual('Visualization',
                         str(result.mantel_scatter_visualization.type))


if __name__ == '__main__':
    unittest.main()
