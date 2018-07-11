# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# TODO Cleanup dependencies
# import io
import unittest

import tempfile
# import biom
import skbio
# import numpy as np
import pandas as pd
# import pandas.util.testing as pdt
import qiime2
from qiime2 import Artifact

from qiime2.plugin.testing import TestPluginBase
# from qiime2 import Artifact, Metadata
# from q2_diversity._beta._beta_correlation import beta_correlation


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
        with tempfile.TemporaryDirectory() as output_dir:
            self.beta_correlation(self.md, self.dm, output_dir)

    def test_outputs(self):
        with tempfile.TemporaryDirectory() as output_dir:
            result = self.beta_correlation(self.md, self.dm, output_dir)
            # correct number of outputs?
            self.assertEqual(2, len(result))
            # correct types?
            self.assertEqual('DistanceMatrix',
                             str(result.metadata_dist_matrix.type))
            self.assertEqual('Visualization',
                             str(result.mantel_scatter_visualization.type))

    def test_conditionals(self):
        # DistanceMatrix missing sample3
        # note: tests covering duplicate sampleIDs and unequal numbers of
        # samples exist elsewhere
        dm = skbio.DistanceMatrix([[0, 1, 2],
                                   [1, 0, 1],
                                   [2, 1, 0]],
                                  ids=['sample1', 'sample2', 'sample3'])
        self.dm = Artifact.import_data('DistanceMatrix', dm)

        self.md = qiime2.NumericMetadataColumn(
            pd.Series([1, 2, 3], name='number',
                      index=pd.Index(['sample1', 'sample2', 'sample4'],
                                     name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(
                ValueError,
                'All samples in distance matrix must be present '
                'and contain data in the sample metadata. The '
                'following samples were present in the distance '
                'matrix, but were missing from the sample metadata '
                    'or had no data: sample3'):
                self.beta_correlation(self.md, self.dm, output_dir)


if __name__ == '__main__':
    unittest.main()
