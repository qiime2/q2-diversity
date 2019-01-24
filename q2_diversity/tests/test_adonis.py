# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
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


class AdonisTests(TestPluginBase):
    package = 'q2_diversity'

    def setUp(self):
        super().setUp()
        self.adonis = self.plugin.actions['adonis']
        dm = skbio.DistanceMatrix([[0, 0.5, 1],
                                   [0.5, 0, 0.75],
                                   [1, 0.75, 0]],
                                  ids=['sample1', 'sample2', 'sample3'])
        self.dm = Artifact.import_data('DistanceMatrix', dm)

    def test_execution(self):
        # does it run?
        md = qiime2.Metadata(
            pd.DataFrame([[1, 'a'], [1, 'b'], [2, 'b']],
                         columns=['number', 'letter'],
                         index=pd.Index(['sample1', 'sample2', 'sample3'],
                                        name='id')))
        self.adonis(self.dm, md, 'letter+number')

    def test_metadata_is_superset(self):
        md = qiime2.Metadata(
            pd.DataFrame([[1, 'a'], [1, 'b'], [2, 'b'], [2, 'a']],
                         columns=['number', 'letter'],
                         index=pd.Index(['sample1', 'sample2', 'sample3', 'F'],
                                        name='id')))
        self.adonis(self.dm, md, 'letter+number')

    def test_metadata_is_subset(self):
        md = qiime2.Metadata(
            pd.DataFrame([[1, 'a'], [1, 'b'], [2, 'b']],
                         columns=['number', 'letter'],
                         index=pd.Index(['sample1', 'sample2', 'peanuts'],
                                        name='id')))
        with self.assertRaisesRegex(ValueError, "Missing samples"):
            self.adonis(self.dm, md, 'letter+number')

    def test_invalid_formula(self):
        md = qiime2.Metadata(
            pd.DataFrame([[1, 'a'], [1, 'b'], [2, 'b']],
                         columns=['number', 'letter'],
                         index=pd.Index(['sample1', 'sample2', 'sample3'],
                                        name='id')))
        with self.assertRaisesRegex(ValueError, "not a column"):
            self.adonis(self.dm, md, 'letter+fakecolumn')


if __name__ == '__main__':
    unittest.main()
