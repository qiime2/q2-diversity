# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import pandas as pd
import qiime2
import skbio

from q2_diversity import filter_distance_matrix


class TestFilterDistanceMatrix(unittest.TestCase):
    def _sorted(self, dm):
        ids = np.array(dm.ids)
        order = np.argsort(ids)
        return skbio.DistanceMatrix(dm[order][:, order], ids[order])

    def test_without_where_no_filtering(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata)

        # There is no guaranteed order to output. Sort by ID before comparing
        # to expected.
        self.assertEqual(self._sorted(filtered), dm)

    def test_without_where_extra_ids(self):
        df = pd.DataFrame(
            {'Subject': ['subject-1', 'subject-1', 'subject-2', 'subject-2'],
             'SampleType': ['gut', 'tongue', 'gut', 'tongue']},
            index=['S1', 'S4', 'S2', 'S5'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata)

        expected = skbio.DistanceMatrix([[0, 1], [1, 0]], ['S1', 'S2'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_without_where_all_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S4', 'S5', 'S6'])

        with self.assertRaisesRegex(ValueError, 'All samples.*filtered'):
            filter_distance_matrix(dm, metadata)

    def test_without_where_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=['S1', 'S2'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata)

        expected = skbio.DistanceMatrix([[0, 1], [1, 0]], ['S1', 'S2'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_where_no_filtering(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(
            dm, metadata, where="SampleType='gut' OR SampleType='tongue'")

        self.assertEqual(self._sorted(filtered), dm)

    def test_with_where_extra_ids(self):
        df = pd.DataFrame(
            {'Subject': ['subject-1', 'subject-1', 'subject-2', 'subject-2'],
             'SampleType': ['gut', 'tongue', 'gut', 'tongue']},
            index=['S1', 'S4', 'S2', 'S5'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(
            dm, metadata, where="SampleType='gut' OR SampleType='tongue'")

        expected = skbio.DistanceMatrix([[0, 1], [1, 0]], ['S1', 'S2'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_where_all_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, 'All samples.*filtered'):
            filter_distance_matrix(dm, metadata, where="SampleType='palm'")

    def test_with_where_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=['S1', 'S2', 'S3'])
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where="Subject='subject-2'")

        expected = skbio.DistanceMatrix([[0]], ['S3'])
        self.assertEqual(filtered, expected)


if __name__ == "__main__":
    unittest.main()
