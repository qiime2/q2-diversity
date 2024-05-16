# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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
from q2_diversity import filter_alpha_diversity_artifact


class TestFilterDistanceMatrix(unittest.TestCase):
    def _sorted(self, dm):
        ids = np.array(dm.ids)
        order = np.argsort(ids)
        return skbio.DistanceMatrix(dm[order][:, order], ids[order])

    def test_without_where_no_filtering(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
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
            index=pd.Index(['S1', 'S4', 'S2', 'S5'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata)

        expected = skbio.DistanceMatrix([[0, 1], [1, 0]], ['S1', 'S2'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_without_where_all_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S4', 'S5', 'S6'])

        with self.assertRaisesRegex(ValueError, 'All samples.*filtered'):
            filter_distance_matrix(dm, metadata)

    def test_without_where_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata)

        expected = skbio.DistanceMatrix([[0, 1], [1, 0]], ['S1', 'S2'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_where_no_filtering(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
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
            index=pd.Index(['S1', 'S4', 'S2', 'S5'], name='id'))
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
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, 'All samples.*filtered'):
            filter_distance_matrix(dm, metadata, where="SampleType='palm'")

    def test_with_where_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where="Subject='subject-2'")

        expected = skbio.DistanceMatrix([[0]], ['S3'])
        self.assertEqual(filtered, expected)

    def test_with_exclude_ids_filter_nothing(self):
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['tongue']},
                          index=pd.Index(['S2000'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where=None,
                                          exclude_ids=True)
        self.assertEqual(self._sorted(filtered), dm)

    def test_with_exclude_ids_filter_one(self):
        df = pd.DataFrame({'Subject': ['subject-1'],
                           'SampleType': ['tongue']},
                          index=pd.Index(['S2'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where=None,
                                          exclude_ids=True)
        expected = skbio.DistanceMatrix([[0, 2], [2, 0]], ['S1', 'S3'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_exclude_ids_filter_two(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where=None,
                                          exclude_ids=True)
        expected = skbio.DistanceMatrix([[0]], ['S3'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_exclude_ids_filter_all(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, "All samples.*filtered"):
            filter_distance_matrix(dm, metadata, where=None, exclude_ids=True)

    def test_with_exclude_ids_where_filter_nothing(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where="SampleType='toe'",
                                          exclude_ids=True)
        self.assertEqual(self._sorted(filtered), dm)

    def test_with_exclude_ids_where_filter_one(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        filtered = filter_distance_matrix(dm, metadata,
                                          where="SampleType='tongue'",
                                          exclude_ids=True)
        expected = skbio.DistanceMatrix([[0, 2], [2, 0]],
                                        ['S1', 'S3'])
        self.assertEqual(self._sorted(filtered), expected)

    def test_with_exclude_ids_where_filter_two(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['elbow', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        where = "SampleType='tongue' OR SampleType='gut'"

        filtered = filter_distance_matrix(dm, metadata,
                                          where,
                                          exclude_ids=True)
        expected = skbio.DistanceMatrix([[0]], ['S1'])
        self.assertEqual(filtered, expected)

    def test_with_exclude_ids_where_filter_all(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        dm = skbio.DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]],
                                  ['S1', 'S2', 'S3'])

        where = "SampleType='tongue' OR SampleType='gut'"

        with self.assertRaisesRegex(ValueError, "All samples.*filtered"):
            filter_distance_matrix(dm, metadata,
                                   where,
                                   exclude_ids=True)


class TestFilterAlphaDiversityArtifact(unittest.TestCase):
    def test_filter_alpha_diversity_artifact(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata)
        self.assertTrue(filtered.equals(alpha_diversity))

    def test_filter_alpha_diversity_artifact_exclude_ids_true_filter_all(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, "All samples.*filtered"):
            filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                            exclude_ids=True)

    def test_filter_alpha_diversity_artifact_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata)

        expected = pd.Series([1.0, 2.0], index=['S1', 'S2'])
        self.assertTrue(filtered.sort_values().equals(expected.sort_values()))

    def test_filter_alpha_diversity_artifact_all_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 't', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S4', 'S5', 'S6'])

        with self.assertRaisesRegex(ValueError, "All samples.*filtered"):
            filter_alpha_diversity_artifact(alpha_diversity, metadata)

    def test_filter_alpha_diversity_artifact_exclude_ids_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                                   exclude_ids=True)

        expected = pd.Series([3.0], index=['S3'])
        self.assertTrue(filtered.sort_values().equals(expected.sort_values()))

    def test_filter_alpha_diversity_artifact_exclude_ids_none_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0], index=['S3', 'S4'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                                   exclude_ids=True)

        self.assertTrue(filtered.equals(alpha_diversity))

    def test_filter_alpha_diversity_artifact_test_where_no_filtering(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(
            alpha_diversity,
            metadata,
            where="SampleType='gut' OR SampleType='tongue'"
            )

        self.assertTrue(filtered.equals(alpha_diversity))

    def test_filter_alpha_diversity_artifact_test_where_some_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                                   where="SampleType='gut'")

        expected = pd.Series([1.0], index=['S1'])
        self.assertTrue(filtered.sort_values().equals(expected.sort_values()))

    def test_filter_alpha_diversity_artifact_test_where_all_filtered(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 't', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        with self.assertRaisesRegex(ValueError, "All samples.*filtered"):
            filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                            where="SampleType='palm'")

    def test_filter_alpha_diversity_artifact_test_where_extra_ids(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2',
                                       'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut', 'tongue']},
                          index=pd.Index(['S1', 'S4', 'S2', 'S5'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0], index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(
            alpha_diversity,
            metadata,
            where="SampleType='gut' OR SampleType='tongue'"
            )

        expected = pd.Series([1.0, 2.0], index=['S1', 'S2'])
        self.assertTrue(filtered.sort_values().equals(expected.sort_values()))

    def test_filter_alpha_diversity_artifact_where_exclude_ids_filter_some(
            self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1'],
                           'SampleType': ['gut', 'tongue']},
                          index=pd.Index(['S1', 'S2'], name='id'))
        metadata = qiime2.Metadata(df)

        alpha_diversity = pd.Series([1.0, 2.0, 3.0],
                                    index=['S1', 'S2', 'S3'])

        filtered = filter_alpha_diversity_artifact(alpha_diversity, metadata,
                                                   where="SampleType='gut'",
                                                   exclude_ids=True)

        expected = pd.Series([2.0, 3.0], index=['S2', 'S3'])

        self.assertTrue(filtered.sort_values().equals(expected.sort_values()))


if __name__ == "__main__":
    unittest.main()
