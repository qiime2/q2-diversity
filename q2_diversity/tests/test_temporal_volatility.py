# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import pandas as pd
import pandas.testing as pdt
import qiime2
import skbio

from q2_diversity import pcoa_centroid_temporal_volatility


class PCoACentroidVolatilityTests(unittest.TestCase):

    def _pcoa(self, samples):
        axes = samples.columns
        return skbio.OrdinationResults(
            short_method_name="PCoA",
            long_method_name="Principal Coordinate Analysis",
            eigvals=pd.Series([3.0, 2.0, 1.0][: len(axes)], index=axes),
            samples=samples,
            proportion_explained=pd.Series([0.5, 0.3, 0.2][: len(axes)], index=axes),
        )

    def _metadata(self, ids, subject, group=None):
        columns = {"subject": subject}
        if group is not None:
            columns["group"] = group
        return qiime2.Metadata(pd.DataFrame(columns, index=pd.Index(ids, name="id")))

    def setUp(self):
        self.samples = pd.DataFrame(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [4.0, 0.0, 0.0],
                [10.0, 10.0, 10.0],
            ],
            index=["a1", "a2", "a3", "a4", "a5", "b1"],
            columns=["Axis 1", "Axis 2", "Axis 3"],
        )
        self.pcoa = self._pcoa(self.samples)
        self.metadata = self._metadata(
            self.samples.index,
            subject=["A", "A", "A", "A", "A", "B"],
            group=["T1", "T1", "T1", "T2", "T2", "T1"],
        )

    def test_default_metric_is_median(self):
        # A:T1 distances to centroid (0,0,0): 0, 1, 1 -> median 1.0
        # A:T2 distances to centroid (2,0,0): 2, 2 -> median 2.0
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.metadata, "subject", "group"
        )

        expected = pd.DataFrame(
            {
                "id": ["A:T1", "A:T2"],
                "measure": [1.0, 2.0],
                "error": [np.std([0, 1, 1], ddof=1) / np.sqrt(3), 0.0],
                "group": ["T1", "T2"],
                "subject": ["A", "A"],
                "n_samples": [3, 2],
            },
            index=pd.Index(["A:T1", "A:T2"], name="id"),
        )

        pdt.assert_frame_equal(observed, expected)

    def test_metric_mean(self):
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.metadata, "subject", "group", metric="mean"
        )
        pdt.assert_series_equal(
            observed["measure"],
            pd.Series(
                [2 / 3, 2.0],
                index=pd.Index(["A:T1", "A:T2"], name="id"),
                name="measure",
            ),
        )

    def test_metric_sum(self):
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.metadata, "subject", "group", metric="sum"
        )
        pdt.assert_series_equal(
            observed["measure"],
            pd.Series(
                [2.0, 4.0], index=pd.Index(["A:T1", "A:T2"], name="id"), name="measure"
            ),
        )

    def test_empty_group_column_groups_by_subject_only(self):
        # all of subject A's samples (across both former groups) are pooled
        # into a single centroid when group_column is left empty.
        metadata = self._metadata(
            self.samples.index, subject=["A", "A", "A", "A", "A", "B"]
        )

        observed = pcoa_centroid_temporal_volatility(self.pcoa, metadata, "subject")

        self.assertEqual(list(observed["id"]), ["A:all"])
        self.assertEqual(list(observed["group"]), ["all"])

    def test_min_group_size_drops_small_groups(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0], [10.0, 10.0]],
            index=["a1", "a2", "b1"],
            columns=["Axis 1", "Axis 2"],
        )
        pcoa = self._pcoa(samples)
        metadata = self._metadata(
            samples.index, subject=["A", "A", "B"], group=["T1", "T1", "T1"]
        )

        observed = pcoa_centroid_temporal_volatility(
            pcoa, metadata, "subject", "group", dimensions=2, min_group_size=2
        )

        self.assertEqual(list(observed["id"]), ["A:T1"])

    def test_n_samples_column(self):
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.metadata, "subject", "group"
        )

        self.assertEqual(list(observed["n_samples"]), [3, 2])

    def test_min_group_size_below_two_raises(self):
        with self.assertRaisesRegex(ValueError, "min_group_size"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.metadata, "subject", "group", min_group_size=1
            )

    def test_dimensions_below_one_raises(self):
        with self.assertRaisesRegex(ValueError, "dimensions"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.metadata, "subject", "group", dimensions=0
            )

    def test_dimensions_exceeds_available_axes_raises(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0]], index=["a1", "a2"], columns=["Axis 1", "Axis 2"]
        )
        pcoa = self._pcoa(samples)
        metadata = self._metadata(samples.index, subject=["A", "A"], group=["T1", "T1"])

        with self.assertRaisesRegex(ValueError, "Cannot compute centroids"):
            pcoa_centroid_temporal_volatility(
                pcoa, metadata, "subject", "group", dimensions=3
            )

    def test_no_overlapping_samples_raises(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0]], index=["a1", "a2"], columns=["Axis 1", "Axis 2"]
        )
        pcoa = self._pcoa(samples)
        metadata = self._metadata(["x1", "x2"], subject=["A", "A"], group=["T1", "T1"])

        with self.assertRaisesRegex(ValueError, "No samples are shared"):
            pcoa_centroid_temporal_volatility(
                pcoa, metadata, "subject", "group", dimensions=2
            )

    def test_unknown_column_raises(self):
        with self.assertRaisesRegex(ValueError, "not found in the metadata"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.metadata, "not_a_column", "group"
            )
