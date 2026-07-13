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

    def _group(self, ids, values):
        return qiime2.CategoricalMetadataColumn(
            pd.Series(values, index=pd.Index(ids, name="id"), name="group")
        )

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
        # a1-a3 -> T1 (centroid (0,0,0)), a4-a5 -> T2 (centroid (2,0,0)),
        # b1 -> T1 (its own centroid, dropped by min_group_size).
        self.group = self._group(
            self.samples.index,
            ["T1", "T1", "T1", "T2", "T2", "T1"],
        )

    def test_returns_metadata(self):
        observed = pcoa_centroid_temporal_volatility(self.pcoa, self.group)
        self.assertIsInstance(observed, qiime2.Metadata)

    def test_default_metric_is_median(self):
        # T1 pools a1-a3 and b1 (grouping is only by the group column now),
        # T2 pools a4-a5.
        observed = pcoa_centroid_temporal_volatility(self.pcoa, self.group)
        observed_df = observed.to_dataframe()

        t1 = self.samples.loc[["a1", "a2", "a3", "b1"]]
        t1_centroid = t1.mean()
        t1_distances = np.sqrt(((t1 - t1_centroid) ** 2).sum(axis=1))

        t2 = self.samples.loc[["a4", "a5"]]
        t2_centroid = t2.mean()
        t2_distances = np.sqrt(((t2 - t2_centroid) ** 2).sum(axis=1))

        expected = pd.DataFrame(
            {
                "measure": [t1_distances.median(), t2_distances.median()],
                "error": [
                    t1_distances.std(ddof=1) / np.sqrt(4),
                    t2_distances.std(ddof=1) / np.sqrt(2),
                ],
                "group": ["T1", "T2"],
                "n_samples": [4.0, 2.0],
            },
            index=pd.Index(["T1", "T2"], name="id"),
        )

        pdt.assert_frame_equal(observed_df, expected)

    def test_metric_mean(self):
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.group, metric="mean"
        )
        observed_df = observed.to_dataframe()

        t1 = self.samples.loc[["a1", "a2", "a3", "b1"]]
        t1_distances = np.sqrt(((t1 - t1.mean()) ** 2).sum(axis=1))
        t2 = self.samples.loc[["a4", "a5"]]
        t2_distances = np.sqrt(((t2 - t2.mean()) ** 2).sum(axis=1))

        pdt.assert_series_equal(
            observed_df["measure"],
            pd.Series(
                [t1_distances.mean(), t2_distances.mean()],
                index=pd.Index(["T1", "T2"], name="id"),
                name="measure",
            ),
        )

    def test_metric_sum(self):
        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, self.group, metric="sum"
        )
        observed_df = observed.to_dataframe()

        t1 = self.samples.loc[["a1", "a2", "a3", "b1"]]
        t1_distances = np.sqrt(((t1 - t1.mean()) ** 2).sum(axis=1))
        t2 = self.samples.loc[["a4", "a5"]]
        t2_distances = np.sqrt(((t2 - t2.mean()) ** 2).sum(axis=1))

        pdt.assert_series_equal(
            observed_df["measure"],
            pd.Series(
                [t1_distances.sum(), t2_distances.sum()],
                index=pd.Index(["T1", "T2"], name="id"),
                name="measure",
            ),
        )

    def test_missing_values_are_dropped(self):
        group = self._group(
            self.samples.index, ["T1", "T1", "T1", "T2", "T2", None]
        )

        observed = pcoa_centroid_temporal_volatility(
            self.pcoa, group, min_group_size=2
        )
        observed_df = observed.to_dataframe()

        self.assertEqual(list(observed_df.index), ["T1", "T2"])
        self.assertEqual(list(observed_df["n_samples"]), [3.0, 2.0])

    def test_min_group_size_drops_small_groups(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0], [10.0, 10.0]],
            index=["a1", "a2", "b1"],
            columns=["Axis 1", "Axis 2"],
        )
        pcoa = self._pcoa(samples)
        group = self._group(samples.index, ["T1", "T1", "T2"])

        observed = pcoa_centroid_temporal_volatility(
            pcoa, group, dimensions=2, min_group_size=2
        )
        observed_df = observed.to_dataframe()

        self.assertEqual(list(observed_df.index), ["T1"])

    def test_n_samples_column(self):
        observed = pcoa_centroid_temporal_volatility(self.pcoa, self.group)
        observed_df = observed.to_dataframe()

        self.assertEqual(list(observed_df["n_samples"]), [4.0, 2.0])

    def test_min_group_size_below_two_raises(self):
        with self.assertRaisesRegex(ValueError, "min_group_size"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.group, min_group_size=1
            )

    def test_dimensions_below_one_raises(self):
        with self.assertRaisesRegex(ValueError, "dimensions"):
            pcoa_centroid_temporal_volatility(self.pcoa, self.group, dimensions=0)

    def test_dimensions_exceeds_available_axes_raises(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0]], index=["a1", "a2"], columns=["Axis 1", "Axis 2"]
        )
        pcoa = self._pcoa(samples)
        group = self._group(samples.index, ["T1", "T1"])

        with self.assertRaisesRegex(ValueError, "Cannot compute centroids"):
            pcoa_centroid_temporal_volatility(pcoa, group, dimensions=3)

    def test_no_overlapping_samples_raises(self):
        samples = pd.DataFrame(
            [[0.0, 0.0], [1.0, 0.0]], index=["a1", "a2"], columns=["Axis 1", "Axis 2"]
        )
        pcoa = self._pcoa(samples)
        group = self._group(["x1", "x2"], ["T1", "T1"])

        with self.assertRaisesRegex(ValueError, "No samples are shared"):
            pcoa_centroid_temporal_volatility(pcoa, group, dimensions=2)

    def test_min_group_size_too_large_raises(self):
        with self.assertRaisesRegex(ValueError, "no volatility scores"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.group, min_group_size=10
            )

    def test_invalid_metric_raises(self):
        with self.assertRaisesRegex(ValueError, "metric must be one of"):
            pcoa_centroid_temporal_volatility(
                self.pcoa, self.group, metric="stdev"
            )
