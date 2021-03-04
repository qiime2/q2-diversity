# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import io
import unittest

import biom
import skbio
import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact, Metadata


class CoreMetricsTests(TestPluginBase):
    package = 'q2_diversity'

    def setUp(self):
        super().setUp()
        self.core_metrics = self.plugin.pipelines['core_metrics']
        self.core_metrics_phylogenetic = self.plugin.pipelines[
            'core_metrics_phylogenetic']

    def test_core_metrics_phylogenetic(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        table = Artifact.import_data('FeatureTable[Frequency]', table)

        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        metadata = Metadata(
            pd.DataFrame({'foo': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        results = self.core_metrics_phylogenetic(table, tree, 13, metadata)

        self.assertEqual(len(results), 17)

        self.assertEqual(repr(results.bray_curtis_distance_matrix.type),
                         'DistanceMatrix')
        self.assertEqual(repr(results.jaccard_emperor.type), 'Visualization')

        # pipelines preserve the output's type, in this case, beta_phylogenetic
        # returns this type, and that is passed through to the final output
        # (as long as the type is a subtype of the signature).
        self.assertEqual(
            repr(results.faith_pd_vector.type),
            "SampleData[AlphaDiversity]")

        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2},
                             name='observed_features')
        pdt.assert_series_equal(results[2].view(pd.Series), expected)

    def test_core_metrics_phylogenetic_multiple_jobs(self):
        table = biom.Table(np.array([[0, 11, 11], [13, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        table = Artifact.import_data('FeatureTable[Frequency]', table)

        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        metadata = Metadata(
            pd.DataFrame({'foo': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        results = self.core_metrics_phylogenetic(table, tree, 13, metadata,
                                                 n_jobs_or_threads=2)

        self.assertEqual(len(results), 17)

        self.assertEqual(repr(results.bray_curtis_distance_matrix.type),
                         'DistanceMatrix')
        self.assertEqual(repr(results.jaccard_emperor.type), 'Visualization')

        # pipelines preserve the output's type, in this case, beta_phylogenetic
        # returns this type, and that is passed through to the final output
        # (as long as the type is a subtype of the signature).
        self.assertEqual(
            repr(results.faith_pd_vector.type),
            "SampleData[AlphaDiversity]")

        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2},
                             name='observed_features')
        pdt.assert_series_equal(results[2].view(pd.Series), expected)

    def test_core_metrics_phylogenetic_rarefy_drops_sample(self):
        table = biom.Table(np.array([[0, 11, 11], [12, 11, 11]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        table = Artifact.import_data('FeatureTable[Frequency]', table)

        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        tree = Artifact.import_data('Phylogeny[Rooted]', tree)

        metadata = Metadata(
            pd.DataFrame({'foo': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        results = self.core_metrics_phylogenetic(table, tree, 13, metadata)

        self.assertEqual(len(results), 17)

        expected = pd.Series({'S2': 2, 'S3': 2},
                             name='observed_features')
        pdt.assert_series_equal(results[2].view(pd.Series), expected)

    def test_core_metrics(self):
        # NOTE: this test uses a table and sampling depth that produce
        # deterministic values
        table = biom.Table(np.array([[150, 100, 100], [50, 100, 100]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        table = Artifact.import_data('FeatureTable[Frequency]', table)

        metadata = Metadata(
            pd.DataFrame({'foo': ['1', '2', '3']},
                         index=pd.Index(['S1', 'S2', 'S3'], name='id')))

        results = self.core_metrics(table=table, sampling_depth=200,
                                    metadata=metadata)

        self.assertEqual(len(results), 10)
        self.assertEqual(repr(results.bray_curtis_distance_matrix.type),
                         'DistanceMatrix')
        self.assertEqual(repr(results.jaccard_emperor.type), 'Visualization')

        obs_feat_exp = pd.Series({'S1': 2, 'S2': 2, 'S3': 2},
                                 name='observed_features')
        shannon_exp = pd.Series({'S1': 0.811278124459, 'S2': 1., 'S3': 1.},
                                name='shannon_entropy')
        pdt.assert_series_equal(results[1].view(pd.Series), obs_feat_exp)
        pdt.assert_series_equal(results[2].view(pd.Series), shannon_exp)


if __name__ == '__main__':
    unittest.main()
