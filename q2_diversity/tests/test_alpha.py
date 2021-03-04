# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import io
import os
import tempfile
import unittest

import biom
import skbio
import qiime2
from qiime2.plugin.testing import TestPluginBase
import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from qiime2 import Artifact
from q2_diversity import (alpha_correlation, alpha_group_significance)


class AlphaTests(TestPluginBase):

    package = 'q2_diversity.tests'

    def setUp(self):
        super().setUp()
        self.alpha = self.plugin.pipelines['alpha']
        self.alpha_phylogenetic = self.plugin.pipelines['alpha_phylogenetic']

        empty_table = self.get_data_path('empty.biom')
        self.empty_table = Artifact.import_data('FeatureTable[Frequency]',
                                                empty_table)

        two_feature_table = self.get_data_path('two_feature_table.biom')
        self.two_feature_table = Artifact.import_data(
                'FeatureTable[Frequency]',
                two_feature_table)

        three_feature_tree = self.get_data_path('three_feature.tree')
        self.three_feature_tree = Artifact.import_data('Phylogeny[Rooted]',
                                                       three_feature_tree)

        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        self.t = Artifact.import_data('FeatureTable[Frequency]', t)
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.tree = Artifact.import_data('Phylogeny[Rooted]', tree)

    def test_alpha(self):
        actual = self.alpha(table=self.t, metric='observed_features')
        actual = actual[0].view(pd.Series)
        # expected computed by hand
        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2},
                             name='observed_features')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_with_passthrough_metric(self):
        actual = self.alpha(table=self.t, metric='singles')
        actual = actual[0].view(pd.Series)
        # expected computed by hand
        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 0},
                             name='singles')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylo_metric(self):
        with self.assertRaisesRegex(TypeError, 'faith_pd.*incompatible'):
            self.alpha(table=self.t, metric='faith_pd')

    def test_alpha_unknown_metric(self):
        with self.assertRaisesRegex(TypeError, 'not-a-metric.*incompatible'):
            self.alpha(table=self.t, metric='not-a-metric')

    def test_alpha_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha(table=self.empty_table, metric='observed_features')

    def test_alpha_phylogenetic(self):
        actual = self.alpha_phylogenetic(table=self.two_feature_table,
                                         phylogeny=self.three_feature_tree,
                                         metric='faith_pd')
        actual = actual[0].view(pd.Series)
        # expected computed with skbio.diversity.alpha_diversity
        expected = pd.Series({'S1': 0.75, 'S2': 1.0, 'S3': 1.0},
                             name='faith_pd')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylogenetic_non_phylo_metric(self):
        with self.assertRaisesRegex(TypeError,
                                    'observed_features.*incompatible'):
            self.alpha_phylogenetic(table=self.two_feature_table,
                                    phylogeny=self.three_feature_tree,
                                    metric='observed_features')

    def test_alpha_phylogenetic_unknown_metric(self):
        with self.assertRaisesRegex(TypeError, 'not-a-metric.*incompatible'):
            self.alpha_phylogenetic(table=self.two_feature_table,
                                    phylogeny=self.three_feature_tree,
                                    metric='not-a-metric')

    def test_alpha_phylogenetic_empty_table(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            self.alpha_phylogenetic(table=self.empty_table,
                                    phylogeny=self.three_feature_tree,
                                    metric='faith_pd')


class AlphaCorrelationTests(unittest.TestCase):

    def test_spearman(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            jsonp_fp = os.path.join(output_dir, 'column-value.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))

            with open(jsonp_fp) as jsonp_fh:
                jsonp_content = jsonp_fh.read()
            self.assertTrue('Spearman' in jsonp_content)
            self.assertTrue('"sampleSize": 3' in jsonp_content)
            self.assertTrue('"data":' in jsonp_content)
            self.assertFalse('filtered' in jsonp_content)

    def test_pearson(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md, method='pearson')
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            jsonp_fp = os.path.join(output_dir, 'column-value.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))

            with open(jsonp_fp) as jsonp_fh:
                jsonp_content = jsonp_fh.read()
            self.assertTrue('Pearson' in jsonp_content)
            self.assertTrue('"sampleSize": 3' in jsonp_content)
            self.assertTrue('"data":' in jsonp_content)
            self.assertFalse('filtered' in jsonp_content)

    def test_bad_method(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaises(ValueError):
                alpha_correlation(output_dir, alpha_div, md, method='bad!')

    def test_non_numeric_metadata(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'col1': [4, 5, 6],
                 'col2': ['a', 'b', 'c']},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-col1.jsonp')))
            self.assertFalse(os.path.exists(
                             os.path.join(output_dir,
                                          'column-col2.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('contain numeric data' in index_content)
            self.assertTrue('<strong>col2' in index_content)

    def test_nan_metadata(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, np.nan]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            jsonp_fp = os.path.join(output_dir, 'column-value.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))
            with open(jsonp_fp) as jsonp_fh:
                jsonp_content = jsonp_fh.read()
            self.assertTrue('"filtered": 2' in jsonp_content)
            self.assertTrue('"initial": 3' in jsonp_content)

    def test_extra_metadata(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0, 4.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3',
                                'sample4'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            jsonp_fp = os.path.join(output_dir, 'column-value.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))

            with open(jsonp_fp) as jsonp_fh:
                self.assertTrue('"sampleSize": 3' in jsonp_fh.read())

    def test_extra_alpha_div_no_intersect(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0, 8.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3',
                                     'sample4'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(ValueError,
                                        'not present.*metadata.*sample4'):
                alpha_correlation(output_dir, alpha_div, md)

    def test_extra_alpha_div_intersect(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0, 8.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3',
                                     'sample4'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'value': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            alpha_correlation(output_dir, alpha_div, md, intersect_ids=True)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            jsonp_fp = os.path.join(output_dir, 'column-value.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))

    def test_all_metadata_columns_filtered(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        # Non-numeric and empty columns are filtered.
        md = qiime2.Metadata(
            pd.DataFrame(
                {'col1': ['a', 'b', 'a'],
                 'col2': [np.nan, np.nan, np.nan]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaisesRegex(
                    ValueError, 'contains only non-numeric or empty columns'):
                alpha_correlation(output_dir, alpha_div, md)


class AlphaGroupSignificanceTests(unittest.TestCase):

    def test_alpha_group_significance(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a or b': ['a', 'b', 'b']},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-a%20or%20b.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('Kruskal-Wallis (all groups)' in index_content)
            self.assertTrue('Kruskal-Wallis (pairwise)' in index_content)

    def test_alpha_group_significance_some_numeric(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a or b': ['a', 'b', 'b'],
                 'bad': [1.0, 2.0, 3.0]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-a%20or%20b.jsonp')))
            self.assertFalse(os.path.exists(
                             os.path.join(output_dir,
                                          'column-bad.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('contain categorical data:' in index_content)
            self.assertTrue('<strong>bad' in index_content)

    def test_alpha_group_significance_one_group_all_unique_values(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a or b': ['a', 'b', 'b'],
                 'bad': ['x', 'y', 'z']},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-a%20or%20b.jsonp')))
            self.assertFalse(os.path.exists(
                             os.path.join(output_dir,
                                          'column-bad.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('number of samples' in index_content)
            self.assertTrue('<strong>bad' in index_content)

    def test_alpha_group_significance_one_group_single_value(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a or b': ['a', 'b', 'b'],
                 'bad': ['x', 'x', 'x']},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-a%20or%20b.jsonp')))
            self.assertFalse(os.path.exists(
                             os.path.join(output_dir,
                                          'column-bad.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('single group' in index_content)
            self.assertTrue('<strong>bad' in index_content)

    def test_alpha_group_significance_KW_value_error(self):
        alpha_div = pd.Series([2.0, 2.0, 3.0, 2.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3',
                                     'sample4'])
        md = qiime2.Metadata(
            pd.DataFrame({'x': ['a', 'b', 'b', 'c']},
                         index=pd.Index(['sample1', 'sample2', 'sample3',
                                         'sample4'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir,
                                         'column-x.jsonp')))
            with open(index_fp) as index_fh:
                index_content = index_fh.read()
            self.assertTrue('pairwise group comparisons have been omitted'
                            in index_content)
            self.assertTrue('x:c (n=1) vs x:a (n=1)' in index_content)

    def test_alpha_group_significance_numeric_only(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'col1': [1, 2, 1],
                 'col2': [4.2, 4.2, 4.3]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            err_msg = ("does not contain any columns that satisfy this "
                       "visualizer's requirements")
            with self.assertRaisesRegex(ValueError, err_msg):
                alpha_group_significance(output_dir, alpha_div, md)

    def test_alpha_group_significance_single_quote(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a or b': ['a', "b'", "b'"]},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            with open(index_fp) as index_fh:
                self.assertTrue("\'" in index_fh.read())

    def test_alpha_group_significance_forward_slash_in_metadata_col(self):
        alpha_div = pd.Series([2.0, 4.0, 6.0], name='alpha-div',
                              index=['sample1', 'sample2', 'sample3'])
        md = qiime2.Metadata(
            pd.DataFrame(
                {'a/b': ['a', 'b', 'b']},
                index=pd.Index(['sample1', 'sample2', 'sample3'], name='id')))

        with tempfile.TemporaryDirectory() as output_dir:
            alpha_group_significance(output_dir, alpha_div, md)
            index_fp = os.path.join(output_dir, 'index.html')
            with open(index_fp) as index_fh:
                self.assertTrue("/" in index_fh.read())
            jsonp_fp = os.path.join(output_dir, 'column-a%2Fb.jsonp')
            self.assertTrue(os.path.exists(jsonp_fp))
            csv_fp = os.path.join(output_dir,
                                  'kruskal-wallis-pairwise-a%2Fb.csv')
            self.assertTrue(os.path.exists(csv_fp))


if __name__ == '__main__':
    unittest.main()
