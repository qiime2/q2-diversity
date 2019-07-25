# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import unittest

import io
import biom
import skbio
import qiime2
from qiime2.plugin.testing import TestPluginBase
import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from q2_diversity import (alpha, alpha_phylogenetic, alpha_correlation,
                          alpha_phylogenetic_alt, alpha_group_significance)


class AlphaTests(TestPluginBase):

    package = 'q2_diversity.tests'

    def test_alpha(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        actual = alpha(table=t, metric='observed_otus')
        # expected computed by hand
        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2},
                             name='observed_otus')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylo_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'Unknown metric'):
            alpha(table=t, metric='faith_pd')

    def test_alpha_unknown_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with self.assertRaisesRegex(ValueError, 'Unknown metric'):
            alpha(table=t, metric='not-a-metric')

    def test_alpha_empty_table(self):
        t = biom.Table(np.array([]), [], [])

        with self.assertRaisesRegex(ValueError, "empty"):
            alpha(table=t, metric='observed_otus')

    def test_alpha_phylogenetic(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        actual = alpha_phylogenetic(table=t, phylogeny=tree, metric='faith_pd')
        # expected computed with skbio.diversity.alpha_diversity
        expected = pd.Series({'S1': 0.75, 'S2': 1.0, 'S3': 1.0},
                             name='faith_pd')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylogenetic_non_phylo_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaisesRegex(ValueError, 'Unknown phylogenetic metric'):
            alpha_phylogenetic(table=t, phylogeny=tree,
                               metric='observed_otus')

    def test_alpha_phylogenetic_unknown_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaisesRegex(ValueError, 'Unknown phylogenetic metric'):
            alpha_phylogenetic(table=t, phylogeny=tree, metric='not-a-metric')

    def test_alpha_phylogenetic_skbio_error_rewriting(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25):0.25, O3:0.75)root;'))
        # Verify through regex that there is a ``feature_ids`` substring
        # followed by a ``phylogeny``
        with self.assertRaisesRegex(skbio.tree.MissingNodeError,
                                    'feature_ids.*phylogeny'):
            alpha_phylogenetic(table=t, phylogeny=tree, metric='faith_pd')

    def test_alpha_phylogenetic_empty_table(self):
        t = biom.Table(np.array([]), [], [])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25):0.25, O3:0.75)root;'))

        with self.assertRaisesRegex(ValueError, "empty"):
            alpha_phylogenetic(table=t, phylogeny=tree, metric='faith_pd')

    def test_alpha_phylogenetic_alt(self):
        table = self.get_data_path('two_feature_table.biom')
        tree = self.get_data_path('three_feature.tree')
        actual = alpha_phylogenetic_alt(table=table,
                                        phylogeny=tree,
                                        metric='faith_pd')
        # expected computed with skbio.diversity.alpha_diversity
        expected = pd.Series({'S1': 0.75, 'S2': 1.0, 'S3': 1.0},
                             name='faith_pd')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylogenetic_alt_non_phylo_metric(self):
        table = self.get_data_path('two_feature_table.biom')
        tree = self.get_data_path('three_feature.tree')
        with self.assertRaisesRegex(ValueError, 'Unknown phylogenetic metric'):
            alpha_phylogenetic_alt(table=table,
                                   phylogeny=tree,
                                   metric='observed_otus')

    def test_alpha_phylogenetic_alt_unknown_metric(self):
        table = self.get_data_path('two_feature_table.biom')
        tree = self.get_data_path('three_feature.tree')
        with self.assertRaisesRegex(ValueError, 'Unknown phylogenetic metric'):
            alpha_phylogenetic_alt(table=table,
                                   phylogeny=tree,
                                   metric='not-a-metric')

    def test_alpha_phylogenetic_alt_skbio_error_rewriting(self):
        table = self.get_data_path('two_feature_table.biom')
        tree = self.get_data_path('vaw.nwk')
        with self.assertRaisesRegex(ValueError, "The table does not "
                                    "appear to be completely represented "
                                    "by the phylogeny."):
            alpha_phylogenetic_alt(table=table,
                                   phylogeny=tree,
                                   metric='faith_pd')

    def test_alpha_phylogenetic_alt_empty_table(self):
        table = self.get_data_path('empty.biom')
        tree = self.get_data_path('three_feature.tree')

        with self.assertRaisesRegex(ValueError, "empty"):
            alpha_phylogenetic_alt(table=table,
                                   phylogeny=tree,
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

    def test_extra_alpha_div(self):
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
