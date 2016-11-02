# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-diversity",
    # TODO stop duplicating version string
    version="0.0.6.dev",
    packages=find_packages(),
    install_requires=['qiime >= 2.0.6', 'q2-feature-table >= 0.0.6',
                      'q2-types >= 0.0.6', 'scikit-bio', 'seaborn',
                      'statsmodels', 'scipy', 'numpy', 'pandas',
                      'biom-format >= 2.1.5, < 2.2.0', 'q2templates >= 0.0.6'],
    package_data={'q2_diversity': ['markdown/*md'],
                  'q2_diversity._alpha': [
                      'alpha_group_significance_assets/index.html',
                      'alpha_group_significance_assets/dst/*',
                      'alpha_correlation_assets/index.html',
                      'alpha_correlation_assets/dst/*'],
                  'q2_diversity._beta': [
                      'beta_group_significance_assets/index.html',
                      'bioenv_assets/index.html',
                   ]},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Core diversity analyses.",
    license='BSD-3-Clause',
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugins': ['q2-diversity=q2_diversity.plugin_setup:plugin']
    }
)
