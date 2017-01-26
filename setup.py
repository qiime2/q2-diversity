# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-diversity",
    version="2017.2.0.dev0",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*', 'q2-feature-table == 2017.2.*',
                      'q2-types == 2017.2.*', 'q2templates == 2017.2.*',
                      'scikit-bio', 'seaborn', 'statsmodels', 'scipy', 'numpy',
                      'pandas', 'biom-format >= 2.1.5, < 2.2.0',
                      # `ipywidgets` included to avoid ShimWarning from
                      # `seaborn` imports:
                      #  https://github.com/mwaskom/seaborn/issues/874
                      'ipywidgets'],
    package_data={'q2_diversity._alpha': [
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
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins': ['q2-diversity=q2_diversity.plugin_setup:plugin']
    }
)
