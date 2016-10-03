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
    version="0.0.2",
    packages=find_packages(),
    install_requires=['q2-feature-table', 'scikit-bio', 'qiime >= 2.0.2',
                      'q2-types >= 0.0.2', 'seaborn', 'statsmodels', 'scipy',
                      'numpy', 'pandas', 'biom-format >= 2.1.5, < 2.2.0',
                      'trender'],
    package_data={'q2_diversity': ['markdown/*md'],
                  'q2_diversity._alpha': ['assets/index.template',
                                          'assets/dst/*']},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Core diversity analyses.",
    license='BSD-3-Clause',
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugins': ['q2-diversity=q2_diversity.plugin_setup:plugin']
    }
)
