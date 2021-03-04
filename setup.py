# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-diversity",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'q2_diversity._alpha': [
                      'alpha_group_significance_assets/index.html',
                      'alpha_group_significance_assets/dist/*',
                      'alpha_correlation_assets/index.html',
                      'alpha_correlation_assets/dist/*',
                      'alpha_rarefaction_assets/index.html',
                      'alpha_rarefaction_assets/dist/*',
                  ],
                  'q2_diversity': ['citations.bib'],
                  'q2_diversity._beta': [
                      'beta_group_significance_assets/index.html',
                      'mantel_assets/index.html',
                      'beta_rarefaction_assets/*',
                      'adonis_assets/*',
                      'bioenv_assets/index.html',
                  ],
                  'q2_diversity.tests': [
                      'data/*'
                  ]},
    scripts=['q2_diversity/_beta/adonis_assets/run_adonis.R'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Core diversity analyses.",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins': ['q2-diversity=q2_diversity.plugin_setup:plugin']
    },
    zip_safe=False,
)
