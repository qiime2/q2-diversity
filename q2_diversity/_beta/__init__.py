# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._method import (beta_phylogenetic, beta,
                      phylogenetic_metrics, non_phylogenetic_metrics,
                      all_metrics)
from ._visualizer import (bioenv, beta_group_significance, mantel, adonis)
from ._beta_rarefaction import beta_rarefaction
from ._beta_correlation import beta_correlation

__all__ = ['beta_phylogenetic', 'beta', 'bioenv',
           'beta_group_significance', 'phylogenetic_metrics',
           'non_phylogenetic_metrics', 'all_metrics', 'mantel',
           'beta_rarefaction', 'beta_correlation', 'adonis']
