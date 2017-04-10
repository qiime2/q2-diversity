# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources

from ._alpha import (alpha, alpha_phylogenetic, alpha_group_significance,
                     alpha_correlation)
from ._beta import (beta, beta_phylogenetic, bioenv, beta_group_significance,
                    beta_correlation)
from ._ordination import pcoa
from ._core_metrics import core_metrics
from ._filter import filter_distance_matrix


__version__ = pkg_resources.get_distribution('q2-diversity').version

__all__ = ['beta', 'beta_phylogenetic', 'alpha', 'alpha_phylogenetic', 'pcoa',
           'alpha_group_significance', 'bioenv', 'beta_group_significance',
           'alpha_correlation', 'core_metrics', 'filter_distance_matrix',
           'beta_correlation']
