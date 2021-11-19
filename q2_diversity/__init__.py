# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._alpha import (alpha, alpha_phylogenetic, alpha_group_significance,
                     alpha_correlation, alpha_rarefaction)
from ._beta import (beta, beta_phylogenetic, bioenv,
                    beta_group_significance, mantel, beta_rarefaction,
                    beta_correlation, adonis)
from ._ordination import pcoa, pcoa_biplot, tsne, umap
from ._procrustes import procrustes_analysis
from ._core_metrics import core_metrics_phylogenetic, core_metrics
from ._filter import filter_distance_matrix
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions


__all__ = ['beta', 'beta_phylogenetic', 'alpha', 'alpha_phylogenetic',
           'pcoa', 'tsne', 'umap', 'pcoa_biplot', 'alpha_group_significance',
           'bioenv', 'beta_group_significance', 'alpha_correlation',
           'core_metrics_phylogenetic', 'core_metrics',
           'filter_distance_matrix', 'mantel', 'alpha_rarefaction',
           'beta_rarefaction', 'procrustes_analysis', 'beta_correlation',
           'adonis',
           ]
