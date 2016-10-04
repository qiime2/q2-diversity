# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import (Plugin, Str, Properties, MetadataCategory, Choices,
                          Metadata)

import q2_diversity
import q2_diversity._alpha as alpha
import q2_diversity._beta as beta
from q2_types import (FeatureTable, Frequency, DistanceMatrix, AlphaDiversity,
                      PCoAResults, SampleData, Phylogeny, Rooted)


plugin = Plugin(
    name='diversity',
    version=q2_diversity.__version__,
    website='https://github.com/qiime2/q2-diversity',
    package='q2_diversity'
)

plugin.methods.register_function(
    function=q2_diversity.beta_phylogenetic,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta.phylogenetic_metrics())},
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.beta,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling')},
    parameters={'metric': Str % Choices(beta.non_phylogenetic_metrics())},
    outputs=[('distance_matrix', DistanceMatrix)],
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.alpha_phylogenetic,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling'),
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(alpha.phylogenetic_metrics())},
    outputs=[('alpha_diversity',
              SampleData[AlphaDiversity] % Properties('phylogenetic'))],
    name='Alpha diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic alpha diversity "
                 "metric for all samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.alpha,
    inputs={'table': FeatureTable[Frequency] % Properties('uniform-sampling')},
    parameters={'metric': Str % Choices(alpha.non_phylogenetic_metrics())},
    outputs=[('alpha_diversity', SampleData[AlphaDiversity])],
    name='Alpha diversity',
    description=("Computes a user-specified alpha diversity metric for all "
                 "samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.pcoa,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={},
    outputs=[('pcoa', PCoAResults)],
    name='Principal Coordinate Analysis',
    description=("Apply principal coordinate analysis.")
)

plugin.visualizers.register_function(
    function=q2_diversity.alpha_compare,
    inputs={'alpha_diversity': SampleData[AlphaDiversity]},
    parameters={'metadata': MetadataCategory},
    name='Alpha diversity comparisons',
    description=("Visually and statistically compare groups of alpha diversity"
                 " values.")
)

plugin.visualizers.register_function(
    function=q2_diversity.bioenv,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'metadata': Metadata},
    name='bioenv',
    description=("Find the subsets of variables in metadata whose Euclidean "
                 "distances are maximally rank-correlated with distance "
                 "matrix. All numeric variables in metadata will be "
                 "considered, and samples which are missing data will be "
                 "dropped. The output visualization will indicate how many "
                 "samples were dropped due to missing data, if any were "
                 "dropped.")
)

plugin.methods.register_markdown('markdown/core_metrics.md')
