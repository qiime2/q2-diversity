# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Properties, MetadataCategory, Choices,
                           Metadata, Int, Bool, Set)

import q2_diversity
from q2_diversity import _alpha as alpha
from q2_diversity import _beta as beta
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.tree import Phylogeny, Rooted
from q2_types.ordination import PCoAResults


plugin = Plugin(
    name='diversity',
    version=q2_diversity.__version__,
    website='https://github.com/qiime2/q2-diversity',
    package='q2_diversity',
    description=('This QIIME 2 plugin supports metrics for calculating '
                 'and exploring community alpha and beta diversity through '
                 'statistics and visualizations in the context of sample '
                 'metadata.'),
    short_description='Plugin for exploring community diversity.'
)

plugin.methods.register_function(
    function=q2_diversity.beta_phylogenetic,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta.phylogenetic_metrics())},
    outputs=[('distance_matrix', DistanceMatrix % Properties('phylogenetic'))],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.'),
        'phylogeny': ('Phylogenetic tree containing tip identifiers that '
                      'correspond to the feature identifiers in the table. '
                      'This tree can contain tip ids that are not present in '
                      'the table, but all feature ids in the table must be '
                      'present in this tree.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.'
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.beta,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'metric': Str % Choices(beta.non_phylogenetic_metrics())},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.'
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.alpha_phylogenetic,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(alpha.phylogenetic_metrics())},
    outputs=[('alpha_diversity',
              SampleData[AlphaDiversity] % Properties('phylogenetic'))],
    input_descriptions={
        'table': ('The feature table containing the samples for which alpha '
                  'diversity should be computed.'),
        'phylogeny': ('Phylogenetic tree containing tip identifiers that '
                      'correspond to the feature identifiers in the table. '
                      'This tree can contain tip ids that are not present in '
                      'the table, but all feature ids in the table must be '
                      'present in this tree.')
    },
    parameter_descriptions={
        'metric': 'The alpha diversity metric to be computed.'
    },
    output_descriptions={
        'alpha_diversity': 'Vector containing per-sample alpha diversities.'
    },
    name='Alpha diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic alpha diversity "
                 "metric for all samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.alpha,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'metric': Str % Choices(alpha.non_phylogenetic_metrics())},
    outputs=[('alpha_diversity', SampleData[AlphaDiversity])],
    input_descriptions={
        'table': ('The feature table containing the samples for which alpha '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The alpha diversity metric to be computed.'
    },
    output_descriptions={
        'alpha_diversity': 'Vector containing per-sample alpha diversities.'
    },
    name='Alpha diversity',
    description=("Computes a user-specified alpha diversity metric for all "
                 "samples in a feature table.")
)

plugin.methods.register_function(
    function=q2_diversity.pcoa,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={},
    outputs=[('pcoa', PCoAResults)],
    input_descriptions={
        'distance_matrix': ('The distance matrix on which PCoA should be '
                            'computed.')
    },
    parameter_descriptions={},
    output_descriptions={'pcoa': 'The resulting PCoA matrix.'},
    name='Principal Coordinate Analysis',
    description=("Apply principal coordinate analysis.")
)

plugin.methods.register_function(
    function=q2_diversity.core_metrics,
    inputs={
        'table': FeatureTable[Frequency],
        'phylogeny': Phylogeny[Rooted]
    },
    parameters={
        'sampling_depth': Int
    },
    outputs=[
        ('faith_pd_vector', SampleData[AlphaDiversity]),
        ('observed_otus_vector', SampleData[AlphaDiversity]),
        ('shannon_vector', SampleData[AlphaDiversity]),
        ('evenness_vector', SampleData[AlphaDiversity]),
        ('unweighted_unifrac_distance_matrix', DistanceMatrix),
        ('weighted_unifrac_distance_matrix', DistanceMatrix),
        ('jaccard_distance_matrix', DistanceMatrix),
        ('bray_curtis_distance_matrix', DistanceMatrix),
        ('unweighted_unifrac_pcoa_results', PCoAResults),
        ('weighted_unifrac_pcoa_results', PCoAResults),
        ('jaccard_pcoa_results', PCoAResults),
        ('bray_curtis_pcoa_results', PCoAResults)
    ],
    input_descriptions={
        'table': ('The feature table containing the samples over which '
                  'diversity metrics should be computed.'),
        'phylogeny': ('Phylogenetic tree containing tip identifiers that '
                      'correspond to the feature identifiers in the table. '
                      'This tree can contain tip ids that are not present in '
                      'the table, but all feature ids in the table must be '
                      'present in this tree.')
    },
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'rarefied to prior to computing diversity metrics.')
    },
    output_descriptions={
        'faith_pd_vector': 'Vector of Faith PD values by sample.',
        'observed_otus_vector': 'Vector of Observed OTUs values by sample.',
        'shannon_vector': 'Vector of Shannon diversity values by sample.',
        'evenness_vector': 'Vector of Pielou\'s evenness values by sample.',
        'unweighted_unifrac_distance_matrix':
            'Matrix of unweighted UniFrac distances between pairs of samples.',
        'weighted_unifrac_distance_matrix':
            'Matrix of weighted UniFrac distances between pairs of samples.',
        'jaccard_distance_matrix':
            'Matrix of Jaccard distances between pairs of samples.',
        'bray_curtis_distance_matrix':
            'Matrix of Bray-Curtis distances between pairs of samples.',
        'unweighted_unifrac_pcoa_results':
            ('PCoA matrix computed from unweighted UniFrac distances between '
             'samples.'),
        'weighted_unifrac_pcoa_results':
            ('PCoA matrix computed from weighted UniFrac distances between '
             'samples.'),
        'jaccard_pcoa_results':
            ('PCoA matrix computed from Jaccard distances between '
             'samples.'),
        'bray_curtis_pcoa_results':
            ('PCoA matrix computed from Bray-Curtis distances between '
             'samples.'),
    },
    name='Core diversity metrics',
    description="Applies a collection of diversity metrics to a feature table."
)

plugin.methods.register_function(
    function=q2_diversity.filter_distance_matrix,
    inputs={
        'distance_matrix': DistanceMatrix
    },
    parameters={
        'metadata': Metadata,
        'where': Str,
        'exclude_ids': Bool
    },
    outputs=[
        ('filtered_distance_matrix', DistanceMatrix)
    ],
    name="Filter samples from a distance matrix.",
    description="Filter samples from a distance matrix, retaining only the "
                "samples matching search criteria specified by "
                "`metadata` and `where` parameters (or retaining only the "
                "samples not matching that criteria, if `exclude_ids` is "
                "True). See the filtering tutorial on "
                "https://docs.qiime2.org for additional details.",
    input_descriptions={
        'distance_matrix': 'Distance matrix to filter by sample.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata used with `where` parameter when '
                    'selecting samples to retain, or with `exclude_ids` '
                    'when selecting samples to discard.',
        'where': 'SQLite WHERE clause specifying sample metadata criteria '
                 'that must be met to be included in the filtered distance '
                 'matrix. If not provided, all samples in `metadata` that are '
                 'also in the input distance matrix will be retained.',
        'exclude_ids': 'If `True`, the samples selected by `metadata` or '
                       '`where` parameters will be excluded from the filtered '
                       'distance matrix instead of being retained.'
    },
    output_descriptions={
        'filtered_distance_matrix': 'Distance matrix filtered to include '
                                    'samples matching search criteria'
    }
)

plugin.visualizers.register_function(
    function=q2_diversity.alpha_group_significance,
    inputs={'alpha_diversity': SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata},
    input_descriptions={
        'alpha_diversity': 'Vector of alpha diversity values by sample.'
    },
    parameter_descriptions={
        'metadata': 'The sample metadata.'
    },
    name='Alpha diversity comparisons',
    description=("Visually and statistically compare groups of alpha diversity"
                 " values.")
)

plugin.visualizers.register_function(
    function=q2_diversity.bioenv,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'metadata': Metadata},
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'metadata': 'The sample metadata.'
    },
    name='bioenv',
    description=("Find the subsets of variables in metadata whose Euclidean "
                 "distances are maximally rank-correlated with distance "
                 "matrix. All numeric variables in metadata will be "
                 "considered, and samples which are missing data will be "
                 "dropped. The output visualization will indicate how many "
                 "samples were dropped due to missing data, if any were "
                 "dropped.")
)

beta_group_significance_methods = \
    list(q2_diversity._beta._visualizer._beta_group_significance_fns)

plugin.visualizers.register_function(
    function=q2_diversity.beta_group_significance,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'method': Str % Choices(beta_group_significance_methods),
                'permutations': Int,
                'metadata': MetadataCategory,
                'pairwise': Bool},
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'method': 'The group significance test to be applied.',
        'permutations': ('The number of permutations to be run when computing '
                         'p-values.'),
        'metadata': 'The sample metadata.',
        'pairwise': ('Perform pairwise tests between all pairs of groups '
                     'in addition to the test across all groups. '
                     'This can be very slow if there are a lot of groups '
                     'in the category.')
    },
    name='Beta diversity group significance',
    description=('Determine whether groups of samples are significantly '
                 'different from one another using a permutation-based '
                 'statistical test.')
)

beta_correlation_methods = ['spearman', 'pearson']

plugin.visualizers.register_function(
    function=q2_diversity.beta_correlation,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'metadata': MetadataCategory,
                'permutations': Int,
                'method': Str % Choices(beta_correlation_methods)},
    name=('Beta diversity correlation'),
    description=('Apply a two-sided Mantel test to identify correlation '
                 'between the distance matrix and a numeric sample metadata '
                 'category. Sample metadata pairwise distances are computed '
                 'as the Euclidean distance between each pair of samples in '
                 'the metadata category.'),
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'method': 'The correlation test to be applied in the Mantel test.',
        'permutations': ('The number of permutations to be run when computing '
                         'p-values.'),
        'metadata': 'The sample metadata.',
    },
)


alpha_correlation_methods = \
    list(q2_diversity._alpha._visualizer._alpha_correlation_fns)

plugin.visualizers.register_function(
    function=q2_diversity.alpha_correlation,
    inputs={'alpha_diversity': SampleData[AlphaDiversity]},
    parameters={'method': Str % Choices(alpha_correlation_methods),
                'metadata': Metadata},
    input_descriptions={
        'alpha_diversity': 'Vector of alpha diversity values by sample.'
    },
    parameter_descriptions={
        'method': 'The correlation test to be applied.',
        'metadata': 'The sample metadata.'
    },
    name='Alpha diversity correlation',
    description=('Determine whether numeric sample metadata category is '
                 'correlated with alpha diversity.')
)

plugin.visualizers.register_function(
    function=q2_diversity.alpha_rarefaction,
    inputs={'feature_table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metrics': Set[Str % Choices(
                                alpha.non_phylogenetic_metrics())],
                'min_depth': Int,
                'max_depth': Int,
                'steps': Int,
                'iterations': Int},
    input_descriptions={
        'feature_table': 'Feature table to be rarefied.',
        'phylogeny': 'Optional phylogeny for phylogenetic metrics.'
    },
    parameter_descriptions={
        'metrics': 'The metrics to be measured.',
        'min_depth': 'The minimum sequencing depth to rarefy.',
        'max_depth': 'The maximum sequencing depth to rarefy.',
        'steps': 'The number of sequencing depths to rarefy.',
        'iterations': 'The number of subsamples to rarefy at each step.'
    },
    name='Alpha rarefaction',
    description=('Compute alpha diversity metrics for subsamples of a '
                 'feature table at sequential sequencing depths, and '
                 'create alpha rarefaction plots for the requested metrics')
)
