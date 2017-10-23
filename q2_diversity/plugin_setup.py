# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Properties, MetadataCategory, Choices,
                           Metadata, Int, Bool, Range, Float, Visualization)

import q2_diversity
from q2_diversity import _alpha as alpha
from q2_diversity import _beta as beta
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.distance_matrix import DistanceMatrix
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.tree import Phylogeny, Rooted
from q2_types.ordination import PCoAResults


sklearn_n_jobs_description = (
    'The number of jobs to use for the computation. This works by breaking '
    'down the pairwise matrix into n_jobs even slices and computing them in '
    'parallel. If -1 all CPUs are used. If 1 is given, no parallel computing '
    'code is used at all, which is useful for debugging. For n_jobs below -1, '
    '(n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one '
    'are used. (Description from sklearn.metrics.pairwise_distances)'
)


plugin = Plugin(
    name='diversity',
    version=q2_diversity.__version__,
    website='https://github.com/qiime2/q2-diversity',
    package='q2_diversity',
    description=('This QIIME 2 plugin supports metrics for calculating '
                 'and exploring community alpha and beta diversity through '
                 'statistics and visualizations in the context of sample '
                 'metadata.'),
    short_description='Plugin for exploring community diversity.',
    citation_text=('Unweighted UniFrac: '
                   'Lozupone and Knight 2005 Appl Environ Microbiol; DOI: '
                   '10.1128/AEM.71.12.8228-8235.2005.\n'
                   'Weighted UniFrac: '
                   'Lozupone et al. 2007 Appl Environ Microbiol; DOI: '
                   '10.1128/AEM.01996-06.\n'
                   'Variance adjusted UniFrac: '
                   'Chang et al. BMC Bioinformatics 2011 '
                   'https://doi.org/10.1186/1471-2105-12-118.\n'
                   'Generalized UniFrac: '
                   'Chen et al. 2012 Bioinformatics; DOI: '
                   '10.1093/bioinformatics/bts342')
)

plugin.methods.register_function(
    function=q2_diversity.beta_phylogenetic,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta.phylogenetic_metrics()),
                'n_jobs': Int % Range(1, None)},
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
        'metric': 'The beta diversity metric to be computed.',
        'n_jobs': '[Excluding weighted_unifrac] - %s' %
                  sklearn_n_jobs_description
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table.")
)


plugin.methods.register_function(
    function=q2_diversity.beta_phylogenetic_alt,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta.phylogenetic_metrics_alt_dict()),
                'n_jobs': Int,
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True),
                'bypass_tips': Bool},
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
        'metric': 'The beta diversity metric to be computed.',
        'n_jobs': 'The number of workers to use.',
        'variance_adjusted': ('Perform variance adjustment based on Chang et '
                              'al. BMC Bioinformatics 2011. Weights distances '
                              'based on the proportion of the relative '
                              'abundance represented between the samples at a'
                              ' given node under evaluation.'),
        'alpha': ('This parameter is only used when the choice of metric is '
                  'generalized_unifrac. The value of alpha controls importance'
                  ' of sample proportions. 1.0 is weighted normalized UniFrac.'
                  ' 0.0 is close to unweighted UniFrac, but only if the sample'
                  ' proportions are dichotomized.'),
        'bypass_tips': ('In a bifurcating tree, the tips make up about 50% of '
                        'the nodes in a tree. By ignoring them, specificity '
                        'can be traded for reduced compute time. This has the'
                        ' effect of collapsing the phylogeny, and is analogous'
                        ' (in concept) to moving from 99% to 97% OTUs')
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity (phylogenetic) - High Performance Computation',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table. This "
                 "implementation is recommended for large datasets, otherwise "
                 "the results are identical to beta_phylogenetic.\n\n"
                 "This method is an implementation of the Strided State "
                 "UniFrac algorithm. Multiple variants of the UniFrac metric "
                 "are available, including Generalized UniFrac (Chen et al. "
                 "2012), Variance Adjusted UniFrac (Chang et al. 2011), "
                 "as well as Weighted normalized and unnormalized UniFrac "
                 "(Lozupone et al. 2007) and unweighted UniFrac "
                 "(Lozupone et al. 2005)")
)


plugin.methods.register_function(
    function=q2_diversity.beta,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'metric': Str % Choices(beta.non_phylogenetic_metrics()),
                'n_jobs': Int},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'n_jobs': sklearn_n_jobs_description
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

plugin.pipelines.register_function(
    function=q2_diversity.core_metrics_phylogenetic,
    inputs={
        'table': FeatureTable[Frequency],
        'phylogeny': Phylogeny[Rooted]
    },
    parameters={
        'sampling_depth': Int % Range(1, None),
        'metadata': Metadata,
        'n_jobs': Int % Range(0, None),
    },
    outputs=[
        ('rarefied_table', FeatureTable[Frequency]),
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
        ('bray_curtis_pcoa_results', PCoAResults),
        ('unweighted_unifrac_emperor', Visualization),
        ('weighted_unifrac_emperor', Visualization),
        ('jaccard_emperor', Visualization),
        ('bray_curtis_emperor', Visualization),
    ],
    input_descriptions={
        'table': 'The feature table containing the samples over which '
                 'diversity metrics should be computed.',
        'phylogeny': 'Phylogenetic tree containing tip identifiers that '
                     'correspond to the feature identifiers in the table. '
                     'This tree can contain tip ids that are not present in '
                     'the table, but all feature ids in the table must be '
                     'present in this tree.'
    },
    parameter_descriptions={
        'sampling_depth': 'The total frequency that each sample should be '
                          'rarefied to prior to computing diversity metrics.',
        'metadata': 'The sample metadata to use in the emperor plots.',
        'n_jobs': '[beta/beta-phylogenetic methods only, excluding weighted_'
                  'unifrac] - %s' % sklearn_n_jobs_description
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.',
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
            'PCoA matrix computed from unweighted UniFrac distances between '
            'samples.',
        'weighted_unifrac_pcoa_results':
            'PCoA matrix computed from weighted UniFrac distances between '
            'samples.',
        'jaccard_pcoa_results':
            'PCoA matrix computed from Jaccard distances between '
            'samples.',
        'bray_curtis_pcoa_results':
            'PCoA matrix computed from Bray-Curtis distances between '
            'samples.',
        'unweighted_unifrac_emperor':
            'Emperor plot of the PCoA matrix computed from unweighted'
            ' UniFrac.',
        'weighted_unifrac_emperor':
            'Emperor plot of the PCoA matrix computed from weighted UniFrac.',
        'jaccard_emperor':
            'Emperor plot of the PCoA matrix computed from Jaccard.',
        'bray_curtis_emperor':
            'Emperor plot of the PCoA matrix computed from Bray-Curtis.',
    },
    name='Core diversity metrics (phylogenetic and non-phylogenetic)',
    description="Applies a collection of diversity metrics (both "
                "phylogenetic and non-phylogenetic) to a feature table.",
)

plugin.pipelines.register_function(
    function=q2_diversity.core_metrics,
    inputs={
        'table': FeatureTable[Frequency],
    },
    parameters={
        'sampling_depth': Int % Range(1, None),
        'metadata': Metadata,
        'n_jobs': Int % Range(0, None),
    },
    outputs=[
        ('rarefied_table', FeatureTable[Frequency]),
        ('observed_otus_vector', SampleData[AlphaDiversity]),
        ('shannon_vector', SampleData[AlphaDiversity]),
        ('evenness_vector', SampleData[AlphaDiversity]),
        ('jaccard_distance_matrix', DistanceMatrix),
        ('bray_curtis_distance_matrix', DistanceMatrix),
        ('jaccard_pcoa_results', PCoAResults),
        ('bray_curtis_pcoa_results', PCoAResults),
        ('jaccard_emperor', Visualization),
        ('bray_curtis_emperor', Visualization),
    ],
    input_descriptions={
        'table': 'The feature table containing the samples over which '
                 'diversity metrics should be computed.',
    },
    parameter_descriptions={
        'sampling_depth': 'The total frequency that each sample should be '
                          'rarefied to prior to computing diversity metrics.',
        'metadata': 'The sample metadata to use in the emperor plots.',
        'n_jobs': '[beta methods only] - %s' % sklearn_n_jobs_description
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.',
        'observed_otus_vector': 'Vector of Observed OTUs values by sample.',
        'shannon_vector': 'Vector of Shannon diversity values by sample.',
        'evenness_vector': 'Vector of Pielou\'s evenness values by sample.',
        'jaccard_distance_matrix':
            'Matrix of Jaccard distances between pairs of samples.',
        'bray_curtis_distance_matrix':
            'Matrix of Bray-Curtis distances between pairs of samples.',
        'jaccard_pcoa_results':
            'PCoA matrix computed from Jaccard distances between samples.',
        'bray_curtis_pcoa_results':
            'PCoA matrix computed from Bray-Curtis distances between samples.',
        'jaccard_emperor':
            'Emperor plot of the PCoA matrix computed from Jaccard.',
        'bray_curtis_emperor':
            'Emperor plot of the PCoA matrix computed from Bray-Curtis.',
    },
    name='Core diversity metrics (non-phylogenetic)',
    description=("Applies a collection of diversity metrics "
                 "(non-phylogenetic) to a feature table."),
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
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(
                                    alpha.alpha_rarefaction_supported_metrics),
                'metadata': Metadata,
                'min_depth': Int % Range(1, None),
                'max_depth': Int % Range(1, None),
                'steps': Int % Range(2, None),
                'iterations': Int % Range(1, None)},
    input_descriptions={
        'table': 'Feature table to compute rarefaction curves from.',
        'phylogeny': 'Optional phylogeny for phylogenetic metrics.',
    },
    parameter_descriptions={
        'metric': ('The metric to be measured. By default computes '
                   'observed_otus, shannon, and if phylogeny is '
                   'provided, faith_pd.'),
        'metadata': 'The sample metadata.',
        'min_depth': 'The minimum rarefaction depth.',
        'max_depth': ('The maximum rarefaction depth. '
                      'Must be greater than min_depth.'),
        'steps': ('The number of rarefaction depths to include '
                  'between min_depth and max_depth.'),
        'iterations': ('The number of rarefied feature tables to '
                       'compute at each step.'),
    },
    name='Alpha rarefaction curves',
    description=('Generate interactive alpha rarefaction curves by computing '
                 'rarefactions between `min_depth` and `max_depth`. The '
                 'number of intermediate depths to compute is controlled by '
                 'the `steps` parameter, with n `iterations` being computed '
                 'at each rarefaction depth. If sample metadata is provided, '
                 'samples may be grouped based on distinct values within a '
                 'metadata column.'),
)

_beta_rarefaction_color_schemes = [
    'BrBG', 'BrBG_r', 'PRGn', 'PRGn_r', 'PiYG', 'PiYG_r',
    'PuOr', 'PuOr_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r',
    'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r']

plugin.visualizers.register_function(
    function=q2_diversity._beta._visualizer.beta_rarefaction,
    inputs={
        'table': FeatureTable[Frequency],
        'phylogeny': Phylogeny[Rooted]},
    parameters={
        'metric': Str % Choices(beta.all_metrics()),
        'sampling_depth': Int % Range(1, None),
        # Need at least two iterations to do a comparison.
        'iterations': Int % Range(2, None),
        'correlation_method': Str % Choices({'spearman', 'pearson'}),
        'color_scheme': Str % Choices(_beta_rarefaction_color_schemes)
    },
    input_descriptions={
        'table': 'Feature table upon which to perform beta diversity '
                 'rarefaction analyses.',
        'phylogeny': 'Phylogenetic tree containing tip identifiers that '
                     'correspond to the feature identifiers in the table. '
                     'This tree can contain tip ids that are not present in '
                     'the table, but all feature ids in the table must be '
                     'present in this tree. [required for phylogenetic '
                     'metrics]'
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'sampling_depth': 'The total frequency that each sample should be '
                          'rarefied to prior to computing diversity '
                          'metrics.',
        'iterations': 'Number of times to rarefy the feature table at a given '
                      'sampling depth.',
        'correlation_method': 'The Mantel correlation test to be applied when '
                              'computing correlation between beta diversity '
                              'distance matrices.',
        'color_scheme': 'The matplotlib color scheme to generate the heatmap '
                        'with.',
    },
    name='Beta diversity rarefaction',
    description='Repeatedly rarefy a feature table to compare beta diversity '
                'results across and within rarefaction depths.\n\n'
                'This visualizer provides a heatmap showing the correlation '
                'between beta diversity distance matrices computed by '
                'repeatedly rarefying at a given sampling depth.\n\n'
                'Note: currently only a single sampling depth is supported. A '
                'range of sampling depths (e.g. similar to '
                '`alpha_rarefaction`) will be supported in the future to '
                'allow comparison between sampling depths. Additional '
                'visualizations will be supported, such as Emperor plots and '
                'UPGMA trees.'
)
