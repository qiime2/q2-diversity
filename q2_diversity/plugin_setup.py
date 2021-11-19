# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations)
import q2_diversity
from q2_diversity import _alpha as alpha
from q2_diversity import _beta as beta
from q2_types.feature_table import (FeatureTable, Frequency, RelativeFrequency,
                                    PresenceAbsence)
from q2_types.distance_matrix import DistanceMatrix
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.tree import Phylogeny, Rooted
from q2_types.ordination import PCoAResults, ProcrustesStatistics

citations = Citations.load('citations.bib', package='q2_diversity')

n_jobs_description = (
    'The number of concurrent jobs to use in performing this calculation. '
    'May not exceed the number of available physical cores. If n_jobs = '
    '\'auto\', one job will be launched for each identified CPU core on the '
    'host.'
)

threads_description = (
    'The number of CPU threads to use in performing this calculation. '
    'May not exceed the number of available physical cores. If threads = '
    '\'auto\', one thread will be created for each identified CPU core on the '
    'host.'
)

n_jobs_or_threads_description = (
    'The number of concurrent jobs or CPU threads to use in performing this '
    'calculation. Individual methods will create jobs/threads as implemented '
    'in q2-diversity-lib dependencies. May not exceed the number of available '
    'physical cores. If n_jobs_or_threads = \'auto\', one thread/job will be '
    'created for each identified CPU core on the host.'
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
)

plugin.pipelines.register_function(
    function=q2_diversity.beta_phylogenetic,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(beta.METRICS['PHYLO']['IMPL'] |
                                        beta.METRICS['PHYLO']['UNIMPL']),
                'threads': Int % Range(1, None) | Str % Choices(['auto']),
                'variance_adjusted': Bool,
                'alpha': Float % Range(0, 1, inclusive_end=True),
                'bypass_tips': Bool},
    outputs=[('distance_matrix', DistanceMatrix)],
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
        'threads': threads_description,
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
    name='Beta diversity (phylogenetic)',
    description=("Computes a user-specified phylogenetic beta diversity metric"
                 " for all pairs of samples in a feature table.")
)

plugin.pipelines.register_function(
    function=q2_diversity.beta,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={'metric': Str % Choices(beta.METRICS['NONPHYLO']['IMPL'] |
                                        beta.METRICS['NONPHYLO']['UNIMPL']),
                'pseudocount': Int % Range(1, None),
                'n_jobs': Int % Range(1, None) | Str % Choices(['auto'])},
    outputs=[('distance_matrix', DistanceMatrix)],
    input_descriptions={
        'table': ('The feature table containing the samples over which beta '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The beta diversity metric to be computed.',
        'pseudocount': ('A pseudocount to handle zeros for compositional '
                        'metrics.  This is ignored for other metrics.'),
        'n_jobs': n_jobs_description
    },
    output_descriptions={'distance_matrix': 'The resulting distance matrix.'},
    name='Beta diversity',
    description=("Computes a user-specified beta diversity metric for all "
                 "pairs of samples in a feature table.")
)

plugin.pipelines.register_function(
    function=q2_diversity.alpha_phylogenetic,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metric': Str % Choices(alpha.METRICS['PHYLO']['IMPL'] |
                                        alpha.METRICS['PHYLO']['UNIMPL'])},
    outputs=[('alpha_diversity',
              SampleData[AlphaDiversity])],
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
    description=('Computes a user-specified phylogenetic alpha diversity '
                 'metric for all samples in a feature table.'),
)


plugin.pipelines.register_function(
    function=q2_diversity.alpha,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={'metric': Str % Choices(alpha.METRICS['NONPHYLO']['IMPL'] |
                                        alpha.METRICS['NONPHYLO']['UNIMPL'])},
    outputs=[('alpha_diversity', SampleData[AlphaDiversity])],
    input_descriptions={
        'table': ('The feature table containing the samples for which alpha '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The alpha diversity metric to be computed. Information '
        'about specific metrics is available at '
        'https://data.qiime2.org/a_diversity_metrics'
    },
    output_descriptions={
        'alpha_diversity': 'Vector containing per-sample alpha diversities.'
    },
    name='Alpha diversity',
    description=('Computes a user-specified alpha diversity metric for all '
                 'samples in a feature table.')
)

plugin.methods.register_function(
    function=q2_diversity.pcoa,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={
        'number_of_dimensions': Int % Range(1, None)
    },
    outputs=[('pcoa', PCoAResults)],
    input_descriptions={
        'distance_matrix': ('The distance matrix on which PCoA should be '
                            'computed.')
    },
    parameter_descriptions={
        'number_of_dimensions': "Dimensions to reduce the distance matrix to. "
                                "This number determines how many "
                                "eigenvectors and eigenvalues are returned,"
                                "and influences the choice of algorithm used "
                                "to compute them. "
                                "By default, uses the default "
                                "eigendecomposition method, SciPy's eigh, "
                                "which computes all eigenvectors "
                                "and eigenvalues in an exact manner. For very "
                                "large matrices, this is expected to be slow. "
                                "If a value is specified for this parameter, "
                                "then the fast, heuristic "
                                "eigendecomposition algorithm fsvd "
                                "is used, which only computes and returns the "
                                "number of dimensions specified, but suffers "
                                "some degree of accuracy loss, the magnitude "
                                "of which varies across different datasets."
    },
    output_descriptions={'pcoa': 'The resulting PCoA matrix.'},
    name='Principal Coordinate Analysis',
    description=("Apply principal coordinate analysis."),
    citations=[citations['legendrelegendre'],
               citations['halko2010']]
)

plugin.methods.register_function(
    function=q2_diversity.pcoa_biplot,
    inputs={'pcoa': PCoAResults,
            'features': FeatureTable[RelativeFrequency]},
    parameters={},
    outputs=[('biplot', PCoAResults % Properties('biplot'))],
    input_descriptions={
        'pcoa': 'The PCoA where the features will be projected onto.',
        'features': 'Variables to project onto the PCoA matrix'
    },
    parameter_descriptions={},
    output_descriptions={'biplot': 'The resulting PCoA matrix.'},
    name='Principal Coordinate Analysis Biplot',
    description="Project features into a principal coordinates matrix. The "
                "features used should be the features used to compute the "
                "distance matrix. It is recommended that these variables be"
                " normalized in cases of dimensionally heterogeneous physical"
                " variables.",
    citations=[citations['legendrelegendre']]
)

plugin.methods.register_function(
    function=q2_diversity.tsne,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={
        'number_of_dimensions': Int % Range(2, None),
        'perplexity': Float % Range(1, None),
        'early_exaggeration': Float % Range(0, None),
        'learning_rate': Float % Range(10.0, None),
        'n_iter': Int % Range(1, None)
    },
    outputs=[('tsne', PCoAResults)],
    input_descriptions={
        'distance_matrix': ('The distance matrix on which t-SNE should be '
                            'computed.')
    },
    parameter_descriptions={
        'number_of_dimensions': "Dimensions to reduce the distance matrix to.",
        'perplexity':           "Provide the balance between local and global "
                                "structure. Low values concentrate on local "
                                "structure. Large values sacrifice local "
                                "details for a broader global embedding. "
                                "The default value is 25 to achieve better "
                                "results for small microbiome datasets.",
        'early_exaggeration':   "Affects the tightness of the shown clusters. "
                                "Larger values increase the distance between "
                                "natural clusters in the embedded space.",
        'learning_rate':        "Controls how much the weights are adjusted "
                                "at each update.",
    },
    output_descriptions={'tsne': 'The resulting t-SNE matrix.'},
    name='t-distributed stochastic neighbor embedding',
    description=("Apply t-distributed stochastic neighbor embedding."),

)

plugin.methods.register_function(
    function=q2_diversity.umap,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={
        'number_of_dimensions': Int % Range(2, None),
        'n_neighbors': Int % Range(1, None),
        'min_dist': Float % Range(0, None),
    },
    outputs=[('umap', PCoAResults)],
    input_descriptions={
        'distance_matrix': ('The distance matrix on which UMAP should be '
                            'computed.')
    },
    parameter_descriptions={
        'number_of_dimensions': "Dimensions to reduce the distance matrix to.",
        'n_neighbors':          "Provide the balance between local and global "
                                "structure. Low values prioritize the "
                                "preservation of local structures. Large "
                                "values sacrifice local details for a "
                                "broader global embedding.",
        'min_dist':             "Controls the cluster size. Low values cause "
                                "clumpier clusters. Higher values preserve a "
                                "broad topological structure. To get "
                                "less overlapping data points the "
                                "default value is set to 0.4. For more "
                                "details visit: "
                                "https://umap-learn.readthedocs.io/en/latest/"
                                "parameters.html",
    },
    output_descriptions={'umap': 'The resulting UMAP matrix.'},
    name='Uniform Manifold Approximation and Projection',
    description=("Apply Uniform Manifold Approximation and Projection."),
    )

plugin.methods.register_function(
    function=q2_diversity.procrustes_analysis,
    inputs={'reference': PCoAResults, 'other': PCoAResults},
    parameters={
        'dimensions': Int % Range(1, None),
        'permutations': Int % Range(1, None) | Str % Choices('disable')
    },
    outputs=[
        ('transformed_reference', PCoAResults),
        ('transformed_other', PCoAResults),
        ('disparity_results', ProcrustesStatistics)
    ],
    input_descriptions={
        'reference': ('The ordination matrix to which data is fitted to.'),
        'other': ("The ordination matrix that's fitted to the reference "
                  "ordination."),
    },
    parameter_descriptions={
        'dimensions': ('The number of dimensions to use when fitting the two '
                       'matrices'),
        'permutations': 'The number of permutations to be run when computing '
                        'p-values. Supplying a value of `disable` will disable'
                        ' permutation testing and p-values will not be '
                        'calculated (this results in *much* quicker execution '
                        'time if p-values are not desired).',

    },
    output_descriptions={
        'transformed_reference': 'A normalized version of the "reference" '
                                 'ordination matrix.',
        'transformed_other': 'A normalized and fitted version of the "other" '
                             'ordination matrix.',
        'disparity_results': 'The sum of the squares of the pointwise '
                             'differences between the two input datasets & '
                             'its p value.'},
    name='Procrustes Analysis',
    description='Fit two ordination matrices with Procrustes analysis'
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
        'n_jobs_or_threads': Int % Range(1, None) | Str % Choices(['auto']),
    },
    outputs=[
        ('rarefied_table', FeatureTable[Frequency]),
        ('faith_pd_vector', SampleData[AlphaDiversity]),
        ('observed_features_vector', SampleData[AlphaDiversity]),
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
        'n_jobs_or_threads': '[beta/beta-phylogenetic methods only] - %s'
                          % n_jobs_or_threads_description
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.',
        'faith_pd_vector': 'Vector of Faith PD values by sample.',
        'observed_features_vector': 'Vector of Observed Features values by '
                                    'sample.',
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
                "phylogenetic and non-phylogenetic) to a feature table."
)

plugin.pipelines.register_function(
    function=q2_diversity.core_metrics,
    inputs={
        'table': FeatureTable[Frequency],
    },
    parameters={
        'sampling_depth': Int % Range(1, None),
        'metadata': Metadata,
        'with_replacement': Bool,
        'n_jobs': Int % Range(1, None) | Str % Choices(['auto']),
    },
    outputs=[
        ('rarefied_table', FeatureTable[Frequency]),
        ('observed_features_vector', SampleData[AlphaDiversity]),
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
        'with_replacement': 'Rarefy with replacement by sampling from the '
                            'multinomial distribution instead of rarefying '
                            'without replacement.',
        'n_jobs': '[beta methods only] - %s' % n_jobs_description
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.',
        'observed_features_vector': 'Vector of Observed Features values by '
                                    'sample.',
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
                 "(non-phylogenetic) to a feature table.")
)

plugin.pipelines.register_function(
    function=q2_diversity.beta_correlation,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={
        'metadata': MetadataColumn[Numeric],
        'method': Str % Choices(['spearman', 'pearson']),
        'permutations': Int % Range(0, None),
        'intersect_ids': Bool,
        'label1': Str,
        'label2': Str
    },
    outputs=[('metadata_distance_matrix', DistanceMatrix),
             ('mantel_scatter_visualization', Visualization)],
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'},
    parameter_descriptions={
        'metadata': 'Numeric metadata column from which to compute pairwise '
                    'Euclidean distances',
        'method': 'The correlation test to be applied in the Mantel test.',
        'permutations': 'The number of permutations to be run when computing '
                        'p-values. Supplying a value of zero will disable '
                        'permutation testing and p-values will not be '
                        'calculated (this results in *much* quicker execution '
                        'time if p-values are not desired).',
        'intersect_ids': 'If supplied, IDs that are not found in both '
                         'distance matrices will be discarded before applying '
                         'the Mantel test. Default behavior is to error on '
                         'any mismatched IDs.',
        'label1': 'Label for `distance_matrix` in the output visualization.',
        'label2': 'Label for `metadata_distance_matrix` in the output '
                  'visualization.'
    },
    output_descriptions={
        'metadata_distance_matrix': 'The Distance Matrix produced from the '
                                    'metadata column and used in the mantel '
                                    'test',
        'mantel_scatter_visualization': 'Scatter plot rendering of the mantel'
                                        'test results'},
    name='Beta diversity correlation',
    description=('Create a distance matrix from a numeric metadata column and '
                 'apply a two-sided Mantel test to identify correlation '
                 'between two distance matrices. Actions used internally: '
                 '`distance-matrix` from q2-metadata and `mantel` from '
                 'q2-diversity.')
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
                 " values."),
    citations=[citations['kruskal1952use']]
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
                 "dropped."),
    citations=[citations['clarke1993method']]
)

beta_group_significance_methods = \
    list(q2_diversity._beta._visualizer._beta_group_significance_fns)

plugin.visualizers.register_function(
    function=q2_diversity.beta_group_significance,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'method': Str % Choices(beta_group_significance_methods),
                'permutations': Int,
                'metadata': MetadataColumn[Categorical],
                'pairwise': Bool},
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'method': 'The group significance test to be applied.',
        'permutations': ('The number of permutations to be run when computing '
                         'p-values.'),
        'metadata': 'Categorical sample metadata column.',
        'pairwise': ('Perform pairwise tests between all pairs of groups '
                     'in addition to the test across all groups. '
                     'This can be very slow if there are a lot of groups '
                     'in the metadata column.')
    },
    name='Beta diversity group significance',
    description=('Determine whether groups of samples are significantly '
                 'different from one another using a permutation-based '
                 'statistical test.'),
    citations=[citations['anderson2001new']]
)

plugin.visualizers.register_function(
    function=q2_diversity.mantel,
    inputs={'dm1': DistanceMatrix,
            'dm2': DistanceMatrix},
    parameters={'permutations': Int % Range(0, None),
                'method': Str % Choices(['spearman', 'pearson']),
                'intersect_ids': Bool,
                'label1': Str,
                'label2': Str},
    name='Apply the Mantel test to two distance matrices',
    description='Apply a two-sided Mantel test to identify correlation '
                'between two distance matrices.\n\nNote: the directionality '
                'of the comparison has no bearing on the results. Thus, '
                'comparing distance matrix X to distance matrix Y is '
                'equivalent to comparing Y to X.\n\nNote: the order of '
                'samples within the two distance matrices does not need to be '
                'the same; the distance matrices will be reordered before '
                'applying the Mantel test.\n\nSee the scikit-bio docs for '
                'more details about the Mantel test:\n\n'
                'http://scikit-bio.org/docs/latest/generated/'
                'skbio.stats.distance.mantel',
    input_descriptions={
        'dm1': 'Matrix of distances between pairs of samples.',
        'dm2': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'method': 'The correlation test to be applied in the Mantel test.',
        'permutations': 'The number of permutations to be run when computing '
                        'p-values. Supplying a value of zero will disable '
                        'permutation testing and p-values will not be '
                        'calculated (this results in *much* quicker execution '
                        'time if p-values are not desired).',
        'intersect_ids': 'If supplied, IDs that are not found in both '
                         'distance matrices will be discarded before applying '
                         'the Mantel test. Default behavior is to error on '
                         'any mismatched IDs.',
        'label1': 'Label for `dm1` in the output visualization.',
        'label2': 'Label for `dm2` in the output visualization.'
    },
    citations=[
        citations['mantel1967detection'],
        citations['pearson1895note'],
        citations['spearman1904proof']]
)

alpha_correlation_methods = \
    list(q2_diversity._alpha._visualizer._alpha_correlation_fns)

plugin.visualizers.register_function(
    function=q2_diversity.alpha_correlation,
    inputs={'alpha_diversity': SampleData[AlphaDiversity]},
    parameters={'method': Str % Choices(alpha_correlation_methods),
                'metadata': Metadata,
                'intersect_ids': Bool},
    input_descriptions={
        'alpha_diversity': 'Vector of alpha diversity values by sample.'
    },
    parameter_descriptions={
        'method': 'The correlation test to be applied.',
        'metadata': 'The sample metadata.',
        'intersect_ids': 'If supplied, IDs that are not found in both '
                         'the alpha diversity vector and metadata will '
                         'be discarded before calculating '
                         'the correlation. Default behavior is to error on '
                         'any mismatched IDs.'
    },
    name='Alpha diversity correlation',
    description=('Determine whether numeric sample metadata columns are '
                 'correlated with alpha diversity.'),
    citations=[citations['pearson1895note'], citations['spearman1904proof']]
)

_metric_set = Set[Str % Choices((alpha.METRICS['PHYLO']['IMPL'] |
                                 alpha.METRICS['PHYLO']['UNIMPL'] |
                                 alpha.METRICS['NONPHYLO']['IMPL'] |
                                 alpha.METRICS['NONPHYLO']['UNIMPL']) -
                                alpha.alpha_rarefaction_unsupported_metrics)]
plugin.visualizers.register_function(
    function=q2_diversity.alpha_rarefaction,
    inputs={'table': FeatureTable[Frequency],
            'phylogeny': Phylogeny[Rooted]},
    parameters={'metrics': _metric_set,
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
        'metrics': ('The metrics to be measured. By default computes '
                    'observed_features, shannon, and if phylogeny is '
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
    function=q2_diversity._beta.beta_rarefaction,
    inputs={
        'table': FeatureTable[Frequency],
        'phylogeny': Phylogeny[Rooted]},
    parameters={
        'metric': Str % Choices(beta.METRICS['NONPHYLO']['IMPL'] |
                                beta.METRICS['NONPHYLO']['UNIMPL'] |
                                beta.METRICS['PHYLO']['IMPL'] |
                                beta.METRICS['PHYLO']['UNIMPL']),
        'clustering_method': Str % Choices({'nj', 'upgma'}),
        'metadata': Metadata,
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
                          'rarefied to prior to computing the diversity '
                          'metric.',
        'clustering_method': 'Samples can be clustered with neighbor joining '
                             'or UPGMA. An arbitrary rarefaction trial will '
                             'be used for the tree, and the remaining trials '
                             'are used to calculate the support of the '
                             'internal nodes of that tree.',
        'metadata': 'The sample metadata used for the Emperor jackknifed PCoA '
                    'plot.',
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
                'results within a given rarefaction depth.\n\n'
                'For a given beta diversity metric, this visualizer will '
                'provide: an Emperor jackknifed PCoA plot, samples clustered '
                'by UPGMA or neighbor joining with support calculation, and '
                'a heatmap showing the correlation between rarefaction trials '
                'of that beta diversity metric.',
    citations=[
        citations['mantel1967detection'],
        citations['pearson1895note'],
        citations['spearman1904proof']]
)

plugin.visualizers.register_function(
    function=q2_diversity.adonis,
    inputs={'distance_matrix': DistanceMatrix},
    parameters={'metadata': Metadata,
                'formula': Str,
                'permutations': Int % Range(1, None),
                'n_jobs': Int % Range(1, None)},
    input_descriptions={
        'distance_matrix': 'Matrix of distances between pairs of samples.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata containing formula terms.',
        'formula': 'Model formula containing only independent terms contained '
                   'in the sample metadata. These can be continuous variables '
                   'or factors, and they can have interactions as in a '
                   'typical R formula. E.g., the formula "treatment+block" '
                   'would test whether the input distance matrix partitions '
                   'based on "treatment" and "block" sample metadata. The '
                   'formula "treatment*block" would test both of those '
                   'effects as well as their interaction. Enclose formulae in '
                   'quotes to avoid unpleasant surprises.',
        'permutations': 'The number of permutations to be run when computing '
                        'p-values.',
        'n_jobs': 'Number of parallel processes to run.'
    },
    name='adonis PERMANOVA test for beta group significance',
    description=('Determine whether groups of samples are significantly '
                 'different from one another using the ADONIS permutation-'
                 'based statistical test in vegan-R. The function partitions '
                 'sums of squares of a multivariate data set, and is directly '
                 'analogous to MANOVA (multivariate analysis of variance). '
                 'This action differs from beta_group_significance in that it '
                 'accepts R formulae to perform multi-way ADONIS tests; '
                 'beta_group_signficance only performs one-way tests. For '
                 'more details, consult the reference manual available '
                 'on the CRAN vegan page: '
                 'https://CRAN.R-project.org/package=vegan'),
    citations=[citations['anderson2001new'], citations['Oksanen2018']]
)
