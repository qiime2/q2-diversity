# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2


# PD Mice Data
pd_alpha_div_faith_pd_url = 'https://data.qiime2.org/usage-examples/' \
                            'pd-mice/core-metrics-results/faith_pd_vector.qza'

pd_metadata_url = (f'https://data.qiime2.org/{qiime2.__release__}/tutorials/'
                   'pd-mice/sample_metadata.tsv')

# Moving Pictures Data
mp_beta_div_jaccard_url = 'data.qiime2.org/usage-examples/' \
                          'moving-pictures/core-metrics-results/' \
                          'jaccard_distance_matrix.qza'

mp_metadata_url = (f'https://data.qiime2.org/{qiime2.__release__}/tutorials/'
                   'moving-pictures/sample_metadata.tsv')


def mp_metadata_column(use):
    md = use.init_metadata_from_url('metadata', mp_metadata_url)
    md_column = use.get_metadata_column('metadata_column', 'month', md)
    return md_column


# Alpha Diversity examples
def alpha_group_significance_faith_pd(use):
    alpha_div_faith_pd = use.init_artifact_from_url('alpha_div_faith_pd',
                                                    pd_alpha_div_faith_pd_url)
    metadata = use.init_metadata_from_url('metadata', pd_metadata_url)

    viz, = use.action(
        use.UsageAction('diversity', 'alpha_group_significance'),
        use.UsageInputs(
            alpha_diversity=alpha_div_faith_pd,
            metadata=metadata
        ),
        use.UsageOutputNames(
            visualization='visualization'
        )
    )

    viz.assert_output_type('Visualization')


def alpha_correlation_faith_pd(use):
    alpha_div_faith_pd = use.init_artifact_from_url('alpha_div_faith_pd',
                                                    pd_alpha_div_faith_pd_url)
    metadata = use.init_metadata_from_url('metadata', pd_metadata_url)

    viz, = use.action(
        use.UsageAction('diversity', 'alpha_correlation'),
        use.UsageInputs(
            alpha_diversity=alpha_div_faith_pd,
            metadata=metadata
        ),
        use.UsageOutputNames(
            visualization='visualization'
        )
    )

    viz.assert_output_type('Visualization')


# Beta Diversity examples
def beta_correlation_jaccard(use):
    beta_div_jaccard = use.init_artifact_from_url('beta_div_jaccard',
                                                  mp_beta_div_jaccard_url)
    metadata_column = mp_metadata_column(use)

    md_distance_matrix, mantel_scatter_viz = use.action(
        use.UsageAction('diversity', 'beta_correlation'),
        use.UsageInputs(
            distance_matrix=beta_div_jaccard,
            metadata=metadata_column,
            intersect_ids=True
        ),
        use.UsageOutputNames(
            metadata_distance_matrix='metadata_distance_matrix',
            mantel_scatter_visualization='mantel_scatter_visualization'
        )
    )

    md_distance_matrix.assert_output_type('DistanceMatrix')
    mantel_scatter_viz.assert_output_type('Visualization')
