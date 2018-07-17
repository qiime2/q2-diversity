# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def beta_correlation(ctx, metadata, distance_matrix, output_dir,
                     method="spearman", permutations=999,
                     intersect_ids=False, label1='Metadata',
                     label2='Distance Matrix'):

    convert_to_dist_matrix = ctx.get_action('metadata', 'distance_matrix')
    mantel = ctx.get_action('diversity', 'mantel')

    results = []
    metadata_dist_matrix, = convert_to_dist_matrix(metadata)
    results.append(metadata_dist_matrix)
    mantel_visualization, = mantel(distance_matrix, metadata_dist_matrix,
                                   method, permutations, intersect_ids, label1,
                                   label2)
    results.append(mantel_visualization)

    return tuple(results)
