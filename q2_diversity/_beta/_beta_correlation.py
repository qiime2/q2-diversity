import skbio


def beta_correlation(ctx, metadata, distance_matrix, output_dir,
                     method="spearman", permutations=999,
                     intersect_ids=False, label1='Distance Matrix 1',
                     label2='Distance Matrix 2'):

    dist_matrix = ctx.get_action('metadata', 'distance_matrix')
    mantel = ctx.get_action('diversity', 'mantel')

    # Raise error if any samples in distance matrix are missing from metadata
    dm = distance_matrix.view(skbio.DistanceMatrix)
    ids_with_missing_metadata = set(dm.ids) - set(metadata.ids)
    if len(ids_with_missing_metadata) > 0:
        raise ValueError('All samples in distance matrix must be present '
                         'and contain data in the sample metadata. The '
                         'following samples were present in the distance '
                         'matrix, but were missing from the sample metadata '
                         'or had no data: %s' %
                         ', '.join(ids_with_missing_metadata))

    results = []
    # Convert metadata column into a distance matrix
    metadata_dist_matrix = dist_matrix(metadata)
    # dist_matrix returns a Results object, which CONTAINS a distance_matrix
    # Reach in to get the relevant property
    metadata_dist_matrix = metadata_dist_matrix.distance_matrix
    results.append(metadata_dist_matrix)

    # as above
    mantel_output = mantel(distance_matrix, metadata_dist_matrix,
                           method, permutations, intersect_ids, label1, label2)
    mantel_output = mantel_output.visualization
    results.append(mantel_output)

    return tuple(results)
