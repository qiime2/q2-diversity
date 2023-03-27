# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2

import numpy as np
import pandas as pd

from skbio import OrdinationResults
from scipy.spatial import procrustes
from scipy.linalg import orthogonal_procrustes
from numpy.random import default_rng


def procrustes_analysis(reference: OrdinationResults, other: OrdinationResults,
                        dimensions: int = 5,
                        permutations: int = 999) -> (OrdinationResults,
                                                     OrdinationResults,
                                                     pd.DataFrame):

    if reference.samples.shape != other.samples.shape:
        raise ValueError('The matrices cannot be fitted unless they have the '
                         'same dimensions')

    if reference.samples.shape[1] < dimensions:
        raise ValueError('Cannot fit fewer dimensions than available')

    # fail if there are any elements in the symmetric difference
    diff = reference.samples.index.symmetric_difference(other.samples.index)
    if not diff.empty:
        raise ValueError('The ordinations represent two different sets of '
                         'samples')

    # make the matrices be comparable
    other.samples = other.samples.reindex(index=reference.samples.index)
    mtx1, mtx2, m2 = procrustes(reference.samples.values[:, :dimensions],
                                other.samples.values[:, :dimensions])

    axes = reference.samples.columns[:dimensions]
    samples1 = pd.DataFrame(data=mtx1,
                            index=reference.samples.index.copy(),
                            columns=axes.copy())
    samples2 = pd.DataFrame(data=mtx2,
                            index=reference.samples.index.copy(),
                            columns=axes.copy())

    info = _procrustes_monte_carlo(reference.samples.values[:, :dimensions],
                                   other.samples.values[:, :dimensions],
                                   m2, permutations)

    out1 = OrdinationResults(
            short_method_name=reference.short_method_name,
            long_method_name=reference.long_method_name,
            eigvals=reference.eigvals[:dimensions].copy(),
            samples=samples1,
            features=reference.features,
            biplot_scores=reference.biplot_scores,
            sample_constraints=reference.sample_constraints,
            proportion_explained=reference.proportion_explained[:dimensions]
            .copy())
    out2 = OrdinationResults(
            short_method_name=other.short_method_name,
            long_method_name=other.long_method_name,
            eigvals=other.eigvals[:dimensions].copy(),
            samples=samples2,
            features=other.features,
            biplot_scores=other.biplot_scores,
            sample_constraints=other.sample_constraints,
            proportion_explained=other.proportion_explained[:dimensions]
            .copy())
    return out1, out2, info


def _procrustes_monte_carlo(reference: np.ndarray, other: np.ndarray,
                            true_m2, permutations) -> (pd.DataFrame):
    '''
    Outputs a dataframe containing:
    0: True M^2 value
    1: p-value for true M^2 value
    2: number of Monte Carlo permutations done in simulation
    '''

    rng = default_rng()

    trials_below_m2 = 0

    if permutations == 'disable':
        permutations = 0

    for i in range(permutations):

        # shuffle rows in np array
        rng.shuffle(other)

        # run procrustes analysis
        _, _, m2 = procrustes(reference, other)

        # check m2 value
        if m2 < true_m2:
            trials_below_m2 += 1

    if permutations == 0:
        p_val = np.nan
    else:
        # mimic the behaviour in scikit-bio's permutation-based tests and avoid
        # returning p-values equal to zero
        p_val = (trials_below_m2 + 1) / (permutations + 1)

    df = pd.DataFrame({'true M^2 value': [true_m2],
                       'p-value for true M^2 value': [p_val],
                       'number of Monte Carlo permutations': [permutations]},
                      index=pd.Index(['results'], name='id'))

    return df


def partial_procrustes(reference: OrdinationResults, other: OrdinationResults,
                       pairing: qiime2.CategoricalMetadataColumn,
                       dimensions: int = 5) -> OrdinationResults:
    if reference.samples.shape[1] < dimensions:
        raise ValueError('Cannot fit fewer dimensions than available')

    if other.samples.shape[1] < dimensions:
        raise ValueError('Cannot fit fewer dimensions than available')

    pairing = pairing.to_series()
    pairing = pairing[~pairing.isnull()]

    if len(pairing) == 0:
        raise ValueError('The metadata are lacking paired samples')

    ref_pairs = sorted(set(pairing.index) & set(reference.samples.index))
    other_pairs = sorted(set(pairing.index) & set(other.samples.index))

    if len(ref_pairs) == 0:
        raise ValueError('The reference frame lacks paired samples')

    if len(other_pairs) == 0:
        raise ValueError('The other frame lacks paired samples')

    ref_order = ref_pairs
    other_order = pairing.loc[ref_pairs].values

    ref_df, other_df = _partial_procrustes(reference.samples,
                                           other.samples,
                                           ref_order, other_order)

    out = OrdinationResults(
            short_method_name=reference.short_method_name,
            long_method_name=reference.long_method_name,
            eigvals=reference.eigvals[:dimensions].copy(),
            samples=pd.concat([ref_df, other_df]),
            features=reference.features,
            biplot_scores=reference.biplot_scores,
            sample_constraints=reference.sample_constraints,
            proportion_explained=reference.proportion_explained[:dimensions]
            .copy())
    return out

def _deconstructed_procrustes(mtx1, mtx2):
    # Derived from scipy procrustes
    # https://github.com/scipy/scipy/blob/d541c752246a9e196034957d3e044950eec75907/scipy/spatial/_procrustes.py#L100-L125
    mtx1 = mtx1.copy()
    mtx2 = mtx2.copy()

    # translate all the data to the origin
    mtx1_translate = np.mean(mtx1, 0)
    mtx2_translate = np.mean(mtx2, 0)
    mtx1 -= mtx1_translate
    mtx2 -= mtx2_translate

    # uniform scaling
    norm1 = np.linalg.norm(mtx1)
    norm2 = np.linalg.norm(mtx2)

    if norm1 == 0 or norm2 == 0:
        raise ValueError('Input matrices must contain >1 unique points')

    # change scaling of data (in rows) such that trace(mtx*mtx') = 1
    mtx1 /= norm1
    mtx2 /= norm2

    R, s = orthogonal_procrustes(mtx1, mtx2)

    return mtx1_translate, mtx2_translate, norm1, norm2, R, s


def _partial_procrustes(df_mtx1, df_mtx2, df_mtx1_pair_ids, df_mtx2_pair_ids):
    df_mtx1 = df_mtx1.copy()
    df_mtx2 = df_mtx2.copy()

    # pull paired samples
    paired_mtx1 = df_mtx1.loc[df_mtx1_pair_ids]
    paired_mtx2 = df_mtx2.loc[df_mtx2_pair_ids]

    # compute procrustes on paired data
    results = _deconstructed_procrustes(paired_mtx1, paired_mtx2)
    mtx1_translate, mtx2_translate, norm1, norm2, R, s = results

    # transform both full input matrices
    df_mtx1 -= mtx1_translate
    df_mtx2 -= mtx2_translate
    df_mtx1 /= norm1
    df_mtx2 /= norm2

    # rotate mtx2 to orient relative to mtx1 (derived from scipy procrustes)
    df_mtx2_mat = np.dot(df_mtx2.values, R.T) * s
    df_mtx2 = pd.DataFrame(df_mtx2_mat, columns=df_mtx2.columns,
                           index=df_mtx2.index)

    return df_mtx1, df_mtx2
