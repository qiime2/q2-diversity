# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np 
import pandas as pd

from skbio import OrdinationResults
from scipy.spatial import procrustes
from numpy.random import default_rng

def procrustes_analysis(reference: OrdinationResults, other: OrdinationResults,
                        dimensions: int = 5) -> (OrdinationResults,
                                                 OrdinationResults,
                                                 pd.DataFrame):

    if reference.samples.shape != other.samples.shape:
        raise ValueError('The matrices cannot be fitted unless they have the '
                         'same dimensions')

    if reference.samples.shape[1] < dimensions:
        raise ValueError('Cannot fit fewer dimensions than available')

    # fail if there are any elements in the symmetric difference
    if not (reference.samples.index ^ other.samples.index).empty:
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
                                       m2,
                                       1000)

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

def _procrustes_monte_carlo(reference: np.ndarray,
                            other: np.ndarray,
                            true_m2, 
                            trials=1000) -> (pd.DataFrame):
    '''Adapted from QIIME 1's transform_coordinate_matrices
    1. set up axes for resulting dataframe 
    2. iterating through total # of permutations 
    2a. shuffle rows of the other matrix 
    2b. run procrustes analysis & store the resulting values 
    2c. store M^2 results into final dataframe 

    separate dataframe 
    true m^2 value 
    p value (# of m^2 trials < true m^2 / n trials)

    remove dimensions from function signature & require that dimensions for ref
    & other are the same 

    '''
    
    df = pd.DataFrame()
    rng = default_rng()

    results = []
    trials_below_m2 = 0
    
    for i in range(trials):

        # shuffle rows in np array 
        rng.shuffle(other)

        # run procrustes analysis
        _, _, m2 = procrustes(reference, other)
        
        # check m2 value 
        if m2 < true_m2: 
            trials_below_m2 += 1
    
    p_val = trials_below_m2 / trials
    results.append(true_m2) 
    results.append(p_val)
    df['Procrustes Results'] = results 

    return df 
