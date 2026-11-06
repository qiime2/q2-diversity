# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import importlib

import numpy as np
import pandas as pd
import skbio
from patsy import ModelDesc

import qiime2
from qiime2.plugin import get_available_cores
import q2templates

from ._visualizer import _run_command


TEMPLATES = importlib.resources.files('q2_diversity') / '_beta'


def adonis(output_dir: str,
           distance_matrix: skbio.DistanceMatrix,
           metadata: qiime2.Metadata,
           formula: str,
           permutations: int = 999,
           n_jobs: int = 1) -> None:
    if n_jobs == 0:
        n_jobs = get_available_cores()

    metadata = _validate_and_filter_metadata(metadata, distance_matrix)

    _validate_formula(formula, metadata)

    results = _run_adonis(
        distance_matrix, metadata, formula, permutations, n_jobs
    )

    results.to_csv(
        os.path.join(output_dir, 'adonis.tsv'),
        sep='\t',
        index_label=False
    )
    results = q2templates.df_to_html(results)
    index = os.path.join(TEMPLATES, 'adonis_assets', 'index.html')
    q2templates.render(index, output_dir, context={'results': results})


def blocked_adonis(
    output_dir: str,
    distance_matrix: skbio.DistanceMatrix,
    metadata: qiime2.Metadata,
    formula: str,
    permutation_group_column: str,
    permutations: int = 999,
    n_jobs: int = 1
) -> None:
    '''
    '''
    if n_jobs == 0:
        n_jobs = get_available_cores()

    metadata = _validate_and_filter_metadata(metadata, distance_matrix)
    _validate_formula(formula, metadata, permutation_group_column)

    terms = ModelDesc.from_formula(formula)
    formula_columns = [i.name() for t in terms.rhs_termlist for i in t.factors]
    target = formula_columns[0]

    observed = _run_adonis(
        distance_matrix, metadata, formula, 1, n_jobs
    )
    observed_f_stat = _extract_f_statistic(observed, target)

    f_statistics = []
    for _ in range(permutations):
        shuffled_md = _permute_target_column(
            metadata, formula_columns[0], permutation_group_column
        )
        results = _run_adonis(
            distance_matrix, shuffled_md, formula, 1, n_jobs
        )

        f_stat = _extract_f_statistic(results, target)
        f_statistics.append(f_stat)

    # Count numerical ties using the same tolerance as vegan.
    tolerance = np.sqrt(np.finfo(float).eps)
    num_sig_fs = (
        pd.Series(f_statistics) >= observed_f_stat - tolerance
    ).sum()
    new_p = (num_sig_fs + 1) / (permutations + 1)

    observed.loc[target, 'Pr(>F)'] = new_p

    results = observed

    results.to_csv(
        os.path.join(output_dir, 'adonis.tsv'),
        sep='\t',
        index_label=False
    )
    results = q2templates.df_to_html(results)
    index = os.path.join(TEMPLATES, 'adonis_assets', 'index.html')
    q2templates.render(index, output_dir, context={'results': results})


def _run_adonis(
    distance_matrix: skbio.DistanceMatrix,
    metadata: qiime2.Metadata,
    formula: str,
    permutations: int,
    n_jobs: int,
) -> pd.DataFrame:
    with tempfile.TemporaryDirectory() as temp_dir:
        dm_fp = os.path.join(temp_dir, 'dm.tsv')
        distance_matrix.write(dm_fp)
        md_fp = os.path.join(temp_dir, 'md.tsv')
        metadata.save(md_fp)

        results_fp = os.path.join(temp_dir, 'adonis.tsv')

        cmd = ['run_adonis.R', dm_fp, md_fp, formula, str(permutations),
               str(n_jobs), results_fp]

        _run_command(cmd)

        return pd.read_csv(results_fp, sep='\t')


def _validate_and_filter_metadata(
    metadata: qiime2.Metadata, distance_matrix: skbio.DistanceMatrix
) -> qiime2.Metadata:
    '''
    Validate that distance matrix IDs are a subset of the metadata IDs. Filter
    metadata to contain only IDs present in the distance matrix.

    Parameters
    ----------
    metadata : qiime2.Metadata
        The sample metadata.
    distance_matrix : skbio.DistanceMatrix
        The beta diversity distance matrix.

    Returns
    -------
    qiime2.Metadata
        The sample metadata filtered to contain only IDs found in the distance
        matrix.

    Raises
    ------
    ValueError
        When there are sample IDs in the distance matrix that are not found
        in the metadata.
    '''
    metadata_ids = set(metadata.ids)
    dm_ids = set(distance_matrix.ids)
    missing_ids = dm_ids.difference(metadata_ids)
    if len(missing_ids) > 0:
        raise ValueError('Missing samples in metadata: %r' % missing_ids)

    filtered_md = metadata.to_dataframe().reindex(distance_matrix.ids)
    filtered_md.index.name = 'sample-id'
    return qiime2.Metadata(filtered_md)


def _validate_formula(
    formula: str,
    metadata: qiime2.Metadata,
    permutation_group_column: str = None
) -> None:
    '''
    Ensure all columns in the formula are in the metadata and that each has
    no missing values. If the permutation group column is provided also ensure
    that there is only one term in the formula.

    Parameters
    ----------
    formula : str
        The model formula.
    metadata : qiime2.Metadata
        The sample metadata filtered to the distance matrix.
    permuatation_group_column : str | None
        The grouping column if using the blocked permutation approach.

    Raises
    ------
    ValueError
        If `permutation_group_column` is not None and there is more than
        one term in the formula.
    ValueError
        If a term in the formula is not a column in the metadata.
    ValueError
        If a metadata column has missing values.
    '''
    terms = ModelDesc.from_formula(formula)
    formula_columns = [i.name() for t in terms.rhs_termlist for i in t.factors]

    if permutation_group_column is not None:
        if len(set(formula_columns)) != 1:
            raise ValueError(
                'Unit-level label permutations require a formula with exactly '
                'one metadata column.'
            )

        formula_columns.append(permutation_group_column)

    for column_name in formula_columns:
        column = metadata.get_column(column_name)
        if column.has_missing_values():
            raise ValueError(
                'adonis requires metadata columns with no NaN values (missing '
                f'values in column {column.name}).'
            )


def _permute_target_column(
    metadata: qiime2.Metadata, target: str, group: str
) -> qiime2.Metadata:
    '''
    Shuffle groups of values in the `target` column defined by `group`.

    Parameters
    ----------
    metadata : qiime2.Metadata
        The sample metadata that contains the group and target columns.
    target : str
        The variable of interest the values of which will be shuffled.
    group : str
        The column specifying the groupings of the `target` values.

    Returns
    -------
    qiime2.Metadata
        The shuffled metadata.

    Raises
    ------
    ValueError
        If a group of target values contains more than one unique value.
    '''
    md = metadata.to_dataframe()
    grouped = md.groupby(group, sort=False)[target]

    if (grouped.nunique(dropna=False) != 1).any():
        raise ValueError(
            f'The target column {target} must be constant within each group.'
        )

    shuffled_target_values = pd.Series(
        grouped.first().sample(frac=1).tolist(), index=grouped.first().index
    )
    md[target] = md[group].map(shuffled_target_values)

    return qiime2.Metadata(md)


def _extract_f_statistic(results: pd.DataFrame, target: str) -> float:
    '''

    '''
    if target not in results.index or pd.isna(results.loc[target, 'F']):
        raise ValueError('An undefined F statistic was created.')

    return results.loc[target, 'F']
