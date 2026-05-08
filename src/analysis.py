"""
Statistical analysis utilities for NRR surface science data.

Includes correlation analysis, linear regression, and feature normalization
used to explore structure-property relationships between electronic descriptors
(Bader charge, electronegativity difference) and adsorption/SOE energies.
"""

import numpy as np
import pandas as pd
import scipy.stats


def linear_regression(x: list | np.ndarray, y: list | np.ndarray) -> dict:
    """
    Perform ordinary least-squares linear regression.

    Parameters
    ----------
    x, y : array-like
        Independent and dependent variables.

    Returns
    -------
    dict with keys: slope, intercept, r2, p_value, std_err, equation.
    """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    r2 = r_value ** 2
    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "p_value": p_value,
        "std_err": std_err,
        "equation": f"y = {slope:.2f}x + {intercept:.2f}  (R² = {r2:.2f})",
    }


def compute_pearson_correlation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute Pearson correlation matrix for numeric columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input data (numeric columns only; non-numeric columns are ignored).

    Returns
    -------
    pd.DataFrame : Symmetric correlation matrix.
    """
    return df.select_dtypes(include="number").corr()


def normalize_minmax(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply min-max normalization: x' = (x - min) / (max - min).

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame to normalize.

    Returns
    -------
    pd.DataFrame with values scaled to [0, 1].
    """
    numeric = df.select_dtypes(include="number")
    return (numeric - numeric.min()) / (numeric.max() - numeric.min())


def normalize_zscore(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply z-score (mean) normalization: x' = (x - mean) / std.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame to normalize.

    Returns
    -------
    pd.DataFrame with zero mean and unit variance per column.
    """
    numeric = df.select_dtypes(include="number")
    return (numeric - numeric.mean()) / numeric.std()


def build_feature_dataframe(
    df: pd.DataFrame,
    systems: list[str],
    dopant_charges: dict[str, float],
    electroneg_diffs: dict[str, float],
    system_col: str = "system",
) -> pd.DataFrame:
    """
    Attach electronic descriptor features to an adsorption energy DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a system identifier column and energy columns.
    systems : list[str]
        System labels to include (in the display-name format, e.g., 'Ag@Cu').
    dopant_charges : dict
        Mapping from system label → Bader charge of the dopant atom (e).
    electroneg_diffs : dict
        Mapping from system label → electronegativity difference X(dopant) - X(host).
    system_col : str
        Name of the system identifier column.

    Returns
    -------
    pd.DataFrame with 'dopant_charge' and 'electroneg_diff' columns appended.
    """
    df = df.copy()
    df["dopant_charge"] = df[system_col].map(dopant_charges)
    df["electroneg_diff"] = df[system_col].map(electroneg_diffs)
    return df


def pathway_selectivity_score(barriers_df: pd.DataFrame) -> pd.Series:
    """
    Compute a simple selectivity proxy: the ratio of the NRR limiting step
    barrier to the HER limiting step barrier for each system.

    A score > 1 suggests NRR is kinetically penalised relative to HER.

    Parameters
    ----------
    barriers_df : pd.DataFrame
        DataFrame with systems as columns and 'nrr_barrier' / 'her_barrier'
        as row index.

    Returns
    -------
    pd.Series with selectivity scores per system.
    """
    if "nrr_barrier" not in barriers_df.index or "her_barrier" not in barriers_df.index:
        raise ValueError("barriers_df must have 'nrr_barrier' and 'her_barrier' rows.")
    return barriers_df.loc["nrr_barrier"] / barriers_df.loc["her_barrier"]
