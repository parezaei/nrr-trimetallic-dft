"""
Energetics calculations for NRR on bimetallic surfaces.

Covers adsorption energies, spillover energies (SOEs),
free energy diagrams for N2 reduction pathways, and HER pathway.
All energies are in eV unless otherwise noted.
"""

import pandas as pd
import numpy as np

RY_TO_EV = 13.605703976

# DFT-computed gas-phase molecule energies (Ry), converted to eV
MOLECULE_ENERGIES_EV = {
    "0_N2":    -40.27971791  * RY_TO_EV,
    "2_N":     -19.5429571333 * RY_TO_EV,
    "1_H2":    -2.3302438141 * RY_TO_EV,
    "3_N2H":   -41.2696631674 * RY_TO_EV,
    "4_NHNH":  -42.5083642397 * RY_TO_EV,
    "5_NNH2":  -42.4405833567 * RY_TO_EV,
    "6_NNH3":  -43.4810886190 * RY_TO_EV,
    "7_NH2NH": -43.6855706118 * RY_TO_EV,
    "8_N2H4":  -44.9391323026 * RY_TO_EV,
    "9_NH2":   -22.3040791262 * RY_TO_EV,
    "10_NH":   -20.9151150555 * RY_TO_EV,
    "11_NH3":  -23.7121681317 * RY_TO_EV,
    "12_NH2NH3": 0.0,
    "13_H":    -0.9233303917 * RY_TO_EV,
}

SYSTEMS_LABEL_MAP = {
    "0_CuNi":     "CuNi",
    "1_Ag_subCu": "Ag@Cu",
    "2_Au_subCu": "Au@Cu",
    "6_Pd_subNi": "Pd@Ni",
    "7_Pt_subNi": "Pt@Ni",
    "9_Fe_subNi": "Fe@Ni",
    "11_Pd_subCu": "Pd@Cu",
    "13_Ag_subNi": "Ag@Ni",
    "14_Au_subNi": "Au@Ni",
}

SYSTEM_SITE_COLUMNS = {
    1: {
        "0_CuNi": 0,
        "1_Ag_subCu": 2,
        "2_Au_subCu": 4,
        "6_Pd_subNi": 7,
        "7_Pt_subNi": 8,
        "9_Fe_subNi": 11,
        "11_Pd_subCu": 13,
        "13_Ag_subNi": 14,
        "14_Au_subNi": 17,
    },
    2: {
        "0_CuNi": 1,
        "1_Ag_subCu": 3,
        "2_Au_subCu": 5,
        "6_Pd_subNi": 6,
        "7_Pt_subNi": 9,
        "9_Fe_subNi": 10,
        "11_Pd_subCu": 12,
        "13_Ag_subNi": 15,
        "14_Au_subNi": 16,
    },
}


def parse_qe_final_energy(output_path: str, unit: str = "ev") -> float:
    """
    Parse the final total energy from a Quantum ESPRESSO output file.

    Parameters
    ----------
    output_path : str
        Path to the QE .out file.
    unit : str
        'ev' (default) or 'ry' to return energy in eV or Rydberg.

    Returns
    -------
    float : Final total energy.
    """
    final_e = None
    with open(output_path) as f:
        for line in f:
            if "Final energy" in line:
                final_e = float(line.split()[-2])
    if final_e is None:
        raise ValueError(f"No 'Final energy' line found in {output_path}")
    return final_e * RY_TO_EV if unit == "ev" else final_e


def compute_adsorption_energy(
    slab_ads_energy: float,
    slab_energy: float,
    molecule_energy: float,
) -> float:
    """
    Compute adsorption energy: E_ads = E(slab+ads) - E(slab) - E(molecule).

    Parameters
    ----------
    slab_ads_energy : float
        DFT energy of the slab with adsorbed molecule (eV).
    slab_energy : float
        DFT energy of the clean slab (eV).
    molecule_energy : float
        DFT energy of the free molecule in gas phase (eV).

    Returns
    -------
    float : Adsorption energy in eV.
    """
    return slab_ads_energy - slab_energy - molecule_energy


def compute_soe(
    df: pd.DataFrame,
    reference_system: str,
    mol_cols: list[str],
    system_col: str = "system",
) -> pd.DataFrame:
    """
    Compute Spillover Energies (SOEs) relative to a reference system.

    SOE = E_ads(dopant) - E_ads(reference)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing adsorption energies per system and molecule.
    reference_system : str
        System label to use as the SOE reference (e.g., '0_CuNi').
    mol_cols : list[str]
        Column names for each NxHy intermediate.
    system_col : str
        Column name identifying the catalyst system.

    Returns
    -------
    pd.DataFrame with SOE values (reference subtracted).
    """
    ref_row = df[df[system_col] == reference_system][mol_cols].values
    if len(ref_row) == 0:
        raise ValueError(f"Reference system '{reference_system}' not found in DataFrame.")
    soe_df = df.copy()
    soe_df[mol_cols] = df[mol_cols].values - ref_row
    return soe_df


def _reference_energy(df_row: pd.Series, mols: dict, n_h2: float) -> float:
    """Compute the common reference energy: E(slab) + E(N2) + n_h2 * E(H2)."""
    return df_row["surf_energy"] + mols["0_N2"] + mols["1_H2"] * n_h2


def compute_reaction_path_assoc_alternating(
    df: pd.DataFrame,
    mols: dict = None,
    site_number: int = 1,
) -> pd.DataFrame:
    """
    Compute free energy levels for the associative alternating NRR pathway.

    Pathway: * + N2 + 3H2 → *N2 → *N2H → *NHNH → *NH2NH → *N2H4 →
             *NH2NH3 → *NH3 + NH3 → 2NH3 + *

    Returns a DataFrame with energy levels per system (columns) and
    reaction step labels in a 'steps' column.
    """
    if mols is None:
        mols = MOLECULE_ENERGIES_EV

    sub = df[df["ad_site_number"] == f"ad_site_{site_number}"].copy()
    ref = _reference_energy(sub, mols, 3.0)

    sub["lev1"] = (sub["surf_energy"] + mols["0_N2"] + mols["1_H2"] * 3) - ref
    sub["lev2"] = sub["0_N2"] + mols["1_H2"] * 3 - ref
    sub["lev3"] = sub["3_N2H"] + mols["1_H2"] * 2 + mols["1_H2"] * 0.5 - ref
    sub["lev4"] = sub["4_NHNH"] + mols["1_H2"] * 2 - ref
    sub["lev5"] = sub["7_NH2NH"] + mols["1_H2"] + mols["1_H2"] * 0.5 - ref
    sub["lev6"] = sub["8_N2H4"] + mols["1_H2"] - ref
    sub["lev7"] = sub["12_NH2NH3"] + mols["1_H2"] * 0.5 - ref
    sub["lev8"] = sub["11_NH3"] + mols["11_NH3"] - ref
    sub["lev9"] = mols["11_NH3"] * 2 + sub["surf_energy"] - ref

    path = sub[["system", "lev1", "lev2", "lev3", "lev4", "lev5",
                "lev6", "lev7", "lev8", "lev9"]].transpose()
    col_map = SYSTEM_SITE_COLUMNS[site_number]
    path = path.rename(columns=col_map)
    path.drop("system", axis=0, inplace=True)
    path["steps"] = [
        "* + N2 + 3H2", "*N2 + 3H2", "*N2H + 2H2 + *H", "*NHNH + 2H2",
        "*NH2NH + H2 + H", "*N2H4 + H2", "*NH2NH3 + H", "*NH3 + NH3", "2NH3 + *",
    ]
    path = path.reset_index(drop=True)
    path.index = path.index + 1
    return path


def compute_reaction_path_assoc_distal(
    df: pd.DataFrame,
    mols: dict = None,
    site_number: int = 1,
) -> pd.DataFrame:
    """
    Compute free energy levels for the associative distal NRR pathway.

    Pathway: * + N2 + 3H2 → *N2 → *N2H → *NNH2 → *NNH3 → *N + NH3 →
             *NH + NH3 → *NH2 + NH3 → *NH3 + NH3 → 2NH3 + *
    """
    if mols is None:
        mols = MOLECULE_ENERGIES_EV

    sub = df[df["ad_site_number"] == f"ad_site_{site_number}"].copy()
    ref = _reference_energy(sub, mols, 3.0)

    sub["lev1"] = (sub["surf_energy"] + mols["0_N2"] + mols["1_H2"] * 3) - ref
    sub["lev2"] = sub["0_N2"] + mols["1_H2"] * 3 - ref
    sub["lev3"] = sub["3_N2H"] + mols["1_H2"] * 2 + mols["1_H2"] * 0.5 - ref
    sub["lev4"] = sub["5_NNH2"] + mols["1_H2"] * 2 - ref
    sub["lev5"] = sub["6_NNH3"] + mols["1_H2"] + mols["1_H2"] * 0.5 - ref
    sub["lev6"] = sub["2_N"] + mols["1_H2"] + mols["1_H2"] * 0.5 + mols["11_NH3"] - ref
    sub["lev7"] = sub["10_NH"] + mols["1_H2"] + mols["11_NH3"] - ref
    sub["lev8"] = sub["9_NH2"] + mols["1_H2"] * 0.5 + mols["11_NH3"] - ref
    sub["lev9"] = sub["11_NH3"] + mols["11_NH3"] - ref
    sub["lev10"] = mols["11_NH3"] * 2 + sub["surf_energy"] - ref

    path = sub[["system", "lev1", "lev2", "lev3", "lev4", "lev5",
                "lev6", "lev7", "lev8", "lev9", "lev10"]].transpose()
    col_map = SYSTEM_SITE_COLUMNS[site_number]
    path = path.rename(columns=col_map)
    path.drop("system", axis=0, inplace=True)
    path["steps"] = [
        "* + N2 + 3H2", "*N2 + 3H2", "*N2H + 2H2 + H", "*NNH2 + 2H2",
        "*NNH3 + H2 + H", "*N + H2 + H + NH3", "*NH + H2 + NH3",
        "*NH2 + H + NH3", "*NH3 + NH3", "2NH3 + *",
    ]
    path = path.reset_index(drop=True)
    path.index = path.index + 1
    return path


def compute_reaction_path_dissociative(
    df: pd.DataFrame,
    mols: dict = None,
    site_number: int = 1,
) -> pd.DataFrame:
    """
    Compute free energy levels for the dissociative NRR pathway.

    Pathway: 1/2N2 + 3/2H2 + * → N* → NH* → NH2* → NH3* → NH3 + *
    """
    if mols is None:
        mols = MOLECULE_ENERGIES_EV

    sub = df[df["ad_site_number"] == f"ad_site_{site_number}"].copy()
    ref = sub["surf_energy"] + mols["0_N2"] * 0.5 + mols["1_H2"] * 1.5

    sub["lev1"] = (sub["surf_energy"] + mols["0_N2"] * 0.5 + mols["1_H2"] * 1.5) - ref
    sub["lev2"] = sub["2_N"] + mols["1_H2"] * 1.5 - ref
    sub["lev3"] = sub["10_NH"] + mols["1_H2"] - ref
    sub["lev4"] = sub["9_NH2"] + mols["1_H2"] * 0.5 - ref
    sub["lev5"] = sub["11_NH3"] - ref
    sub["lev6"] = mols["11_NH3"] + sub["surf_energy"] - ref

    path = sub[["system", "lev1", "lev2", "lev3", "lev4", "lev5", "lev6"]].transpose()
    col_map = SYSTEM_SITE_COLUMNS[site_number]
    path = path.rename(columns=col_map)
    path.drop("system", axis=0, inplace=True)
    path["steps"] = [
        "1/2N2 + 3/2H2 + *", "N* + 3/2H2", "NH* + H2",
        "NH2* + 1/2H2", "NH3*", "NH3 + *",
    ]
    path = path.reset_index(drop=True)
    path.index = path.index + 1
    return path


def compute_her_path(
    df: pd.DataFrame,
    mols: dict = None,
    site_number: int = 1,
) -> pd.DataFrame:
    """
    Compute the Hydrogen Evolution Reaction (HER) pathway free energy levels.

    Pathway: * + H+ + e- → *H → * + 1/2 H2
    """
    if mols is None:
        mols = MOLECULE_ENERGIES_EV

    sub = df[df["ad_site_number"] == f"ad_site_{site_number}"].copy()
    ref = sub["surf_energy"] + mols["1_H2"] * 0.5

    sub["lev1"] = (sub["surf_energy"] + mols["1_H2"] * 0.5) - ref
    sub["lev2"] = sub["13_H"] - ref
    sub["lev3"] = (sub["surf_energy"] + mols["1_H2"] * 0.5) - ref

    path = sub[["system", "lev1", "lev2", "lev3"]].transpose()
    col_map = SYSTEM_SITE_COLUMNS[site_number]
    path = path.rename(columns=col_map)
    path.drop("system", axis=0, inplace=True)
    path["steps"] = ["* + H⁺ + e⁻", "*H", "* + ½H₂"]
    path = path.reset_index(drop=True)
    path.index = path.index + 1
    return path


def compute_step_barriers(path_df: pd.DataFrame, system_cols: list[str]) -> pd.DataFrame:
    """
    Compute the energy barrier at each elementary step for all systems.

    Parameters
    ----------
    path_df : pd.DataFrame
        Output of any compute_reaction_path_* function.
    system_cols : list[str]
        Column names corresponding to the catalyst systems.

    Returns
    -------
    pd.DataFrame with step names as index and one column per system.
    """
    numeric = path_df[system_cols]
    barriers = numeric.diff().dropna()
    barriers.index = path_df["steps"].iloc[1:].values
    return barriers


def find_top_barriers(
    barriers_df: pd.DataFrame, n: int = 3
) -> pd.DataFrame:
    """
    For each catalyst system, identify the top-N largest energy barriers.

    Parameters
    ----------
    barriers_df : pd.DataFrame
        Output of compute_step_barriers (rows = steps, columns = systems).
    n : int
        Number of top barriers to report (default: 3).

    Returns
    -------
    pd.DataFrame with columns: system, barrier_{1..n}_step, barrier_{1..n}_val.
    """
    records = []
    for system in barriers_df.columns:
        sorted_steps = barriers_df[system].sort_values(ascending=False)
        row = {"system": system}
        for rank, (step, val) in enumerate(sorted_steps.head(n).items(), start=1):
            row[f"barrier_{rank}_step"] = step
            row[f"barrier_{rank}_val"] = val
        records.append(row)
    return pd.DataFrame(records)
