"""Tests for thermochemistry.py — adsorption energies, SOEs, reaction paths."""

import numpy as np
import pandas as pd
import pytest

from src.thermochemistry import (
    RY_TO_EV,
    MOLECULE_ENERGIES_EV,
    compute_adsorption_energy,
    compute_soe,
    compute_step_barriers,
    find_top_barriers,
    compute_reaction_path_assoc_alternating,
    compute_reaction_path_assoc_distal,
    compute_reaction_path_dissociative,
    compute_her_path,
    parse_qe_final_energy,
)


# ── Constants ──────────────────────────────────────────────────────────────────

def test_ry_to_ev_value():
    assert abs(RY_TO_EV - 13.605703976) < 1e-9


def test_molecule_energies_all_negative():
    # All DFT molecules (except the placeholder NH2NH3 = 0.0) should be large negative
    for key, val in MOLECULE_ENERGIES_EV.items():
        if key != "12_NH2NH3":
            assert val < 0, f"{key} energy should be negative, got {val}"


def test_molecule_energies_count():
    assert len(MOLECULE_ENERGIES_EV) == 14


# ── Adsorption energy ─────────────────────────────────────────────────────────

def test_compute_adsorption_energy_simple():
    # E_ads = -5.0 - (-4.0) - (-2.0) = 1.0
    result = compute_adsorption_energy(-5.0, -4.0, -2.0)
    assert abs(result - 1.0) < 1e-10


def test_compute_adsorption_energy_negative():
    # Favourable adsorption → negative E_ads
    result = compute_adsorption_energy(-7.0, -4.0, -2.0)
    assert result < 0


# ── Spillover energies ────────────────────────────────────────────────────────

def _make_ads_df():
    """Minimal DataFrame mimicking adsorption energy data."""
    return pd.DataFrame({
        "system": ["0_CuNi", "1_Ag_subCu", "2_Au_subCu"],
        "N2": [-0.50, -0.30, -0.70],
        "NH3": [-0.80, -0.60, -0.90],
    })


def test_compute_soe_reference_row_is_zero():
    df = _make_ads_df()
    soe = compute_soe(df, "0_CuNi", ["N2", "NH3"])
    ref_row = soe[soe["system"] == "0_CuNi"]
    assert abs(ref_row["N2"].values[0]) < 1e-12
    assert abs(ref_row["NH3"].values[0]) < 1e-12


def test_compute_soe_correct_difference():
    df = _make_ads_df()
    soe = compute_soe(df, "0_CuNi", ["N2", "NH3"])
    ag_row = soe[soe["system"] == "1_Ag_subCu"]
    assert abs(ag_row["N2"].values[0] - 0.20) < 1e-10   # -0.30 - (-0.50)
    assert abs(ag_row["NH3"].values[0] - 0.20) < 1e-10  # -0.60 - (-0.80)


def test_compute_soe_missing_reference_raises():
    df = _make_ads_df()
    with pytest.raises(ValueError, match="not found"):
        compute_soe(df, "nonexistent_system", ["N2"])


# ── Step barriers ─────────────────────────────────────────────────────────────

def _make_path_df():
    """Minimal 4-level, 2-system path DataFrame."""
    return pd.DataFrame({
        "sys_A": [0.0, -0.5, 0.3, -0.2],
        "sys_B": [0.0, 0.2, -0.1, 0.5],
        "steps": ["step1", "step2", "step3", "step4"],
    })


def test_compute_step_barriers_shape():
    df = _make_path_df()
    barriers = compute_step_barriers(df, ["sys_A", "sys_B"])
    assert barriers.shape == (3, 2)   # 4 levels → 3 diffs


def test_compute_step_barriers_values():
    df = _make_path_df()
    barriers = compute_step_barriers(df, ["sys_A", "sys_B"])
    assert abs(barriers["sys_A"].iloc[0] - (-0.5)) < 1e-10  # -0.5 - 0.0
    assert abs(barriers["sys_A"].iloc[1] - 0.8) < 1e-10     # 0.3 - (-0.5)


def test_find_top_barriers_returns_n_per_system():
    df = _make_path_df()
    barriers = compute_step_barriers(df, ["sys_A", "sys_B"])
    top = find_top_barriers(barriers, n=2)
    assert "barrier_1_step" in top.columns
    assert "barrier_2_step" in top.columns
    assert len(top) == 2  # one row per system


# ── Reaction path shapes ──────────────────────────────────────────────────────

def _make_reaction_df():
    """Minimal DataFrame that mimics corrected_df.csv for path functions."""
    systems = ["0_CuNi", "1_Ag_subCu"]
    rows = []
    mol_cols = [
        "0_N2", "2_N", "1_H2", "3_N2H", "4_NHNH", "5_NNH2",
        "6_NNH3", "7_NH2NH", "8_N2H4", "9_NH2", "10_NH",
        "11_NH3", "12_NH2NH3", "13_H",
    ]
    rng = np.random.default_rng(42)
    for site in [1, 2]:
        for sys in systems:
            row = {"system": sys, "ad_site_number": f"ad_site_{site}",
                   "surf_energy": rng.uniform(-4000, -3000)}
            for col in mol_cols:
                row[col] = rng.uniform(-800, -100)
            rows.append(row)
    return pd.DataFrame(rows)


def test_assoc_alternating_has_9_steps():
    df = _make_reaction_df()
    path = compute_reaction_path_assoc_alternating(df, site_number=1)
    assert len(path) == 9


def test_assoc_distal_has_10_steps():
    df = _make_reaction_df()
    path = compute_reaction_path_assoc_distal(df, site_number=1)
    assert len(path) == 10


def test_dissociative_has_6_steps():
    df = _make_reaction_df()
    path = compute_reaction_path_dissociative(df, site_number=1)
    assert len(path) == 6


def test_her_path_has_3_steps():
    df = _make_reaction_df()
    path = compute_her_path(df, site_number=1)
    assert len(path) == 3


def test_path_has_steps_column():
    df = _make_reaction_df()
    path = compute_reaction_path_dissociative(df, site_number=1)
    assert "steps" in path.columns


# ── QE parser ─────────────────────────────────────────────────────────────────

def test_parse_qe_final_energy_ev(tmp_path):
    qe_out = tmp_path / "calc.out"
    qe_out.write_text("Some header\n  Final energy  =    -2952.3456789  Ry\nOther text\n")
    energy_ev = parse_qe_final_energy(str(qe_out), unit="ev")
    expected = -2952.3456789 * RY_TO_EV
    assert abs(energy_ev - expected) < 1e-6


def test_parse_qe_final_energy_ry(tmp_path):
    qe_out = tmp_path / "calc.out"
    qe_out.write_text("Final energy  =    -100.0  Ry\n")
    energy_ry = parse_qe_final_energy(str(qe_out), unit="ry")
    assert abs(energy_ry - (-100.0)) < 1e-10


def test_parse_qe_final_energy_missing_raises(tmp_path):
    qe_out = tmp_path / "calc.out"
    qe_out.write_text("No energy here\n")
    with pytest.raises(ValueError, match="No 'Final energy'"):
        parse_qe_final_energy(str(qe_out))


def test_parse_qe_final_energy_takes_last_occurrence(tmp_path):
    qe_out = tmp_path / "calc.out"
    qe_out.write_text(
        "Final energy  =    -50.0  Ry\n"
        "Final energy  =    -75.0  Ry\n"
    )
    energy_ry = parse_qe_final_energy(str(qe_out), unit="ry")
    assert abs(energy_ry - (-75.0)) < 1e-10
