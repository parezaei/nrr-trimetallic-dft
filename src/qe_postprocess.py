"""
Quantum ESPRESSO output post-processing utilities.

Parses QE output files to check convergence, extract final relaxed coordinates,
and generate restart input files for continued or restarted calculations.
"""

import re
import os
import shutil
from pathlib import Path


def check_convergence(output_lines: list[str]) -> bool:
    """Return True if a QE relaxation calculation converged."""
    return any("Begin final coordinates" in line for line in output_lines)


def get_nat(input_lines: list[str]) -> int:
    """Extract the number of atoms (nat) from a QE input file."""
    for line in input_lines:
        if "nat" in line:
            return int(re.findall(r"\d+", line)[0])
    raise ValueError("Could not find 'nat' in input file.")


def extract_input_header(input_lines: list[str]) -> list[str]:
    """Return all input lines up to and including the ATOMIC_POSITIONS card."""
    for i, line in enumerate(input_lines):
        if "ATOMIC_POSITIONS" in line:
            return input_lines[: i + 1]
    raise ValueError("ATOMIC_POSITIONS card not found in input file.")


def extract_kpoints_cellpar(input_lines: list[str]) -> list[str]:
    """Return input lines from K_POINTS onwards (k-points and cell parameters)."""
    for i, line in enumerate(input_lines):
        if "K_POINTS" in line:
            return input_lines[i:]
    raise ValueError("K_POINTS card not found in input file.")


def extract_final_coordinates(output_lines: list[str], nat: int) -> list[str]:
    """
    Return the final relaxed atomic coordinates from a QE output file.

    Finds the last occurrence of 'ATOMIC_POSITIONS (crystal)' and returns
    the following `nat` lines.
    """
    coord_indices = [
        i
        for i, line in enumerate(output_lines)
        if "ATOMIC_POSITIONS (crystal)" in line
    ]
    if not coord_indices:
        raise ValueError("No 'ATOMIC_POSITIONS (crystal)' block found in output.")
    start = max(coord_indices) + 1
    return output_lines[start : start + nat]


def generate_restart_input(
    input_lines: list[str], output_lines: list[str]
) -> list[str]:
    """
    Build a restart QE input by combining the original namelist/tags with
    the final relaxed coordinates from the output file.

    Sets restart_mode = 'restart' in the &CONTROL namelist.
    """
    nat = get_nat(input_lines)
    header = extract_input_header(input_lines)
    final_coords = extract_final_coordinates(output_lines, nat)
    kpoints_cellpar = extract_kpoints_cellpar(input_lines)

    restart_header = []
    for line in header:
        if "restart_mode" in line:
            line = re.sub(r"restart_mode.*", "restart_mode = 'restart',", line)
        restart_header.append(line)

    return restart_header + final_coords + ["\n\n"] + kpoints_cellpar


def process_site_range(
    main_dir: str,
    site_min: int,
    site_max: int,
    file_format: str = "ad_site_",
    overwrite_original: bool = True,
) -> dict[str, bool]:
    """
    Generate restart input files for a range of adsorption site calculations.

    For each site index in [site_min, site_max], reads `{file_format}{n}.in`
    and `{file_format}{n}.out`, checks convergence, and writes a restart file.

    Parameters
    ----------
    main_dir : str
        Directory containing the .in and .out files.
    site_min : int
        First site index (inclusive).
    site_max : int
        Last site index (inclusive).
    file_format : str
        Filename prefix (default: 'ad_site_').
    overwrite_original : bool
        If True, copy the restart file back to main_dir and remove the temp copy.

    Returns
    -------
    dict mapping site label → convergence status (True/False).
    """
    main_dir = Path(main_dir)
    final_coord_dir = main_dir / "final_coord"
    final_coord_dir.mkdir(exist_ok=True)

    results = {}
    for n in range(site_min, site_max + 1):
        label = f"{file_format}{n}"
        in_path = main_dir / f"{label}.in"
        out_path = main_dir / f"{label}.out"

        with open(in_path) as f:
            input_lines = f.readlines()
        with open(out_path) as f:
            output_lines = f.readlines()

        converged = check_convergence(output_lines)
        results[label] = converged

        if converged:
            print(f"{label}: Converged — writing restart file.")
            restart_lines = generate_restart_input(input_lines, output_lines)
            restart_path = final_coord_dir / f"RP_{label}.in"
            with open(restart_path, "w") as f:
                f.writelines(restart_lines)
            if overwrite_original:
                shutil.copy(restart_path, main_dir / f"RP_{label}.in")
                restart_path.unlink()
        else:
            print(f"{label}: NOT converged — skipping.")

    return results


def process_molecules_wholesale(
    base_dir: str,
    mol_dirs: list[str],
    site_min: int = 1,
    site_max: int = 2,
    file_format: str = "ad_site_",
) -> dict[str, dict[str, bool]]:
    """
    Run process_site_range across multiple molecule subdirectories.

    Parameters
    ----------
    base_dir : str
        Parent directory containing the molecule subdirectories.
    mol_dirs : list[str]
        List of molecule subdirectory names (e.g., ['2_N', '4_NHNH', ...]).
    site_min, site_max : int
        Site range to process in each molecule directory.
    file_format : str
        Input/output filename prefix.

    Returns
    -------
    Nested dict: mol_dir → {site_label → converged (bool)}.
    """
    base_dir = Path(base_dir)
    all_results = {}
    for mol in mol_dirs:
        mol_path = base_dir / mol
        print(f"\n--- Processing molecule: {mol} ---")
        all_results[mol] = process_site_range(
            mol_path, site_min, site_max, file_format
        )
    return all_results


def generate_rerun_input(
    main_dir: str,
    site_min: int,
    site_max: int,
    file_format: str = "ad_site_",
) -> None:
    """
    Prepare inputs for re-running unconverged QE calculations.

    Reads the latest partial output, extracts the last available coordinates,
    and writes a new input file to `final_coord/` for re-submission.
    """
    main_dir = Path(main_dir)
    final_coord_dir = main_dir / "final_coord"
    final_coord_dir.mkdir(exist_ok=True)

    for n in range(site_min, site_max + 1):
        label = f"{file_format}{n}"
        in_path = main_dir / f"{label}.in"
        out_path = main_dir / f"{label}.out"

        with open(in_path) as f:
            input_lines = f.readlines()
        with open(out_path) as f:
            output_lines = f.readlines()

        restart_lines = generate_restart_input(input_lines, output_lines)
        restart_path = final_coord_dir / f"{label}.in"
        with open(restart_path, "w") as f:
            f.writelines(restart_lines)

        shutil.copy(restart_path, main_dir / f"{label}.in")
        restart_path.unlink()
        print(f"{label}: Rerun input written.")
