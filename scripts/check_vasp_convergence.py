#!/usr/bin/env python3
"""
CLI script: check VASP calculation convergence and prepare rerun inputs.

Scans OUTCAR files under a given material/molecule directory tree,
reports converged/not-converged status, and optionally updates INCAR +
POSCAR for resubmission.

Usage
-----
    python check_vasp_convergence.py --material_dir /scratch/ML_project/Satvik/0_Al13Fe4 \
        --mols 2_N2 3_H2 --cluster narval [--rerun]
"""

import argparse
import os
import re
import shutil
import filecmp
from pathlib import Path


CONVERGENCE_MARKER = "reached required accuracy - stopping structural energy minimisation"


def check_vasp_outcar(outcar_path: str) -> bool:
    """Return True if VASP OUTCAR indicates convergence."""
    with open(outcar_path) as f:
        return any(CONVERGENCE_MARKER in line for line in f)


def prepare_vasp_rerun(site_dir: Path, cluster_name: str, rerun: bool = True) -> None:
    """
    Update INCAR for restart (ISTART=1, ICHARG=1), copy CONTCAR → POSCAR,
    and optionally resubmit the SLURM job.
    """
    incar_path = site_dir / "INCAR"
    contcar_path = site_dir / "CONTCAR"
    outcar_path = site_dir / "OUTCAR"

    with open(incar_path) as f:
        lines = f.readlines()

    updated = []
    for line in lines:
        if "ISTART" in line:
            line = re.sub(r"ISTART.*", "ISTART = 1", line)
        if "ICHARG" in line:
            line = re.sub(r"ICHARG.*", "ICHARG = 1", line)
        updated.append(line)

    incar_path.unlink()
    with open(incar_path, "w") as f:
        f.writelines(updated)

    shutil.copy(contcar_path, site_dir / "POSCAR")

    if rerun:
        outcars = sorted(
            [f for f in os.listdir(site_dir) if re.match(r"OUTCAR\.\d+", f)],
            key=lambda x: int(re.findall(r"\d+", x)[0]),
        )
        if not outcars:
            shutil.copy(outcar_path, site_dir / "OUTCAR.1")
            outcars = ["OUTCAR.1"]
        last_n = int(re.findall(r"\d+", outcars[-1])[0])
        backup = site_dir / f"OUTCAR.{last_n + 1}"
        if not filecmp.cmp(outcar_path, site_dir / outcars[-1]):
            shutil.copy(outcar_path, backup)
        os.chdir(site_dir)
        os.system(f"sbatch submitVASP_{cluster_name}.sh")


def main():
    parser = argparse.ArgumentParser(
        description="Check VASP convergence and prepare rerun inputs."
    )
    parser.add_argument(
        "--material_dir", required=True,
        help="Root directory for a single material (contains molecule subdirs).",
    )
    parser.add_argument(
        "--mols", nargs="+", default=["2_N2", "3_H2"],
        help="Molecule subdirectory names to check.",
    )
    parser.add_argument(
        "--cluster", default="narval",
        help="Cluster name used in SLURM batch script filename (default: narval).",
    )
    parser.add_argument(
        "--rerun", action="store_true",
        help="If set, resubmit not-converged calculations.",
    )
    args = parser.parse_args()

    material_dir = Path(args.material_dir)
    not_converged = []

    for mol in args.mols:
        mol_dir = material_dir / mol
        if not mol_dir.exists():
            print(f"  Skipping {mol} — directory not found.")
            continue
        sites = sorted(
            [d for d in mol_dir.iterdir() if d.is_dir() and "ad_site" in d.name],
            key=lambda d: int(re.findall(r"\d+", d.name)[0]),
        )
        for site_dir in sites:
            outcar = site_dir / "OUTCAR"
            if not outcar.exists():
                print(f"  {mol}/{site_dir.name}: NO OUTCAR")
                continue
            converged = check_vasp_outcar(outcar)
            status = "CONVERGED" if converged else "NOT CONVERGED"
            print(f"  {mol}/{site_dir.name}: {status}")
            if not converged:
                not_converged.append(f"{mol}/{site_dir.name}")
                prepare_vasp_rerun(site_dir, args.cluster, rerun=args.rerun)

    print(f"\nNot converged: {not_converged}")


if __name__ == "__main__":
    main()
