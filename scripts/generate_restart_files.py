#!/usr/bin/env python3
"""
CLI script: generate restart input files for converged QE calculations.

Usage
-----
Single directory (interactive site range):
    python generate_restart_files.py --dir ./0_CuNi/2_N --min 1 --max 2

All molecule subdirectories in a catalyst folder:
    python generate_restart_files.py --dir ./0_CuNi --wholesale --min 1 --max 2
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.qe_postprocess import process_site_range, process_molecules_wholesale

MOLECULE_DIRS = [
    "2_N", "3_N2H", "4_NHNH", "5_NNH2", "6_NNH3",
    "7_NH2NH", "8_N2H4", "9_NH2", "10_NH", "11_NH3", "12_NH2NH3",
]


def main():
    parser = argparse.ArgumentParser(
        description="Generate QE restart input files from converged calculations."
    )
    parser.add_argument(
        "--dir", required=True,
        help="Target directory (molecule dir for single mode, catalyst dir for --wholesale).",
    )
    parser.add_argument(
        "--min", type=int, default=1, dest="site_min",
        help="First adsorption site index (default: 1).",
    )
    parser.add_argument(
        "--max", type=int, default=2, dest="site_max",
        help="Last adsorption site index (default: 2).",
    )
    parser.add_argument(
        "--format", default="ad_site_", dest="file_format",
        help="Input/output filename prefix (default: 'ad_site_').",
    )
    parser.add_argument(
        "--wholesale", action="store_true",
        help="Process all molecule subdirectories under --dir.",
    )
    args = parser.parse_args()

    if args.wholesale:
        results = process_molecules_wholesale(
            base_dir=args.dir,
            mol_dirs=MOLECULE_DIRS,
            site_min=args.site_min,
            site_max=args.site_max,
            file_format=args.file_format,
        )
        total = sum(len(v) for v in results.values())
        converged = sum(sum(v.values()) for v in results.values())
        print(f"\nDone. {converged}/{total} calculations converged.")
    else:
        results = process_site_range(
            main_dir=args.dir,
            site_min=args.site_min,
            site_max=args.site_max,
            file_format=args.file_format,
        )
        converged = sum(results.values())
        print(f"\nDone. {converged}/{len(results)} calculations converged.")


if __name__ == "__main__":
    main()
