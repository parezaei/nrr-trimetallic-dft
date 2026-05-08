# Designing Trimetallic Single-Doped Alloy Catalysts for Sustainable Ammonia Production

**Code for:** *Designing Trimetallic Single-Doped Alloy Catalysts for Sustainable Ammonia Production: The Role of Dopants in Active Site Engineering*
**Published in:** ACS Catalysis — [https://doi.org/10.1021/acscatal.5c01049](https://doi.org/10.1021/acscatal.5c01049)

**Author:** Parastoo Agharezaei
**Advisor:** Prof. Kulbir K. Ghuman
**Laboratory:** [Insilicomatters Laboratory](https://github.com/insilicomatters), INRS-EMT

This repository provides the Python analysis toolkit and post-processing scripts
used in our computational study of nitrogen reduction reaction (NRR) activity on
trimetallic single-doped alloy surfaces. We screened nine Cu/Ni-based catalysts — including
single-atom dopant systems (Ag@Cu, Au@Cu, Pd@Ni, Pt@Ni, Fe@Ni, Pd@Cu, Ag@Ni,
Au@Ni) — using density functional theory (DFT) and compared three NRR pathways:
associative alternating, associative distal, and dissociative.

---

## Key Findings

- Identified electronic descriptors (Bader charge, electronegativity difference)
  that correlate with surface occupation energies (SOEs) for NxHy intermediates
- Compared free energy barriers across three N2 reduction pathways and the
  competing HER pathway for all nine catalysts
- Showed that dopant charge is a strong predictor of N2 and NH3 adsorption
  strength across the bimetallic series

---

## Repository Structure

```
├── src/                        # Core Python library
│   ├── thermochemistry.py           # Adsorption energies, SOEs, reaction path calculations
│   ├── analysis.py             # Correlation analysis, linear regression, normalization
│   ├── visualization.py        # Publication-quality plotting functions
│   └── qe_postprocess.py       # Quantum ESPRESSO output parsing & restart generation
│
├── notebooks/                  # Numbered analysis notebooks (outputs included)
│   ├── 01_bader_charge_n2.ipynb
│   ├── 02_bader_charge_orientations.ipynb
│   ├── 03_n2_adsorption_descriptors.ipynb
│   ├── 04_n2_adsorption_energy.ipynb
│   ├── 05_h2_adsorption_energy.ipynb
│   ├── 06_reaction_path_analysis.ipynb
│   ├── 07_soe_analysis.ipynb
│   └── 08_final_comparison.ipynb
│
├── scripts/                    # Command-line automation tools
│   ├── generate_restart_files.py   # Post-process QE outputs → restart inputs
│   └── check_vasp_convergence.py   # Check VASP convergence & prepare reruns
│
├── data/                       # Processed CSV data files (see data/README.md)
└── requirements.txt
```

---

## Computational Setup

- **DFT code:** Quantum ESPRESSO (PW)
- **Exchange-correlation:** PBE-GGA
- **Surfaces:** 4-layer p(3×3) Cu(111) and Ni(111) slabs with single-atom dopant substitutions
- **Adsorption sites:** Top, bridge, FCC hollow, HCP hollow (2 symmetry-distinct sites reported)
- **Charge analysis:** Bader charge decomposition (VASP + Henkelman group code)

---

## Installation

```bash
git clone https://github.com/parezaei/nrr-trimetallic-dft.git
cd nrr-trimetallic-dft
pip install -r requirements.txt
```

---

## Quick Start

```python
import pandas as pd
from src import (
    compute_reaction_path_assoc_alternating,
    compute_step_barriers,
    find_top_barriers,
    plot_reaction_energy_diagram,
)

# Load DFT energy data
df = pd.read_csv("data/corrected_df.csv")

# Compute associative alternating pathway energies for adsorption site 1
path = compute_reaction_path_assoc_alternating(df, site_number=1)

# Find the top 3 rate-limiting steps for each catalyst
barriers = compute_step_barriers(path, system_cols=list(path.columns[:-1]))
top_barriers = find_top_barriers(barriers, n=3)
print(top_barriers)
```

---

## Notebooks

| Notebook | Description |
|----------|-------------|
| `01_bader_charge_n2.ipynb` | Bader charge analysis of N2 adsorption on all catalyst systems |
| `02_bader_charge_orientations.ipynb` | Bader charges for N2 vertical, horizontal, and N2H orientations |
| `03_n2_adsorption_descriptors.ipynb` | N2 adsorption descriptors: charge transfer, bond length, site preference |
| `04_n2_adsorption_energy.ipynb` | N2 adsorption energy vs. N–N bond length correlations |
| `05_h2_adsorption_energy.ipynb` | H2 adsorption energy analysis (HER competing pathway) |
| `06_reaction_path_analysis.ipynb` | Free energy diagrams for all three NRR pathways and HER; barrier analysis |
| `07_soe_analysis.ipynb` | Surface occupation energies; correlation with Bader charge and ΔX |
| `08_final_comparison.ipynb` | Final ranking of all catalysts; NRR vs HER selectivity |

---

## Scripts

### `generate_restart_files.py`

Parses Quantum ESPRESSO `.out` files, checks convergence, and generates
restart `.in` files with updated atomic coordinates.

```bash
# Single molecule directory
python scripts/generate_restart_files.py --dir ./0_CuNi/2_N --min 1 --max 2

# All molecules in a catalyst directory
python scripts/generate_restart_files.py --dir ./0_CuNi --wholesale --min 1 --max 2
```

### `check_vasp_convergence.py`

Checks VASP OUTCAR convergence, updates INCAR/POSCAR for restart, and
optionally resubmits SLURM jobs.

```bash
python scripts/check_vasp_convergence.py \
    --material_dir /scratch/ML_project/Satvik/0_Al13Fe4 \
    --mols 2_N2 3_H2 --cluster narval --rerun
```

---

## Citation

If you use this code, please cite:

```bibtex
@article{agharezaei2025trimetallic,
  title   = {Designing Trimetallic Single-Doped Alloy Catalysts for Sustainable
             Ammonia Production: The Role of Dopants in Active Site Engineering},
  author  = {Agharezaei, Parastoo and Ghuman, Kulbir K.},
  journal = {ACS Catalysis},
  year    = {2025},
  doi     = {10.1021/acscatal.5c01049},
  url     = {https://doi.org/10.1021/acscatal.5c01049}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
