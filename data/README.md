# Data

Processed CSV files derived from Quantum ESPRESSO and VASP DFT calculations.
Raw `.in`/`.out`/`OUTCAR` files are not included due to size.

## Directory Structure

```
data/
├── reaction_paths/     # Energetics for NRR and HER free energy diagrams
├── soe/                # Spillover Energies relative to CuNi
└── bader_charge/
    ├── N2/             # Bader charge analysis for N2 adsorption
    ├── H2/             # Bader charge analysis for H2 adsorption
    └── orientations/   # N2 vertical, horizontal, and N2H orientations
```

## reaction_paths/

| File | Description |
|------|-------------|
| `corrected_df.csv` | DFT total energies for all NxHy adsorbates on each catalyst system and adsorption site |
| `df_adsE.csv` | Computed adsorption energies for all intermediates |
| `raw_df.csv` | Raw uncorrected DFT energies |
| `HER_path.csv` | Free energy levels for the HER pathway |
| `path1_s1.csv`, `path1_s2.csv` | Associative alternating pathway energy levels (sites 1 & 2) |
| `path2_s1.csv`, `path2_s2.csv` | Associative distal pathway energy levels (sites 1 & 2) |
| `path3_s1.csv`, `path3_s2.csv` | Dissociative pathway energy levels (sites 1 & 2) |
| `multip_path_*_site*.csv` | Interpolated path data for plotting plateau connectors |

## soe/

| File | Description |
|------|-------------|
| `df_dopant_SOEs_cunisite{1-4}.csv` | SOEs relative to pristine CuNi at each of the 4 adsorption sites |
| `df_dopant_adsE.csv` | Raw adsorption energies per dopant system |
| `SOEs_dopant_host.csv` | SOE broken down by dopant and host atom contributions |
| `raw_df_dopant_sites.csv` | Raw DFT energies for dopant-site calculations |
| `raw_df_host_sites.csv` | Raw DFT energies for host-site calculations |
| `corrected_df_moststable.csv` | Energies for the most stable adsorption site per system |

## bader_charge/N2/

| File | Description |
|------|-------------|
| `N2_final_comparison.csv` | Summary: N2 adsorption energy and bond length per system |
| `comparison_energy_blength.csv` | Correlation between E_ads and N–N bond length |
| `dopant_charge_diff.csv` | Bader charge on dopant atom before/after N2 adsorption |
| `medians_N2.csv`, `min_N2.csv` | Statistical summaries across adsorption sites |
| `merged_df.csv` | Combined adsorption energy + Bader charge data |

## bader_charge/H2/

| File | Description |
|------|-------------|
| `H2_final_comparison.csv` | Summary: H2 adsorption energy and bond length per system |
| `median_H2.csv`, `min_H2.csv` | Statistical summaries across adsorption sites |

## bader_charge/orientations/

| File | Description |
|------|-------------|
| `N2_vert_charge.csv` | Bader charges for vertically adsorbed N2 (all systems, all neighbours) |
| `N2_horz_charge.csv` | Bader charges for horizontally adsorbed N2 |
| `N2H_charge.csv` | Bader charges for adsorbed N2H intermediate |

## Units

- All energies: **eV**
- Charges: **electrons (e)**
- Bond lengths: **Å**
