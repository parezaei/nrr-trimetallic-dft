"""
nrr-bimetallic-dft: Analysis toolkit for nitrogen reduction reaction (NRR)
on bimetallic Cu/Ni-based surfaces from DFT calculations.
"""

from .energetics import (
    MOLECULE_ENERGIES_EV,
    SYSTEMS_LABEL_MAP,
    RY_TO_EV,
    parse_qe_final_energy,
    compute_adsorption_energy,
    compute_soe,
    compute_reaction_path_assoc_alternating,
    compute_reaction_path_assoc_distal,
    compute_reaction_path_dissociative,
    compute_her_path,
    compute_step_barriers,
    find_top_barriers,
)
from .analysis import (
    linear_regression,
    compute_pearson_correlation,
    normalize_minmax,
    normalize_zscore,
    build_feature_dataframe,
)
from .visualization import (
    SYSTEM_COLORS,
    MOLECULE_LABELS,
    plot_reaction_energy_diagram,
    plot_soe_scatter_vs_feature,
    plot_correlation_heatmap,
    plot_spider_chart,
    plot_adsorption_scatter_by_pathway,
)
from .qe_postprocess import (
    check_convergence,
    generate_restart_input,
    process_site_range,
    process_molecules_wholesale,
)

__all__ = [
    "MOLECULE_ENERGIES_EV", "SYSTEMS_LABEL_MAP", "RY_TO_EV",
    "parse_qe_final_energy", "compute_adsorption_energy", "compute_soe",
    "compute_reaction_path_assoc_alternating", "compute_reaction_path_assoc_distal",
    "compute_reaction_path_dissociative", "compute_her_path",
    "compute_step_barriers", "find_top_barriers",
    "linear_regression", "compute_pearson_correlation",
    "normalize_minmax", "normalize_zscore", "build_feature_dataframe",
    "SYSTEM_COLORS", "MOLECULE_LABELS",
    "plot_reaction_energy_diagram", "plot_soe_scatter_vs_feature",
    "plot_correlation_heatmap", "plot_spider_chart",
    "plot_adsorption_scatter_by_pathway",
    "check_convergence", "generate_restart_input",
    "process_site_range", "process_molecules_wholesale",
]
