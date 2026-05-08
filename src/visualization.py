"""
Visualization utilities for NRR on bimetallic surfaces.

Provides publication-quality plotting functions for:
- Free energy reaction diagrams
- Adsorption energy scatter plots
- Correlation heatmaps
- Radar (spider) charts
- SOE vs electronic descriptor scatter plots with linear fits
"""

from math import pi

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import pandas as pd
import seaborn as sns

from .analysis import linear_regression

SYSTEM_COLORS: dict[str, str] = {
    "CuNi":  "#1f77b4",
    "Ag@Cu": "#ff7f0e",
    "Au@Cu": "#2ca02c",
    "Pd@Ni": "#d62728",
    "Pt@Ni": "#9467bd",
    "Fe@Ni": "#8c564b",
    "Pd@Cu": "#e377c2",
    "Ag@Ni": "#7f7f7f",
    "Au@Ni": "#bcbd22",
    "Pt@Cu": "palevioletred",
    "Fe@Cu": "indianred",
    "Co@Ni": "green",
    "Co@Cu": "sandybrown",
    "Ru@Ni": "blue",
    "Ru@Cu": "orange",
}

MOLECULE_LABELS: dict[str, str] = {
    "0_N2":    "N₂",
    "2_N":     "N",
    "3_N2H":   "N₂H",
    "4_NHNH":  "NHNH",
    "5_NNH2":  "NNH₂",
    "6_NNH3":  "NNH₃",
    "7_NH2NH": "NH₂NH",
    "8_N2H4":  "N₂H₄",
    "9_NH2":   "NH₂",
    "10_NH":   "NH",
    "11_NH3":  "NH₃",
    "12_NH2NH3": "NH₂NH₃",
    "13_H":    "H",
}


def plot_reaction_energy_diagram(
    path_df: pd.DataFrame,
    test_path: pd.DataFrame,
    xticks: list[str],
    systems: list[str],
    legends: list[str],
    colors: list[str] | None = None,
    title: str = "",
    figsize: tuple = (35, 20),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot a free energy reaction path diagram with step plateaus and connectors.

    Parameters
    ----------
    path_df : pd.DataFrame
        Energy levels DataFrame (output of compute_reaction_path_*).
        Rows = energy levels, columns include system labels.
    test_path : pd.DataFrame
        Interpolated path for dashed connectors between plateaus.
    xticks : list[str]
        Labels for each reaction step along the x-axis.
    systems : list[str]
        System column labels in path_df.
    legends : list[str]
        Display names for the legend (same order as systems).
    colors : list[str], optional
        One colour per system. Defaults to SYSTEM_COLORS.
    title : str
        Plot title.
    figsize : tuple
        Figure size in inches.
    save_path : str, optional
        If provided, saves the figure to this path (PDF recommended).

    Returns
    -------
    matplotlib Figure.
    """
    if colors is None:
        colors = [SYSTEM_COLORS.get(leg, "#333333") for leg in legends]

    plt.rcParams["font.family"] = "Arial"
    fig, ax = plt.subplots(figsize=figsize)

    step_width = 4
    gap = 1
    for i, system in enumerate(systems):
        for step_idx, level in enumerate(path_df.index):
            xmin = step_idx * (step_width + gap)
            xmax = xmin + step_width
            ax.hlines(
                y=path_df[system][level],
                xmin=xmin,
                xmax=xmax,
                color=colors[i],
                lw=15,
                alpha=1,
            )
        ax.plot(
            test_path.index,
            test_path[system],
            ls="--",
            lw=5,
            color=colors[i],
            alpha=1,
        )

    plt.xticks(test_path.index, xticks, rotation=30)
    ax.tick_params(axis="both", which="both", length=0, labelsize=60, pad=20)
    leg = plt.legend(legends, prop={"size": 42}, loc="lower left")
    for line in leg.get_lines():
        line.set_linestyle("-")
    for handle in leg.legend_handles:
        handle.set_linewidth(12)
    plt.ylabel("Energy (eV)", fontsize=80, weight="bold", labelpad=20)
    plt.xlabel("Reaction Coordinate", fontsize=80, weight="bold", labelpad=20)
    if title:
        plt.title(title, fontsize=40, y=1.05)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_soe_scatter_vs_feature(
    feature_vals: list[float],
    soe_vals: list[float],
    systems: list[str],
    x_range: list[float] | None = None,
    xlabel: str = "Descriptor",
    ylabel: str = "SOE (eV)",
    title: str = "",
    colors: list[str] | None = None,
    figsize: tuple = (12, 8),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Scatter plot of a SOE value against an electronic descriptor with a
    linear regression fit line.

    Parameters
    ----------
    feature_vals : list[float]
        x-axis descriptor values (e.g., Bader charges).
    soe_vals : list[float]
        y-axis SOE values.
    systems : list[str]
        System labels corresponding to each data point.
    x_range : list[float], optional
        [x_min, x_max] for drawing the fit line. Auto-detected if None.
    xlabel, ylabel : str
        Axis labels.
    title : str
        Plot title.
    colors : list[str], optional
        Colours for each data point. Uses SYSTEM_COLORS if None.
    figsize : tuple
        Figure size.
    save_path : str, optional
        Save path for the figure.

    Returns
    -------
    matplotlib Figure.
    """
    if colors is None:
        colors = [SYSTEM_COLORS.get(s, "#333333") for s in systems]

    reg = linear_regression(feature_vals, soe_vals)
    if x_range is None:
        x_range = [min(feature_vals) - 0.05, max(feature_vals) + 0.05]
    x_line = np.linspace(x_range[0], x_range[1], 100)
    y_line = reg["slope"] * x_line + reg["intercept"]

    pylab.rcParams["xtick.major.pad"] = "10"
    pylab.rcParams["ytick.major.pad"] = "10"

    fig, ax = plt.subplots(figsize=figsize)
    for x, y, color, label in zip(feature_vals, soe_vals, colors, systems):
        ax.scatter(x, y, s=500, color=color, edgecolors="black", zorder=3, label=label)
    ax.plot(x_line, y_line, color="black", lw=2, linestyle="--", label=reg["equation"])
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.legend(fontsize=12, bbox_to_anchor=(1.05, 1), loc="upper left")
    if title:
        ax.set_title(title, fontsize=16)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_correlation_heatmap(
    df: pd.DataFrame,
    title: str = "Pearson Correlation Matrix",
    figsize: tuple = (20, 18),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot a lower-triangular Pearson correlation heatmap.

    Parameters
    ----------
    df : pd.DataFrame
        Numeric DataFrame whose correlations will be plotted.
    title : str
        Plot title.
    figsize : tuple
        Figure size.
    save_path : str, optional
        Save path.

    Returns
    -------
    matplotlib Figure.
    """
    corr = df.select_dtypes(include="number").corr()
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True

    sns.set(font_scale=2.0)
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        corr,
        center=0,
        cmap="coolwarm",
        annot=True,
        linewidths=0.3,
        linecolor="white",
        fmt=".2f",
        annot_kws={"size": 14},
        mask=mask,
        square=True,
        ax=ax,
    )
    ax.set_title(title, fontsize=20, pad=15)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_spider_chart(
    df: pd.DataFrame,
    group_col: str = "group",
    col_labels: dict | None = None,
    figsize: tuple = (14, 14),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Plot radar (spider/polar) charts — one subplot per catalyst system.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with a `group_col` column and numeric feature columns.
    group_col : str
        Name of the column that identifies the catalyst system.
    col_labels : dict, optional
        Mapping from column name → display label. Uses column names if None.
    figsize : tuple
        Figure size.
    save_path : str, optional
        Save path.

    Returns
    -------
    matplotlib Figure.
    """
    categories = [c for c in df.columns if c != group_col]
    display_labels = (
        [col_labels.get(c, c) for c in categories] if col_labels else categories
    )
    n = len(categories)
    angles = [i / float(n) * 2 * pi for i in range(n)]
    angles += angles[:1]

    n_rows = int(np.ceil(len(df) / 3))
    fig = plt.figure(figsize=figsize)
    palette = plt.cm.get_cmap("Set2", len(df))

    for idx, (_, row) in enumerate(df.iterrows()):
        ax = fig.add_subplot(n_rows, 3, idx + 1, polar=True)
        ax.set_theta_offset(pi / 2)
        ax.set_theta_direction(-1)
        plt.xticks(angles[:-1], display_labels, color="grey", size=9)
        ax.set_rlabel_position(0)

        values = row.drop(group_col).values.flatten().tolist()
        values += values[:1]
        color = palette(idx)
        ax.plot(angles, values, color=color, linewidth=2, linestyle="solid")
        ax.fill(angles, values, color=color, alpha=0.4)
        plt.title(row[group_col], size=11, color=color, y=1.1)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_adsorption_scatter_by_pathway(
    df: pd.DataFrame,
    intermediates: list[str],
    systems: list[str],
    system_col: str = "system",
    ylabel: str = "Adsorption Energy (eV)",
    figsize: tuple = (14, 8),
    save_path: str | None = None,
) -> plt.Figure:
    """
    Scatter plot comparing adsorption energies across a pathway for all systems.

    Each point represents one (system, intermediate) combination.

    Parameters
    ----------
    df : pd.DataFrame
        Adsorption energy DataFrame.
    intermediates : list[str]
        Column names for the NxHy intermediates to plot.
    systems : list[str]
        System display names (matched against `system_col`).
    system_col : str
        Column identifying each catalyst system.
    ylabel : str
        y-axis label.
    figsize : tuple
        Figure size.
    save_path : str, optional
        Save path.

    Returns
    -------
    matplotlib Figure.
    """
    x_labels = [MOLECULE_LABELS.get(m, m) for m in intermediates]
    pylab.rcParams["xtick.major.pad"] = "10"
    pylab.rcParams["ytick.major.pad"] = "10"

    fig, ax = plt.subplots(figsize=figsize)
    for i, system in enumerate(systems):
        row = df[df[system_col] == system]
        if row.empty:
            continue
        y_vals = row[intermediates].values.flatten()
        ax.scatter(
            x_labels,
            y_vals,
            s=300,
            color=SYSTEM_COLORS.get(system, "#333333"),
            edgecolors="black",
            label=system,
            zorder=3,
        )
    ax.set_ylabel(ylabel, fontsize=16)
    ax.tick_params(axis="x", labelsize=14, rotation=30)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=12)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig
