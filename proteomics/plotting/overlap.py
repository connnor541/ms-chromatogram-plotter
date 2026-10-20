"""UpSet plot of peptide overlap across fractions."""
import numpy as np
import streamlit as st

from matplotlib.ticker import MaxNLocator
from proteomics.analysis.peptide_stats import compute_membership_combos
from proteomics.plotting.base import fmt_fraction, locked, new_figure, with_method


@locked
def plot_peptide_overlap(peptide_sets, top_n=25, method=None):
    """
    Visualize peptide overlap across fractions as an UpSet plot: a bar chart of
    exclusive intersection sizes on top, with a dot-matrix below showing which
    fractions make up each intersection. Scales cleanly to any number of fractions.
    top_n specifies how many columns are shown in the upSet plot. 
    """
    fractions = sorted(peptide_sets.keys())
    n = len(fractions)
 
    if n < 2:
        st.warning("Need at least 2 fractions to show peptide overlap.")
        return None
 
    combo_counts = compute_membership_combos(peptide_sets)
    if len(combo_counts) == 0:
        st.warning("No peptides found to plot overlap.")
        return None
 
    truncated = len(combo_counts) > top_n
    total_combos = len(combo_counts)
    combo_counts = combo_counts.iloc[:top_n]
    n_combos = len(combo_counts)
 
    fig = new_figure(figsize=(max(8, n_combos * 0.45), 7), layout='constrained')
    gs = fig.add_gridspec(2, 1, height_ratios=[2.2, max(1.0, n * 0.35)], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_dots = fig.add_subplot(gs[1], sharex=ax_bar)
 
    x = np.arange(n_combos)
    ax_bar.bar(x, combo_counts.values, color='#333333', width=0.6)
    for xi, val in zip(x, combo_counts.values):
        ax_bar.text(xi, val, f'{val}', ha='center', va='bottom', fontsize=8)
    ax_bar.set_ylabel('Peptide count', fontsize=11)
    ax_bar.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_bar.set_title(with_method("Peptide Overlap Across Fractions (UpSet plot)", method), fontsize=14, fontweight='bold')
    ax_bar.spines[['top', 'right']].set_visible(False)
    ax_bar.tick_params(axis='x', bottom=False, labelbottom=False)
 
    for yi, fraction in enumerate(fractions):
        ax_dots.axhline(yi, color='#dddddd', linewidth=0.8, zorder=0)
 
    for xi, combo in enumerate(combo_counts.index):
        y_positions = [fractions.index(f) for f in combo]
        ax_dots.plot([xi] * len(y_positions), y_positions, color='#333333',
                     linewidth=1.5, zorder=1)
        for yi in range(n):
            in_combo = yi in y_positions
            ax_dots.scatter(xi, yi,
                             color='#333333' if in_combo else '#e0e0e0',
                             s=60, zorder=2)
 
    ax_dots.set_yticks(range(n))
    ax_dots.set_yticklabels([f"Fraction {fmt_fraction(f)}" for f in fractions], fontsize=10)
    ax_dots.set_ylim(-0.5, n - 0.5)
    ax_dots.invert_yaxis()
    ax_dots.set_xlim(-0.5, n_combos - 0.5)
    ax_dots.set_xticks([])
    ax_dots.spines[['top', 'right', 'bottom']].set_visible(False)
 
    if truncated:
        fig.text(0.5, 0.005, f"Showing top {top_n} of {total_combos} combinations",
                  ha='center', fontsize=9, style='italic', color='gray')
 
    return fig
