"""Boxplot grids (one subplot per file, one box per fraction) for peptide properties and intensities."""
import math
import streamlit as st

from proteomics.plotting.base import fmt_fraction, locked, new_figure

@locked
def _boxplot_grid(panels, ylabel, suptitle, color, log_scale=False, zero_line=False, max_cols=3):
    """
    One big figure with one subplot per file. Each subplot is titled with the
    file's fractionation method and shows one box per fraction.
    panels: {method_label: (fractions_list, [array_per_fraction])}
    All subplots share the y-axis so the methods can be compared directly.
    """
    n = len(panels)
    ncols = min(n, max_cols)
    nrows = math.ceil(n / ncols)
    max_frac = max(len(fr) for fr, _ in panels.values())
    panel_w = max(4.0, max_frac * 0.8)

    fig = new_figure(figsize=(panel_w * ncols, 4.3 * nrows), layout='constrained')
    axes = fig.subplots(nrows, ncols, sharey=True, squeeze=False)
    flat = axes.ravel()

    for ax, (label, (fractions, data)) in zip(flat, panels.items()):
        ax.set_title(label, fontsize=12, fontweight='bold')
        ax.set_xlabel('Fraction', fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.3, axis='y')

        if not any(len(d) for d in data):
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes, color='gray')
            continue

        bp = ax.boxplot(data, patch_artist=True)
        # set ticks manually
        ax.set_xticks(range(1, len(fractions) + 1))
        ax.set_xticklabels([fmt_fraction(f) for f in fractions])
        for patch in bp['boxes']:
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        for median in bp['medians']:
            median.set_color('black')
        if zero_line:
            ax.axhline(0, color='gray', linestyle=':', linewidth=1)

    for ax in flat[n:]:
        ax.set_visible(False)
    if log_scale:
        for ax in flat[:n]:
            ax.set_yscale('log')
    for r in range(nrows):
        axes[r, 0].set_ylabel(ylabel, fontsize=10)

    fig.suptitle(suptitle, fontsize=14, fontweight='bold')
    return fig


def plot_property_boxplot_grid(props_by_method, property_col, ylabel, title, color='#6baed6'):
    """
    Plot boxplots for PI, GRAVY and Length
    props_by_method: {method_label: df_props}
    """
    panels = {}
    for label, df_props in props_by_method.items():
        if property_col not in df_props.columns:
            st.warning(f"{label}: column '{property_col}' not found - skipped in boxplot.")
            continue
        fractions = sorted(df_props['Fraction'].dropna().unique())
        dedup = df_props.drop_duplicates(subset=['Fraction', 'Sequence'])
        data = [dedup.loc[dedup['Fraction'] == f, property_col].dropna().values for f in fractions]
        panels[label] = (fractions, data)

    if not panels:
        st.warning(f"No data available for {property_col} boxplot.")
        return None
    return _boxplot_grid(panels, ylabel, title, color, zero_line=(property_col == 'GRAVY'))


def plot_intensity_boxplot_grid(dfs_by_method, ylabel='Intensity (a.u.)',
                                title='Peptide Intensity Distribution by Fraction',
                                color='#74c476', log_scale=True):
    """
    Plot Intensity boxplot
    dfs_by_method: {method_label: df_clean}
    """
    panels = {}
    for label, df_clean in dfs_by_method.items():
        if 'Intensity' not in df_clean.columns:
            st.warning(f"{label}: column 'Intensity' not found - skipped in boxplot.")
            continue
        fractions = sorted(df_clean['Fraction'].dropna().unique())
        data = [df_clean.loc[df_clean['Fraction'] == f, 'Intensity'].dropna().values for f in fractions]
        panels[label] = (fractions, data)

    if not panels:
        st.warning("No data available for the intensity boxplot.")
        return None
    return _boxplot_grid(panels, ylabel, title, color, log_scale=log_scale)
