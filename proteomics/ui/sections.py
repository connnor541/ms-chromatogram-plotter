"""Main-page plot sections: chromatograms, combined visualizations, boxplots."""
import streamlit as st

from proteomics.analysis.peptide_stats import compute_peptide_fraction_sets
from proteomics.plotting.base import fmt_fraction
from proteomics.plotting.boxplots import plot_intensity_boxplot_grid, plot_property_boxplot_grid
from proteomics.plotting.chromatograms import plot_chromatogram_with_ma, plot_stacked_fractions
from proteomics.plotting.overlap import plot_peptide_overlap
from proteomics.plotting.waterfall import default_waterfall_colors, plot_waterfall_3d
from proteomics.ui.downloads import slug


# ------------------------------------------------------------------
# 5. Chromatograms (one section per selected file)
# ------------------------------------------------------------------
def render_chromatograms(datasets, settings, traces, display):
    chromo_labels = [l for l in datasets if l in settings.show_chromo]
    if chromo_labels:
        st.header("Chromatograms")
        for label in chromo_labels:
            ds = datasets[label]
            all_exact, all_binned = traces.get(label)
            st.subheader(label)
            for fraction in ds['fractions']:
                stats = ds['fraction_stats'][fraction]
                fig = plot_chromatogram_with_ma(
                    all_exact[fraction], all_binned[fraction], fraction,
                    stats['n_peptides'], stats['n_unique'], stats['cumulative_intensity'],
                    settings.bar_width, settings.smooth_mode, settings.show_filter, ds['x_min'], ds['x_max'],
                    unique_intensity=stats['unique_intensity'], method=label)
                if fig is not None:
                    display.show_fig(fig, f"{slug(label)}_FRACTION_{fmt_fraction(fraction)}")


# ------------------------------------------------------------------
# 6. Stacked / 3D waterfall / overlap (one row per selected file)
# ------------------------------------------------------------------
PLOT_TITLES = {"stacked": ("Stacked Fractions", "STACKED"),
               "waterfall": ("3D Waterfall", "WATERFALL"),
               "overlap": ("Peptide Overlap (UpSet)", "PEPTIDE_OVERLAP")}

def _pick_waterfall_colors(label, fractions):
    """One color picker per fraction (inside a collapsed expander). Starts at the default colors,
    so the plot only changes when the user picks something. Returns {fraction: hex color}."""
    defaults = default_waterfall_colors(fractions)
    colors = {}
    with st.expander("Waterfall colors"):
        cols = st.columns(3)
        for i, (fraction, default) in enumerate(zip(fractions, defaults)):
            # the number of fractions is part of the key, so the defaults reset when a different file is loaded under the same label
            colors[fraction] = cols[i % 3].color_picker(
                f"Fraction {fmt_fraction(fraction)}", default,
                key=f"wf_color_{label}_{len(fractions)}_{fmt_fraction(fraction)}")
    return colors

def render_combined_visualizations(datasets, settings, traces, display):
    combined_labels = [l for l in datasets
                       if l in settings.show_stacked or l in settings.show_waterfall or l in settings.show_overlap]
    if combined_labels:
        st.header("Combined Visualizations")
        for label in combined_labels:
            ds = datasets[label]
            kinds = [k for k, sel in (("stacked", settings.show_stacked),
                                      ("waterfall", settings.show_waterfall),
                                      ("overlap", settings.show_overlap)) if label in sel]
            st.subheader(label)
            cols = st.columns(len(kinds))
            for col, kind in zip(cols, kinds):
                with col:
                    title, key = PLOT_TITLES[kind]
                    st.markdown(f"**{title}**")
                    if kind == "stacked":
                        all_exact, all_binned = traces.get(label)
                        fig = plot_stacked_fractions(all_exact, all_binned, settings.show_filter,
                                                     ds['x_min'], ds['x_max'], settings.bar_width, method=label)
                    elif kind == "waterfall":
                        _, all_binned = traces.get(label)
                        colors = _pick_waterfall_colors(label, ds['fractions'])
                        fig = plot_waterfall_3d(all_binned, settings.smooth_mode,
                                                elev=settings.elev, azim=settings.azim, method=label,
                                                custom_colors=colors)
                    else:
                        peptide_sets = compute_peptide_fraction_sets(ds['df_clean'])
                        fig = plot_peptide_overlap(peptide_sets, method=label)
                    if fig is not None:
                        display.show_fig(fig, f"{slug(label)}_{key}")


# ------------------------------------------------------------------
# 7. Boxplots: one figure per property, one subplot per file
# ------------------------------------------------------------------
def render_boxplots(datasets, display):
    st.header("Boxplots")
    st.caption("Tick a box to display the boxplot (one subplot per file, one box per fraction).")

    c1, c2, c3, c4 = st.columns(4)
    show_pi = c1.checkbox("pI", value=False)
    show_gravy = c2.checkbox("GRAVY", value=False)
    show_length = c3.checkbox("Peptide length", value=False)
    show_intensity = c4.checkbox("Intensity", value=False)

    props_by_method = {label: ds['df_props'] for label, ds in datasets.items()}

    if show_pi:
        fig_pi = plot_property_boxplot_grid(
            props_by_method, 'pI', ylabel='Isoelectric Point (pI)',
            title='Theoretical pI Distribution by Fraction', color='#6baed6')
        if fig_pi is not None:
            display.show_fig(fig_pi, "PI_BOXPLOT")

    if show_gravy:
        fig_gravy = plot_property_boxplot_grid(
            props_by_method, 'GRAVY', ylabel='GRAVY Index',
            title='GRAVY Hydrophobicity Distribution by Fraction', color='#fd8d3c')
        if fig_gravy is not None:
            display.show_fig(fig_gravy, "GRAVY_BOXPLOT")

    if show_length:
        fig_length = plot_property_boxplot_grid(
            props_by_method, 'Length', ylabel='Peptide length (amino acids)',
            title='Peptide Length Distribution by Fraction', color='#9e9ac8')
        if fig_length is not None:
            display.show_fig(fig_length, "LENGTH_BOXPLOT")

    if show_intensity:
        log_scale_intensity = st.checkbox("Log scale (Intensity)", value=True)
        fig_intensity = plot_intensity_boxplot_grid(
            {label: ds['df_clean'] for label, ds in datasets.items()}, log_scale=log_scale_intensity)
        if fig_intensity is not None:
            display.show_fig(fig_intensity, "INTENSITY_BOXPLOT")
