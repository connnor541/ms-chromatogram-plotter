"""Orthogonality section: results table and joint entropy / mutual information plot."""
import streamlit as st

from proteomics.analysis.orthogonality import compute_orthogonality
from proteomics.plotting.orthogonality_plots import plot_joint_entropy_mi, plot_orthogonality_heatmaps


def render_orthogonality(datasets, display):
    st.header("Fractionation Orthogonality (Joint Entropy & Mutual Information)")

    if len(datasets) < 2:
        st.info("Load at least two files to calculate orthogonality.")
        return

    only_sharp = st.checkbox(
        "Only use peptides with a sharp elution peak in every file",
        value=False,
        help="A peptide is 'sharp' when at least 80% of its signal lies within one fraction of its "
             "apex fraction. Diffuse peptides are excluded when this is ticked.")

    result = compute_orthogonality({label: ds['df_clean'] for label, ds in datasets.items()},
                                   only_sharp=only_sharp)
    if result is None:
        st.warning("No peptides are available for the calculation (the loaded files share no peptides, "
                   "or none of them has a sharp elution peak in every file).")
        return
    pairs_df, n_used, n_fractions, summary = result

    st.write(f"**{n_used} peptides** are used for the calculation ({n_fractions} fractions). "
             f"{summary['n_shared']} peptides are present in all {len(datasets)} files, "
             f"of which {summary['n_sharp']} are sharp in every file.")
    st.dataframe(pairs_df.round(3), hide_index=True, width='stretch')

    fig = plot_joint_entropy_mi(pairs_df, summary['marginals'])
    if fig is not None:
        display.show_fig(fig, "ENTROPY_MI")
    
    st.subheader("Peptide distribution per method pair")
    fig_hm = plot_orthogonality_heatmaps(summary['fractions'], n_fractions)
    if fig_hm is not None:
        display.show_fig(fig_hm, "ORTHOGONALITY_HEATMAPS")