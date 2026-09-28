import streamlit as st

import visualization_logic as vl
import method_detection as md
import pandas as pd
import io
import re
import zipfile

st.set_page_config(page_title="Proteomics Visualizer", layout="wide")

st.title("Proteomics Data Visualization Dashboard")


def slug(label):
    """'RP high pH' -> 'RP_high_pH' (safe for file names)."""
    return re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_')


# ------------------------------------------------------------------
# 1. Sidebar: files + fractionation methods
# ------------------------------------------------------------------
st.sidebar.header("Data & Filter Settings")
uploaded_files = st.sidebar.file_uploader(
    "Upload Proteomics CSV file(s)", type=['csv'], accept_multiple_files=True
)

active = []  # list of (uploaded_file, method)
if uploaded_files:
    auto_methods, notes = md.auto_assign_methods([f.name for f in uploaded_files])
    options = md.METHOD_OPTIONS + [md.SKIP]

    with st.sidebar.expander("Fractionation methods", expanded=True):
        st.caption("Detected from the file name. Change it if it is wrong, "
                   "or choose 'Skip file' to leave a file out.")
        for note in notes:
            st.warning(note)
        for i, (f, detected) in enumerate(zip(uploaded_files, auto_methods)):
            default = detected if detected in md.METHOD_OPTIONS else md.SKIP
            # the detected default is part of the key so the dropdown resets when detection changes
            choice = st.selectbox(f.name, options, index=options.index(default),
                                  key=f"method_{i}_{f.name}_{default}")
            if choice != md.SKIP:
                active.append((f, choice))

labels = md.make_unique_labels([m for _, m in active])

# ------------------------------------------------------------------
# 2. Sidebar: processing settings
# ------------------------------------------------------------------
smooth_mode = st.sidebar.selectbox("Smoothing Mode", ['none', 'savgol', 'gaussian', 'decimate'])
bin_width = st.sidebar.slider("Bin Width (min)", 0.01, 0.5, 0.1, 0.01)
window = st.sidebar.slider("Smoothing Window (min)", 0.1, 5.0, 1.0, 0.1)
bar_width = st.sidebar.slider("Bar Width", 0.01, 0.2, 0.05, 0.01)
show_filter = st.sidebar.checkbox("Show Filter Overlay", value=False)
combine_mode = st.sidebar.radio("Aggregation Method", ['max', 'sum'])
filter_quant = st.sidebar.checkbox("Enforce 'Use in quantitation' Filter", value=True)

# ------------------------------------------------------------------
# 3. Sidebar: which files to show per plot type (all OFF by default)
# ------------------------------------------------------------------
st.sidebar.header("Per-file Plots")
st.sidebar.caption("Pick the files (fractionation methods) to display. Nothing is shown by default.")
show_chromo = st.sidebar.multiselect("Chromatograms", labels, default=[])
show_stacked = st.sidebar.multiselect("Stacked Fractions", labels, default=[])
show_waterfall = st.sidebar.multiselect("3D Waterfall", labels, default=[])
show_overlap = st.sidebar.multiselect("Peptide Overlap (UpSet)", labels, default=[])

st.sidebar.header("3D Waterfall View")
elev = st.sidebar.slider("Elevation", 0, 90, 23)
azim = st.sidebar.slider("Azimuth", -180, 180, -83)

# ------------------------------------------------------------------
# 4. Load + clean every file (needed for the combined boxplots)
# ------------------------------------------------------------------
if not active:
    st.info("Upload one or more proteomics CSV files in the sidebar. "
            f"The file name must contain the fractionation method ({', '.join(md.ALLOWED_METHODS)}).")
    st.stop()

datasets = {}
for (f, method), label in zip(active, labels):
    with st.expander(f"Processing log - {label}  ({f.name})", expanded=False):
        try:
            df = vl.load_proteomics_data(f)
            df_clean, intensity_col, invalid_stats = vl.clean_data(
                df, filter_quant=filter_quant, combine_mode=combine_mode)
        except Exception as e:
            st.error(f"'{f.name}' could not be processed and is skipped: {e}")
            continue

        if invalid_stats["TOTAL"] > 0:
            st.warning(
                f"**{invalid_stats['TOTAL']} invalid rows removed:** "
                f"({invalid_stats['RT']} RT, {invalid_stats['Intensity']} Intensity, "
                f"{invalid_stats['Accession']} Accession)")
        else:
            st.success("Data clean! No invalid rows detected.")

        if df_clean.empty:
            st.error(f"'{f.name}' has no usable rows after cleaning and is skipped.")
            continue

        datasets[label] = {
            'file': f.name,
            'df_clean': df_clean,
            'df_props': vl.compute_aa_properties(df_clean),
            'fractions': sorted(df_clean['Fraction'].dropna().unique()),
            'x_min': df_clean['Retention_time'].min() - 1.0,
            'x_max': df_clean['Retention_time'].max() + 1.0,
            'fraction_stats': vl.compute_fraction_peptide_stats(df_clean),
        }

if not datasets:
    st.stop()

# Overview of what was loaded
st.subheader("Loaded files")
st.dataframe(pd.DataFrame([{
    'Method': label,
    'File': ds['file'],
    'Fractions': len(ds['fractions']),
    'Unique peptides': ds['df_clean']['Sequence'].nunique(),
    'Rows after cleaning': len(ds['df_clean']),
} for label, ds in datasets.items()]), hide_index=True, width='stretch')


# Traces are only built for files that need them (chromatogram / stacked / waterfall)
_trace_cache = {}
def get_traces(label):
    if label not in _trace_cache:
        ds = datasets[label]
        all_exact, all_binned = {}, {}
        with st.expander(f"Trace & smoothing log - {label}", expanded=False):
            for fraction in ds['fractions']:
                exact = vl.build_exact_trace(ds['df_clean'], fraction, combine_mode)
                binned_raw = vl.build_binned_trace(ds['df_clean'], fraction, bin_width,
                                                   ds['x_min'], ds['x_max'], combine_mode)
                all_exact[fraction] = exact
                all_binned[fraction] = vl.apply_smoothing_pipeline(binned_raw, smooth_mode, window, bin_width)
        _trace_cache[label] = (all_exact, all_binned)
    return _trace_cache[label]


download_pngs = {}   # file name (no extension) -> PNG bytes


def show_fig(fig, name):
    """Render the figure once (thread-safe) and use the same PNG for display and download."""
    png = vl.fig_to_png(fig)
    st.image(png, width='stretch')
    download_pngs[name] = png

# ------------------------------------------------------------------
# 5. Chromatograms (one section per selected file)
# ------------------------------------------------------------------
chromo_labels = [l for l in datasets if l in show_chromo]
if chromo_labels:
    st.header("Chromatograms")
    for label in chromo_labels:
        ds = datasets[label]
        all_exact, all_binned = get_traces(label)
        st.subheader(label)
        for fraction in ds['fractions']:
            stats = ds['fraction_stats'][fraction]
            fig = vl.plot_chromatogram_with_ma(
                all_exact[fraction], all_binned[fraction], fraction,
                stats['n_peptides'], stats['n_unique'], stats['cumulative_intensity'],
                bar_width, smooth_mode, show_filter, ds['x_min'], ds['x_max'],
                unique_intensity=stats['unique_intensity'], method=label)
            if fig is not None:
                show_fig(fig, f"{slug(label)}_FRACTION_{vl.fmt_fraction(fraction)}")

# ------------------------------------------------------------------
# 6. Stacked / 3D waterfall / overlap (one row per selected file)
# ------------------------------------------------------------------
combined_labels = [l for l in datasets
                   if l in show_stacked or l in show_waterfall or l in show_overlap]
if combined_labels:
    st.header("Combined Visualizations")
    for label in combined_labels:
        ds = datasets[label]
        kinds = [k for k, sel in (("stacked", show_stacked),
                                  ("waterfall", show_waterfall),
                                  ("overlap", show_overlap)) if label in sel]
        st.subheader(label)
        cols = st.columns(len(kinds))
        for col, kind in zip(cols, kinds):
            with col:
                if kind == "stacked":
                    all_exact, all_binned = get_traces(label)
                    fig = vl.plot_stacked_fractions(all_exact, all_binned, show_filter,
                                                    ds['x_min'], ds['x_max'], bar_width, method=label)
                    title, key = "Stacked Fractions", "STACKED"
                elif kind == "waterfall":
                    _, all_binned = get_traces(label)
                    fig = vl.plot_waterfall_3d(all_binned, smooth_mode, elev=elev, azim=azim, method=label)
                    title, key = "3D Waterfall", "WATERFALL"
                else:
                    peptide_sets = vl.compute_peptide_fraction_sets(ds['df_clean'])
                    fig = vl.plot_peptide_overlap(peptide_sets, method=label)
                    title, key = "Peptide Overlap (UpSet)", "PEPTIDE_OVERLAP"
                st.markdown(f"**{title}**")
                if fig is not None:
                    show_fig(fig, f"{slug(label)}_{key}")

# ------------------------------------------------------------------
# 7. Boxplots: always shown, one figure per property, one subplot per file
# ------------------------------------------------------------------
st.header("Peptide Biophysical Properties")
props_by_method = {label: ds['df_props'] for label, ds in datasets.items()}

fig_pi = vl.plot_property_boxplot_grid(
    props_by_method, 'pI', ylabel='Isoelectric Point (pI)',
    title='Theoretical pI Distribution by Fraction', color='#6baed6')
if fig_pi is not None:
    show_fig(fig_pi, "PI_BOXPLOT")

fig_gravy = vl.plot_property_boxplot_grid(
    props_by_method, 'GRAVY', ylabel='GRAVY Index',
    title='GRAVY Hydrophobicity Distribution by Fraction', color='#fd8d3c')
if fig_gravy is not None:
    show_fig(fig_gravy, "GRAVY_BOXPLOT")

st.header("Peptide Intensity by Fraction")
log_scale_intensity = st.checkbox("Log scale (Intensity)", value=True)
fig_intensity = vl.plot_intensity_boxplot_grid(
    {label: ds['df_clean'] for label, ds in datasets.items()}, log_scale=log_scale_intensity)
if fig_intensity is not None:
    show_fig(fig_intensity, "INTENSITY_BOXPLOT")

# ------------------------------------------------------------------
# 8. Download Center (only the plots that are currently displayed)
# ------------------------------------------------------------------
st.divider()
st.subheader("Download Center")

if download_pngs:
    cols = st.columns(4)
    for i, (name, png) in enumerate(download_pngs.items()):
        cols[i % 4].download_button(
            label=f"Download {name}",
            data=png,
            file_name=f"{name.lower()}.png",
            mime="image/png"
        )

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zf:
        for name, png in download_pngs.items():
            zf.writestr(f"{name.lower()}.png", png)

    st.download_button(
        label="📦 Download ALL Displayed Plots (ZIP)",
        data=zip_buffer.getvalue(),
        file_name="all_plots.zip",
        mime="application/zip",
        width='stretch'
    )
else:
    st.caption("No plots to download yet.")