import matplotlib 
matplotlib.use('Agg')  
import matplotlib.pyplot as plt
import streamlit as st

import visualization_logic as vl
import pandas as pd
import io
import zipfile

st.set_page_config(page_title="Proteomics Visualizer", layout="wide")

st.title("Proteomics Data Visualization Dashboard")

# 1. Sidebar: User Settings
st.sidebar.header("Data & Filter Settings")
uploaded_file = st.sidebar.file_uploader("Upload Proteomics CSV", type=['csv'])

smooth_mode = st.sidebar.selectbox("Smoothing Mode", ['none', 'savgol', 'gaussian', 'decimate'])
bin_width = st.sidebar.slider("Bin Width (min)", 0.01, 0.5, 0.1, 0.01)
window = st.sidebar.slider("Smoothing Window (min)", 0.1, 5.0, 1.0, 0.1)
bar_width = st.sidebar.slider("Bar Width", 0.01, 0.2, 0.05, 0.01)
show_filter = st.sidebar.checkbox("Show Filter Overlay", value=False)
combine_mode = st.sidebar.radio("Aggregation Method", ['max', 'sum'])
filter_quant = st.sidebar.checkbox("Enforce 'Use in quantitation' Filter", value=True)

st.sidebar.header("3D Waterfall View")
elev = st.sidebar.slider("Elevation", 0, 90, 23)
azim = st.sidebar.slider("Azimuth", -180, 180, -83)

# 2. Execution Logic
if uploaded_file:
    df = vl.load_proteomics_data(uploaded_file)
    df_clean, intensity_col, invalid_stats = vl.clean_data(df, filter_quant=filter_quant)
    if invalid_stats["TOTAL"] > 0:
        st.warning(
            f"**{invalid_stats['TOTAL']} invalid rows removed:** "
            f"({invalid_stats['RT']} RT, {invalid_stats['Intensity']} Intensity, "
            f"{invalid_stats['Accession']} Accession)")
    else: 
        st.success("Data clean! No invalid rows detected.")
    
    fractions = sorted(df_clean['Fraction'].dropna().unique())
    x_min = df_clean['Retention_time'].min() - 1.0
    x_max = df_clean['Retention_time'].max() + 1.0

    all_exact = {}
    all_binned = {}
    fraction_figs = {}

    for fraction in fractions:
        exact = vl.build_exact_trace(df_clean, fraction, combine_mode)
        binned_raw = vl.build_binned_trace(df_clean, fraction, bin_width, x_min, x_max, combine_mode)
        binned_proc = vl.apply_smoothing_pipeline(binned_raw, smooth_mode, window, bin_width)
        
        all_exact[fraction] = exact
        all_binned[fraction] = binned_proc
        
        n_peptides = df_clean[df_clean['Fraction'] == fraction]['Accession'].nunique()
        fig = vl.plot_chromatogram_with_ma(exact, binned_proc, fraction, n_peptides, 
                                           bar_width, smooth_mode, show_filter, x_min, x_max)
        fraction_figs[fraction] = fig
        st.pyplot(fig)

    # Combined Plots
    st.header("Combined Visualizations")
    col1, col2 = st.columns(2)
    
    fig_stacked = vl.plot_stacked_fractions(all_exact, all_binned, show_filter, x_min, x_max, bar_width)
    fig_3d = vl.plot_waterfall_3d(all_binned, smooth_mode, elev=elev, azim=azim)
    
    with col1:
        st.subheader("Stacked Fractions")
        st.pyplot(fig_stacked)
    with col2:
        st.subheader("3D Waterfall")
        st.pyplot(fig_3d)
    
    # 3. Download Center
    st.divider()
    st.subheader("Download Center")

    # 1. Prepare all figures
    download_figs = {
        "STACKED": fig_stacked,
        "WATERFALL": fig_3d
    }
    for f in fractions:
        download_figs[f"FRACTION_{f}"] = fraction_figs[f]

    # 2. Display Individual Download Buttons in a Grid
    cols = st.columns(4)
    for i, (name, fig) in enumerate(download_figs.items()):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight')
        cols[i % 4].download_button(
            label=f"Download {name}",
            data=buf.getvalue(),
            file_name=f"{name.lower()}.png",
            mime="image/png"
        )

    # 3. "Download All" as a ZIP file
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zf:
        for name, fig in download_figs.items():
            img_buf = io.BytesIO()
            fig.savefig(img_buf, format="png", bbox_inches='tight')
            zf.writestr(f"{name.lower()}.png", img_buf.getvalue())
    
    for fig in download_figs.values():
        plt.close(fig)

    st.download_button(
        label="📦 Download ALL Plots (ZIP)",
        data=zip_buffer.getvalue(),
        file_name="all_plots.zip",
        mime="application/zip",
        use_container_width=True
    )