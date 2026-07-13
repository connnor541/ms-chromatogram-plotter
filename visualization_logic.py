import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.signal import savgol_filter
from mpl_toolkits.mplot3d import proj3d
import streamlit as st

def load_proteomics_data(file_input):
    try:
        df = pd.read_csv(file_input, sep=',', encoding='utf-8')
        st.info(f"Data loaded succesfully! Shape: {df.shape[0]} rows, {df.shape[1]} columns.")    

        return df
    except Exception as e:
        st.error(f"Error loading file: {e}")
        raise

def clean_data(df, filter_quant=True):
    df.columns = df.columns.str.strip() 

    if '#' in df.columns:
        df = df.rename(columns={'#': 'Index'})

    # 1. Retention time 
    rt_col = next((c for c in df.columns if 'Retention time' in c), None) #find next column with 'Retention time' in its name
    if rt_col is None:
        st.error("CRITICAL ERROR: No 'Retention time' column found!")
        raise ValueError("No retention time column found!")
    df['Retention_time'] = pd.to_numeric(df[rt_col], errors='coerce')

    # 2. Fraction 
    frac_col = next((c for c in df.columns if 'Fraction' in c), None)
    if frac_col:
        df['Fraction'] = pd.to_numeric(df[frac_col], errors='coerce')
    else:
        st.warning("No fraction column found, using all data")
        df['Fraction'] = 1

    # 3. Quantitation filter
    quant_col = next((c for c in df.columns if 'Use in quantitation' in c), None)
    if quant_col and filter_quant:
        mask = df[quant_col].astype(str).str.upper().str.strip() == 'TRUE'
        df_clean = df[mask].copy()
        st.success(f"Quantitation Filter applied ('{quant_col}'): {len(df)} -> {len(df_clean)} rows.")
    else:
        df_clean = df.copy()
        st.info("Skipping 'Use in quantitation' filtration")

    # 4. Intensity 
    intensity_col = next((c for c in df_clean.columns if 'Sample' in c), None)
    if intensity_col is None:
        numeric_cols = df_clean.select_dtypes(include=[np.number]).columns 
        exclude_cols = ['Retention_time', 'Fraction', 'Charge', 'Score', 'Mass error', 'm/z']
        numeric_cols = [c for c in numeric_cols if c not in exclude_cols]
        if not numeric_cols:
            st.error("No intensity column found!")
            raise ValueError("No intensity column found!")
        intensity_col = numeric_cols[0]

    df_clean['Intensity'] = pd.to_numeric(df_clean[intensity_col], errors='coerce')

    # Accession & Cleaning 
    accession_col = next((c for c in df_clean.columns if 'Accession' in c), None)
    if accession_col is None:
        st.error("CRITICAL ERROR: No 'Accession' column found!")
        raise ValueError("No 'Accession' column found.")
    
    df_clean['Accession'] = df_clean[accession_col].astype(str).str.strip()
    df_clean = df_clean[~((df_clean['Accession'] == '') | (df_clean['Accession'].str.upper() == 'NAN'))]
    
    df_clean = df_clean.dropna(subset=['Retention_time', 'Intensity'])

    # 6. Collapse
    rows_before = len(df_clean)
    df_clean = df_clean.groupby(['Fraction', 'Retention_time', 'Accession'], as_index=False)['Intensity'].sum()
    
    st.write(f"Collapse complete: {rows_before} rows -> {len(df_clean)} rows.")
    
    return df_clean, intensity_col


def build_exact_trace(df, fraction, combine_mode):
    # Filter by fraction if applicable
    if fraction is not None and 'Fraction' in df.columns:
        df_frac = df[df['Fraction'] == fraction].copy()
    else:
        df_frac = df.copy()

    # Handle empty data
    if len(df_frac) == 0:
        st.warning(f"Fraction {fraction}: No data found.")
        return pd.DataFrame(columns=['Retention_time', 'Intensity'])

    # AGGREGATION: Combine distinct accession signals at shared retention timestamps
    exact = df_frac.groupby('Retention_time', as_index=False)['Intensity'].agg(combine_mode)
    exact = exact.sort_values('Retention_time').reset_index(drop=True)
    
    # UI Output - Updated to show all three metrics
    n_start = len(df_frac)
    n_end = len(exact)
    n_collapsed = n_start - n_end
    
    st.write(f"**[EXACT TRACE] Fraction {fraction}** - Combine Mode: {combine_mode.upper()}")
    st.write(f"Rows: Start: {n_start} | End: {n_end} | Collapsed: {n_collapsed}")
    
    return exact

def build_binned_trace(df, fraction, bin_width_min, x_min, x_max, combine_mode):
    if fraction is not None and 'Fraction' in df.columns:
        df_frac = df[df['Fraction'] == fraction].copy()
    else:
        df_frac = df.copy()

    if len(df_frac) == 0:
        return pd.DataFrame(columns=['bin_center', 'Intensity'])

    # Create the bins
    bin_edges = np.arange(x_min, x_max + bin_width_min, bin_width_min)
    df_frac['bin'] = pd.cut(df_frac['Retention_time'], bins=bin_edges, include_lowest=True)

    # Aggregate data within corresponding bin
    binned = df_frac.groupby('bin', observed=False)['Intensity'].agg(combine_mode).fillna(0).reset_index()
    binned['bin_center'] = binned['bin'].apply(lambda b: b.mid).astype(float)
    binned = binned[['bin_center', 'Intensity']].sort_values('bin_center').reset_index(drop=True)

    # UI Output
    st.write(f"**[BINNING] Fraction {fraction}** (Width: {bin_width_min} min)")
    st.write(f"- Total Bins: {len(binned)} | Non-Zero Bins: {np.count_nonzero(binned['Intensity'])}")
    
    return binned

def apply_smoothing_pipeline(binned_df, mode, window_minutes, bin_width_min):
    if len(binned_df) == 0:
        return binned_df.assign(Processed=[])

    out = binned_df.copy()
    window_bins = max(1, round(window_minutes / bin_width_min))
    if window_bins % 2 == 0: 
        window_bins += 1

    # --------------------NoFilter-------------------
    if mode == 'none' or mode is None:
        out['Processed'] = out['Intensity']

    # --------------Savitzky-GolayFilter-------------
    elif mode == 'savgol':
        savgol_window = max(5, window_bins)
        valid_points = np.count_nonzero(out['Intensity'])
        if len(out['Intensity']) < savgol_window or valid_points < 3:
            st.warning(f"Insufficient data for SavGol window ({savgol_window}). Falling back.")
            out['Processed'] = out['Intensity']
        else:
            try:
                smoothed = savgol_filter(out['Intensity'], window_length=savgol_window, polyorder=2)
                out['Processed'] = np.clip(smoothed, 0, None)
            except ValueError as e:
                st.warning(f"Savgol filter failed: {e}. Falling back.")
                out['Processed'] = out['Intensity']

    # ----------RollingGaussianWindowFilter---------
    elif mode == 'gaussian':
        gaussian_window = max(3, window_bins)
        if len(out['Intensity']) < 2:
            out['Processed'] = out['Intensity']
        else:
            smoothed = out['Intensity'].rolling(
                window=gaussian_window, center=True, min_periods=1, win_type='gaussian'
            ).mean(std=gaussian_window / 4)
            out['Processed'] = smoothed.fillna(0)

    # ----------DecimatingMovingAverageFilter---------
    elif mode == 'decimate':
        n = 4
        if len(out) <= n:
            n = max(1, len(out) // 2) if len(out) > 2 else 1
        if n > 1:
            out['block_id'] = np.arange(len(out)) // n
            out = out.groupby('block_id', observed=False).agg({
                'bin_center': 'mean',
                'Intensity': 'mean'
            }).reset_index(drop=True)
            out['Processed'] = out['Intensity']
        else:
            out['Processed'] = out['Intensity']
    else:
        st.error(f"Unknown smoothing mode: {mode}")
        raise ValueError(f"Unknown smoothing mode: {mode}")

    # UI status report
    st.write(f"**[SMOOTHING]** Mode: {str(mode).upper()}")
    
    return out

def plot_chromatogram_with_ma(exact_df, binned_df, fraction, n_peptides, bar_width,
                               smooth_mode, show_filter, x_min, x_max,
                               overlay_mode='twin_axis'):
    if len(exact_df) == 0:
        st.warning("No data to plot")
        return None

    fig, ax = plt.subplots(figsize=(10, 5)) # Reduced slightly for better web viewing
    label_text = f"Profile (filter = {smooth_mode})"
    pep_label = f'Peptide identifications (n={n_peptides})'

    if overlay_mode == 'pct_of_max':
        bar_vals = exact_df['Intensity'] / exact_df['Intensity'].max() * 100
        ax.bar(exact_df['Retention_time'], bar_vals, width=bar_width,
               color='black', edgecolor='black', linewidth=0.3, alpha=0.9,
               label=pep_label)
        if show_filter:
            line_vals = binned_df['Processed'] / binned_df['Processed'].max() * 100 if binned_df['Processed'].max() > 0 else binned_df['Processed']
            ax.plot(binned_df['bin_center'], line_vals, color='red', alpha=0.55, linewidth=2,
                    label=f'{label_text} (% of max)')
        ax.set_ylabel('Relative Intensity (% of max)', fontsize=12)
        ax.legend(loc='upper right')
    else:
        ax.bar(exact_df['Retention_time'], exact_df['Intensity'], width=bar_width,
               color='black', edgecolor='black', linewidth=0.3, alpha=0.9,
               label=pep_label)
        ax.set_ylabel('Intensity (a. u.)', fontsize=12)

        if show_filter:
            ax2 = ax.twinx()
            ax2.plot(binned_df['bin_center'], binned_df['Processed'], color='red', alpha=0.55, linewidth=1.8,
                     label=label_text)
            ax2.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)

            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        else:
            ax.legend(loc='upper right')

    if x_min is not None and x_max is not None:
        ax.set_xlim(x_min, x_max)

    title = f'MS Chromatogram - Fraction {fraction}' if fraction is not None else 'MS Chromatogram'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Retention Time (min)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.3)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.tight_layout()
    
    # Return the figure object instead of saving it
    return fig


def compute_zlabel_rotation(ax, flip=False):
    """
    Calculate the on-screen rotation for the z-axis label. We calculate the angle of rotation
    so that the z-label is perpendicular with the axis. Normally this angle is calculated internally by Matplotlib. 
    However, Matplotlib will always put the text in top to bottom orientation, bottom to top orientation looks better graphically. 
    So we need to calculate the angle rotation ourselves and then set the rotation of the z-label to that angle and then 
    flip the text upside down. 
    """
    renderer = ax.figure.canvas.get_renderer()
    ax.draw(renderer) 

    zaxis = ax.zaxis
    mins, maxs, tc, highs = zaxis._get_coord_info()
    minmax = np.where(highs, maxs, mins)
    maxmin = np.where(~highs, maxs, mins)

    edgep1s, edgep2s, _ = zaxis._get_all_axis_line_edge_points(
        minmax, maxmin, zaxis._label_position)
    edgep1, edgep2 = edgep1s[0], edgep2s[0]

    # ax.M will exist because we called ax.draw()
    pep = np.asarray(proj3d._proj_trans_points([edgep1, edgep2], ax.M))
    
    dx, dy = (ax.transAxes.transform([pep[0:2, 1]]) -
              ax.transAxes.transform([pep[0:2, 0]]))[0]

    angle = np.degrees(np.arctan2(dy, dx))
    if flip:
        angle += 180
    return angle

def place_zaxis_multiplier_text(ax, text, fontsize=10, pad_points=10):
    """
    Places a small multiplier label (e.g. r'$\\times 10^6$') right next to the
    TOP of the z-axis, following it correctly for ANY elev/azim.
    """

    fig = ax.get_figure()
    fig.canvas.draw()  # finalize the current projection matrix (ax.M) and layout

    zaxis = ax.zaxis
    mins, maxs, tc, highs = zaxis._get_coord_info()
    minmax = np.where(highs, maxs, mins)
    maxmin = np.where(~highs, maxs, mins)

    edgep1s, edgep2s, _ = zaxis._get_all_axis_line_edge_points(
        minmax, maxmin, zaxis._label_position)
    edgep_top = np.asarray(edgep2s[0], dtype=float)  # 3D point: top of the z-axis line

    # Outward direction in DATA space: away from the box's center, in the
    # (x, y) plane only (the corner's x,y is fixed along the whole vertical
    # z-axis edge, so this is well defined regardless of azimuth/elevation).
    x0, x1 = ax.get_xlim3d()
    y0, y1 = ax.get_ylim3d()
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0

    out_x, out_y = edgep_top[0] - cx, edgep_top[1] - cy
    norm_xy = np.hypot(out_x, out_y)
    if norm_xy < 1e-9:
        out_x, out_y = 1.0, 0.0  # degenerate fallback
    else:
        out_x, out_y = out_x / norm_xy, out_y / norm_xy

    scale = 0.08 * max(x1 - x0, y1 - y0)
    outward_point_3d = np.array([edgep_top[0] + out_x * scale,
                                  edgep_top[1] + out_y * scale,
                                  edgep_top[2]])

    # Project the axis-top point and the outward point into 2D display space
    proj = np.asarray(proj3d._proj_trans_points([edgep_top, outward_point_3d], ax.M))
    p_top_disp = ax.transData.transform(proj[0:2, 0])
    p_out_disp = ax.transData.transform(proj[0:2, 1])

    direction = p_out_disp - p_top_disp
    norm = np.hypot(*direction)
    direction = direction / norm if norm > 1e-6 else np.array([0.0, 1.0])

    pad_px = pad_points * fig.dpi / 72.0  # points -> pixels, AT THE CURRENT dpi
    text_disp = p_top_disp + direction * pad_px

    # Convert to DPI-independent figure-fraction coordinates before placing,
    # so a later savefig(dpi=...) at a different resolution can't shift it.
    text_fig_frac = fig.transFigure.inverted().transform(text_disp)

    fig.text(text_fig_frac[0], text_fig_frac[1], text, fontsize=fontsize,
              ha='center', va='center', transform=fig.transFigure)
            
def plot_waterfall_3d(all_binned, smooth_mode, elev=23, azim=-83):
    fractions = sorted(list(all_binned.keys()))
    if len(fractions) == 0:
        st.warning("No data available for the waterfall plot.")
        return None

    # Find max intensity across all fractions for scaling
    all_intensities = pd.concat([df['Processed'] for df in all_binned.values()])
    max_val = all_intensities.max()
    
    # Calculate scale factor: nearest power of 1000 (e.g., 10^3, 10^6, 10^9)
    # If max_val is 0, default to 1 to avoid division errors
    exponent = int(np.log10(max_val) // 3 * 3) if max_val > 0 else 0
    scale_factor = 10**exponent

    colors = plt.cm.brg(np.linspace(0, 1, len(fractions)))
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    for i, fraction in enumerate(fractions):
        binned = all_binned[fraction]
        x = binned['bin_center'].values
        y_val = i * 2.5
        depth = np.full_like(x, y_val, dtype=float)
        
        # 2. Apply Dynamic Scale
        z = binned['Processed'].values / scale_factor 
        
        ax.plot(x, depth, z, color=colors[i], linewidth=1.8, alpha=0.8, label=str(fraction))

    ax.zaxis._axinfo['juggled'] = (1, 2, 0)
    ax.set_xlabel('Retention time (min)', fontsize=12, labelpad=10)
    ax.set_zlabel('Intensity (a.u.)', fontsize=12, labelpad=6)
    ax.zaxis.set_rotate_label(False)
    ax.set_yticks([])
    ax.set_ylabel('')

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 1))
        axis.pane.set_edgecolor('lightgray')
        axis._axinfo['grid']['color'] = (0.85, 0.85, 0.85, 1)

    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((1.5, 2.0, 0.8))

    # Apply the rotation function
    zlabel_angle = compute_zlabel_rotation(ax, flip=False)
    ax.zaxis.get_label().set_rotation(zlabel_angle)    
    
    ax.legend(loc='upper right', frameon=True, fontsize=10)
    plt.tight_layout()

    # 3. Dynamic Scaling Label (e.g., x 10^6)
    place_zaxis_multiplier_text(ax, rf'$\times 10^{{{exponent}}}$', fontsize=10, pad_points=48)   

    return fig

def plot_stacked_fractions(all_exact, all_binned, show_filter, x_min, x_max, bar_width):
    n_fractions = len(all_exact)
    if n_fractions == 0:
        st.warning("No data available for stacked plotting.")
        return None

    # Use constrained_layout=True for better spacing between stacked subplots
    fig, axes = plt.subplots(n_fractions, 1, figsize=(10, 3 * n_fractions), 
                             sharex=True, constrained_layout=True)
    
    if n_fractions == 1:
        axes = [axes]

    for idx, fraction in enumerate(all_exact.keys()):
        exact_df = all_exact[fraction]
        binned_df = all_binned[fraction]
        ax = axes[idx]
        
        ax.bar(exact_df['Retention_time'], exact_df['Intensity'], width=bar_width,
               color='black', edgecolor='black', linewidth=0.3, alpha=0.8)

        if show_filter:
            ax2 = ax.twinx()
            ax2.plot(binned_df['bin_center'], binned_df['Processed'], color='red', 
                     alpha=0.55, linewidth=1.8)
            ax2.yaxis.set_visible(False)

        if x_min is not None and x_max is not None:
            ax.set_xlim(x_min, x_max)

        ax.set_ylabel(f'Frac {fraction}', fontsize=11, rotation=0, labelpad=30)
        ax.grid(True, linestyle='--', alpha=0.3)

        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axes[-1].set_xlabel('Retention Time (min)', fontsize=12)
    fig.suptitle('MS Chromatograms - All Fractions (Stacked)', fontsize=14, fontweight='bold')

    return fig
