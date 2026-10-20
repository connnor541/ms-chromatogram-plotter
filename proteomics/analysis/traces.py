"""Chromatographic traces: exact trace, binned trace and smoothing."""
import numpy as np
import pandas as pd
import streamlit as st
from scipy.signal import savgol_filter


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
    #calculate how many bins you need
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
