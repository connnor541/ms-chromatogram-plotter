"""Chromatogram plots: single fraction (exact bars + smoothed profile) and stacked fractions."""
import streamlit as st
from matplotlib.ticker import ScalarFormatter

from proteomics.plotting.base import fmt_fraction, locked, new_figure, with_method


@locked
def plot_chromatogram_with_ma(exact_df, binned_df, fraction, n_peptides, n_unique, cumulative_intensity, bar_width,
                               smooth_mode, show_filter, x_min, x_max,
                               overlay_mode='twin_axis', unique_intensity=None, method=None):
    if len(exact_df) == 0:
        st.warning("No data to plot")
        return None

    fig = new_figure(figsize=(10, 5))  
    ax = fig.subplots()
    label_text = f"Profile (filter = {smooth_mode})"
    #percentage of unique peptides
    pct_unique = (n_unique / n_peptides * 100) if n_peptides > 0 else 0.0
     # Legend will look as follows -> number of bars drawn on the plot, how many peptides there are, how many exclusive peptides there are, 
    # cumulative intensity, exclusive peptide intensity 
    pep_label = (f"Peaks: {len(exact_df)}\n"
                 f"Unique peptides: {n_peptides} \n"
                 f"Exclusive peptides: {n_unique} ({pct_unique:.0f}%)\n"
                 f"Cumulative intensity: {cumulative_intensity:.2e}"
            )    
                  
    if unique_intensity is not None:
        pct_unique_intensity = (unique_intensity / cumulative_intensity * 100) if cumulative_intensity > 0 else 0.0
        pep_label += f"\nExclusive peptide intensity: {unique_intensity:.2e} ({pct_unique_intensity:.0f}%)"
    
    #if you want intensity values displayed as a % of max change the overlay mode to 'pct_of_max', twin axis mode will be used by default
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
        ax.set_ylim(0, exact_df['Intensity'].max() * 1.5)
        ax.set_ylabel('Intensity (a. u.)', fontsize=12)

        if show_filter:
            ax2 = ax.twinx()
            ax2.plot(binned_df['bin_center'], binned_df['Processed'], color='red', alpha=0.55, linewidth=1.8,
                     label=label_text)
            ax2.set_ylim(0, binned_df['Processed'].max() * 1.5)
            ax2.tick_params(axis='y', left=False, right=False, labelleft=False, labelright=False)

            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
        else:
            ax.legend(loc='upper right')

    if x_min is not None and x_max is not None:
        ax.set_xlim(x_min, x_max)

    title = f'MS Chromatogram - Fraction {fmt_fraction(fraction)}' if fraction is not None else 'MS Chromatogram'
    title = with_method(title, method)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel('Retention Time (min)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.3)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    fig.tight_layout()
    
    # Return the figure object instead of saving it
    return fig


@locked
def plot_stacked_fractions(all_exact, all_binned, show_filter, x_min, x_max, bar_width, method=None):
    n_fractions = len(all_exact)
    if n_fractions == 0:
        st.warning("No data available for stacked plotting.")
        return None

    # Use constrained_layout=True for better spacing between stacked subplots
    fig = new_figure(figsize=(10, 3 * n_fractions), layout='constrained')
    axes = fig.subplots(n_fractions, 1, sharex=True)
    
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

        ax.set_ylabel(f'Frac {fmt_fraction(fraction)}', fontsize=11, rotation=0, labelpad=30)
        ax.grid(True, linestyle='--', alpha=0.3)

        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axes[-1].set_xlabel('Retention Time (min)', fontsize=12)
    fig.suptitle(with_method('MS Chromatograms - All Fractions (Stacked)', method), fontsize=14, fontweight='bold')

    return fig
