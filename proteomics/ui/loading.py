"""Loading and cleaning of the uploaded files into one dataset per fractionation method."""
import pandas as pd
import streamlit as st

import proteomics.analysis.method_detection as md
from proteomics.analysis.peptide_stats import compute_aa_properties, compute_fraction_peptide_stats
from proteomics.data.cleaning import clean_data, clean_plgs_data
from proteomics.data.loaders import load_plgs_data, load_proteomics_data


# ------------------------------------------------------------------
# 4. Load + clean files
# ------------------------------------------------------------------
def stop_if_no_files(settings):
    #if no files are uploaded or if all files are set to skip, give an info message and stop the rest of the code from running
    if not settings.active:
        kind = "PLGS Excel (.xlsx)" if settings.is_plgs else "proteomics CSV"
        st.info(f"Upload one or more {kind} files in the sidebar. "
                f"The file name must contain the fractionation method ({', '.join(md.ALLOWED_METHODS)}).")
        st.stop()

def load_datasets(settings):
    """Read + clean every active file. Returns {label: dataset dict}; files that fail are skipped."""
    datasets = {}
    #loop over every file in format ((filename, method), unique_label)
    for (f, method), label in zip(settings.active, settings.labels):
        #create collapsable container
        with st.expander(f"Processing log - {label}  ({f.name})", expanded=False):
            try:
                if settings.is_plgs:
                    df = load_plgs_data(f)
                    df_clean, intensity_col, invalid_stats = clean_plgs_data(
                        df, separate_modified=settings.separate_modified, rt_tolerance=settings.rt_tolerance)
                else:
                    df = load_proteomics_data(f)
                    df_clean, intensity_col, invalid_stats = clean_data(
                        df, filter_quant=settings.filter_quant, combine_mode=settings.combine_mode,
                        separate_modified=settings.separate_modified, rt_tolerance=settings.rt_tolerance)
            except Exception as e:
                st.error(f"'{f.name}' could not be processed and is skipped: {e}")
                continue

            if invalid_stats["TOTAL"] > 0:
                detail = ", ".join(f"{n} {k}" for k, n in invalid_stats.items() if k != "TOTAL")
                st.warning(f"**{invalid_stats['TOTAL']} invalid rows removed:** ({detail})")
            else:
                st.success("Data clean! No invalid rows detected.")

            if df_clean.empty:
                st.error(f"'{f.name}' has no usable rows after cleaning and is skipped.")
                continue

            datasets[label] = {
                'file': f.name,
                'df_clean': df_clean,
                'df_props': compute_aa_properties(df_clean),
                'fractions': sorted(df_clean['Fraction'].dropna().unique()),
                'x_min': df_clean['Retention_time'].min() - 1.0,
                'x_max': df_clean['Retention_time'].max() + 1.0,
                'fraction_stats': compute_fraction_peptide_stats(df_clean),
            }
    return datasets


def show_loaded_overview(datasets):
    # Overview of what was loaded
    st.subheader("Loaded files")
    st.dataframe(pd.DataFrame([{
        'Method': label,
        'File': ds['file'],
        'Fractions': len(ds['fractions']),
        'Unique peptides': ds['df_clean']['Peptide'].nunique(),
        'Rows after cleaning': len(ds['df_clean']),
    } for label, ds in datasets.items()]), hide_index=True, width='stretch')
