"""Per-fraction peptide statistics, biophysical properties and fraction overlap sets."""
import numpy as np
import pandas as pd
import streamlit as st

from proteomics.analysis.aa_properties import compute_peptide_properties


def compute_fraction_peptide_stats(df_clean):
    #counts in many the fractions each sequence occurs
    seq_fraction_counts = df_clean.groupby('Peptide')['Fraction'].nunique()
    #calculates what sequence are exclusive (only in 1 fraction)
    exclusive_peptides = set(seq_fraction_counts[seq_fraction_counts == 1].index)

    stats = {}
    for fraction, group in df_clean.groupby('Fraction'):
        sequences = set(group['Peptide'])
        #pct_unique is the percentage of exclusive sequences in the fraction
        n_peptides = len(sequences)
        n_unique = len(sequences & exclusive_peptides)
        pct_unique = (n_unique / n_peptides * 100) if n_peptides > 0 else 0.0

        #pct_unique_intensity is percentage of cumulative intensity for the exclusive sequences
        cumulative_intensity = group['Intensity'].sum()
        unique_mask = group['Peptide'].isin(exclusive_peptides)
        unique_intensity = group.loc[unique_mask, 'Intensity'].sum()
        pct_unique_intensity = (unique_intensity / cumulative_intensity * 100) if cumulative_intensity > 0 else 0.0

        stats[fraction] = {
            'n_peptides': n_peptides,
            'n_unique': n_unique,
            'pct_unique': pct_unique,
            'cumulative_intensity': cumulative_intensity,
            'unique_intensity': unique_intensity,
            'pct_unique_intensity': pct_unique_intensity,
        }
    return stats 


def compute_aa_properties(df_clean):
    #Computes theoretical pI, GRAVY, and molecular weight for every unique
    #peptide Sequence in df_clean, then adds them to the dataframe

    unique_seqs = df_clean['Sequence'].dropna().unique()

    records = []
    n_failed = 0
    for seq in unique_seqs:
        try:
            mw, pi, gravy = compute_peptide_properties(seq)
        except Exception:
            mw, pi, gravy = np.nan, np.nan, np.nan
            n_failed += 1
        records.append({'Sequence': seq, 'Length': len(str(seq.strip())), 'MW': mw, 'pI': pi, 'GRAVY': gravy})

    props_df = pd.DataFrame(records)
    df_props = df_clean.merge(props_df, on='Sequence', how='left')

    if n_failed > 0:
        st.warning(
            f"Could not compute biophysical properties for {n_failed} sequence(s) "
            f"(likely non-standard amino acid codes)."
        )

    return df_props


def compute_peptide_fraction_sets(df_clean):
    #Map each fraction to the set of peptide Sequences identified in it. Used for UpSet plot
    return {fraction: set(group['Peptide']) for fraction, group in df_clean.groupby('Fraction')}


def compute_membership_combos(peptide_sets):
    #For each peptide, determine the exact combination of fractions it appears in,
    #then count how many peptides share each combination (used in UpSet plot)
    #Returns a pandas Series indexed by tuple-of-fractions, sorted by count descending.
    membership = {}
    for fraction, seqs in peptide_sets.items():
        for seq in seqs:
            membership.setdefault(seq, set()).add(fraction)
 
    combo_counts = {}
    for combo in membership.values():
        key = tuple(sorted(combo))
        combo_counts[key] = combo_counts.get(key, 0) + 1
 
    # NOTE: combo tuples have varying lengths (e.g. (1,) vs (1, 2, 3)). Passing a
    # dict with mixed-length tuple keys straight to pd.Series silently builds a
    # MultiIndex and pads shorter tuples with NaN, corrupting the combos. Building
    # the Index explicitly with tupleize_cols=False keeps each tuple as one scalar key.
    keys = list(combo_counts.keys())
    values = list(combo_counts.values())
    idx = pd.Index(keys, tupleize_cols=False)
    series = pd.Series(values, index=idx).sort_values(ascending=False)
    return series
