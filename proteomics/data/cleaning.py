#Data cleaning: raw DataFrames -> tidy (Fraction, Retention_time, Sequence, Intensity) tables.
import re
import numpy as np
import pandas as pd
import streamlit as st

#----------------------------------------------------
# Shared helpers: modification column + RT-tolerant collapse
#----------------------------------------------------
def _clean_modifications(df):
    #Returns the modification string of every row ('' = unmodified). 
    mod_col = next((c for c in df.columns if c.strip().lower() == 'modifications'), None)
    if mod_col is None:
        return pd.Series('', index=df.index)
    return df[mod_col].fillna('').astype(str).str.strip()

def _add_peptide_column(df, separate_modified):
    # 'Peptide' is the identity used for exclusivity / overlap statistics.
    #Modified forms only count as distinct peptides when the user keeps them separate.
    df = df.copy()
    if separate_modified:
        mods = df['Modifications'].fillna('').astype(str)
        df['Peptide'] = df['Sequence'] + mods.map(lambda m: f' [{m}]' if m else '')
    else:
        df['Peptide'] = df['Sequence']
    return df


def _join_modifications(mods):
    return '; '.join(sorted({m for m in mods if m}))


def _collapse_by_rt(df, keys, rt_tolerance, agg):
    
#    Collapses rows that share all `keys` and whose retention times lie within `rt_tolerance` (min) of the first
#    retention time of their cluster. Intensity is combined with `agg` ('sum' / 'max'), the retention time of the
#    most intense row is kept. A tolerance of 0 only merges identical retention times.
#    If 'Modifications' is not one of the keys, the distinct modifications of the merged rows are joined together.
    
    df = df.sort_values(keys + ['Retention_time']).copy()

    #for each group of identical keys, _cluster_ids walks through the RTs in order. start is the RT that opened the cluster
    #If a row is more than the rt_tolerance from the start, it opens a new cluster. Otherwise it joins the current one.
    def _cluster_ids(rts):
        ids, start, current = [], None, -1
        for rt in rts.values:
            if start is None or rt - start > rt_tolerance:
                current += 1
                start = rt
            ids.append(current)
        return pd.Series(ids, index=rts.index)

    #assign each row a cluster number
    df['_cluster'] = df.groupby(keys, sort=False)['Retention_time'].transform(_cluster_ids)

    #sort by intensity so that 'last' = most intense row of the cluster
    df = df.sort_values('Intensity', kind='stable')
    #collapse each cluster to 1 row, this only happens when the table actually has a modifications column 
    #and when modification is not one of the grouping columns. this is the case when "Keep modified forms separate" is off. 
    named_aggs = {'Retention_time': ('Retention_time', 'last'), 'Intensity': ('Intensity', agg)}
    if 'Modifications' in df.columns and 'Modifications' not in keys:
        #merged row keeps track if it contains a modified sequence
        named_aggs['Modifications'] = ('Modifications', _join_modifications)
    out = df.groupby(keys + ['_cluster'], as_index=False, sort=False).agg(**named_aggs).drop(columns='_cluster')

    first = ['Fraction', 'Retention_time']
    last = ['Intensity']
    middle = [c for c in out.columns if c not in first + last]
    return out[first + middle + last].sort_values(['Fraction', 'Retention_time']).reset_index(drop=True)


#----------------------------------------------------
#Progenesis input 
#----------------------------------------------------
def clean_data(df, filter_quant=True, combine_mode='max', separate_modified=False, rt_tolerance=0.0):
    """
    Cleans and aggregates raw data from progenesis input. Same retention here is exactly same retention time unless a RT tolerance is 
    specified by the user. User can specify to treat modified sequence as separate sequence to its original even if they
    have the same fraction and retention time. 
    Returns:
        -clean dataframe
        -name of the column that was detected as the intensity column
        -statistics about faulty data
    """
    df.columns = df.columns.str.strip() 

    if '#' in df.columns:
        df = df.rename(columns={'#': 'Index'})

    # 1. Retention time 
    #find next column with 'Retention time or RT' in its name
    rt_col = next((c for c in df.columns if 'Retention time' in c or 'RT' in c), None) 
    if rt_col is None:
        st.error("CRITICAL ERROR: No 'Retention time' column found!")
        raise ValueError("No retention time column found!")
    df['Retention_time'] = pd.to_numeric(df[rt_col], errors='coerce')    

    # 2. Fraction 
    frac_col = next((c for c in df.columns if 'Fraction' in c), None)
    if frac_col:
        df['Fraction'] = pd.to_numeric(df[frac_col], errors='coerce')
    else:
        #if no fraction column is found, all data will have fraction 1
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
    #First search for an column with Intensity column with Sample is used as a backup, last fallback first numeric column that is 
    #not on the exclude list
    intensity_col = _find_col(df_clean, lambda c: c.lower() == 'intensity')
    if intensity_col is None:
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
    #count rows with invalid values
    invalid_rt_count = df_clean['Retention_time'].isna().sum()
    invalid_intensity_count = df_clean['Intensity'].isna().sum()    

    # Accession & Cleaning 
    accession_col = next((c for c in df_clean.columns if 'Accession' in c), None)
    if accession_col is None:
        st.error("CRITICAL ERROR: No 'Accession' column found!")
        raise ValueError("No 'Accession' column found.")
    
    
    df_clean['Accession'] = df_clean[accession_col].astype(str).str.strip()
    #count rows with invalid accession number    
    is_invalid = df_clean['Accession'].isna() | \
             df_clean['Accession'].isin(['', 'nan', 'NaN', 'None', 'NULL'])

    invalid_accession_count = is_invalid.sum()
    #drop rows with invalid values
    df_clean = df_clean[~is_invalid].copy()    
    df_clean = df_clean.dropna(subset=['Retention_time', 'Intensity'])

    seq_col = next((c for c in df_clean.columns if 'Sequence' in c), None)
    if seq_col is None: 
        st.error("CRITICAL ERROR: No 'Sequence' column found!")
        raise ValueError("No 'Sequence' column found.")
    
    df_clean['Sequence'] = df_clean[seq_col].astype(str).str.strip()
    df_clean['Modifications'] = _clean_modifications(df_clean)
    bad_seq = df_clean['Sequence'].isin(['', 'nan', 'NaN', 'None', 'NULL'])
    invalid_sequence_count = bad_seq.sum()
    df_clean = df_clean[~bad_seq].copy()

    invalid_stats = {
        "RT": invalid_rt_count,
        "Intensity": invalid_intensity_count,
        "Accession": invalid_accession_count,
        "Sequence": invalid_sequence_count,
        "TOTAL": invalid_rt_count + invalid_intensity_count + invalid_accession_count + invalid_sequence_count
    }

    # 6. Collapse
    rows_before = len(df_clean)
    #modified forms only stay separate if the user asks for it, otherwise they are merged into their unmodified sequence
    mod_key = ['Modifications'] if separate_modified else []
    #peptides with same fraction, retention time, accession number and sequence should be summed. 
    #Retention times within the RT tolerance count as the same, modification is part of the key only if modified forms are kept separate
    df_clean = _collapse_by_rt(df_clean, ['Fraction', 'Accession', 'Sequence'] + mod_key, rt_tolerance, 'sum')
    rows_after_accession_collapse = len(df_clean)
    
    st.write(f"**[COLLAPSE 1/2 - by Accession]** {rows_before} rows -> {rows_after_accession_collapse} rows "
             f"(merged charge-state duplicates of the same peptide/accession, RT tolerance {rt_tolerance} min).")
    
    #peptides with same fraction retention time sequence, but different accession number should be aggregated based on user-specified mode
    #same RT tolerance concept applies here
    df_clean = _collapse_by_rt(df_clean.drop(columns='Accession'), ['Fraction', 'Sequence'] + mod_key, rt_tolerance, combine_mode)
    st.write(
        f"**[COLLAPSE 2/2 - by Sequence, mode={combine_mode.upper()}]** "
        f"{rows_after_accession_collapse} rows -> {len(df_clean)} rows "
        f"(merged peptides shared across multiple protein accessions"
        f"{'; modified forms kept separate' if separate_modified else '; modified forms merged'}).")
    
    df_clean = _add_peptide_column(df_clean, separate_modified)
    return df_clean, intensity_col, invalid_stats

# ------------------------------------------------------------------
# PLGS input (Excel workbook, one sheet per fraction)
# ------------------------------------------------------------------
def strip_flanking_residues(seq):
    #Remove the bracketed aa residues at beginning and end.
    # (K)HLVDEPQNLIK(Q)' -> 'HLVDEPQNLIK',  '(K)LVVSTQTALA(-)' -> 'LVVSTQTALA'
    if pd.isna(seq):
        return ''
    s = str(seq).strip()
    s = re.sub(r'^\([^)]*\)', '', s)
    s = re.sub(r'\([^)]*\)$', '', s)
    return s.strip()


def _find_col(df, predicate):
    # Find first column whose (stripped) name satisfies predicate, else None.
    return next((c for c in df.columns if predicate(c)), None)


def clean_plgs_data(df, separate_modified=False, rt_tolerance=0.05):
    """
    Cleans and aggregates raw data from PLGS input. Less stringent than Progenesis, same retention is considered 0.05 or less 
    difference in RT. User can specify to treat modified sequence as separate sequence to its original even if they
    have the same fraction and retention time. 
    Returns:
        -clean dataframe
        -name of the column that was detected as the intensity column
        -statistics about faulty data
    """
    df = df.copy()
    df.columns = df.columns.astype(str).str.strip()

    # Case-insensitive lookups for column names. Can't use RT here since it is present in another column (Products RMS RT Error (min))
    rt_col = _find_col(df, lambda c: c.lower().startswith('retention time'))
    int_col = _find_col(df, lambda c: c.lower() == 'intensity')
    seq_col = _find_col(df, lambda c: c.lower() == 'sequence')
    for col, name in ((rt_col, 'Retention Time'), (int_col, 'Intensity'), (seq_col, 'Sequence')):
        if col is None:
            st.error(f"CRITICAL ERROR: No '{name}' column found in the PLGS data!")
            raise ValueError(f"No '{name}' column found.")

    out = pd.DataFrame({
        'Fraction': df['Fraction'],
        'Retention_time': pd.to_numeric(df[rt_col], errors='coerce'),
        'Sequence': df[seq_col].map(strip_flanking_residues),
        'Intensity': pd.to_numeric(df[int_col], errors='coerce'),
        'Modifications': _clean_modifications(df),
    })

    bad_rt = out['Retention_time'].isna()
    bad_int = out['Intensity'].isna()
    bad_seq = out['Sequence'] == ''
    invalid = bad_rt | bad_int | bad_seq
    invalid_stats = {
        "RT": int(bad_rt.sum()),
        "Intensity": int(bad_int.sum()),
        "Sequence": int(bad_seq.sum()),
        "TOTAL": int(invalid.sum()),
    }
    df_clean = out[~invalid].reset_index(drop=True)

    n_raw_seq = df[seq_col].nunique()
    st.write(f"**[PLGS - SEQUENCES]** bracketed residues removed: {n_raw_seq} distinct sequences "
             f"-> {df_clean['Sequence'].nunique()}.")
    rows_before = len(df_clean)
    #modified forms only stay separate if the user asks for it, otherwise they are merged into their unmodified sequence
    mod_key = ['Modifications'] if separate_modified else []
    df_clean = _collapse_by_rt(df_clean, ['Fraction', 'Sequence'] + mod_key, rt_tolerance, 'sum')
    st.write(f"**[PLGS - COLLAPSE]** {rows_before} rows -> {len(df_clean)} rows "
             f"(rows of the same sequence within {rt_tolerance} min are merged; "
             f"{'modified forms stay separate bars' if separate_modified else 'modified forms are merged'}).")
    
    df_clean = _add_peptide_column(df_clean, separate_modified)
    return df_clean, int_col, invalid_stats
