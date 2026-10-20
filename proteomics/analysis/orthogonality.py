"""
Orthogonality of fractionation methods via joint entropy and mutual information.

Pipeline
--------
1. keep only the peptides present in ALL methods
2. per method: TIC-normalise, log2-transform, fit a normal distribution to each peptide's
   profile over the fractions, round the fitted mean -> one fraction number per peptide
3. per pair of methods: count matrix -> joint probability matrix -> entropies -> mutual information
"""
import itertools

import numpy as np
import pandas as pd

WINDOW = 1             # centroid is computed over apex +/- WINDOW fractions
SHARP_THRESHOLD = 0.8  # "sharp" = at least 80% of the signal lies inside that window

def shared_peptides(dfs_clean):
    #calculate set of peptides present in every method. 
    return set.intersection(*[set(df['Peptide']) for df in dfs_clean.values()])


def assign_fractions(df_clean, n_fractions, window=WINDOW, sharp_threshold=SHARP_THRESHOLD):
    """Returns a DataFrame indexed by Peptide with columns Fraction, Confidence, Sharp."""
    d = df_clean.dropna(subset=['Fraction', 'Peptide', 'Intensity'])
    fraction = d['Fraction'].round().astype(int)

    # one row per peptide, one column per fraction, summed intensity
    profile = (d.groupby(['Peptide', fraction])['Intensity'].sum()
                 .unstack('Fraction', fill_value=0.0)
                 .reindex(columns=range(1, n_fractions + 1), fill_value=0.0))

    # TIC normalisation -> divide by total intensity of each fraction
    tic = profile.sum(axis=0).replace(0, np.nan)
    profile = profile.div(tic, axis=1).fillna(0.0)

    # remove empty peptides, peptides without any signal cannot be placed
    profile = profile[profile.sum(axis=1) > 0]
    w = profile.to_numpy()
    fractions = np.arange(1, n_fractions + 1)

    # find the apex = fraction with the highest intensity
    apex = fractions[np.argmax(w, axis=1)]

    # keep only apex +/- window, everything else gets weight 0
    in_window = np.abs(fractions[None, :] - apex[:, None]) <= window
    w_win = np.where(in_window, w, 0.0)

    # calculate weighted average of the intensities inside the centroid window
    mu = (w_win * fractions).sum(axis=1) / w_win.sum(axis=1)

    # confidence = share of the peptide's total signal inside the window
    confidence = w_win.sum(axis=1) / w.sum(axis=1)

    assigned = np.clip(np.floor(mu + 0.5), 1, n_fractions).astype(int)
    return pd.DataFrame({'Fraction': assigned, 'Confidence': confidence,
                         'Sharp': confidence >= sharp_threshold}, index=profile.index)


def entropy(p):
    #convert input to flat array of floats
    p = np.asarray(p, dtype=float).ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def joint_probability(assign_a, assign_b, n_fractions):
    """n x n matrix: cell (i, j) = share of peptides in fraction i of method A and fraction j of method B."""
    counts = np.zeros((n_fractions, n_fractions))
    np.add.at(counts, (assign_a.to_numpy() - 1, assign_b.loc[assign_a.index].to_numpy() - 1), 1)
    return counts / counts.sum()


def compute_orthogonality(dfs_clean, only_sharp=False):
    """Returns (pairs_df, n_used, n_fractions, summary), or None when no peptides can be used."""
    shared = shared_peptides(dfs_clean)
    if not shared:
        return None

    dfs_clean = {label: df[df['Peptide'].isin(shared)] for label, df in dfs_clean.items()}
    n_fractions = int(max(df['Fraction'].round().max() for df in dfs_clean.values()))
    assignments = {label: assign_fractions(df, n_fractions) for label, df in dfs_clean.items()}

    # peptides placed in every method, and of those the ones that are sharp in every method
    usable = set.intersection(*[set(a.index) for a in assignments.values()])
    sharp_all = {pep for pep in usable if all(a.at[pep, 'Sharp'] for a in assignments.values())}
    used = sharp_all if only_sharp else usable
    if not used:
        return None
    used_idx = sorted(used)
    fractions = {label: a.loc[used_idx, 'Fraction'] for label, a in assignments.items()}

    rows = []
    for a, b in itertools.combinations(fractions, 2):
        P = joint_probability(fractions[a], fractions[b], n_fractions)
        h_a, h_b, h_joint = entropy(P.sum(axis=1)), entropy(P.sum(axis=0)), entropy(P)
        rows.append({'Pair': f"{a} - {b}", 'H(M1)': h_a, 'H(M2)': h_b, 'Joint entropy': h_joint,
                     'Mutual information': max(0.0, h_a + h_b - h_joint)})
    marginals = {}
    for label, fr in fractions.items():
        counts = np.bincount(fr.to_numpy() - 1, minlength=n_fractions)
        marginals[label] = entropy(counts / counts.sum())
    summary = {'n_shared': len(shared), 'n_sharp': len(sharp_all), 'n_used': len(used), 'marginals': marginals, 'fractions': fractions}
    return pd.DataFrame(rows), len(used), n_fractions, summary