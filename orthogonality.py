"""
Orthogonality of fractionation methods via Shannon entropy and mutual information,
following Biba et al., J. Chromatogr. Open 8 (2025) 100262 (sections 2.6 and 2.7).

Pipeline
--------
1. peptide_fraction_profiles()  df_clean  ->  table  Sequence x fraction  (TIC-normalised intensity)
2. assign_fractions()           profile   ->  one fraction number per peptide  (paper: fitted mean, rounded)
3. compute_orthogonality()      assignments of all methods -> 1D entropies, pair metrics, 6x6 count matrices
4. plot_*()                     Fig. 6A / 6B / 7 equivalents

Everything up to step 3 is plain pandas/numpy (no Streamlit), so it can be unit-tested on its own.
"""
import itertools
import math

import numpy as np
import pandas as pd

from visualization_logic import MPL_LOCK, _locked, _new_figure, fmt_fraction

ASSIGN_MEAN = "Intensity-weighted mean (paper)"
ASSIGN_APEX = "Apex fraction (most intense)"
WEIGHT_LINEAR = "linear"
WEIGHT_LOG2 = "log2"


# ------------------------------------------------------------------
# 1. Peptide x fraction intensity profiles
# ------------------------------------------------------------------
def peptide_fraction_profiles(df_clean, n_fractions=None, aggregate="sum", normalize_tic=True):
    """
    One row per peptide Sequence, one column per fraction (1..n_fractions), 0 where the
    peptide was not seen.

    df_clean       : output of clean_data / clean_plgs_data (Fraction, Retention_time, Sequence, Intensity).
                     A sequence can have several rows in one fraction (charge states, modified forms,
                     several retention times); they are merged here with `aggregate` ('sum' or 'max').
    normalize_tic  : divide every fraction by its total intensity (sum over ALL rows of that fraction).
                     This is only a proxy for the real TIC, because the exported tables only contain
                     identified peptides.
    """
    d = df_clean.dropna(subset=["Fraction", "Sequence", "Intensity"]).copy()
    d["Sequence"] = d["Sequence"].astype(str).str.strip().str.upper()
    d["Fraction"] = d["Fraction"].round().astype(int)

    if n_fractions is None:
        n_fractions = int(d["Fraction"].max())

    profile = (d.groupby(["Sequence", "Fraction"])["Intensity"].agg(aggregate)
                 .unstack("Fraction", fill_value=0.0))
    profile = profile.reindex(columns=range(1, n_fractions + 1), fill_value=0.0)

    if normalize_tic:
        tic = d.groupby("Fraction")["Intensity"].sum().reindex(profile.columns)
        profile = profile.div(tic.replace(0, np.nan), axis=1).fillna(0.0)
    return profile


# ------------------------------------------------------------------
# 2. Fraction assignment per peptide
# ------------------------------------------------------------------
def assign_fractions(profile, method=ASSIGN_MEAN, weight=WEIGHT_LINEAR):
    """
    Returns a Series  Sequence -> fraction number (int).

    ASSIGN_MEAN : the intensity profile is treated as a distribution over the fraction index; its mean
                  mu is rounded to the nearest whole fraction. A tie (mu = x.5) goes to the neighbouring
                  fraction with the higher intensity (paper, 2.6.2).
    ASSIGN_APEX : the fraction holding the highest intensity (robust for peptides that elute in two
                  far-apart fractions, where a mean would land in between).

    weight (only for the mean): 'linear' uses the normalised intensities as weights, 'log2' uses
    log2(1 + intensity in ppm of the fraction total), a flattened version that follows the paper's
    "log2-transformed" wording.
    """
    fractions = np.asarray(profile.columns, dtype=int)
    values = profile.to_numpy(dtype=float)

    if method == ASSIGN_APEX:
        idx = values.argmax(axis=1)
        return pd.Series(fractions[idx], index=profile.index, name="Fraction")

    w = np.log2(1.0 + values * 1e6) if weight == WEIGHT_LOG2 else values
    totals = w.sum(axis=1)
    mu = np.divide((w * fractions).sum(axis=1), totals, out=np.full(len(w), np.nan), where=totals > 0)

    assigned = np.zeros(len(mu), dtype=int)
    for i, m in enumerate(mu):
        if np.isnan(m):
            assigned[i] = 0          # peptide without any signal; dropped below
            continue
        lo = math.floor(m)
        if math.isclose(m - lo, 0.5, abs_tol=1e-9):
            hi = lo + 1
            lo_i = np.where(fractions == lo)[0]
            hi_i = np.where(fractions == hi)[0]
            i_lo = values[i, lo_i[0]] if len(lo_i) else -1
            i_hi = values[i, hi_i[0]] if len(hi_i) else -1
            assigned[i] = hi if i_hi > i_lo else lo
        else:
            assigned[i] = int(math.floor(m + 0.5))   # not Python round(): that rounds .5 to even
    out = pd.Series(assigned, index=profile.index, name="Fraction")
    out = out[out > 0]
    return out.clip(lower=int(fractions.min()), upper=int(fractions.max()))


# ------------------------------------------------------------------
# 3. Entropy and mutual information
# ------------------------------------------------------------------
def _entropy(p):
    p = np.asarray(p, dtype=float).ravel()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def contingency_matrix(assign_a, assign_b, n_fractions):
    """n x n counts: cell (i, j) = peptides in fraction i of method A and fraction j of method B."""
    common = assign_a.index.intersection(assign_b.index)
    m = np.zeros((n_fractions, n_fractions), dtype=int)
    for a, b in zip(assign_a.loc[common], assign_b.loc[common]):
        m[a - 1, b - 1] += 1
    return m


def pair_metrics(counts):
    """Joint / marginal entropies, mutual information, and sum of conditional entropies (bits)."""
    total = counts.sum()
    if total == 0:
        return dict(n_peptides=0, H1=np.nan, H2=np.nan, joint=np.nan, MI=np.nan, cond=np.nan)
    P = counts / total
    h_joint = _entropy(P)
    h1 = _entropy(P.sum(axis=1))
    h2 = _entropy(P.sum(axis=0))
    mi = max(0.0, h1 + h2 - h_joint)          # Eq. 4; clip tiny negative float error
    return dict(n_peptides=int(total), H1=h1, H2=h2, joint=h_joint, MI=mi,
                cond=h_joint - mi)            # = H(M1|M2) + H(M2|M1), the blue part of Fig. 6A


def compute_orthogonality(assignments, n_fractions, scope="all"):
    """
    assignments : {label: Series Sequence -> fraction}
    scope       : 'all'  -> only peptides present in EVERY method are used (paper, 2.6.1)
                  'pair' -> each pair uses the peptides the two methods share (keeps more peptides)

    Returns (singles_df, pairs_df, matrices)
      singles_df : label, n_peptides, entropy
      pairs_df   : pair, A, B, n_peptides, joint, MI, cond, H_A, H_B  (sorted by joint entropy, high to low)
      matrices   : {(A, B): count matrix}
    """
    labels = list(assignments)
    if scope == "all" and len(labels) > 1:
        universe = set.intersection(*[set(a.index) for a in assignments.values()])
    else:
        universe = None

    def _restrict(a):
        return a if universe is None else a[a.index.isin(universe)]

    singles = []
    for label in labels:
        a = _restrict(assignments[label])
        counts = np.bincount(a.to_numpy(), minlength=n_fractions + 1)[1:]
        singles.append({"label": label, "n_peptides": int(counts.sum()),
                        "entropy": _entropy(counts / counts.sum()) if counts.sum() else np.nan})

    pairs, matrices = [], {}
    for A, B in itertools.combinations(labels, 2):
        counts = contingency_matrix(_restrict(assignments[A]), _restrict(assignments[B]), n_fractions)
        met = pair_metrics(counts)
        matrices[(A, B)] = counts
        pairs.append({"pair": f"{A} - {B}", "A": A, "B": B, "n_peptides": met["n_peptides"],
                      "joint": met["joint"], "MI": met["MI"], "cond": met["cond"],
                      "H_A": met["H1"], "H_B": met["H2"]})

    pairs_df = pd.DataFrame(pairs)
    if not pairs_df.empty:
        pairs_df = pairs_df.sort_values("joint", ascending=False).reset_index(drop=True)
    return pd.DataFrame(singles), pairs_df, matrices


def assignment_table(assignments):
    """Sequence x method table of assigned fractions (blank = not detected in that method)."""
    return pd.DataFrame(assignments).sort_index().astype("Int64")


# ------------------------------------------------------------------
# 4. Plots (Fig. 6A, 6B and 7 of the paper)
# ------------------------------------------------------------------
@_locked
def plot_pair_entropy(pairs_df, singles_df, n_fractions):
    """Left: stacked bars per method pair (blue = conditional entropies, purple = mutual information,
    total height = joint entropy). Right: 1D entropy per method."""
    if pairs_df.empty and singles_df.empty:
        return None
    fig = _new_figure(figsize=(max(9, 0.75 * len(pairs_df) + 5), 5.5), layout="constrained")
    ax_a, ax_b = fig.subplots(1, 2, gridspec_kw={"width_ratios": [max(2, len(pairs_df)), max(1, len(singles_df))]})

    if not pairs_df.empty:
        x = np.arange(len(pairs_df))
        ax_a.bar(x, pairs_df["cond"], color="#4db3ff", edgecolor="black", linewidth=0.6,
                 label="H(M1|M2) + H(M2|M1)")
        ax_a.bar(x, pairs_df["MI"], bottom=pairs_df["cond"], color="#b44cf0", edgecolor="black",
                 linewidth=0.6, label="Mutual information I(M1;M2)")
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(pairs_df["pair"], rotation=90, fontsize=9)
        ax_a.set_ylabel("Entropy [bits]")
        ax_a.set_title("Joint entropy per method pair", fontweight="bold")
        ax_a.legend(loc="lower left", fontsize=8)
        ax_a.grid(True, axis="y", linestyle="--", alpha=0.3)
    else:
        ax_a.axis("off")

    x = np.arange(len(singles_df))
    ax_b.bar(x, singles_df["entropy"], color="#b44cf0", edgecolor="black", linewidth=0.6)
    ax_b.axhline(math.log2(n_fractions), color="gray", linestyle=":", linewidth=1)
    ax_b.text(len(singles_df) - 0.5, math.log2(n_fractions), f" max {math.log2(n_fractions):.2f}",
              va="bottom", ha="right", fontsize=8, color="gray")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(singles_df["label"], rotation=90, fontsize=9)
    ax_b.set_title("1D entropy per method", fontweight="bold")
    ax_b.grid(True, axis="y", linestyle="--", alpha=0.3)
    return fig


@_locked
def plot_contingency_grid(matrices, n_fractions, max_cols=3):
    """One heat map per method pair: peptide count per (fraction in A, fraction in B)."""
    if not matrices:
        return None
    n = len(matrices)
    ncols = min(n, max_cols)
    nrows = math.ceil(n / ncols)
    fig = _new_figure(figsize=(4.2 * ncols, 4.2 * nrows), layout="constrained")
    axes = fig.subplots(nrows, ncols, squeeze=False).ravel()
    vmax = max(1, max(int(m.max()) for m in matrices.values()))

    im = None
    for ax, ((A, B), m) in zip(axes, matrices.items()):
        # rows = method A (y axis), columns = method B (x axis); fraction 1 at the bottom like Fig. 7
        im = ax.imshow(m, origin="lower", cmap="cool", vmin=0, vmax=vmax)
        for i in range(n_fractions):
            for j in range(n_fractions):
                ax.text(j, i, str(m[i, j]), ha="center", va="center", fontsize=9)
        ticks = range(n_fractions)
        ax.set_xticks(ticks)
        ax.set_xticklabels([fmt_fraction(t + 1) for t in ticks])
        ax.set_yticks(ticks)
        ax.set_yticklabels([fmt_fraction(t + 1) for t in ticks])
        ax.set_xlabel(f"{B} fraction")
        ax.set_ylabel(f"{A} fraction")
    for ax in axes[n:]:
        ax.set_visible(False)
    fig.colorbar(im, ax=list(axes[:n]), shrink=0.8, label="Number of peptides")
    fig.suptitle("Peptide distribution over fraction pairs", fontsize=14, fontweight="bold")
    return fig
