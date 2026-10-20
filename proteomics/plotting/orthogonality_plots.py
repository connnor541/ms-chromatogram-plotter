"""Orthogonality plot: joint entropy + mutual information per method pair, marginal entropy per method."""
import itertools
import math
import numpy as np

from matplotlib.patches import Patch
from proteomics.plotting.base import locked, new_figure

BLUE, PURPLE = "#4db3ff", "#b44cf0"
EDGE = dict(edgecolor="black", linewidth=0.5)

@locked
def plot_joint_entropy_mi(pairs_df, marginals):
    """Panel A: per method pair, mutual information stacked on top of the joint entropy (bits).
    Panel B: marginal entropy of every method (bits). `marginals` is {method label: entropy}."""
    if pairs_df.empty:
        return None
    df = pairs_df.sort_values("Joint entropy", ascending=False).reset_index(drop=True)
    je, mi = df["Joint entropy"].to_numpy(), df["Mutual information"].to_numpy()
    xa = np.arange(len(df))
    labels_b = list(marginals)
    xb = np.arange(len(labels_b))
    width = 0.8    
    
    fig = new_figure(figsize=(max(8, 0.55 * len(df) + 0.55 * len(labels_b) + 3), 5.5),
                     layout="constrained")
    ax_a, ax_b = fig.subplots(1, 2, width_ratios=[len(df) + 1, len(labels_b) + 1])

    # --- A: stacked bars ---
    ax_a.bar(xa, je, width, color=BLUE, label="Joint entropy", **EDGE)
    ax_a.bar(xa, mi, width, bottom=je, color=PURPLE, label="Mutual information", **EDGE)
    ax_a.set_xticks(xa)
    ax_a.set_xticklabels([p.replace(" - ", "-") for p in df["Pair"]], rotation=90, ha="center", fontsize=9)
    ax_a.set_ylim(0, (je + mi).max() * 1.25)
    ax_a.legend(fontsize=9, loc="upper right")

    # --- B: marginal entropies ---
    ax_b.bar(xb, list(marginals.values()), width, color=PURPLE, label="Marginal entropy", **EDGE)
    ax_b.set_xticks(xb)
    ax_b.set_xticklabels([f"{m}-{m}" for m in labels_b], rotation=90, ha="center", fontsize=9)
    ax_b.set_ylim(0, max(marginals.values()) * 1.25)
    ax_b.legend(fontsize=9, loc="upper right")

    for ax, letter in ((ax_a, "A"), (ax_b, "B")):
        ax.set_xlabel("Method Pairs")
        ax.set_ylabel("Entropy [Bits]")
        ax.text(0, 1.02, letter, transform=ax.transAxes, fontsize=14, fontweight="bold")
    return fig

@locked
def plot_orthogonality_heatmaps(fractions, n_fractions, ncols=5):
    """One heat map per method pair: number of peptides in each (fraction of method 1, fraction of
    method 2) cell. `fractions` is {method label: Series of assigned fraction per peptide}."""
    pairs = list(itertools.combinations(fractions, 2))
    if not pairs:
        return None
    ncols = min(ncols, len(pairs))
    nrows = math.ceil(len(pairs) / ncols)

    # count matrices first, so every map shares one colour scale
    mats = []
    for a, b in pairs:
        counts = np.zeros((n_fractions, n_fractions), dtype=int)
        # rows = fraction of method b (y-axis), columns = fraction of method a (x-axis)
        np.add.at(counts, (fractions[b].loc[fractions[a].index].to_numpy() - 1,
                           fractions[a].to_numpy() - 1), 1)
        mats.append(counts)
    vmax = max(m.max() for m in mats)

    fig = new_figure(figsize=(3.0 * ncols + 1.2, 3.0 * nrows), layout="constrained")
    axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
    ticks = np.arange(n_fractions)
    for ax, (a, b), counts in zip(axes, pairs, mats):
        im = ax.imshow(counts, origin="lower", cmap="cool", vmin=0, vmax=vmax, aspect="equal")
        for i, j in np.ndindex(counts.shape):
            ax.text(j, i, counts[i, j], ha="center", va="center", fontsize=8)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks + 1, fontsize=8)
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks + 1, fontsize=8)
        ax.set_xlabel(f"{a} fraction", fontsize=9)
        ax.set_ylabel(f"{b} fraction", fontsize=9)
    for ax in axes[len(pairs):]:
        ax.set_visible(False)
    fig.colorbar(im, ax=list(axes[:len(pairs)]), label="Number of Peptides", shrink=0.8)
    return fig