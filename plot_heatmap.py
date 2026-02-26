"""
plot_heatmap.py
---------------
Produces heatmaps of mean shape values across all 147 positions for
each class (nucleosomal / linker), plus a difference heatmap.

Layout per figure:
    Row 1 — nucleosomal mean
    Row 2 — linker mean
    Row 3 — difference (nucleosomal − linker)

One combined heatmap showing all 5 features stacked.

Output: figures/heatmap/
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\figures\heatmap"
BASE_NAME     = "nucleosomes_vs_linkers_sapiens_{feature}.csv"

FEATURES  = ["EP", "HelT", "MGW", "ProT", "Roll"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]
POSITIONS = np.arange(1, N_POS + 1)
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def mean_profile(df, label):
    sub = df[df["label"] == label][POS_COLS].replace(0.0, np.nan)
    return sub.mean(axis=0).values


def load_means(feature):
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)
    nuc  = mean_profile(df, "nucleosomal")
    lnk  = mean_profile(df, "linker")
    return nuc, lnk


def plot_single_feature_heatmap(feature):
    nuc, lnk = load_means(feature)
    diff = nuc - lnk

    # Stack as 3-row matrix
    mat = np.vstack([nuc, lnk, diff])

    fig, ax = plt.subplots(figsize=(16, 3))

    # Use a diverging colormap for the difference row and sequential for means
    # We'll plot as a single heatmap with a shared scale for simplicity
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
                   extent=[0.5, N_POS + 0.5, 2.5, -0.5])

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["Nucleosomal", "Linker", "Difference\n(Nuc − Lnk)"],
                       fontsize=9)
    ax.set_xlabel("Position (bp)", fontsize=10)
    ax.set_title(f"Mean {feature} Heatmap across 147 positions", fontsize=11)

    # x-ticks every 10 bp
    xticks = np.arange(10, N_POS + 1, 10)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=8)

    plt.colorbar(im, ax=ax, shrink=0.8, label=feature)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"heatmap_{feature}.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Saved → {out}")


def plot_combined_heatmap():
    """
    One large heatmap: rows = features × classes (nucleosomal, linker, diff),
    cols = positions.
    """
    n_features = len(FEATURES)
    n_rows     = n_features * 3          # 3 rows per feature
    mat        = np.zeros((n_rows, N_POS))
    row_labels = []

    for i, feat in enumerate(FEATURES):
        nuc, lnk = load_means(feat)
        # Normalise each feature to z-score so different scales are comparable
        all_vals = np.concatenate([nuc[~np.isnan(nuc)], lnk[~np.isnan(lnk)]])
        mu, sd   = np.nanmean(all_vals), np.nanstd(all_vals)
        sd       = sd if sd > 0 else 1.0
        nuc_z    = (nuc - mu) / sd
        lnk_z    = (lnk - mu) / sd
        diff_z   = nuc_z - lnk_z

        base = i * 3
        mat[base]     = np.nan_to_num(nuc_z)
        mat[base + 1] = np.nan_to_num(lnk_z)
        mat[base + 2] = np.nan_to_num(diff_z)

        row_labels += [f"{feat} Nuc", f"{feat} Lnk", f"{feat} Diff"]

    fig, ax = plt.subplots(figsize=(18, n_rows * 0.55 + 1.5))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3,
                   extent=[0.5, N_POS + 0.5, n_rows - 0.5, -0.5])

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_xlabel("Position (bp)", fontsize=11)
    ax.set_title("Combined DNA Shape Heatmap (z-scored per feature)\n"
                 "Nucleosomal | Linker | Difference", fontsize=12)

    xticks = np.arange(10, N_POS + 1, 10)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=8)

    # Horizontal lines separating features
    for i in range(1, n_features):
        ax.axhline(i * 3 - 0.5, color="black", linewidth=1.2)

    plt.colorbar(im, ax=ax, shrink=0.6, label="z-score")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "heatmap_combined_all_features.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out}")


def main():
    print(f"\nOutput dir: {OUTPUT_DIR}\n")
    for feat in FEATURES:
        print(f"Plotting heatmap for {feat} …")
        plot_single_feature_heatmap(feat)
    print("\nPlotting combined heatmap …")
    plot_combined_heatmap()
    print("\nDone.\n")


if __name__ == "__main__":
    main()