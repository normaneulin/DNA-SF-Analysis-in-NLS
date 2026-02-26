"""
plot_distribution.py
--------------------
For each of the 5 DNA shape features, produces a violin + box plot
comparing the overall value distributions of nucleosomal vs linker
sequences (all positions pooled together).

Useful for seeing global shifts in the feature ranges between classes.

Output: one PNG per feature + one combined 5-panel PNG.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\figures\distribution"
BASE_NAME     = "nucleosomes_vs_linkers_sapiens_{feature}.csv"

FEATURES  = ["EP", "HelT", "MGW", "ProT", "Roll"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]

PALETTE   = {"nucleosomal": "#2166AC", "linker": "#D6604D"}
LABELS    = ["nucleosomal", "linker"]
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_flat(feature: str):
    """Return dict {label: 1-D array of all non-zero values across all positions}."""
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)
    result = {}
    for lbl in LABELS:
        vals = df[df["label"] == lbl][POS_COLS].values.flatten()
        vals = vals[vals != 0.0]        # remove boundary padding zeros
        vals = vals[~np.isnan(vals)]
        result[lbl] = vals
    return result


def violin_plot(ax, data_dict, feature):
    """Draw side-by-side violin + overlaid box for two classes."""
    positions = [1, 2]
    parts = ax.violinplot(
        [data_dict[lbl] for lbl in LABELS],
        positions=positions,
        showmedians=False,
        showextrema=False,
        widths=0.7,
    )
    for body, lbl in zip(parts["bodies"], LABELS):
        body.set_facecolor(PALETTE[lbl])
        body.set_edgecolor("black")
        body.set_linewidth(0.6)
        body.set_alpha(0.7)

    # Overlay box plot
    bp = ax.boxplot(
        [data_dict[lbl] for lbl in LABELS],
        positions=positions,
        widths=0.12,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="white", linewidth=2),
        boxprops=dict(facecolor="black", linewidth=0.8),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )
    _ = bp  # suppress unused warning

    ax.set_xticks(positions)
    ax.set_xticklabels(LABELS, fontsize=9)
    ax.set_ylabel(feature, fontsize=10)
    ax.set_title(feature, fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)


def main():
    print(f"\nOutput dir: {OUTPUT_DIR}\n")

    # ── individual plots ──────────────────────────────────────────────────────
    for feat in FEATURES:
        print(f"Plotting {feat} …")
        data = load_flat(feat)
        fig, ax = plt.subplots(figsize=(5, 5))
        violin_plot(ax, data, feat)
        ax.set_title(
            f"Value Distribution — {feat}\n(all positions pooled)",
            fontsize=11)
        plt.tight_layout()
        out = os.path.join(OUTPUT_DIR, f"distribution_{feat}.png")
        plt.savefig(out, dpi=150)
        plt.close()
        print(f"  Saved → {out}")

    # ── combined 5-panel plot ─────────────────────────────────────────────────
    print("\nPlotting combined figure …")
    fig, axes = plt.subplots(1, 5, figsize=(18, 5))
    for ax, feat in zip(axes, FEATURES):
        data = load_flat(feat)
        violin_plot(ax, data, feat)

    legend_patches = [
        mpatches.Patch(color=PALETTE[lbl], label=lbl) for lbl in LABELS
    ]
    fig.legend(handles=legend_patches, loc="upper right",
               fontsize=10, framealpha=0.9)
    fig.suptitle("DNA Shape Feature Distributions: Nucleosomal vs Linker",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "distribution_all_features.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out}")

    print("\nDone.\n")


if __name__ == "__main__":
    main()