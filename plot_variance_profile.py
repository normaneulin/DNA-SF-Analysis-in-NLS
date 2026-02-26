"""
plot_variance_profile.py
------------------------
Your explicit goal: analyse variance of shape features between
nucleosomal and linker sequences.

This script plots the per-position variance (and standard deviation)
for each class across all 147 positions, and also plots the
variance ratio (nucleosomal / linker) to highlight positions where
one class is markedly more variable than the other.

Output: figures/variance/
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\figures\variance"
BASE_NAME     = "nucleosomes_vs_linkers_sapiens_{feature}.csv"

FEATURES  = ["EP", "HelT", "MGW", "ProT", "Roll"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]
POSITIONS = np.arange(1, N_POS + 1)

COLORS    = {"nucleosomal": "#2166AC", "linker": "#D6604D"}
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def variance_profile(df, label):
    sub = df[df["label"] == label][POS_COLS].replace(0.0, np.nan)
    return sub.var(axis=0, ddof=1).values


def load_variances(feature):
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)
    nuc  = variance_profile(df, "nucleosomal")
    lnk  = variance_profile(df, "linker")
    return nuc, lnk


def plot_feature_variance(feature):
    nuc_var, lnk_var = load_variances(feature)

    fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)

    # ── Panel 1: variance per class ───────────────────────────────────────────
    ax = axes[0]
    for label, var in [("nucleosomal", nuc_var), ("linker", lnk_var)]:
        ax.plot(POSITIONS, var, color=COLORS[label], linewidth=1.3,
                label=label)
    ax.set_ylabel("Variance", fontsize=10)
    ax.set_title(f"{feature} — Per-position Variance", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.5)

    # ── Panel 2: standard deviation per class ─────────────────────────────────
    ax = axes[1]
    for label, var in [("nucleosomal", nuc_var), ("linker", lnk_var)]:
        ax.plot(POSITIONS, np.sqrt(var), color=COLORS[label], linewidth=1.3,
                label=label)
    ax.set_ylabel("Std Dev", fontsize=10)
    ax.set_title(f"{feature} — Per-position Standard Deviation", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.5)

    # ── Panel 3: variance ratio (Nuc / Lnk) ───────────────────────────────────
    ax = axes[2]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(lnk_var > 0, nuc_var / lnk_var, np.nan)
    ax.plot(POSITIONS, ratio, color="#5E4FA2", linewidth=1.3)
    ax.axhline(1.0, color="grey", linestyle="--", linewidth=1.0,
               label="ratio = 1 (equal variance)")
    ax.set_ylabel("Variance Ratio\n(Nuc / Lnk)", fontsize=10)
    ax.set_xlabel("Position (bp)", fontsize=10)
    ax.set_title(f"{feature} — Variance Ratio per Position", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.5)

    for ax in axes:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"variance_profile_{feature}.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Saved → {out}")


def plot_combined_variance():
    """
    Single summary figure: one subplot per feature showing variance
    for both classes side by side.
    """
    fig, axes = plt.subplots(len(FEATURES), 1,
                             figsize=(16, 3.5 * len(FEATURES)),
                             sharex=True)

    for ax, feat in zip(axes, FEATURES):
        nuc_var, lnk_var = load_variances(feat)
        for label, var in [("nucleosomal", nuc_var), ("linker", lnk_var)]:
            ax.plot(POSITIONS, var, color=COLORS[label], linewidth=1.2,
                    label=label)
        ax.set_ylabel(f"{feat}\nVariance", fontsize=9)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.4)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

    axes[-1].set_xlabel("Position (bp)", fontsize=11)
    fig.suptitle("Per-position Variance: Nucleosomal vs Linker\n(all 5 features)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "variance_profile_all_features.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out}")


def save_variance_csv():
    """Export a summary CSV: feature, label, position, variance, std."""
    rows = []
    for feat in FEATURES:
        path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feat))
        df   = pd.read_csv(path)
        for label in ["nucleosomal", "linker"]:
            sub = df[df["label"] == label][POS_COLS].replace(0.0, np.nan)
            var = sub.var(axis=0, ddof=1).values
            std = np.sqrt(var)
            for pos_i, (v, s) in enumerate(zip(var, std), start=1):
                rows.append({
                    "feature":  feat,
                    "label":    label,
                    "position": pos_i,
                    "variance": round(v, 6) if not np.isnan(v) else None,
                    "std_dev":  round(s, 6) if not np.isnan(s) else None,
                })
    out_csv = os.path.join(OUTPUT_DIR, "variance_summary.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"  Variance summary CSV → {out_csv}")


def main():
    print(f"\nOutput dir: {OUTPUT_DIR}\n")
    for feat in FEATURES:
        print(f"Plotting variance for {feat} …")
        plot_feature_variance(feat)
    print("\nPlotting combined variance figure …")
    plot_combined_variance()
    print("\nSaving variance summary CSV …")
    save_variance_csv()
    print("\nDone.\n")


if __name__ == "__main__":
    main()