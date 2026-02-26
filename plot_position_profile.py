"""
plot_position_profile.py
------------------------
For each of the 5 DNA shape features, plots the mean value at every
position (1–147) separately for nucleosomal and linker sequences.

This is the most biologically informative view: nucleosomal DNA often
shows ~10 bp periodicity due to helical winding around the histone core.

Output: one PNG per feature in the figures/ subfolder.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\figures\position_profile"
BASE_NAME     = "nucleosomes_vs_linkers_sapiens_{feature}.csv"

FEATURES  = ["EP", "HelT", "MGW", "ProT", "Roll"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]

COLORS    = {"nucleosomal": "#2166AC", "linker": "#D6604D"}
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_shape(feature: str) -> pd.DataFrame:
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)
    # Replace padded zeros at boundary positions with NaN so they don't
    # drag the mean down (pos_1, pos_2, pos_146, pos_147 for Roll/HelT)
    return df


def get_mean_profile(df: pd.DataFrame, label: str):
    sub  = df[df["label"] == label][POS_COLS]
    # Treat 0.0 at boundary positions as NaN for mean calculation
    sub  = sub.replace(0.0, np.nan)
    return sub.mean(axis=0).values   # shape (147,)


def plot_feature(feature: str):
    df = load_shape(feature)

    positions = np.arange(1, N_POS + 1)

    fig, ax = plt.subplots(figsize=(14, 4))

    for label, color in COLORS.items():
        mean_vals = get_mean_profile(df, label)
        n         = (df["label"] == label).sum()
        ax.plot(positions, mean_vals, color=color, linewidth=1.4,
                label=f"{label} (n={n:,})")
        # ± 1 SD shading
        sub  = df[df["label"] == label][POS_COLS].replace(0.0, np.nan)
        std  = sub.std(axis=0).values
        ax.fill_between(positions,
                        mean_vals - std, mean_vals + std,
                        color=color, alpha=0.15)

    ax.set_xlabel("Position along sequence (bp)", fontsize=11)
    ax.set_ylabel(f"Mean {feature}", fontsize=11)
    ax.set_title(f"Positional Profile — {feature}\n"
                 f"Mean ± 1 SD across nucleosomal vs linker sequences",
                 fontsize=12)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.legend(fontsize=10)
    ax.grid(axis="x", which="major", linestyle="--", linewidth=0.4, alpha=0.5)
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.4, alpha=0.5)

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"position_profile_{feature}.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Saved → {out}")


def main():
    print(f"\nOutput dir: {OUTPUT_DIR}\n")
    for feat in FEATURES:
        print(f"Plotting {feat} …")
        plot_feature(feat)
    print("\nDone.\n")


if __name__ == "__main__":
    main()