"""
plot_shape_profile_centered.py
-------------------------------
Produces a centered positional mean profile plot for each of the 5 DNA
shape features, styled after the classic DNAshapeR R output:

  - X-axis: bp position relative to the dyad center (-73 ... 0 ... +73)
  - Y-axis: mean shape value
  - Scatter dots (translucent) + smoothing spline line for each class
  - Blue  = nucleosomal  |  Red = linker
  - Vertical dashed line at the center (position 0 / dyad)
  - One PNG per feature + one combined 5-panel PNG

Dependencies: numpy, pandas, matplotlib, scipy
    pip install numpy pandas matplotlib scipy

Usage
-----
    python plot_shape_profile_centered.py
    python plot_shape_profile_centered.py --shape EP
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_smoothing_spline

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Drosophila_melanogaster\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Drosophila_melanogaster\figures\shape_profile_centered"
BASE_NAME     = "nucleosomes_vs_linkers_melanogaster_{feature}.csv"

FEATURES  = ["HelT", "Rise", "Roll", "Shift", "Slide", "Tilt", "Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "EP", "MGW"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]

# X-axis: center 147 positions so position 74 becomes 0  (-73 ... +73)
CENTER_IDX = 74
X_AXIS     = np.arange(1, N_POS + 1) - CENTER_IDX

# Smoothing lambda: None = scipy auto-selects (recommended starting point)
# Increase (e.g. 1e2, 1e4) for a flatter/smoother curve
# Decrease (e.g. 1e-2) to follow the data more closely
SMOOTH_LAM = None

COLORS = {
    "nucleosomal": "#2166AC",   # blue
    "linker":      "#D6604D",   # red
}

Y_LABELS = {
    "HelT":    "Mean HelT (°)",     # Helix Twist — degrees
    "Rise":    "Mean Rise (Å)",     # Rise — Angstroms
    "Roll":    "Mean Roll (°)",     # Roll — degrees
    "Shift":   "Mean Shift (Å)",    # Shift — Angstroms
    "Slide":   "Mean Slide (Å)",    # Slide — Angstroms
    "Tilt":    "Mean Tilt (°)",     # Tilt — degrees
    "Buckle":  "Mean Buckle (°)",   # Buckle — degrees
    "Opening": "Mean Opening (°)",  # Opening — degrees
    "ProT":    "Mean ProT (°)",     # Propeller Twist — degrees
    "Shear":   "Mean Shear (Å)",    # Shear — Angstroms
    "Stagger": "Mean Stagger (Å)",  # Stagger — Angstroms
    "Stretch": "Mean Stretch (Å)",  # Stretch — Angstroms
    "EP":      "Mean EP (kT/e)",    # Electrostatic Potential — thermal energy per elementary charge
    "MGW":     "Mean MGW (Å)",      # Minor Groove Width — Angstroms
}
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def trim_edge_zeros(row: np.ndarray) -> np.ndarray:
    """
    Replace leading and trailing runs of 0.0 with NaN on a single row.

    Only contiguous zeros touching each end are masked — interior zeros
    (e.g. a genuine Roll ~ 0) are left untouched.

    [0, 0, 1.2, 0.0, 3.4, 0, 0]  ->  [NaN, NaN, 1.2, 0.0, 3.4, NaN, NaN]
    """
    result = row.astype(float).copy()
    n = len(result)

    i = 0
    while i < n and result[i] == 0.0:
        result[i] = np.nan
        i += 1

    j = n - 1
    while j >= 0 and result[j] == 0.0:
        result[j] = np.nan
        j -= 1

    return result


def load_mean_profile(feature: str):
    """
    Load shape CSV and return per-position mean for each class.
    Edge zeros are trimmed per-sequence before averaging.
    """
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)

    out    = {}
    counts = {}
    for label in ("nucleosomal", "linker"):
        raw     = df[df["label"] == label][POS_COLS].values.astype(float)
        trimmed = np.vstack([trim_edge_zeros(row) for row in raw])
        out[label]    = np.nanmean(trimmed, axis=0)   # (147,)
        counts[label] = raw.shape[0]

    return out, counts


def smooth_spline(x: np.ndarray, y: np.ndarray, lam=SMOOTH_LAM):
    """
    Fit a smoothing spline to (x, y), skipping NaN positions.
    Returns fitted y values of the same length; NaN where input was NaN.

    make_smoothing_spline (scipy >= 1.10) gives the same R-style
    smooth.spline / loess feel without requiring statsmodels.

    Tune SMOOTH_LAM at the top of the file:
      None        -> scipy auto-selects (good default)
      small value -> wiggly, closely follows data
      large value -> very smooth / flattened
    """
    valid  = ~np.isnan(y)
    result = np.full_like(y, np.nan, dtype=float)

    if valid.sum() < 6:
        return result

    spline        = make_smoothing_spline(x[valid], y[valid], lam=lam)
    result[valid] = spline(x[valid])
    return result


def plot_single(feature: str, ax=None, show_legend=True):
    """
    Draw the centered profile on `ax`.
    If ax is None a standalone figure is created and saved.
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(7, 5))

    means, counts = load_mean_profile(feature)

    for label in ("nucleosomal", "linker"):
        mean_vals = means[label]
        smoothed  = smooth_spline(X_AXIS.astype(float), mean_vals)
        color     = COLORS[label]
        n         = counts[label]

        # Translucent scatter dots — raw per-position means
        valid_pts = ~np.isnan(mean_vals)
        ax.scatter(X_AXIS[valid_pts], mean_vals[valid_pts],
                   color=color, alpha=0.25, s=13, linewidths=0, zorder=2)

        # Smooth spline line — drawn only over valid positions
        valid_fit = ~np.isnan(smoothed)
        ax.plot(X_AXIS[valid_fit], smoothed[valid_fit],
                color=color, linewidth=1.9,
                label=f"{label} (n={n:,})", zorder=3)

    # Center dyad dashed line
    ax.axvline(0, color="black", linestyle="--", linewidth=1.0,
               zorder=4, label="Center (dyad)")

    ax.set_xlabel("Position relative to dyad (bp)", fontsize=10)
    ax.set_ylabel(Y_LABELS.get(feature, f"Mean {feature}"), fontsize=10)
    ax.set_title(feature, fontsize=11, fontweight="bold")
    ax.set_xlim(X_AXIS[0] - 1, X_AXIS[-1] + 1)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.5)

    # X-tick labels: numeric with "Center" at 0
    major_ticks = np.arange(-60, 61, 20)
    ax.set_xticks(major_ticks)
    ax.set_xticklabels(
        ["Center" if t == 0 else (f"+{t}" if t > 0 else str(t))
         for t in major_ticks],
        fontsize=9,
    )
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))

    if show_legend:
        ax.legend(fontsize=9, framealpha=0.9)

    if standalone:
        plt.tight_layout()
        out = os.path.join(OUTPUT_DIR, f"shape_profile_{feature}.png")
        plt.savefig(out, dpi=150)
        plt.close()
        print(f"  Saved -> {out}")


def plot_combined():
    """One figure with all 5 features side by side."""
    fig, axes = plt.subplots(1, 5, figsize=(22, 5), sharey=False)

    for ax, feat in zip(axes, FEATURES):
        plot_single(feat, ax=ax, show_legend=(feat == FEATURES[0]))

    handles, labels_ = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_,
               loc="upper center", ncol=3, fontsize=10,
               bbox_to_anchor=(0.5, 1.03), framealpha=0.9)

    fig.suptitle(
        "DNA Shape Feature Profiles: Nucleosomal vs Linker\n"
        "(mean per position, centered at dyad)",
        fontsize=13, fontweight="bold", y=1.07,
    )
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "shape_profile_all_features.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved -> {out}")


def main():
    parser = argparse.ArgumentParser(
        description="Centered positional mean profile plot for DNA shape features.")
    parser.add_argument("--shape", choices=FEATURES, default=None,
                        help="Plot a single feature (default: all five)")
    args = parser.parse_args()

    features = [args.shape] if args.shape else FEATURES

    print(f"\nOutput dir: {OUTPUT_DIR}\n")
    for feat in features:
        print(f"Plotting {feat} ...")
        plot_single(feat)

    if not args.shape:
        print("Plotting combined figure ...")
        plot_combined()

    print("\nDone.\n")


if __name__ == "__main__":
    main()