"""
plot_statistical_test.py
------------------------
For each of the 147 positions and all 5 shape features, runs a
Mann-Whitney U test comparing nucleosomal vs linker distributions.

Mann-Whitney U is used instead of t-test because:
  - We do not assume normality of shape value distributions.
  - It is robust with large sample sizes (n ≈ 2,273 / 2,300).

Produces:
  1. Per-feature plot of -log10(p-value) across positions, with a
     Bonferroni-corrected significance threshold line.
  2. A combined 5-panel figure.
  3. A CSV of all test results (feature, position, U-stat, p-value,
     significant flag, effect size r = Z / sqrt(N)).

Requires: scipy
    pip install scipy
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import mannwhitneyu

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_CSV_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\shape"
OUTPUT_DIR    = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\figures\statistical_test"
BASE_NAME     = "nucleosomes_vs_linkers_sapiens_{feature}.csv"

FEATURES  = ["EP", "HelT", "MGW", "ProT", "Roll"]
N_POS     = 147
POS_COLS  = [f"pos_{i}" for i in range(1, N_POS + 1)]
POSITIONS = np.arange(1, N_POS + 1)
ALPHA     = 0.05
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(OUTPUT_DIR, exist_ok=True)


def run_tests(feature: str):
    path = os.path.join(SHAPE_CSV_DIR, BASE_NAME.format(feature=feature))
    df   = pd.read_csv(path)

    nuc_df = df[df["label"] == "nucleosomal"][POS_COLS].replace(0.0, np.nan)
    lnk_df = df[df["label"] == "linker"][POS_COLS].replace(0.0, np.nan)

    results = []
    bonferroni_n = N_POS   # tests per feature

    for pos_i, col in enumerate(POS_COLS, start=1):
        nuc_vals = nuc_df[col].dropna().values
        lnk_vals = lnk_df[col].dropna().values

        if len(nuc_vals) < 5 or len(lnk_vals) < 5:
            # Boundary position — skip
            results.append({
                "feature": feature, "position": pos_i,
                "U_stat": np.nan, "p_value": np.nan,
                "neg_log10_p": np.nan, "effect_r": np.nan,
                "significant_bonferroni": False,
            })
            continue

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stat, p = mannwhitneyu(nuc_vals, lnk_vals, alternative="two-sided")

        # Effect size r = Z / sqrt(N)
        n_total = len(nuc_vals) + len(lnk_vals)
        # Approximate Z from U
        mu_u  = len(nuc_vals) * len(lnk_vals) / 2
        sigma = np.sqrt(len(nuc_vals) * len(lnk_vals) *
                        (len(nuc_vals) + len(lnk_vals) + 1) / 12)
        z     = (stat - mu_u) / sigma if sigma > 0 else 0.0
        r     = abs(z) / np.sqrt(n_total)

        results.append({
            "feature": feature,
            "position": pos_i,
            "U_stat": stat,
            "p_value": p,
            "neg_log10_p": -np.log10(p) if p > 0 else np.nan,
            "effect_r": r,
            "significant_bonferroni": p < (ALPHA / bonferroni_n),
        })

    return pd.DataFrame(results)


def plot_feature_tests(res_df, feature):
    positions = res_df["position"].values
    nlp       = res_df["neg_log10_p"].values
    eff       = res_df["effect_r"].values

    bonf_threshold = -np.log10(ALPHA / N_POS)
    nominal_thresh = -np.log10(ALPHA)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    # ── Panel 1: -log10(p) ────────────────────────────────────────────────────
    ax1.plot(positions, nlp, color="#4D9DE0", linewidth=1.2)
    ax1.axhline(bonf_threshold, color="red", linestyle="--", linewidth=1.0,
                label=f"Bonferroni (α/147 = {ALPHA/N_POS:.2e})")
    ax1.axhline(nominal_thresh, color="orange", linestyle=":", linewidth=1.0,
                label=f"Nominal α = {ALPHA}")
    ax1.set_ylabel("−log₁₀(p-value)", fontsize=10)
    ax1.set_title(f"{feature} — Mann-Whitney U per position", fontsize=11)
    ax1.legend(fontsize=9)
    ax1.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.4)

    # Shade significant positions
    sig = res_df["significant_bonferroni"].values
    for i, (pos, s) in enumerate(zip(positions, sig)):
        if s:
            ax1.axvspan(pos - 0.5, pos + 0.5, color="red", alpha=0.15,
                        linewidth=0)

    # ── Panel 2: effect size r ────────────────────────────────────────────────
    ax2.plot(positions, eff, color="#E15759", linewidth=1.2)
    ax2.axhline(0.1, color="grey", linestyle="--", linewidth=0.8,
                label="small (r=0.1)")
    ax2.axhline(0.3, color="grey", linestyle="-.", linewidth=0.8,
                label="medium (r=0.3)")
    ax2.set_ylabel("Effect size r", fontsize=10)
    ax2.set_xlabel("Position (bp)", fontsize=10)
    ax2.legend(fontsize=9)
    ax2.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.4)

    for ax in (ax1, ax2):
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"stat_test_{feature}.png")
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"  Saved → {out}")


def plot_combined(all_results):
    fig, axes = plt.subplots(len(FEATURES), 1,
                             figsize=(16, 3.5 * len(FEATURES)),
                             sharex=True)
    bonf_threshold = -np.log10(ALPHA / N_POS)

    for ax, feat in zip(axes, FEATURES):
        res = all_results[feat]
        ax.plot(res["position"], res["neg_log10_p"],
                color="#4D9DE0", linewidth=1.1)
        ax.axhline(bonf_threshold, color="red", linestyle="--",
                   linewidth=0.9, label="Bonferroni")
        ax.set_ylabel(f"{feat}\n−log₁₀(p)", fontsize=9)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="x", linestyle="--", linewidth=0.4, alpha=0.4)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

    axes[-1].set_xlabel("Position (bp)", fontsize=11)
    fig.suptitle("Mann-Whitney U: Nucleosomal vs Linker — all 5 features",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "stat_test_all_features.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out}")


def main():
    print(f"\nOutput dir: {OUTPUT_DIR}\n")
    all_results = {}
    all_rows    = []

    for feat in FEATURES:
        print(f"Testing {feat} …")
        res = run_tests(feat)
        all_results[feat] = res
        all_rows.append(res)

        n_sig = res["significant_bonferroni"].sum()
        print(f"  Significant positions (Bonferroni): {n_sig} / {N_POS}")
        plot_feature_tests(res, feat)

    print("\nPlotting combined figure …")
    plot_combined(all_results)

    print("\nSaving results CSV …")
    out_csv = os.path.join(OUTPUT_DIR, "mannwhitney_results.csv")
    pd.concat(all_rows, ignore_index=True).to_csv(out_csv, index=False)
    print(f"  Saved → {out_csv}")

    print("\nDone.\n")


if __name__ == "__main__":
    main()