"""
shape_fasta_to_csv.py
---------------------
Parses DNAshapeR output files (.EP, .HelT, .MGW, .ProT, .Roll) into
per-position CSVs that mirror the column layout of the sequence CSV:

    sequence_id | label | length | pos_1 | pos_2 | … | pos_147

Rules applied
-------------
- "NA" values are replaced with 0.
- Roll and HelT naturally produce 146 values per sequence; a 0 is
  appended at pos_147 so every shape CSV has the same 147-column width.
- Label (nucleosomal / linker) is inferred from the sequence header.

Usage
-----
    python shape_fasta_to_csv.py          # processes all 5 shape files
    python shape_fasta_to_csv.py --shape MGW   # single feature
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter

# ── paths ─────────────────────────────────────────────────────────────────────
SHAPE_DIR  = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Caenorhabditis_elegans\fasta\shape"
OUTPUT_DIR = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Caenorhabditis_elegans\csv_processed\shape"
BASE_NAME  = "nucleosomes_vs_linkers_elegans.fas"

# Features and their canonical lengths AFTER padding
FEATURES = {
    "MGW":  147,
    "ProT": 147,
    "Roll": 147,   # native 146 → padded to 147
    "HelT": 147,   # native 146 → padded to 147
    "EP":   147,
}

# Native (pre-padding) lengths for Roll / HelT
NATIVE_LENGTHS = {
    "Roll": 146,
    "HelT": 146,
}
# ──────────────────────────────────────────────────────────────────────────────


def classify_label(header: str) -> str:
    h = header.lower()
    if "nucleosom" in h:
        return "nucleosomal"
    if "linker" in h:
        return "linker"
    return "unknown"


def parse_shape_fasta(filepath: str):
    """
    Generator → (header, [float, ...])

    The shape files look like normal FASTA but values are comma-separated
    floats that may span multiple continuation lines.  Whitespace (including
    newlines) between values is ignored.
    """
    header    = None
    raw_parts = []          # accumulate raw text until next header

    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, _parse_values("".join(raw_parts))
                # values may follow the header on the same line
                parts = line[1:].split(None, 1)
                header    = parts[0]
                raw_parts = [parts[1]] if len(parts) > 1 else []
            else:
                raw_parts.append(" " + line)

    if header is not None:
        yield header, _parse_values("".join(raw_parts))


def _parse_values(text: str):
    """Split comma/space-delimited text into floats; NA → 0.0."""
    tokens = re.split(r"[,\s]+", text.strip())
    values = []
    for t in tokens:
        if not t:
            continue
        if t.upper() == "NA":
            values.append(0.0)
        else:
            try:
                values.append(float(t))
            except ValueError:
                print(f"  [WARNING] Could not parse token '{t}', substituting 0.0",
                      file=sys.stderr)
                values.append(0.0)
    return values


def build_rows(filepath: str, feature: str, target_len: int):
    """Return list of row dicts and a summary Counter."""
    native_len = NATIVE_LENGTHS.get(feature, target_len)
    rows       = []
    issues     = Counter()

    for header, values in parse_shape_fasta(filepath):
        n = len(values)

        # ── length validation / correction ────────────────────────────────
        if n == native_len and target_len > native_len:
            # Expected short feature (Roll / HelT): pad with one 0
            values.append(0.0)
            n += 1

        if n < target_len:
            print(f"  [WARNING] {header}: got {n} values, expected {target_len}. "
                  f"Right-padding with zeros.", file=sys.stderr)
            values.extend([0.0] * (target_len - n))
            issues["short"] += 1
        elif n > target_len:
            print(f"  [WARNING] {header}: got {n} values, expected {target_len}. "
                  f"Truncating to {target_len}.", file=sys.stderr)
            values = values[:target_len]
            issues["long"] += 1

        row = {
            "sequence_id": header,
            "label":       classify_label(header),
            "length":      target_len,
        }
        for i, v in enumerate(values, start=1):
            row[f"pos_{i}"] = v

        rows.append(row)

    return rows, issues


def write_csv(rows, target_len: int, output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    pos_cols   = [f"pos_{i}" for i in range(1, target_len + 1)]
    fieldnames = ["sequence_id", "label", "length"] + pos_cols

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            for col in pos_cols:
                row.setdefault(col, 0.0)
            writer.writerow(row)

    print(f"  Wrote {len(rows):,} rows × {len(fieldnames):,} columns → {output_path}")


def summarise(rows, feature: str, issues: Counter):
    label_counts = Counter(r["label"] for r in rows)
    print(f"\n── {feature} Summary ─────────────────────────────────────")
    print(f"  Total sequences : {len(rows):,}")
    for lbl, cnt in sorted(label_counts.items()):
        print(f"  {lbl:>15s}  : {cnt:,}")
    if issues:
        for k, v in issues.items():
            print(f"  [!] length {k}  : {v} sequence(s) corrected")
    print("─────────────────────────────────────────────────────\n")


def process_feature(feature: str):
    target_len = FEATURES[feature]
    input_path = os.path.join(SHAPE_DIR, f"{BASE_NAME}.{feature}")
    output_path = os.path.join(OUTPUT_DIR,
                               f"nucleosomes_vs_linkers_elegans_{feature}.csv")

    print(f"\nProcessing {feature} …")
    print(f"  Input  : {input_path}")

    if not os.path.isfile(input_path):
        print(f"  [ERROR] File not found, skipping.", file=sys.stderr)
        return

    rows, issues = build_rows(input_path, feature, target_len)

    if not rows:
        print(f"  [ERROR] No sequences parsed for {feature}.", file=sys.stderr)
        return

    summarise(rows, feature, issues)
    write_csv(rows, target_len, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Convert DNAshapeR FASTA shape files to per-position CSVs.")
    parser.add_argument(
        "--shape", choices=list(FEATURES.keys()), default=None,
        help="Process a single shape feature (default: all five)")
    args = parser.parse_args()

    features = [args.shape] if args.shape else list(FEATURES.keys())

    print(f"\nShape dir  : {SHAPE_DIR}")
    print(f"Output dir : {OUTPUT_DIR}")
    print(f"Features   : {', '.join(features)}")

    for feat in features:
        process_feature(feat)

    print("\nAll done.\n")


if __name__ == "__main__":
    main()