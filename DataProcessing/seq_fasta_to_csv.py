"""
seq_fasta_to_csv.py
---------------
Parses a FASTA file containing nucleosomal and linker sequences into a
per-base CSV where each row is one sequence and each base occupies its
own column (base_1, base_2, ..., base_N).

This layout matches the per-position output of DNAshapeR so that shape
feature values (MGW, ProT, Roll, HelT, EP) can be merged directly.

Usage
-----
Run as-is (paths are hard-coded below) or pass arguments on the command line:
    python seq_fasta_to_csv.py
    python seq_fasta_to_csv.py --input path/to/file.fas --output path/to/file.csv
"""

import argparse
import csv
import os
import re
import sys


# ── hard-coded defaults ────────────────────────────────────────────────────────
DEFAULT_INPUT  = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\fasta\sequence\nucleosomes_vs_linkers_sapiens.fas"
DEFAULT_OUTPUT = r"C:\Software Projects\Bioinformatics\DNA-SF-Analysis-in-NLS\data\DatasetNup_1\Homo_sapiens\csv_processed\sequence\nucleosomes_vs_linkers_sapiens.csv"
# ──────────────────────────────────────────────────────────────────────────────


def classify_label(header: str) -> str:
    """Return 'nucleosomal' or 'linker' based on the FASTA header."""
    h = header.lower()
    if "nucleosom" in h:
        return "nucleosomal"
    if "linker" in h:
        return "linker"
    return "unknown"


def parse_fasta(filepath: str):
    """
    Generator that yields (header, sequence) tuples.
    Handles multi-line sequences and sequences whose bases are on the
    same line as the header (some non-standard FASTA files).
    """
    header  = None
    seq_parts = []

    with open(filepath, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Flush previous record
                if header is not None:
                    yield header, "".join(seq_parts)

                # Some files put sequence on the same line as the header,
                # separated by whitespace: ">header_1 ATCG..."
                parts = line[1:].split(None, 1)          # split on first whitespace
                header = parts[0]
                seq_parts = [re.sub(r"[^ACGTNacgtn]", "", parts[1])] if len(parts) > 1 else []
            else:
                # Normal continuation line — strip non-IUPAC characters
                seq_parts.append(re.sub(r"[^ACGTNacgtn]", "", line))

    # Flush last record
    if header is not None:
        yield header, "".join(seq_parts)


def build_rows(filepath: str):
    """
    Return:
        rows      – list of dicts, one per sequence
        max_len   – maximum sequence length seen (used to build column headers)
    """
    rows    = []
    max_len = 0

    for header, seq in parse_fasta(filepath):
        seq_upper = seq.upper()
        length    = len(seq_upper)
        if length == 0:
            print(f"  [WARNING] Empty sequence skipped: {header}", file=sys.stderr)
            continue

        max_len = max(max_len, length)

        row = {
            "sequence_id": header,
            "label":       classify_label(header),
            "length":      length,
            # "sequence":    seq_upper,   # full sequence as a single string (handy reference)
        }

        # Per-base columns: base_1 … base_N
        for i, base in enumerate(seq_upper, start=1):
            row[f"base_{i}"] = base

        rows.append(row)

    return rows, max_len


def write_csv(rows, max_len: int, output_path: str):
    """Write the list of row dicts to a CSV file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Column order: metadata first, then per-base columns
    base_cols   = [f"base_{i}" for i in range(1, max_len + 1)]
    fieldnames  = ["sequence_id", "label", "length"] + base_cols

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()

        for row in rows:
            # Fill missing base columns with empty string (shorter sequences)
            for col in base_cols:
                row.setdefault(col, "")
            writer.writerow(row)

    print(f"  Wrote {len(rows):,} rows × {len(fieldnames):,} columns → {output_path}")


def summarise(rows):
    """Print a quick sanity-check summary."""
    from collections import Counter
    label_counts = Counter(r["label"] for r in rows)
    lengths      = [r["length"] for r in rows]

    print("\n── Summary ──────────────────────────────────────────")
    print(f"  Total sequences : {len(rows):,}")
    for lbl, cnt in sorted(label_counts.items()):
        print(f"  {lbl:>15s}  : {cnt:,}")
    print(f"  Min seq length  : {min(lengths)}")
    print(f"  Max seq length  : {max(lengths)}")
    print(f"  Mean seq length : {sum(lengths)/len(lengths):.1f}")
    print("─────────────────────────────────────────────────────\n")


def main():
    parser = argparse.ArgumentParser(description="Convert a FASTA file to a per-base CSV.")
    parser.add_argument("--input",  default=DEFAULT_INPUT,  help="Path to the input .fas / .fasta file")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Path for the output .csv file")
    args = parser.parse_args()

    print(f"\nInput  : {args.input}")
    print(f"Output : {args.output}\n")

    if not os.path.isfile(args.input):
        sys.exit(f"[ERROR] Input file not found: {args.input}")

    print("Parsing FASTA …")
    rows, max_len = build_rows(args.input)

    if not rows:
        sys.exit("[ERROR] No sequences were parsed. Check the input file format.")

    summarise(rows)

    print("Writing CSV …")
    write_csv(rows, max_len, args.output)
    print("Done.\n")


if __name__ == "__main__":
    main()