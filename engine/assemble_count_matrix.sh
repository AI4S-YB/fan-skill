#!/bin/bash
# assemble_count_matrix.sh — merge per-sample count files into count matrix
# Usage: bash assemble_count_matrix.sh <count_dir> [output_prefix]
set -euo pipefail

COUNT_DIR="${1:?Usage: assemble_count_matrix.sh <count_dir> [output_prefix]}"
PREFIX="${2:-count_matrix}"

echo "=== Assemble Count Matrix ==="
echo "Directory: $COUNT_DIR"

export COUNT_DIR PREFIX
python3 << 'PYEOF'
import os, sys, glob
import pandas as pd

count_dir = os.environ["COUNT_DIR"]
prefix = os.environ["PREFIX"]

files = sorted(glob.glob(f"{count_dir}/*.count") + glob.glob(f"{count_dir}/*.txt"))
if not files:
    print("ERROR: No .count or .txt files found")
    sys.exit(1)

print(f"Sample files: {len(files)}")

first = pd.read_csv(files[0], sep="\t")

# Find gene ID column
id_col = first.columns[0]
for col in first.columns:
    if "gene" in col.lower() or "id" in col.lower():
        id_col = col
        break

# Find count column
count_col = "counts"
if count_col not in first.columns:
    for col in first.columns:
        if "count" in col.lower() or "num" in col.lower():
            count_col = col
            break

print(f"Gene ID column: {id_col}, Count column: {count_col}")

# Merge all files
matrix = first[[id_col, count_col]].copy()
matrix.columns = [id_col, os.path.basename(files[0]).replace(".count","").replace(".txt","")]

for f in files[1:]:
    df = pd.read_csv(f, sep="\t")
    sample_name = os.path.basename(f).replace(".count","").replace(".txt","")
    sub = df[[id_col, count_col]].copy()
    sub.columns = [id_col, sample_name]
    matrix = matrix.merge(sub, on=id_col, how="outer")

matrix = matrix.set_index(id_col).fillna(0).astype(int)

outfile = f"{prefix}.csv"
matrix.to_csv(outfile)
print(f"Saved: {outfile} ({matrix.shape[0]} genes x {matrix.shape[1]} samples)")
PYEOF
