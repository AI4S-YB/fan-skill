#!/bin/bash
# check_env.sh — environment self-check, discover tool providers by priority
# Usage: bash check_env.sh

echo "=== Environment Check ==="
echo "Host: $(hostname 2>/dev/null || echo unknown)"
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "CPU: $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo unknown) cores"

# Detect tool provider
check_tool() {
    local name=$1
    local found=""

    # Priority: pixi > conda > singularity > docker > system
    if command -v pixi &>/dev/null && pixi run "$name" --version &>/dev/null 2>&1; then
        found="pixi"
    elif command -v conda &>/dev/null && conda run -n base "$name" --version &>/dev/null 2>&1; then
        found="conda:base"
    elif ls ./*.sif &>/dev/null 2>&1; then
        found="singularity"
    elif command -v docker &>/dev/null && docker run --rm "$name" --version &>/dev/null 2>&1; then
        found="docker"
    elif command -v "$name" &>/dev/null; then
        found="system"
    fi

    if [ -n "$found" ]; then
        echo "  [OK] $name → $found"
    else
        echo "  [MISSING] $name — no available provider found"
    fi
}

echo "--- CLI Tools ---"
check_tool "plink"
check_tool "plink2"
check_tool "R"

# Check R packages
echo "--- R Package Check ---"
Rscript -e '
pkgs <- c("ggplot2", "data.table", "DESeq2", "edgeR", "limma", "WGCNA",
          "GAPIT3", "rrBLUP", "BGLR", "sommer", "lme4", "qqman",
          "clusterProfiler", "pheatmap", "ape", "phytools")
for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  [OK] %s\n", pkg))
    } else {
        cat(sprintf("  [MISSING] %s\n", pkg))
    }
}
' 2>/dev/null || echo "  [WARN] R not available, cannot check R packages"

echo "=== Environment Check Complete ==="
