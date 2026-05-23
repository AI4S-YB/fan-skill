#!/bin/bash
# install_deps.sh — auto-install missing dependencies
set -euo pipefail
echo "=== Fan-Skill Dependency Installer ==="

# Priority: pixi > conda > system
install_with_pixi() {
    if command -v pixi &>/dev/null; then
        echo "[pixi] Installing packages..."
        pixi install 2>/dev/null && return 0
    fi
    return 1
}

install_with_conda() {
    if command -v conda &>/dev/null; then
        echo "[conda] Creating environment from environment.yaml..."
        conda env create -f environment.yaml 2>/dev/null && return 0
    fi
    return 1
}

check_r_packages() {
    echo "--- R package check ---"
    Rscript -e '
    pkgs <- c("ggplot2", "data.table", "DESeq2", "edgeR", "limma", "WGCNA",
              "GAPIT3", "rrBLUP", "BGLR", "sommer", "lme4", "qqman",
              "clusterProfiler", "pheatmap", "ape", "phytools")
    for (p in pkgs) {
        cat(if (requireNamespace(p, quietly=TRUE)) paste0("  [OK] ", p, "\n")
            else paste0("  [MISSING] ", p, " — install.packages(\"", p, "\") or BiocManager::install(\"", p, "\")\n"))
    }'
}

install_with_pixi || install_with_conda || echo "[WARN] No pixi or conda found. Install dependencies manually."
check_r_packages
echo "=== Done ==="
