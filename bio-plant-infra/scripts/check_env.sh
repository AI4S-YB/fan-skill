#!/bin/bash
# check_env.sh — 环境自检，按优先级发现工具提供方式
# 用法: bash check_env.sh

echo "=== 环境检测 ==="
echo "主机: $(hostname 2>/dev/null || echo unknown)"
echo "时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "CPU: $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo unknown) cores"

# 检测工具提供方式
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
        echo "  [MISSING] $name — 未找到可用提供方式"
    fi
}

echo "--- CLI 工具 ---"
check_tool "plink"
check_tool "plink2"
check_tool "R"

# 检查 R 包
echo "--- R 包检测 ---"
Rscript -e '
pkgs <- c("ggplot2", "qqman", "data.table", "GAPIT3")
for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  [OK] %s\n", pkg))
    } else {
        cat(sprintf("  [MISSING] %s\n", pkg))
    }
}
' 2>/dev/null || echo "  [WARN] R 不可用，无法检测 R 包"

echo "=== 环境检测完成 ==="
