#!/bin/bash
# install.sh — 安装 plant-bioinfo-skills 到 Claude Code
set -e

SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Plant-Bioinfo-Skills Install ==="
echo "安装目标: $SKILLS_DIR"

mkdir -p "$SKILLS_DIR"

for skill in bio-plant-infra bio-plant-gwas bio-plant-population bio-plant-genomic-selection bio-plant-rnaseq bio-plant-comparative bio-plant-marker; do
    if [ -d "$SCRIPT_DIR/$skill" ]; then
        echo "  安装: $skill"
        cp -r "$SCRIPT_DIR/$skill" "$SKILLS_DIR/$skill"
    fi
done

echo "=== 安装完成 ==="
echo "Skills 已安装到: $SKILLS_DIR"
echo ""
echo "使用方式: 在 Claude Code 中引用 Skill 名称即可"
echo "  bio-plant-gwas"
echo "  bio-plant-population"
echo "  ..."
