#!/bin/bash
# install.sh — 安装 fan-skill 到 Claude Code
set -euo pipefail

SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== Fan-Skill Install ==="
echo "安装目标: $SKILLS_DIR"

mkdir -p "$SKILLS_DIR"

for skill_dir in "$SCRIPT_DIR"/bio-plant-*/; do
    skill=$(basename "$skill_dir")
    if [ -d "$SKILLS_DIR/$skill" ]; then
        echo "  覆盖: $skill (已存在)"
    fi
    echo "  安装: $skill"
    cp -rp "$SCRIPT_DIR/$skill" "$SKILLS_DIR/$skill"
done

echo "=== 安装完成 ==="
echo "Skills 已安装到: $SKILLS_DIR"
echo ""
echo "使用方式: 在 Claude Code 中引用 Skill 名称即可"
echo "  bio-plant-gwas"
echo "  bio-plant-population"
echo "  ..."
