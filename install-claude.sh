#!/bin/bash
# install-claude.sh — install fan-skill for Claude Code
set -euo pipefail

SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALL_DIR="$SKILLS_DIR/fan-skill"

echo "=== Fan-Skill Install (Claude Code) ==="
echo "Install target: $INSTALL_DIR"

# Create parent directory first
mkdir -p "$SKILLS_DIR"

# Backup existing installation if user has modifications
if [ -d "$INSTALL_DIR" ]; then
    echo "  Existing installation found, removing..."
    rm -rf "$INSTALL_DIR"
fi
mkdir -p "$INSTALL_DIR"

# Copy all required components (including skill.yaml and docs)
for component in SKILL.md skill.yaml engine knowledge-base tool-registry references templates theme docs; do
    if [ -e "$SCRIPT_DIR/$component" ]; then
        echo "  Installing: $component"
        cp -rp "$SCRIPT_DIR/$component" "$INSTALL_DIR/$component"
    fi
done

# Verify installation
if [ ! -f "$INSTALL_DIR/skill.yaml" ]; then
    echo "Warning: skill.yaml not found, installation may be incomplete"
fi

echo "=== Installation complete ==="
echo "fan-skill installed to: $INSTALL_DIR"
echo "To use: just describe your biological question in Claude Code."
