#!/bin/bash
# uninstall.sh — uninstall fan-skill from all AI coding agents
set -euo pipefail

echo "=== Fan-Skill Uninstall ==="

# Claude Code
CLAUDE_SKILLS_DIR="${CLAUDE_SKILLS_DIR:-$HOME/.claude/skills}"
if [ -d "$CLAUDE_SKILLS_DIR/fan-skill" ]; then
    echo "  Removing from Claude Code: $CLAUDE_SKILLS_DIR/fan-skill"
    rm -rf "$CLAUDE_SKILLS_DIR/fan-skill"
fi

# Codex CLI
CODEX_SKILLS_DIR="${CODEX_SKILLS_DIR:-$HOME/.codex/skills}"
if [ -d "$CODEX_SKILLS_DIR/fan-skill" ]; then
    echo "  Removing from Codex CLI: $CODEX_SKILLS_DIR/fan-skill"
    rm -rf "$CODEX_SKILLS_DIR/fan-skill"
fi

# Gemini CLI
GEMINI_SKILLS_DIR="${GEMINI_SKILLS_DIR:-$HOME/.gemini/skills}"
if [ -d "$GEMINI_SKILLS_DIR/fan-skill" ]; then
    echo "  Removing from Gemini CLI: $GEMINI_SKILLS_DIR/fan-skill"
    rm -rf "$GEMINI_SKILLS_DIR/fan-skill"
fi

# OpenCode (project-level)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OPENCODE_SKILL_DIR="$SCRIPT_DIR/.opencode/skill/fan-skill"
if [ -d "$OPENCODE_SKILL_DIR" ]; then
    echo "  Removing from OpenCode: $OPENCODE_SKILL_DIR"
    rm -rf "$OPENCODE_SKILL_DIR"
fi

echo "=== Uninstall complete ==="
echo "fan-skill has been removed from all detected locations."
