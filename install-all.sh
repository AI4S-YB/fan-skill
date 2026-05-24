#!/bin/bash
# install-all.sh — install fan-skill for all detected AI coding agents
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
INSTALLED=0

echo "=== Fan-Skill Install (All Agents) ==="

# Claude Code
if [ -d "$HOME/.claude" ] || command -v claude &>/dev/null 2>&1; then
    echo "  Detected: Claude Code"
    bash "$SCRIPT_DIR/install-claude.sh"
    INSTALLED=$((INSTALLED + 1))
fi

# Codex CLI
if [ -d "$HOME/.codex" ] || command -v codex &>/dev/null 2>&1; then
    echo "  Detected: Codex CLI"
    bash "$SCRIPT_DIR/install-codex.sh"
    INSTALLED=$((INSTALLED + 1))
fi

# Gemini CLI
if [ -d "$HOME/.gemini" ] || command -v gemini &>/dev/null 2>&1; then
    echo "  Detected: Gemini CLI"
    bash "$SCRIPT_DIR/install-gemini.sh"
    INSTALLED=$((INSTALLED + 1))
fi

if [ $INSTALLED -eq 0 ]; then
    echo "No supported agents detected. Install manually:"
    echo "  bash install-claude.sh    # Claude Code"
    echo "  bash install-codex.sh     # Codex CLI"
    echo "  bash install-gemini.sh    # Gemini CLI"
    echo "  bash install-opencode.sh  # OpenCode"
fi

echo "=== Done ($INSTALLED agents) ==="
