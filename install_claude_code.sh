#!/usr/bin/env bash
set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
err()   { echo -e "${RED}[ERR]${NC}   $*"; }
header(){ echo -e "\n${CYAN}=== $* ===${NC}"; }

export NVM_DIR="${NVM_DIR:-$HOME/.nvm}"
NVM_VERSION="v0.40.3"
NVM_GITEE="https://gitee.com/mirrors/nvm.git"
NODE_MIRROR="https://npmmirror.com/mirrors/node"
NPM_REGISTRY="https://registry.npmmirror.com"
CLAUDE_PKG="@anthropic-ai/claude-code"

# ── Platform detection ──────────────────────────────────────────────

OS="$(uname -s)"
case "$OS" in
    Darwin)  PKG_MGR="brew" ;;
    Linux)   PKG_MGR="apt-get" ;;
    *)       err "Unsupported OS: $OS"; exit 1 ;;
esac

header "检测系统环境"
info "OS: $OS | 包管理: $PKG_MGR"

# ── Network detection ─────────────────────────────────────────────────

header "检测网络环境"

USE_FOREIGN=false

check_site() {
    local name="$1"
    local url="$2"
    local code
    code=$(curl -sL -o /dev/null -w "%{http_code}" --connect-timeout 5 --max-time 10 "$url" 2>/dev/null || echo "000")
    case "$code" in
        200|301|302)
            info "$name → 可达"
            return 0
            ;;
        *)
            warn "$name → 不可达"
            return 1
            ;;
    esac
}

if check_site "GitHub"  "https://github.com"; then
    USE_FOREIGN=true
    info "GitHub 可达，将使用国外源安装"
else
    info "GitHub 不可达，将使用国内镜像源安装"
fi

# ── Install target selection ────────────────────────────────────────

header "选择安装目标"

echo ""
echo "  1) 完整安装 (nvm + Node.js + Claude Code)"
echo "  2) 仅安装 nvm 和 Node.js"
echo "  3) 跳过 nvm/Node.js，仅安装 Claude Code"
echo ""

while true; do
    read -r -p "请输入选项 [1-3] (默认: 1): " choice
    choice="${choice:-1}"
    case "$choice" in
        1) INSTALL_NVM=true;  INSTALL_CLAUDE=true;  break ;;
        2) INSTALL_NVM=true;  INSTALL_CLAUDE=false; break ;;
        3) INSTALL_NVM=false; INSTALL_CLAUDE=true;  break ;;
        *) warn "无效选项，请输入 1, 2 或 3" ;;
    esac
done

info "nvm/Node.js: ${INSTALL_NVM} | Claude Code: ${INSTALL_CLAUDE}"

header "安装系统依赖 (curl, git)"

if ! command -v curl &>/dev/null; then
    warn "curl 未安装，正在安装..."
    case "$PKG_MGR" in
        brew)    brew install curl ;;
        apt-get) sudo apt-get update -qq && sudo apt-get install -y -qq curl ;;
    esac
fi

if ! command -v git &>/dev/null; then
    warn "git 未安装，正在安装..."
    case "$PKG_MGR" in
        brew)    brew install git ;;
        apt-get) sudo apt-get install -y -qq git ;;
    esac
fi

info "curl: $(curl --version 2>/dev/null | head -1 || echo 'N/A')"
info "git:  $(git --version 2>/dev/null || echo 'N/A')"

if [ "$INSTALL_NVM" = true ]; then

# ── Install nvm ─────────────────────────────────────────────────────

header "安装 nvm"

if [ -d "$NVM_DIR" ] && [ -s "$NVM_DIR/nvm.sh" ]; then
    info "nvm 已安装: $NVM_DIR"
else
    if [ "$USE_FOREIGN" = true ]; then
        info "从 GitHub 安装 nvm..."
        if curl -fsSL "https://raw.githubusercontent.com/nvm-sh/nvm/$NVM_VERSION/install.sh" | bash; then
            info "nvm 安装成功 (GitHub)"
        else
            err "nvm 安装失败，请检查网络"
            exit 1
        fi
    else
        info "从 Gitee 镜像克隆 nvm..."
        git clone --depth 1 -b "$NVM_VERSION" "$NVM_GITEE" "$NVM_DIR" 2>/dev/null
        info "nvm 克隆完成: $NVM_DIR"
    fi
fi

# ── Load nvm ────────────────────────────────────────────────────────

export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && source "$NVM_DIR/nvm.sh"

if ! command -v nvm &>/dev/null; then
    err "nvm 加载失败，请检查安装"
    exit 1
fi

info "nvm: $(nvm --version)"

# ── Configure nvm node mirror ───────────────────────────────────────

header "配置 Node.js 镜像"

if [ "$USE_FOREIGN" = false ]; then
    configure_rc() {
        local rc="$1"
        local marker="# nvm-mirror (added by install.sh)"

        [ ! -f "$rc" ] && touch "$rc"

        if grep -qF "$marker" "$rc" 2>/dev/null; then
            info "镜像已配置: $rc"
            return
        fi

        cat <<'BLOCK' >> "$rc"

# nvm-mirror (added by install.sh)
export NVM_NODEJS_ORG_MIRROR="https://npmmirror.com/mirrors/node"
BLOCK
        info "已写入镜像配置: $rc"
    }

    for rc in "$HOME/.bashrc" "$HOME/.zshrc" "$HOME/.bash_profile"; do
        configure_rc "$rc"
    done

    export NVM_NODEJS_ORG_MIRROR="$NODE_MIRROR"
    info "Node.js 下载使用 npmmirror 镜像"
else
    info "Node.js 下载使用官方源 (nodejs.org)"
fi

# ── Install Node.js LTS ─────────────────────────────────────────────

header "安装 Node.js LTS"

CURRENT_NODE="$(node --version 2>/dev/null || true)"

if [ -n "$CURRENT_NODE" ]; then
    info "Node.js 已安装: $CURRENT_NODE"
    LTS_LATEST="$(nvm version-remote --lts 2>/dev/null | tail -1 || true)"
    if [ "$CURRENT_NODE" = "$LTS_LATEST" ]; then
        info "已是最新 LTS，跳过安装"
    else
        info "当前: $CURRENT_NODE → 最新 LTS: $LTS_LATEST"
        nvm install --lts
        nvm alias default lts/*
    fi
else
    info "正在安装 Node.js LTS..."
    nvm install --lts
    nvm alias default lts/*
    nvm use default
fi

info "node: $(node --version)"
info "npm:  $(npm --version)"

# ── Set npm registry ────────────────────────────────────────────────

header "配置 npm registry"

if [ "$USE_FOREIGN" = false ]; then
    CURRENT_REGISTRY="$(npm config get registry 2>/dev/null || true)"

    if [ "$CURRENT_REGISTRY" = "$NPM_REGISTRY" ]; then
        info "npm registry 已设为国内源: $NPM_REGISTRY"
    else
        npm config set registry "$NPM_REGISTRY"
        info "npm registry 已切换: $CURRENT_REGISTRY → $NPM_REGISTRY"
    fi
else
    info "npm registry 使用官方源 (registry.npmjs.org)"
fi

fi  # INSTALL_NVM

if [ "$INSTALL_NVM" = false ] && [ "$INSTALL_CLAUDE" = true ]; then
    header "加载已有 nvm"
    [ -s "$NVM_DIR/nvm.sh" ] && source "$NVM_DIR/nvm.sh"
    if ! command -v nvm &>/dev/null; then
        err "未检测到 nvm，请先安装 nvm 和 Node.js (选择选项 1 或 2)"
        exit 1
    fi
    info "已加载已有 nvm: $(nvm --version)"
    info "node: $(node --version 2>/dev/null || echo '未检测到')"
fi

# ── Install Claude Code CLI ─────────────────────────────────────────

if [ "$INSTALL_CLAUDE" = true ]; then

header "安装 Claude Code CLI"

if command -v claude &>/dev/null; then
    info "Claude Code CLI 已安装: $(claude --version 2>/dev/null || echo 'ok')"
    info "尝试更新到最新版本..."
    npm install -g "$CLAUDE_PKG" 2>&1 | tail -1 || warn "更新失败，当前版本仍可用"
else
    info "正在安装 Claude Code CLI..."
    npm install -g "$CLAUDE_PKG"
fi

fi  # INSTALL_CLAUDE

# ── Verify ──────────────────────────────────────────────────────────

header "安装结果"

echo ""
if [ "$INSTALL_NVM" = true ]; then
    echo -e "  ${GREEN}nvm${NC}    $(nvm --version 2>/dev/null || echo 'N/A')"
    echo -e "  ${GREEN}node${NC}   $(node --version 2>/dev/null || echo 'N/A')"
    echo -e "  ${GREEN}npm${NC}    $(npm --version 2>/dev/null || echo 'N/A')    registry: $(npm config get registry 2>/dev/null || echo 'N/A')"
fi
if [ "$INSTALL_CLAUDE" = true ]; then
    echo -e "  ${GREEN}claude${NC} $(claude --version 2>/dev/null || echo 'N/A')"
fi
echo ""

info "全部安装完成！新开终端后生效，或执行: source ~/.bashrc (或 ~/.zshrc)"
