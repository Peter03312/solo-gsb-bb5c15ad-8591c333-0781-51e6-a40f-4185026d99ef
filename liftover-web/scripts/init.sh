#!/usr/bin/env bash
# 依赖初始化：在宿主机上创建虚拟环境并安装依赖（含开发/测试依赖）。
# 仅用于本地开发；Docker 部署不需要执行本脚本。
set -euo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON:-python3}"

if [ ! -d ".venv" ]; then
  echo "[init] 创建虚拟环境 .venv"
  "$PYTHON_BIN" -m venv .venv
else
  echo "[init] 虚拟环境 .venv 已存在，跳过创建"
fi

# shellcheck disable=SC1091
source .venv/bin/activate

echo "[init] 升级 pip 并安装依赖（含测试依赖）"
pip install --upgrade pip
pip install -r requirements-dev.txt

echo "[init] 完成。启动开发服务器："
echo "       .venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8080"
