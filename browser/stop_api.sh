#!/usr/bin/env bash
# Stop the baseline genome browser API on PORT (default 8765)
set -euo pipefail
PORT="${PORT:-8765}"
pids="$(pgrep -f "uvicorn browser.api.main:app.*--port ${PORT}" 2>/dev/null || true)"
if [ -z "${pids}" ]; then
  if command -v lsof >/dev/null 2>&1; then
    pids="$(lsof -ti :"${PORT}" 2>/dev/null || true)"
  fi
fi
if [ -z "${pids}" ]; then
  echo "[browser] no process listening on port ${PORT}"
  exit 0
fi
echo "[browser] stopping: ${pids}"
# shellcheck disable=SC2086
kill ${pids} 2>/dev/null || true
sleep 1
# shellcheck disable=SC2086
kill -9 ${pids} 2>/dev/null || true
echo "[browser] stopped"
