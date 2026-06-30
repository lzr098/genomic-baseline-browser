#!/usr/bin/env bash
# Start the baseline genome browser API (default http://127.0.0.1:8765)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PORT="${PORT:-8765}"
cd "$ROOT"

if ! python3 -c "import fastapi, uvicorn" 2>/dev/null; then
  echo "[browser] installing dependencies from browser/requirements.txt ..."
  python3 -m pip install -r browser/requirements.txt
fi

stop_port() {
  local pids=""
  if command -v lsof >/dev/null 2>&1; then
    pids="$(lsof -ti :"${PORT}" 2>/dev/null || true)"
  elif command -v fuser >/dev/null 2>&1; then
    pids="$(fuser "${PORT}/tcp" 2>/dev/null | tr -s ' ' || true)"
  fi
  if [ -n "${pids}" ]; then
    echo "[browser] stopping existing listener(s) on port ${PORT}: ${pids}"
    # shellcheck disable=SC2086
    kill ${pids} 2>/dev/null || true
    sleep 1
    # shellcheck disable=SC2086
    kill -9 ${pids} 2>/dev/null || true
  fi
}

stop_port

echo "[browser] root: ${ROOT}"
echo "[browser] open: http://127.0.0.1:${PORT}/"
exec python3 -m uvicorn browser.api.main:app --host 127.0.0.1 --port "${PORT}" --reload
