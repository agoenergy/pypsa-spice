#!/usr/bin/env bash
set -euo pipefail

BACKEND_HOST="${BACKEND_HOST:-127.0.0.1}"
BACKEND_PORT="${BACKEND_PORT:-8000}"
FRONTEND_HOST="${FRONTEND_HOST:-127.0.0.1}"
FRONTEND_PORT="${FRONTEND_PORT:-5173}"
FRONTEND_DIR="pypsa_spice_web/frontend"

backend_pid=""

cleanup() {
  if [ -n "$backend_pid" ] && kill -0 "$backend_pid" >/dev/null 2>&1; then
    kill "$backend_pid" >/dev/null 2>&1 || true
    wait "$backend_pid" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT INT TERM

if ! command -v npm >/dev/null 2>&1; then
  echo "Node.js and npm are required to run the React development server." >&2
  exit 1
fi

if [ ! -d "$FRONTEND_DIR/node_modules/react" ] || [ ! -d "$FRONTEND_DIR/node_modules/vite" ]; then
  echo "Frontend dependencies are missing. Run once after checkout or dependency changes:" >&2
  echo "  npm install --prefix $FRONTEND_DIR" >&2
  exit 1
fi

python_command=()
if python -c "import fastapi, uvicorn, plotly, yaml, ruamel.yaml" >/dev/null 2>&1; then
  python_command=(python)
elif command -v conda >/dev/null 2>&1 && conda run -n hotpot python -c "import fastapi, uvicorn, plotly, yaml, ruamel.yaml" >/dev/null 2>&1; then
  python_command=(conda run --no-capture-output -n hotpot python)
else
  echo "Web Python dependencies are missing. Run once after checkout or dependency changes:" >&2
  echo "  python -m pip install -r requirements-web.txt" >&2
  exit 1
fi

echo "Starting FastAPI backend on http://$BACKEND_HOST:$BACKEND_PORT"
"${python_command[@]}" -m uvicorn pypsa_spice_web.app:app --host "$BACKEND_HOST" --port "$BACKEND_PORT" &
backend_pid="$!"

echo "Starting React development server on http://$FRONTEND_HOST:$FRONTEND_PORT/ui/"
echo "Open http://$FRONTEND_HOST:$FRONTEND_PORT/ui/ for the local React app."
echo "FastAPI remains available at http://$BACKEND_HOST:$BACKEND_PORT for API requests."

npm --prefix "$FRONTEND_DIR" run dev -- --host "$FRONTEND_HOST" --port "$FRONTEND_PORT"
