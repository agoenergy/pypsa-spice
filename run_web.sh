#!/usr/bin/env bash
set -euo pipefail

HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-8000}"

if ! command -v npm >/dev/null 2>&1; then
  echo "Node.js and npm are required to build the React interface." >&2
  exit 1
fi

if [ ! -d pypsa_spice_web/frontend/node_modules/react ] || [ ! -d pypsa_spice_web/frontend/node_modules/vite ]; then
  echo "Frontend dependencies are missing. Run: npm install --prefix pypsa_spice_web/frontend" >&2
  exit 1
fi

npm --prefix pypsa_spice_web/frontend run build

if python -c "import fastapi, uvicorn, jinja2, plotly, yaml, ruamel.yaml" >/dev/null 2>&1; then
  exec python -m uvicorn pypsa_spice_web.app:app --host "$HOST" --port "$PORT"
fi

if command -v conda >/dev/null 2>&1 && conda run -n hotpot python -c "import fastapi, uvicorn, jinja2, plotly, yaml, ruamel.yaml" >/dev/null 2>&1; then
  exec conda run --no-capture-output -n hotpot python -m uvicorn pypsa_spice_web.app:app --host "$HOST" --port "$PORT"
fi

echo "Web dependencies are missing. Install them with:" >&2
echo "  python -m pip install -r requirements-web.txt" >&2
exit 1
