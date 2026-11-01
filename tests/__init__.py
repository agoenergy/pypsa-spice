"""Tests for the PyPSA-SPICE web workspace."""

import sys
from pathlib import Path

# The web backend lives in "pypsa-spice-web/", which is not an importable package name,
# so its modules (app, input_editor, ...) are imported from that folder directly.
WEB_DIR = Path(__file__).resolve().parents[1] / "pypsa-spice-web"
if str(WEB_DIR) not in sys.path:
    sys.path.insert(0, str(WEB_DIR))
