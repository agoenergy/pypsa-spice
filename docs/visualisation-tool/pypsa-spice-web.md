# PyPSA-SPICE output explorer

The output explorer is the web-based replacement for the output area of the Streamlit visualisation interface. It reads the existing result CSV structure directly and does not modify model outputs. The current release covers Power, Industry, Transport, Emissions, and Costs. Comparison mode places both scenarios side by side for each indicator; the optional difference view calculates `comparison − primary` for every aligned series and timestamp or model year.

## Development boundary

The web interface is an independent React, TypeScript, CSS, and FastAPI application. When developing `pypsa-spice-vis-ui` or the code in `pypsa_spice_web/`, do not use the legacy Streamlit application in `pypsa-spice-vis/` as UI or UX guidance. Its layouts, controls, component choices, styles, and interaction patterns are not a design reference for the web interface.

Use the existing React frontend, this web documentation, and explicit product requirements as the sources of truth for web UI decisions. The Streamlit implementation may be consulted only when necessary to understand data semantics, calculations, or feature behaviour; implement those requirements with web-native UI patterns rather than copying the Streamlit interface.

## Run locally

From the repository root:

```bash
python -m pip install -r requirements-web.txt
npm install --prefix pypsa_spice_web/frontend
./run_web.sh
```

Open <http://127.0.0.1:8000>. The launcher uses the active Python environment when the web dependencies are installed, and falls back to the existing `hotpot` Conda environment when available.

Set `HOST` or `PORT` to override the defaults:

```bash
HOST=0.0.0.0 PORT=8080 ./run_web.sh
```

## Result layout

The app discovers data in the established structure:

```text
data/<dataset>/<project>/results/<scenario>/csvs/<sector>/<year>/<table>.csv
```

Yearly visualisations combine the tables found across modelled-year folders. Hourly visualisations use the selected year. Large hourly results are sampled for display while CSV downloads always return the full source file.

## Current scope

This release provides read-only output visualisation and scenario comparison for:

- Power
- Industry
- Transport
- Emissions
- Costs

Model configuration, input editing, and model execution remain intentionally unavailable. The output explorer never writes to result or input files.

The current data-management strategy is local discovery of the checked-out `data/` tree. Uploads, object storage, authentication, and hosted data persistence are deliberately deferred until a deployment target is selected.
