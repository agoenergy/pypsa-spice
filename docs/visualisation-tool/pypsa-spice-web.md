# PyPSA-SPICE model workspace

The model workspace is the React and FastAPI interface for working with a local PyPSA-SPICE project. It provides technology-first input editing, scenario configuration, and long-form result analysis while continuing to use the checked-out CSV and YAML files as the source of truth. The Results area covers Power, Industry, Transport, Emissions, and Costs. Comparison mode places both result runs side by side for each indicator; the optional difference view calculates `comparison − primary` for every aligned series and timestamp or model year.

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

## Workspace structure

The persistent workspace context separates concepts that have different lifecycles:

- **Project** selects the local model project and normally changes infrequently.
- **Scenario** selects the editable input/configuration scenario.
- **Result run** selects an available result folder and is intentionally distinct from the input scenario.
- **Country** is shared across Results, Inputs, and country-aware configuration sections.

The main workspace pages are:

- **Inputs**, organised primarily by technology with a table/file view for advanced editing.
- **Scenario configuration**, covering scenario settings, CO₂ management, and custom constraints.
- **Results**, providing one long analytical page for each energy-system section.

The Results pages provide visualisation and scenario comparison for:

- Power
- Industry
- Transport
- Emissions
- Costs

Result files remain read-only. Input and scenario-configuration changes are written directly to the selected local CSV or YAML file only after the user explicitly saves them. Large timeseries input tables remain read-only in the browser.

The current data-management strategy is local discovery of the checked-out `data/` tree. This local-file workflow remains the primary operating mode even when the interface is deployed alongside a model checkout. Uploads, object storage, authentication, and hosted data persistence are separate future concerns.
