# PyPSA-SPICE model workspace

The model workspace is the React and FastAPI interface for working with a local PyPSA-SPICE project. It provides technology-first input editing, scenario configuration, Snakemake-backed model runs, and long-form result analysis while continuing to use the checked-out CSV and YAML files as the source of truth. The Results area covers Power, Industry, Transport, Emissions, and Costs. Comparison mode places both result runs side by side for each indicator; the optional difference view calculates `comparison − primary` for every aligned series and timestamp or model year.

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
- **Country** is chart-specific in Results so figures can be inspected side by side for different countries. In Inputs it appears only on tables with country-specific rows; Configuration keeps the shared country context used by country-aware sections.
- **Compare with** is a persistent Results control in the top workspace bar.

The main workspace pages are:

- **Inputs**, organised primarily by technology with a table/file view for advanced editing.
- **Configure & run**, covering scenario settings, CO₂ management, custom constraints, and a final review-and-run stage.
- **Results**, providing one long analytical page for each energy-system section.

The Results pages provide visualisation and scenario comparison for:

- Power
- Industry
- Transport
- Emissions
- Costs

Result files remain read-only. Input and scenario-configuration changes are written directly to the selected local CSV or YAML file only after the user explicitly saves them. Large timeseries input tables remain read-only in the browser.

## Configure and run workflow

The Configure & run page uses the repository-root `base_config.yaml` as its fixed configuration template. Project and input-scenario context come from the shared workspace bar. A new scenario can be created locally by atomically duplicating a complete existing scenario folder; it uses the normal model-visible path and adds no web-only metadata file.

The Review & run stage shows scenario year and resolution, model regions, sectors and currency, the output folder, core count, workflow target, and the dataset Git branch, commit, and local-change status when the dataset is its own worktree. It never commits, pulls, or pushes data. While a web-launched run is active, the backend rejects web edits and scenario creation; command-line access remains unchanged.

Starting a run creates `logs/web_runs/<run-id>/base_config.yaml` and replaces only these values in that run-local copy:

- `path_configs.data_folder_name`
- `path_configs.project_name`
- `path_configs.input_scenario_name`
- `path_configs.output_scenario_name`

The checked-in root `base_config.yaml` is not mutated. Each run also writes `run_manifest.json` with the selection, target, hotpot environment, data-repository state, and SHA-256 hashes of the selected scenario and shared global inputs. FastAPI then starts the root `Snakefile` with the `solve_all_networks` target. Only one web-launched run can be active at a time. Current rule, logo-chilli progress, final state, and the Snakemake log tail are available in Review & run; complete logs and persisted state remain under `logs/web_runs/<run-id>/`.

Model execution always uses the `hotpot` Conda environment. When the web server is already running in `hotpot`, the runner uses that environment's Snakemake installation. Otherwise it launches `conda run --no-capture-output -n hotpot snakemake`. `SNAKEMAKE_COMMAND` remains available for an explicit local override.

```bash
SNAKEMAKE_COMMAND="conda run -n hotpot snakemake" ./run_web.sh
```

## Result-file access

Result chart and download requests do not accept arbitrary CSV paths. The FastAPI layer:

- accepts only table names present in the configured chart definitions;
- requires the requested hourly/yearly mode to match the chart definition;
- accepts hourly years only when the corresponding numeric year directory exists; and
- resolves every result path and confirms that it remains inside the selected sector directory.

This is a filesystem-boundary check for the local application, not an authentication or hosted-access system. It prevents malformed URL parameters from selecting CSV files elsewhere in the model checkout.

## Input filtering and pagination

Input table filtering is performed by FastAPI before pagination. The table endpoint supports technology, technology carrier, country, configured filter value, free-text query, offset, and limit parameters. Responses include:

- the current page of typed rows;
- original CSV row IDs, which remain stable when filtered rows are saved;
- total source and matching-row counts;
- country and configured-filter options derived from the complete matching dataset; and
- offset, limit, and truncation metadata.

The React table editor requests 100 rows per page. Technology matching therefore covers the complete source table rather than only a browser-loaded prefix. Unsaved edits are retained while moving between pages and are cleared when a confirmed project or scenario change remounts the editor. Per-table loading errors remain visible instead of being treated as empty technology matches.

Input and scenario-configuration requests use an `AbortController` owned by the corresponding React effect. Changing context or unmounting an editor aborts its active request so an older response cannot overwrite the newly selected table, project, or scenario.

## Implementation record

### 2026-07-15: result paths and input-editor data flow

- Confined result chart/download paths to configured tables, valid modes, existing years, and the selected sector directory.
- Corrected React request cleanup in the CSV and scenario-configuration editors.
- Moved technology, country, configured-value, and text filtering to the backend before pagination.
- Replaced the misleading fixed-row-limit message with matching-row and page counts.
- Added an input-table API version marker so a stale backend produces an explicit restart message instead of breaking the Inputs page.
- Added regression coverage for result-path confinement and filtering-before-pagination.

Verification commands:

```bash
npm --prefix pypsa_spice_web/frontend run build
conda run -n hotpot python -m unittest discover -s tests -p 'test*.py' -q
git diff --check -- pypsa_spice_web tests docs/visualisation-tool/pypsa-spice-web.md
```

The current data-management strategy is local discovery of the checked-out `data/` tree. This local-file workflow remains the primary operating mode even when the interface is deployed alongside a model checkout. Uploads, object storage, authentication, and hosted data persistence are separate future concerns.
