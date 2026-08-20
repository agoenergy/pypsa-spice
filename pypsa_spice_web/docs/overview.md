# PyPSA-SPICE model workspace

The model workspace is a React and FastAPI interface for a local PyPSA-SPICE project. Use it to edit inputs grouped by technology, configure scenarios, run Snakemake, and inspect results. The checked-out CSV and YAML files remain the source of truth.

Results cover Power, Industry, Transport, Emissions, and Costs. Comparison mode shows two result runs beside each other for every indicator. The difference view calculates `comparison − primary` for each aligned series and timestamp or model year.

The [web application code structure and logic guide](code-structure.html) maps the modules, request flow, stored data, and safety checks.

The [UI and UX backlog](ui-review-backlog.md) records open interface findings with file references and a suggested order of work.

## Do not copy the Streamlit UI

The web interface is an independent React, TypeScript, CSS, and FastAPI application. When developing `pypsa-spice-vis-ui` or the code in `pypsa_spice_web/`, do not use the legacy Streamlit application in `pypsa-spice-vis/` as UI or UX guidance. Its layouts, controls, component choices, styles, and interaction patterns are not a design reference for the web interface.

Base web UI decisions on the current React frontend, this documentation, and explicit product requirements. Consult the Streamlit implementation only when you need to understand data semantics, calculations, or feature behaviour. Reimplement that behaviour for the web instead of copying the Streamlit interface.

## Local development

### 1. Install dependencies

From the repository root:

```bash
python -m pip install -r requirements-web.txt
npm install --prefix pypsa_spice_web/frontend
```

This installs Python and Node dependencies for initial set-up. Rerun only if dependencies change, or if you switch/create Python or Node environments.

### 2. Run app locally

For convenience, a `run-web-locally.sh` script can be used to run the react app locally. 

From the repository root:

```bash
./run-web-locally.sh
```

This starts FastAPI at http://127.0.0.1:8000, using the active Python environment if it contains the web dependencies, otherwise, it uses the existing `hotpot` Conda environment when available.

The app is then served at Vite's default dev-server port. Open the app in your browser at http://127.0.0.1:5713/ui/.

## Result layout

The app expects this structure:

```text
data/<dataset>/<project>/results/<scenario>/csvs/<sector>/<year>/<table>.csv
```

Yearly visualisations combine the tables in each modelled-year folder. Hourly visualisations use the selected year. The server samples large hourly results for display. CSV downloads always return the full source file.

## Workspace structure

The workspace keeps these selections separate because they change at different times:

- **Project** selects the local model project. It usually changes less often than the other selections.
- **Scenario** selects the editable input/configuration scenario.
- **Result run** selects an available result folder. It is separate from the input scenario and does not have to share its name.
- **Country** belongs to each Results chart, so charts can show different countries side by side. Inputs shows it only for tables with country-specific rows. Configuration uses one shared country selection for country-specific sections.
- **Compare with** is a persistent Results control in the top workspace bar.

The main workspace pages are:

- **Inputs** groups data by technology and includes a table and file view for advanced editing.
- **Configure & run** covers scenario settings, CO₂ management, custom constraints, and the final review and run step.
- **Results** has one analysis page per energy-system section: Power, Industry, Transport, Emissions, and Costs. Each page plots that section and can compare two runs.

## Custom dashboards

The **Dashboards** page collects existing Results charts in the current browser.
Each dashboard belongs to one result dataset and project. Its workspace is an
ordered list of rows. Each row contains either an editable heading or one
full-width chart. You can reorder or remove rows and replace a chart through its
settings.

A chart card can show one or two result runs from its project. With two runs, it
can instead calculate `comparison − reference`. Normal mode puts the two plots
side by side in the chart row. Difference mode uses the full row for one plot.
Each chart keeps its own country, filters, hourly year, and time range. On
import, the app keeps the first two unique scenarios if a chart lists more than
two. It also converts old two-chart rows into two full-width rows.

The browser saves dashboard definitions as you edit them. They contain chart IDs,
result-run selections, filters, titles, and row order. They do not contain
result data. The same browser profile and application address restore the last
open dashboard. Clearing site storage deletes these definitions.

Schema 2 stores each chart configuration inside its `chart` row, alongside its
position in the dashboard. When the app reads schema-1 browser data or an old
export, it converts the definition to schema 2. All later saves and exports use
schema 2.

Use **Export** to download a versioned `pypsa-spice-dashboard` JSON file.
**Import** validates the file, shows a summary, and creates a new local
dashboard without overwriting an existing one. An import keeps references to
charts or scenarios it cannot find. Open the matching result checkout to
repair them. An exported configuration contains no result data and does not
publish the dashboard.

Chart components reach saved dashboards through a storage interface, and the
current implementation is `localStorage`. Moving to a server-backed store later
means writing a second implementation of that interface. The versioned JSON
format and the Results chart registry stay as they are.

Result files remain read-only. The app writes input and scenario configuration changes to the selected local CSV or YAML file only when the user saves. The browser also keeps large time-series input tables read-only.

## Configure and run workflow

The Configure & run page uses the repository-root `base_config.yaml` as its fixed template. The shared workspace bar supplies the project and input scenario. To create a scenario, the app atomically copies a complete existing scenario folder. It writes the copy to the normal model path and adds no web-only metadata.

Review & run shows the scenario year and resolution, model regions, sectors, currency, output folder, core count, and workflow target. If the dataset has its own worktree, it also shows the Git branch, commit, and local changes. The app never commits, pulls, or pushes data. While a web-launched run is active, the backend rejects web edits and new scenarios. Command-line access remains unchanged.

Starting a run creates `logs/web_runs/<run-id>/base_config.yaml` and replaces only these values in that run-local copy:

- `path_configs.data_folder_name`
- `path_configs.project_name`
- `path_configs.input_scenario_name`
- `path_configs.output_scenario_name`

The app does not change the checked-in `base_config.yaml`. Each run also writes `run_manifest.json` with the selection, target, hotpot environment, data repository state, and SHA-256 hashes for the selected scenario and shared global inputs. FastAPI then runs the root `Snakefile` with the `solve_all_networks` target. Only one web-launched run can be active. Review & run shows the current rule, logo-chilli progress, final state, and the tail of the Snakemake log. Complete logs and saved state remain under `logs/web_runs/<run-id>/`.

Model runs always use the `hotpot` Conda environment. If the web server already runs in `hotpot`, the runner uses its Snakemake installation. Otherwise, it launches `conda run --no-capture-output -n hotpot snakemake`. Set `SNAKEMAKE_COMMAND` to override this command locally.

```bash
SNAKEMAKE_COMMAND="conda run -n hotpot snakemake" ./run_web.sh
```

## Result-file access

Chart and download requests cannot supply arbitrary CSV paths. FastAPI:

- accepts only table names present in the configured chart definitions
- requires the requested hourly/yearly mode to match the chart definition
- accepts an hourly year only when its numeric year directory exists
- resolves every result path and confirms it stays inside the selected sector directory

These checks keep a request inside the result folder. They are not authentication or hosted access control. A malformed URL cannot reach a CSV elsewhere in the model checkout.

## Input filtering and pagination

FastAPI filters input tables before pagination. The table endpoint accepts technology, technology carrier, country, configured filter value, free-text query, offset, and limit parameters. It returns:

- the current page of typed rows
- original CSV row IDs, which stay stable when the user saves filtered rows
- total source and matching-row counts
- country and configured-filter options from the complete matching dataset
- offset, limit, and truncation metadata

The React table editor requests 100 rows per page. FastAPI matches technologies against the complete source table instead of limiting the search to rows already loaded in the browser. The editor keeps unsaved changes while the user moves between pages. After the user confirms a project or scenario change, the remounted editor clears them. A table loading error stays visible instead of appearing as an empty technology match.

Each React effect owns the `AbortController` for its input or scenario configuration request. Changing context or unmounting the editor aborts the active request. An older response therefore cannot overwrite a newly selected table, project, or scenario.

## Implementation record

### 2026-08-20: type scale

- Added seven size tokens to `App.css` with a 12px floor and mapped all 158 size declarations across 13 stylesheets onto them.
- Dropped the page title from `clamp(34px, 4vw, 52px)` to 28px and retuned its leading and tracking for the smaller size.
- Added a `chartFont` constant in `Plot.tsx` mirroring the tokens, because Plotly draws its own SVG text. Axis ticks moved from 9px to 12px and plot margins grew to fit the larger labels.
- Widened the chart legend columns, the Home project picker and the count badges, which had started truncating.
- Aligned the sidebar brand border with the workspace bar border, previously 1px out.
- Recorded the open findings from the same review in the [UI and UX backlog](ui-review-backlog.md).

Verification commands:

```bash
npm --prefix pypsa_spice_web/frontend run build
npm --prefix pypsa_spice_web/frontend test
```

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
git diff --check -- pypsa_spice_web tests
```

The app reads the checked-out `data/` tree, including deployments that sit beside a model checkout. Uploads, object storage, authentication, and hosted persistence are out of scope.
