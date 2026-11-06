<!--
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# PyPSA-SPICE model workspace

> **Beta version.** The web workspace is under active development and is not a fully released product. Features, layouts, and the API may change, and some workflows are still incomplete. Mark user-facing text, the app version (`0.2.0-beta`), and new documentation accordingly until the first stable release. User documentation lives in [`docs/visualisation-tool/web-workspace.md`](../../docs/visualisation-tool/web-workspace.md).

The model workspace is a React and FastAPI interface for a local PyPSA-SPICE project. Use it to edit inputs grouped by technology, configure scenarios, run Snakemake, and inspect results. The checked-out CSV and YAML files remain the source of truth.

Results cover Power, Industry, Transport, Emissions, and Costs. Comparison mode shows two result runs beside each other for every indicator. The difference view calculates `comparison − primary` for each aligned series and timestamp or model year.

The [web application code structure and logic guide](code-structure.html) maps the modules, request flow, stored data, and safety checks.

## Do not copy the Streamlit UI

The web interface is an independent React, TypeScript, SCSS, and FastAPI application. When developing `pypsa-spice-vis-ui` or the code in `pypsa-spice-web/`, do not use the legacy Streamlit application in `pypsa-spice-vis/` as UI or UX guidance. Its layouts, controls, component choices, styles, and interaction patterns are not a design reference for the web interface.

Base web UI decisions on the current React frontend, this documentation, and explicit product requirements. Consult the Streamlit implementation only when you need to understand data semantics, calculations, or feature behaviour. Reimplement that behaviour for the web instead of copying the Streamlit interface.

## Local development

### 1. Install dependencies

Install [Node.js](https://nodejs.org/) with `npm`, then from the repository root:

```bash
conda env update -n hotpot -f envs/environment.yaml --prune
```

Rerun it when the Python dependencies change. `app.py` installs the frontend packages itself when
`frontend/node_modules/` is missing. When the frontend dependencies change, `app.py` will update them automatically.

### 2. Run app locally

From the repository root, first activate the `hotpot` environment and then use one of the following commands:

```bash
python pypsa-spice-web/app.py        # use the app
python pypsa-spice-web/app.py --dev  # develop the app
```

Without `--dev`, FastAPI serves the API and the built frontend at <http://127.0.0.1:8000/ui/>.
Before starting, `app.py` builds the frontend into `frontend/dist/` when a file in `frontend/` is newer than the last build.

With `--dev`, Vite serves the frontend at <http://127.0.0.1:5173/ui/> and reloads the page on every frontend change, and FastAPI restarts when backend code changes. Vite forwards `/api` and `/vendor`
requests to FastAPI on port 8000.

`Ctrl+C` stops FastAPI and, with `--dev`, Vite. If port 8000 is already in use, for example by another copy of the app, FastAPI stops with "address already in use". The host and ports are the `HOST`, `PORT`, and `VITE_PORT` constants in `app.py`; the Vite proxy in
`frontend/vite.config.ts` repeats port 8000.

### 3. Format and check frontend changes

From the repository root:

```bash
npm run format --prefix pypsa-spice-web/frontend
npm run format:check --prefix pypsa-spice-web/frontend
npm test --prefix pypsa-spice-web/frontend
npm run build --prefix pypsa-spice-web/frontend
```

Prettier uses a 120-character print width and formats frontend source and configuration
files. It excludes generated files, public assets, dependencies, and the lockfile.
The Web frontend GitHub Actions workflow runs the formatting check, tests, and
production build for frontend changes on pull requests and pushes to main/develop.

Shared controls live in `frontend/src/components/`: `Button.tsx` owns text actions,
`IconButton.tsx` owns icon buttons and native icon links, and `FormControls.tsx`
owns field wrappers, labelled selects, search fields, and toggles. Their base
styles live in collocated SCSS modules. Icon controls require an `aria-label`;
submit buttons specify `type="submit"`. Existing page styles use scoped
`data-control` selectors for layout and `--icon-control-*` variables for compact
icon sizes, without depending on generated module class names.

All frontend styles use SCSS. Component and page classes belong in `.module.scss`
files; `global.scss` contains only fonts, design tokens, resets, focus treatment,
and reduced-motion rules. `App.module.scss` owns the shell. Shared modules such
as `ChartSurface`, `DataTable`, `DialogSurface`, and `EditorPanel` own patterns
used by multiple views. Import their classes in JSX, or use explicit `@value`
imports when a module needs a shared class in a selector. Vite's scoped-name
configuration keeps those imports consistent. `styles.test.ts` checks module
boundaries, JSX class references, and shared class resolution.

The app shares one run-status monitor between its navigation and run page. Active
runs refresh one second after each successful response. A status request times out
after 15 seconds; connection failures retain the last known status and retry with
increasing delays up to 30 seconds. Finishing a run stops polling. Idle workspaces
check once on startup; reload the app to discover runs started in another session.

## Chart code

Chart responsibilities are split across these files under `frontend/src/`:

| File | Responsibility |
| --- | --- |
| `components/ChartCard.tsx` and `components/DashboardChartCard.tsx` | Request result rows and own chart filters, comparison controls, and hidden legend values. |
| `components/Plot.tsx` | Load Plotly, render typed traces with the chart layout, resize expanded plots, and clean up the plot. |
| `components/ChartLegend.tsx` | Render legend buttons and report toggles to the owning card. Its styles live in `ChartLegend.module.scss`. |
| `shared/chartData.ts` | Aggregate rows, derive legend values, align scenario series, and build difference rows for the data view. |
| `shared/plotTraces.ts` | Build series, comparison, difference, and stacked-bar total traces through `buildPlotTraces`. Return `hasColumnTotals` so Plot can reserve label space. |
| `shared/chartPresentation.ts` | Provide mapped labels and fallback colours for legends and traces, plus `chartFont` for Plotly. |
| `shared/types.ts` and `plotly.ts` | Define the Plotly trace contract and load the runtime script from `/vendor/plotly.min.js`. |

Difference calculations align each legend value by model year or hourly snapshot
and calculate `comparison − primary`, treating a missing value as zero. Column
totals exclude hidden and secondary-axis series. Legend labels and colours use
the same helpers in Results & analysis and the Dashboard builder.

Trace builders use the development dependency `@types/plotly.js`. The shared
`ScatterTextTrace` type extends those definitions for per-point text positions,
which the runtime supports but version 3.0.13 of the definitions omits. FastAPI
still serves the Plotly runtime from the Python installation.

The calculation tests live in `shared/chartData.test.ts`, and trace behaviour
is covered by `shared/plotTraces.test.ts`. Run both from the repository root:

```bash
npm --prefix pypsa-spice-web/frontend test -- src/shared/chartData.test.ts src/shared/plotTraces.test.ts
npm --prefix pypsa-spice-web/frontend run check
```

## Result layout

The app expects this structure:

```text
data/<dataset>/<project>/results/<scenario>/csvs/<sector>/<year>/<table>.csv
```

Yearly visualisations read the combined table in `csvs/<sector>/all_years/<table>.csv`, written by the `combine_summaries` rule. If that file is missing, they combine the tables in each modelled-year folder instead. Hourly visualisations use the selected year. The server samples large hourly results for display. CSV downloads return the source file unchanged; when the table comes from several year folders, the download stacks those files as written.

## Workspace structure

The workspace keeps these selections separate because they change at different times:

- **Project** selects the local model project. It usually changes less often than the other selections.
- **Scenario** selects the editable input/configuration scenario.
- **Result run** selects an available result folder. It is separate from the input scenario and does not have to share its name.
- **Country** belongs to each Results chart, so charts can show different countries side by side. Inputs shows it only for tables with country-specific rows. Configuration uses one shared country selection for country-specific sections.
- **Compare with** is a persistent Results control in the top workspace bar.

The sidebar follows the modelling workflow:

- **Home** shows available projects, saved dashboards and the workflow steps.
- **Inputs & configuration** groups the pages that prepare a scenario:
  - **Input data** groups data by technology and includes a table and file view for advanced editing.
  - **Compare inputs** shows configuration and input-data differences between two scenarios.
  - **Scenario configuration** covers general settings, CO₂ management and custom constraints.
- **Review & run** checks the model scope and starts a Snakemake run.
- **Results & analysis** has one analysis page per energy-system section: Power, Industry, Transport, Emissions, and Costs. Each page plots that section and can compare two runs.
- **Dashboard builder** combines Results charts into custom dashboards.

## Custom dashboards

The **Dashboard builder** page collects existing Results charts in the current browser.
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

The Scenario configuration and Review & run pages use the repository-root `base_config.yaml` as its fixed template. The shared workspace bar supplies the project and input scenario. To create a scenario, the app atomically copies a complete existing scenario folder. It writes the copy to the normal model path and adds no web-only metadata.

Review & run shows the scenario year and resolution, model regions, sectors, currency, output folder, core count, and workflow target. If the dataset has its own worktree, it also shows the Git branch, commit, and local changes. The app never commits, pulls, or pushes data. While a web-launched run is active, the backend rejects web edits and new scenarios. Command-line access remains unchanged.

Starting a run creates `logs/web_runs/<run-id>/base_config.yaml` and replaces only these values in that run-local copy:

- `path_configs.data_folder_name`
- `path_configs.project_name`
- `path_configs.input_scenario_name`
- `path_configs.output_scenario_name`

The app does not change the checked-in `base_config.yaml`. Each run also writes `run_manifest.json` with the selection, target, hotpot environment, data repository state, and SHA-256 hashes for the selected scenario and shared global inputs. FastAPI then runs the root `Snakefile` with the `solve_all_networks` target. Only one web-launched run can be active. Review & run shows the current rule, logo-chilli progress, final state, and the tail of the Snakemake log. Complete logs and saved state remain under `logs/web_runs/<run-id>/`.

Model runs always use the `hotpot` Conda environment. If the web server already runs in `hotpot`, the runner uses its Snakemake installation. Otherwise, it launches `conda run --no-capture-output -n hotpot snakemake`. Set `SNAKEMAKE_COMMAND` to override this command locally.

```bash
SNAKEMAKE_COMMAND="conda run -n hotpot snakemake" python pypsa-spice-web/app.py
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

The app reads the checked-out `data/` tree, including deployments that sit beside a model checkout. Uploads, object storage, authentication, and hosted persistence are out of scope.
