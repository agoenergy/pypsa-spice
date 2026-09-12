
### Suggestions batch 2

Implementation update, 2026-09-11: addressed automatic formatting under General 4
and reliable polling under RunModel 1.

Implementation update, 2026-09-12: General 1, 2, and 3 are fully addressed.
All frontend styles now use SCSS, and all component, page, and shared-pattern
classes live in SCSS modules. Global styles contain only fonts, design tokens,
resets, focus treatment, and reduced-motion rules. DataDialog 2 is also fixed.
Commit references below distinguish implemented fixes from partial progress and open suggestions.

#### General

1. Uses a mix of scss and css, probably because it implemented SCSS from the last batch of comments. Mixing is not inherently problematic, but for the sake of cleanliness i'd suggest just using one or the other

    > Reply: Fully addressed on 2026-09-12. Converted all 13 remaining plain CSS files to SCSS modules, including charts, dialogs, dashboards, input editors, scenario configuration, comparison, and run controls. There are no plain CSS files left in frontend/src. The only non-module stylesheet is global.scss, which contains application-wide foundations only.
    >
    > Commit: [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126).

2. Bigger issue: mixes global css and modules (e.g. WorkspaceCard.module.scss): non-module css files have global classes, which can lead to conflicts if two component css files have the same selectors defined. Suggest to scope all component and page css (or scss) files as modules

    > Reply: Fully addressed on 2026-09-12. Every component and page class is now scoped through .module.scss imports. App.module.scss owns the application shell; shared chart frames, tables, dialogs, editor panels, navigation bars, notices, and configuration layouts have explicit shared modules. JSX uses module exports, and cross-module selectors import their shared classes with @value. Vite generates stable scoped names so direct and shared imports resolve identically. Removed the sidebar's global overlay selector, and scoped result-table menu queries to a dialog ref. global.scss contains no component classes or unscoped table rules. Regression checks enforce stylesheet boundaries and verify that JSX references and shared selector imports resolve to real, matching exports. Validation: all 38 frontend tests, formatting, TypeScript, and the production build pass. Browser checks covered home, inputs, configuration, constraints, run review, comparison, results, and the dashboard chart picker in light/dark themes. Measured input and scenario-settings layouts match their pre-migration baselines, and table menus close correctly with Escape and outside clicks.
    >
    > Commit: [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126).

3. Repeated control patterns, e.g., icon-button appears in a ton of components, and they all reuse the same styling concept. This should be its own component and imported throughout. DataDialog and RunModel (and probably some other components) contain the same generic action-button convention, just with different class names like primary/secondary. Warrants a shared button component, probably with variant as a prop. I'd suggest doing a sweep through of the whole codebase and extracting all these repeated control patterns (probably not just limited to buttons) and turning them into shared components

    > Reply: Fixed. Extended the existing Button into its own component with primary, secondary, and danger variants, and reused it for form, run, table-dialog, and dashboard actions. Added IconButton with a required accessible label, native disabled behaviour, and plain, surface, and toolbar variants; IconLink shares its appearance while keeping CSV downloads as native links. Dashboard chart settings now expose their pressed state. Buttons default to type="button", while form submissions explicitly retain type="submit". Button and icon styles live in collocated SCSS modules. Field, SelectField, SearchField, and ToggleField now share FormControls.module.scss, so they no longer depend on chart or scenario-page styles; repeated field wrappers use Field, and SearchField has an accessible name. Page styles retain contextual layout and size overrides through data-control attributes and icon size variables. Specialised navigation, table sorting, and chart legend controls keep their own behaviour. Validation: formatting, TypeScript, all 36 frontend tests, and the production build pass. Browser checks covered light/dark rendering, edit/discard, new-scenario and dashboard dialogs, result-table controls, disabled buttons, and download-link semantics, with no console warnings or errors.
    >
    > Commit: [d434a56](https://github.com/agoenergy/pypsa-spice/commit/d434a56ae5b8b9e0d5fbcf6917027d39079dcbfe).

4. JSX throughout is not at all formatted at all, which makes the whole codebase look like an explosion. Please introduce some normal code formatting (suggest prettier package, but even just a line char limit and auto saving should do the trick! :)

    > Reply: Fixed. Added a pinned Prettier development dependency, a frontend configuration with a 120-character print width, and `npm run format` / `npm run format:check`. Formatted the frontend source and configuration files in a mechanical pass before the polling edits. Generated files, bundled public assets, dependencies, and the lockfile are excluded. The Web frontend GitHub Actions workflow checks formatting, runs tests, and builds the frontend on relevant pull requests and pushes to main/develop. Local formatting checks, all 36 frontend tests, and the production build pass.
    >
    > Commit: [7eac4fc](https://github.com/agoenergy/pypsa-spice/commit/7eac4fc3bee383cdf63b0846d7d00ac2005cc61b).

5. Testing is minimal and a bit random at the moment. The tests that are there are fairly sensible, but I'm not sure why there are only tests for a couple of components. If tests are desired, I'd suggest to actually have a test file per component (or conceptual 'group of behaviour'), living alongside the component, as is more conventional. However, tests are normally run separately in the development process anyway, and at specific points e.g., after fixing a bug, before pushing/building/deploying etc., and the dev needs to know when to run them. For purposes of this repo (assuming it's going to be run by non-devs), it might not make sense to include tests since they would not be so meaningful for people looking at the codebase and they would not know how to interpret them
    > Reply: SKu: I think having tests would be good so agents can also run tests to check changes. I am not sure which approach of grouping test files is better.

    > Follow-up: I agree with keeping tests so developers and agents can catch regressions; non-dev users do not need to run or interpret them, and CI already runs them automatically. I'd favour tests grouped by behaviour, such as run polling, configuration save/discard, or chart calculations, alongside the code that owns that behaviour. Add dedicated component tests where a component has meaningful interactions; a test file for every visual wrapper would add maintenance without much value.
    >
    > Partially addressed: [7eac4fc](https://github.com/agoenergy/pypsa-spice/commit/7eac4fc3bee383cdf63b0846d7d00ac2005cc61b) adds CI and run-monitor regression tests; [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126) adds stylesheet regression checks. The grouping recommendation is recorded in [4bd3582](https://github.com/agoenergy/pypsa-spice/commit/4bd3582f551d8c607b296e0e54786af549bb841b), but broader behaviour coverage remains open.

#### `Plot.tsx`
1. Too much responsibility, seems to be handling a lot of chart utility (see exported functions that are reused in `ChartCard` and `DashboardChartCard` components). Suggest moving the chart data processing functions (`aggregate`, `getLegendValues`, `differenceAggregates`, `buildDifferenceRows`) into a `chartData.ts` or `chartUtils.ts`, which this component, `ChartCard`, and `DashboardChartCard` import from as necessary. Move `ChartLegend` out as a separate component which is also imported by these three components

    > Open. No addressing commit yet; chart data helpers and `ChartLegend` still live in `Plot.tsx`.

2. Also consider moving logic that acts on processed chart data to derive the Plotly trace objects into a separate file like `plotTraces.ts` (`traces`, `differenceTraces`, `stackedBarTotalTrace`, etc.). They are only used in this Plot component but are conceptually separate from the actual Plot lifecycle and orchestration which the component handles. Minor suggestion to move them out, and expose a single main function from `plotTraces.ts` which internally calls all those functions

    > Open. No addressing commit yet; trace-building functions still live in `Plot.tsx`.

3. Typing trace as `Record<string, unknown>` is probably too loose for an object that is supposed to conform to Plotly's trace shape. It would be better to type it as a Plotly trace type to catch invalid property values that may get passed in

    > Open. No addressing commit yet; the trace still uses `Record<string, unknown>`.

#### `DataDialog`
1. Similar issue of too many responsibilities in one component. Consider turning it into a small module instead. `DataDialog` component itself should only handle deciding whether to render hourly or yearly, then `HourlyDataDialog` and `YearlyDataDialog` as components that it imports. Current `YearlyDataDialog` is a bit of a monster and contains several different concepts - suggest splitting the jsx into smaller components `YearlyTableToolbar`, `YearlyResultsTable`, `YearlyTableFooter` to start. Move table specific functions into some table utils file. Suggested structure:
```bash
DataDialog/
├── DataDialog.tsx
├── DataDialog.css
├── HourlyDataDialog.tsx (imported by DataDialog)
├── YearlyDataDialog.tsx (imported by DataDialog)
├── YearlyTableToolbar.tsx (imported by YearlyDataDialog)
├── YearlyResultsTable.tsx (ditto)
├── YearlyTableFooter.tsx (ditto)
└── tableUtils.ts
```

    > Open. No addressing commit yet; both dialog components and table helpers still live in `DataDialog.tsx`.

2. component currently queries the DOM to find out the `<details>` state and decide behaviour. This is a) not very react-y / recommended, and b) fragile because its querying the whole document and searching for a CSS selector, rather than just looking in the yearly dialog that it actually wants to query. Simplest fix without a huge rewrite with state management is to `useRef` to reference yearly data dialog, attach it, then query that instead. Alternatively, use state to manage the menus so react owns it and is able to recognise when menus are open/closed, but this would be much more of a rewrite and may be too much hassle for the gain

    > Reply to DataDialog 2: Fixed on 2026-09-12. YearlyDataDialog now attaches a useRef to its dialog section and queries only that section for open details menus. Escape closes open menus before closing the dialog, and outside clicks dismiss menus. These behaviours were verified in the browser after the style migration.
    >
    > Commit: [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126).

#### `RunModel`
1. `useEffect` poll is issuing a request per second without waiting for the previous `getModelRun` request to complete. I'm not sure how long a run would take to complete, but it might lead to overwriting responses if the API is slow. Suggest using recursive `setTimeout` instead of `setInterval` (though note this needs a decision on whether it should continue indefinitely or not if the polling fails)
```typescript
useEffect(() => {
  if (!run || !activeStatuses.has(run.status)) return;

  let cancelled = false;
  let timeoutId: number;

  const poll = async () => {
    // try getModelRun

      timeoutId = window.setTimeout(poll, 1000);
    }
  };

  timeoutId = window.setTimeout(poll, 1000);

  return () => {
    cancelled = true;
    window.clearTimeout(timeoutId);
  };
}, [run?.id]);
// Essentially this does request -> wait for response -> wait 1s -> request again
```

    > Reply: Fixed. `App` now owns one shared monitor through `useModelRun`, and passes its run state to `RunModel`. This replaces both status intervals and the run page's separate latest-run lookup. Each active-run request finishes before a one-second delay starts. Cleanup and start/cancel responses abort and invalidate pending status requests so late responses cannot overwrite newer state. Requests time out after 15 seconds. Connection failures retain the last known status, show a retry message, and retry after 2, 4, 8, 16, then at most 30 seconds; a successful response restores the normal delay. Polling stops for succeeded, failed, or canceled runs. Idle workspaces perform only the initial lookup; reload the app to discover a run started elsewhere. Fourteen regression cases cover slow requests, active and terminal statuses, idle workspaces, cleanup, cancellation and startup races, timeout, retry recovery, and stop/start lifecycle handling.
    >
    > Commit: [7eac4fc](https://github.com/agoenergy/pypsa-spice/commit/7eac4fc3bee383cdf63b0846d7d00ac2005cc61b).

2. the JSX is giant and the component is doing a lot of different things. Conceptually, i'd suggest separating out `RunSummary`, `RunConfiguration`, and `RunMonitor` as their own components

    > Open. No addressing commit yet; `RunSummary`, `RunConfiguration`, and `RunMonitor` have not been extracted as UI components.

3. encoded repititons again occur inside the JSX, e.g., `ReviewItem` has already been extracted out, but the same structure essentially appears again in run-dimensions and in repository-review. Suggest using a generic presentation primitive of span and b within a div across all these (essentially what `ReviewItem` is already, just use it more consistently). Check for other repeated patterns -- header within section also seems to happen alot

    > Open. No addressing commit yet; the dimension and repository sections still repeat the label/value markup.

4. `useMemo` is unnecessary here as `dimensions` and `scenarioSummary` are tiny computations, and memo-ising adds unnecessary overhead

    > Open. No addressing commit yet; `dimensions` and `scenarioSummary` still use `useMemo`.


### Changes made:

- Added missing Vite Typescript declaration file (vite-env.d.ts). Without this, Typescript tooling sometimes does not know that side-effect imports like `import "./ScenarioConfigEditor.css"` are valid and you get these `Cannot find module or type declarations for side-effect import of ...` errors.

    > Reply: Kept. `vite-env.d.ts` references Vite's client types, so TypeScript recognises stylesheet and asset imports used by the frontend. `npm run check` now completes without the side-effect import errors shown below.
    >
    > Commit: [01fb6e9](https://github.com/agoenergy/pypsa-spice/commit/01fb6e978b08ac7f49daf81a7c07421d891dd8a4).

    ![Image](https://github.com/user-attachments/assets/ddf9e0e7-a48a-487d-adf8-4669314959ab)

- Added run-web-locally.sh. Existing run file builds the app (generates dist/ folder), which is only needed for deployment and not for local dev. Also updated the relevant section in `overview.md` with instructions on how to set-up the app (first time only) and how to run locally. Currently, to access the app locally you visit http://127.0.0.1:5173/ui/. The /ui/ path is specified in your vite.config.ts. You can change it to an empty string if you want to get rid of that from the URL, but you might have to update your frontend API calls as well.

    > Reply: Kept. `run-web-locally.sh` starts Uvicorn and Vite directly, checks that both Python and Node dependencies are available, and stops the backend when the development server exits. It does not build `dist/`. The `/ui/` base remains intentional because FastAPI mounts the production frontend there, so local and production URLs follow the same routing contract.
    >
    > Commit: [6341e59](https://github.com/agoenergy/pypsa-spice/commit/6341e595c3eebac3fecc581da0bb28e970ee0684).

### Further comments/suggestions for improvement:

- Your frontend seems to be polling your backend every 2 seconds (can be seen when you run `run-web-locally.sh`, and you see the GET request logged every 2s). This comes from line 83 in `App.tsx`, which calls `refreshRunStatus` every 2000 ms). I'm not sure if this is a deliberate design choice, but if it's not a technical requirement to actually keep the model-run status fresh, I would suggest just fetching the latest run once, then keep polling only when a run is queued/running/cancelling. Alternatively maybe poll less often like every 10-20s or so. Every 2s is a lot of background activity.

    > Reply: Fixed, with a further update on 2026-09-11. `App.tsx` now owns the shared `useModelRun` monitor described in batch 2. It fetches the latest run on startup, then sequentially fetches that run while queued, running, or canceling. Terminal responses stop polling. Idle use produces no recurring requests, and the run page shares the same state rather than starting another poll.
    >
    > Commits: [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) limits background polling to active runs; [7eac4fc](https://github.com/agoenergy/pypsa-spice/commit/7eac4fc3bee383cdf63b0846d7d00ac2005cc61b) replaces it with the shared sequential monitor.

- Project currently does not declare Node type definitions, so your vite.config.ts file complains. 

    > Reply: Fixed. Added `@types/node` to `frontend/package.json` and the lockfile. This supplies the types used by `vite.config.ts`, including the `node:path` import and Node globals, so the Vite configuration is covered by the normal TypeScript check.
    >
    > Commit: [5e97545](https://github.com/agoenergy/pypsa-spice/commit/5e97545d0267d97c246f82c76c6a2be4dfb17ceb).

    ![Image](https://github.com/user-attachments/assets/a5b4309f-e0f3-4809-a04b-77ad83023b7f)

    Fix is to install node types in your environment:

```bash
npm install --prefix pypsa_spice_web/frontend --save-dev @types/node
```

- The chain to access the Flexo font file right now is quite fragile -- it lives in `pypsa-spice-vis/design/`, and Vite forwards the `/brand` request from the browser to FastAPI to access the font. If the Streamlit directory is renamed/removed/somehow becomes inaccessible to the web app, the font will break. I'd suggest bundling whatever custom font you use with the web app itself, in `public/` or in `assets/` so the app owns it directly.

    > Reply: Fixed. Copied Flexo to `frontend/public/fonts/Flexo-Medium.woff2` and the logo to `frontend/public/pypsa-logo.svg`. Frontend styles and components now reference those bundled files. The `/brand` Vite proxy was removed, so local rendering no longer depends on files inside `pypsa-spice-vis/`.
    >
    > Commits: [5e97545](https://github.com/agoenergy/pypsa-spice/commit/5e97545d0267d97c246f82c76c6a2be4dfb17ceb) bundles the assets; [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) updates their references and removes the `/brand` proxy.

- `utility.tsx` does not seem to contain any tsx, so rename to `utility.ts`. In general .tsx and .jsx extension should only be used when the file contains tsx/jsx

    > Reply: Fixed. Renamed `utility.tsx` to `utility.ts` because it contains no JSX. All imports and the utility test now point to the `.ts` module, and the TypeScript build resolves it without an extension-specific workaround.
    >
    > Commit: [02d5bbc](https://github.com/agoenergy/pypsa-spice/commit/02d5bbca854e12aea48aac0e762ac2c583e925a2).

- Several types and interfaces in `utility.tsx` are actually also exported and used elsewhere, e.g., in particular the ones for the dashboard. Since you have a `types.ts` file for global types/interfaces, these should be moved there. Alternatively, store the dashboard types/interfaces with the Dashboard component, and import in utils. Usually either all types/interfaces live in a types file, or types live close to the code that owns them plus a global types file for all other shared types. Currently it's a mix in the codebase.

    > Reply: Fixed. Moved the dashboard schema, dashboard rows and chart configuration, navigation view type, workspace options, and shared selection models into `types.ts`. `utility.ts` now imports those definitions and contains parsing, storage, selection, and navigation behaviour only. Components use the same shared types instead of importing type declarations from a utility module.
    >
    > Commit: [02d5bbc](https://github.com/agoenergy/pypsa-spice/commit/02d5bbca854e12aea48aac0e762ac2c583e925a2).

- `Homepage`: This file contains a lot of components and non-component helpers and is generally trying to handle a lot of different things. As a first step, I'd suggest:
    - moving WorkspaceCard into `components/` with its own tsx and css - this can include ScenarioList and ProjectDashboardList
    - moving WorkflowStep into `components/` also and importing here. Additionally suggest renaming this to something more semantic like WorkspaceActionCard or WorkflowStepCard

    > Reply: Fixed. `HomePage.tsx` now owns the page-level inventory and selected-project state. `WorkspaceCard.tsx` owns the input scenario, result run, and saved-dashboard lists, while `WorkspaceActionCard.tsx` owns the workflow actions. Both components have collocated SCSS modules, which removes their markup and styling details from the page component.
    >
    > Commit for both extractions: [01dfc26](https://github.com/agoenergy/pypsa-spice/commit/01dfc269a60dc5f41cc86f45cb358fb8f68319a0).

- `InputEditor`:
    1. There are many repeated control patterns in here - extract Dropdown (labelled select elements), SearchField, and Button out as components. Ideally, these could also be reused in other pages where you have similar components encoded in the page (e.g., there are many dropdowns in your scenario settings sub page), assuming you don't need different appearance/UI in the other pages.
    2. The page currently bundles the 'By technology' and 'By table' subpages, but these have substantially different child views and would warrant extracting them into their own files. I'd recommend:
        - keep InputEditor as the overall page container component
        - extract TechnologyEditor and TableEditor as their own separate components
        - let CellEditor live with TableEditor, since it is currently only used there. Since it's (relatively) small, you could even consider just integrating it into TableEditor's return, rather than having it live outside

    > Reply: Fixed. Added typed `SelectField`, `SearchField`, `Button`, and `ToggleField` components in `FormControls.tsx` and reused them across the editor flows. `InputEditor.tsx` now handles the page mode and shared navigation only. `TechnologyEditor.tsx` owns the technology view, while `TableEditor.tsx` owns table loading, filtering, pagination, and its private `CellEditor`.
    >
    > Commits: [01dfc26](https://github.com/agoenergy/pypsa-spice/commit/01dfc269a60dc5f41cc86f45cb358fb8f68319a0) addresses issues 1 and 2 with shared controls and separate editor components; [d434a56](https://github.com/agoenergy/pypsa-spice/commit/d434a56ae5b8b9e0d5fbcf6917027d39079dcbfe) further consolidates controls for issue 1.

- `ScenarioConfigEditor`:
    1. Similar to above, extract repeated control patterns into components (dropdowns, toggles, buttons)
    2. Please consider taking out ScenarioSettings and CO2Editor out as separate components too - this already seems to have been done for the third subpage (RunModel), so currently the organisation is rather confusing, and the file is trying to do too much and too many different things
    3. There are many helper functions living in this file. Consider a dedicated `ScenarioConfigEditorUtils.ts` for things like `validateSection`, `validateFuelLimits` and so on
    4. Loading config, creating original draft, dirty calculation, validation, saving and discarding are all currently being handled in several different places across the code and interleaved with other concerns. I'd suggest creating a custom hook to just handle such stateful workflow

    > Reply: Fixed. Split the two large subpages into `ScenarioSettings.tsx` and `Co2Editor.tsx`, alongside the existing `RunModel.tsx`. Repeated mapping-table and form controls now come from `ScenarioConfigControls.tsx` and `FormControls.tsx`. Pure conversion and validation functions live in `ScenarioConfigEditorUtils.ts`. The `useScenarioConfigEditor` hook now owns the request lifecycle, original and draft values, dirty checks, validation, saving, and discarding, leaving `ScenarioConfigEditor.tsx` to compose the page.
    >
    > Commits: [9211eeb](https://github.com/agoenergy/pypsa-spice/commit/9211eeb4b742631b02e07b746f0bafab09006756) addresses issues 1 through 4 with shared controls, separate subpages, utilities, and the editor hook; [d434a56](https://github.com/agoenergy/pypsa-spice/commit/d434a56ae5b8b9e0d5fbcf6917027d39079dcbfe) further consolidates controls for issue 1.

- `Sidebar`:
    1. Sidebar css classes are currently located in `App.css` -- these should be moved to a `Sidebar.css`
    2. In general the Sidebar component is not well-structured and contains multiple components within its file (SidebarSection, HomeSidebarSection, InputSidebarSection etc.) -- not unacceptable, but currently each sidebar section component is a tiny wrapper function that just passes constants into SidebarSection. I would suggest simplifying to four components:
        - Sidebar: the main exported sidebar component that handles layout
        - SidebarSection: a single generic sidebar menu item (reuse for Home, Inputs, Scenario differences, Configure and run, Dashboards)
        - ResultsSidebarSection: seems different enough to warrant its own component for custom behaviour
        - SidebarFooter
    
        They could all still live in the same Sidebar tsx file, but Sidebar is the only exported component, and it comprises the other components. This would make for a cleaner and more readable return statement.

    > Reply: Fixed. Moved all sidebar rules out of the application stylesheet and into `Sidebar.module.scss`. `Sidebar.tsx` defines one generic item renderer for Home, Inputs, Scenario differences, Configure and run, and Dashboards. Results keeps its own nested section because it displays chart counts and an active subsection. The footer remains separate, and only the top-level `Sidebar` component is exported.
    >
    > Commits: [01dfc26](https://github.com/agoenergy/pypsa-spice/commit/01dfc269a60dc5f41cc86f45cb358fb8f68319a0) adds the stylesheet for issue 1; [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) connects it, removes the old global rules, and simplifies the components for issue 2.

- App.css:
    1. Besides sidebar classes, this file contains some other classes that are not used in `App.tsx` at all, and used across other components, e.g., dialog-backdrop, country-config, amongst others. Move these classes to the relevant component's css file, or if they are used across more than one component (e.g., dialog-backdrop), I suggest creating a `global.css` file and storing shared classes there
    2. The first 78 lines of App.css are also application global and conceptually belongs in a global style file instead. 

    > Reply: Fixed, with the full module migration completed on 2026-09-12. main.tsx imports global.scss once for fonts, tokens, resets, focus treatment, and reduced-motion rules. App.module.scss owns the application shell. DialogSurface.module.scss owns shared dialog defaults, and scenario layouts use explicit shared and page modules. All component classes are scoped; global.scss contains none.
    >
    > Commits: [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) extracts global foundations for issue 2; [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126) completes issue 1 by moving component and shared styles into scoped modules.

- Currently media query breakpoints are variable across components. I'd suggest a standardised breakpoint scale (typical breakpoints are something like 640px, 768px, 1024px, and 1280px). You *could* additionally consider moving to scss, which allows you to define breakpoint variables like this:

    > Reply: Fixed. Replaced the one-off 700, 900, 1100, and similar breakpoints with 640, 768, 1024, and 1280 pixels according to the nearest intended layout change. All frontend styles now use SCSS, with component and page rules in modules. They retain the same numeric breakpoint scale.
    >
    > Commits: [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) standardises the numeric breakpoints; [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126) completes the SCSS module migration. Shared breakpoint variables were not introduced.

```scss
$breakpoint-s: 640px;
$breakpoint-m: 768px;
$breakpoint-l: 1024px;
```

```scss
// Then in component's scss file
@media (max-width: variables.$breakpoint-s) {...}
```

- In general the pages and components are all rather messy and are trying to handle a lot of different concerns at once. Perhaps just additionally run the whole project through AI and ask it to generally refactor the codebase around clear component responsibilities, like extract repeated UI/control patterns, split large components into meaningful feature-level components, move shared styles into appropriate global/shared stylesheets, and remove duplicated markup, styling, and logic.

    > Reply: Addressed in the reviewed areas. The homepage, input editor, scenario configuration editor, controls, sidebar, shared types, and global styles now each have a clear owner. Repeated controls and state workflows were extracted only where more than one view benefits from them. Chart rendering and dashboard behaviour stayed outside this refactor, which kept the change set reviewable and avoided mixing structural cleanup with feature changes. TypeScript checks, frontend tests, and the production build pass after the split.
    >
    > Commits: [02d5bbc](https://github.com/agoenergy/pypsa-spice/commit/02d5bbca854e12aea48aac0e762ac2c583e925a2) for shared types, [01dfc26](https://github.com/agoenergy/pypsa-spice/commit/01dfc269a60dc5f41cc86f45cb358fb8f68319a0) for homepage and input editors, [9211eeb](https://github.com/agoenergy/pypsa-spice/commit/9211eeb4b742631b02e07b746f0bafab09006756) for scenario configuration, [a7c157f](https://github.com/agoenergy/pypsa-spice/commit/a7c157f4cda6810815ced76b6faac949ed1d41f2) for the shell and sidebar, [d434a56](https://github.com/agoenergy/pypsa-spice/commit/d434a56ae5b8b9e0d5fbcf6917027d39079dcbfe) for controls, and [0c41834](https://github.com/agoenergy/pypsa-spice/commit/0c41834477e09641c837e067b092f15467792126) for scoped styles. The remaining batch-2 refactors are still open.
