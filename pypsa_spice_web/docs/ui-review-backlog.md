# UI and UX backlog

Findings from a review of the web workspace on 2026-08-20, taken by walking Home,
Inputs, Configure & run, Results, comparison mode and Dashboards at 1680x1050 in
light and dark mode, then reading the source.

One item is closed. The rest are open and listed here so the next person does not
have to rediscover them. Line references were accurate on 2026-08-20; selector and
symbol names are the more durable pointers.

## Closed

**Type scale.** 97% of text-bearing elements on the Power results page rendered
below 12px (328 at 9px, 371 at 10px, 56 at 11px) while the page title sat at 52px.
`global.scss` now defines seven size tokens with a 12px floor, all 158 size
declarations across 13 stylesheets map onto them, and `Plot.tsx` carries a matching
`chartFont` constant because Plotly draws its own SVG text and cannot read CSS.
Measured after the change: no element renders below 12px on any page, in either
theme.

## Correctness

These two change what the user concludes from the screen, so they come first.

### 1. Deep links drop the result run and the comparison

Opening `?run=scenario_03_run` lands on `scenario_01-1_run`. Adding
`&compare=scenario_03_run` drops the comparison entirely. Reproduced twice.

`App.tsx` starts `loadCatalog()` and `loadInputs()` in parallel, and the effect that
writes state back into the URL is guarded by
`if (!(selection.project || inputSelection.project)) return;` (`App.tsx:122`). When
the input catalog resolves first, that guard opens while `selection` is still empty,
so the effect rewrites the URL without `run` or `compare`. The results catalog then
calls `resolveOutputSelection`, which reads a URL that no longer holds them.

Read `locationParams()` once at mount into a ref and resolve from that, or hold the
URL writer until both catalogs have resolved. Until this is fixed nobody can share
or bookmark a comparison, which is the workspace's main analytical output.

### 2. Side-by-side comparison plots autoscale independently

Each plot in comparison mode computes its own axis range. Observed on one page:

- Interconnector capacity factor: left 0 to 1.5 p.u., right 0 to 2 p.u.
- Nodal flow between regions: left -4000 to 2000 MW, right -2000 to 2000 MW
- Hourly elec price: left x-axis ends Jan 2026, right ends Nov 2025

Two bars of equal pixel height therefore stand for different values, in the one view
built for visual comparison. Compute a shared `[min, max]` across both responses and
pass an explicit `yaxis.range` to both plots, plus `xaxis.range` for hourly charts.

## Layout and density

### 3. The two-column chart grid is dead code

`.chart-grid` sets `repeat(2, minmax(0, 1fr))` at `ChartCard.css:80` and a later rule
at `ChartCard.css:230` overrides it to a single column. The
`@media (max-width: 1150px)` fallback below is now a no-op. Every Results chart is
full width, so a six-bar stacked chart stretches to roughly 1300x340 and the Power
page runs about 8,000px tall for 16 charts. Restore the two-up grid outside
comparison mode, keep full width when comparing, and raise `.chart-body` above 340px
(`ChartCard.css:189`) so plots sit nearer 4:3.

### 4. All 16 charts fetch on page load

Every `ChartCard` fires its request on mount. Wrap in an `IntersectionObserver` so
off-screen cards hold a placeholder until they are scrolled near.

### 5. Two competing toolbars per chart

The card's own actions (view data, download CSV, expand) sit about 25px above
Plotly's modebar (download PNG, zoom, pan, autoscale, reset). Comparison mode draws
two modebars per card, 32 on the page. Trim the modebar to `['resetScale2d']` and
fold PNG export into the card's action row, or drop the modebar and expose zoom and
reset as card actions.

### 6. Empty charts cost a full card

"Battery's E/P ratio" draws a full-height card, a country picker and a Difference
toggle around the words "No values in this result table". Collapse charts with no
data to one compact row, or group them under a "No data in this run (2)" disclosure.

### 7. Per-chart controls have no page-level default

16 identical "All countries" selects. Per-chart scope is a deliberate choice and a
good one, but add a page-level default that seeds every card, with per-card override.

### 8. Result runs are unidentifiable

The picker offers `scenario_01_tag1`, `scenario_01_tag1 copy`, `scenario_03-98_ru2`
and `scenario_03-25_2` with nothing to tell them apart. Add each run's mtime and
modelled-year span to the option label and sort newest first.

### 9. Page identity repeats three and four times

Home shows "Workspace overview" in the top bar and again as the h1. Configure stacks
sidebar "Configure & run", sub-item "Scenario settings", h1 "Configure scenario_01"
and panel h2 "Scenario settings". Keep the h1, drop the top-bar copy and the panel h2
that repeats the active sub-nav item. This got more visible once the h1 shrank to
28px.

### 10. Home spends the fold on documentation

"How to use the workspace" takes six cards and about 250px on every visit, and its
first card is a disabled "Build input data skeleton / NOT IMPLEMENTED", so the
workflow opens on a dead end. Make the panel collapsible and remember the state.
Either drop the unimplemented step or move it last. Below it, "Result runs" spans two
grid columns and is usually empty.

## Forms

### 11. Configure & run keeps a disabled save bar on screen

The floating Save and Discard pair is always mounted with both buttons disabled, and
the pink Save at 42% opacity reads as an enabled primary action at a glance. Show the
bar only when the form is dirty. It is also 36px tall against the 44px control height
used everywhere else.

### 12. Inputs has no global save

Each table panel carries its own Save and Discard pair, three or more on a technology
page, so editing across tables means saving each one separately with no running total
of unsaved changes. Add one sticky bar with a count that commits every dirty table.

### 13. Number fields mix decimal conventions

`type="number"` fields render `0,05` and `0,1` from the OS locale directly under help
text that reads "0.05 means 5%". Switch to `type="text"` with
`inputMode="decimal"` and explicit parsing and formatting so the field and the
documentation agree.

### 14. The config form validates nothing

`data/example/project_01/input/scenario_01/scenario_config.yaml` holds
`start: 20266-01-01`, which is bad source data rather than a UI bug, and the form
shows year 20266 with no warning. Add range checks and an inline error, since this
form writes to disk.

### 15. Absolute paths leak on Configure

Configure prints `/Users/<name>/Desktop/.../scenario_config.yaml` while Inputs prints
the relative `example/project_01/input/...`. Use relative paths everywhere and put
the absolute path in a `title` attribute.

### 16. Pagination shows when there is one page

The footer still reads "Page 1 of 1" with a row count and both arrows disabled.
Hide the control when `pageCount === 1`.

### 17. Field widths ignore content

"Model year" gets a 622px input for four digits. Size numeric fields to their content
and let text fields expand.

## Dashboards

### 18. Delete sits in a row of unlabelled icons

New, Duplicate, Import, Export and Delete are five 40px icon buttons, with the red
trash 44px from Export. `window.confirm` catches the mistake, but move Delete into an
overflow menu or separate it with a divider.

### 19. The header is spatially disconnected

Selects end near x=714, the action cluster starts near x=1333, refresh sits at
x=1588. Two large gaps. Group the actions beside the selects or right-align them with
refresh.

### 20. Two titles compete

The h1 "Custom dashboards" sits above the editable "Untitled dashboard". The editable
one is the real title, so make it the visual h1 and drop the generic heading.

### 21. The description textarea is one clipped line

18px tall at 850px wide with a visible resize handle. Give it three rows and
`resize: vertical`.

### 22. Rows are full width only

Documented as deliberate, but it means heavy scrolling and no at-a-glance comparison.
A two-up option per row would make dashboards far more useful. This is a product
question, not a defect.

## Accessibility

### 23. The expanded chart is a modal that is not one

`.chart-card.expanded` (`ChartCard.css:94`) is `position: fixed; inset: 18px;
z-index: 100` with no `role="dialog"`, no `aria-modal`, no Escape handler, no focus
trap, and the page behind still scrolls. Escape to close is the reflex here.

### 24. Dialog behaviour is inconsistent

`DataDialog` and `NewScenarioDialog` do this well, with `role="dialog"`,
`aria-modal` and Escape. The Dashboards add-chart and import dialogs
(`DashboardPage.tsx:415` and `:449`) close on backdrop click only. No dialog traps
focus or locks body scroll. Extract one `Modal` that owns Escape, focus trap,
focus restore and scroll lock, then use it for all five.

### 25. Icon-only buttons carry no visible label

Refresh and the five dashboard actions rely on `title` and `aria-label`. The
`aria-label`s are present and correct throughout, which is genuinely well done, but a
first-time user hovers five icons to find Export.

## CSS structure

### 26. Page chrome lives in a component stylesheet

`.main-column`, `main`, `.button`, `.menu`, `.search`, `.notice` and `.page-title`
are all defined in `components/ChartCard.css` (lines 1 to 74 and 233). Move them to
`global.scss` or a new `shell.css`.

### 27. Corrections pile up instead of edits

`global.scss` ends in several blocks that undo earlier declarations. One zeroes
`border-left` on `.notice` and `.editor-warning` under the comment "Left-edge accent
bars are intentionally not part of this app's visual language". `.sidebar nav`
re-sets `flex: 1` after `.sidebar .workspace-nav` already did. `.sidebar-foot button`
re-sets `margin-left: 0` after `margin-left: auto`. The dead `.chart-grid` override
in item 3 is the same habit. Edit the original declarations and delete the
corrections, and drop dead selectors such as `.sidebar nav .disabled`.

Five `font-size` declarations stayed hardcoded during the type-scale work:
`.menu`, `.icon-button`, `.search span`, `.chart-actions a` and
`.sidebar-foot button`. They size glyphs on elements that now hold only SVGs with
explicit width and height, so they render nothing. Clean them up here.

## Suggested order

1. Shared axis ranges in comparison mode, item 2. Correctness, small diff.
2. Restore the two-up chart grid and delete the dead override, item 3.
3. Fix the deep-link race, item 1.
4. Lazy-render off-screen charts, item 4, and tame the Plotly modebar, item 5.
5. Form work: dirty-only save bar, global Inputs save, decimal handling, items 11
   to 13.
6. Shared `Modal` and Escape on the expanded chart, items 23 and 24.
7. CSS consolidation, items 26 and 27.
