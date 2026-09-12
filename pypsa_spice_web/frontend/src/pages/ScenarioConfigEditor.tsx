import styles from "./ScenarioConfigEditor.module.scss";
import editorPanelStyles from "../components/EditorPanel.module.scss";
import workspaceFeedbackStyles from "../components/WorkspaceFeedback.module.scss";
import scenarioConfigLayoutStyles from "./ScenarioConfigLayout.module.scss";
import scenarioConfigControlsStyles from "./ScenarioConfigControls.module.scss";
import dataTableStyles from "../components/DataTable.module.scss";
import workspaceUtilitiesStyles from "../components/WorkspaceUtilities.module.scss";
import Button from "../components/Button";
import IconButton from "../components/IconButton";
import { Field, ToggleField } from "../components/FormControls";
import { useEffect, useState } from "react";
import { createPortal } from "react-dom";
import { ClipboardCheck, Code2, Plus, Settings2, Trash2 } from "lucide-react";

import Co2Editor from "./Co2Editor";
import PageHeader from "../components/PageHeader";
import ResultsToc from "../components/ResultsToc";
import RunModel from "../components/RunModel";
import SaveDiscardActions from "../components/SaveDiscardActions";
import sidebarStyles from "../components/Sidebar.module.scss";
import type { InputRow, InputSelection, ModelRunMonitor } from "../types";
import { confirmDiscardChanges } from "../utility";
import { MappingTable } from "./ScenarioConfigControls";
import {
  COMBINED_SECTION,
  CONSTRAINT_DEFAULTS,
  SECTION_LABELS,
  YEAR_FIELDS,
  constraintAnchor,
  inputNumber,
  isScalar,
  isYearKey,
  isYearMatrix,
  normaliseSection,
  objectValue,
  prettyConfigLabel,
} from "./ScenarioConfigEditorUtils";
import ScenarioSettings from "./ScenarioSettings";
import useScenarioConfigEditor from "./useScenarioConfigEditor";

const icons = { scenario_configs: Settings2, [COMBINED_SECTION]: Code2, review_run: ClipboardCheck };

export default function ScenarioConfigEditor({
  selection,
  country,
  runMonitor,
  onOpenResults,
}: {
  selection: InputSelection;
  country: string;
  runMonitor: ModelRunMonitor;
  onOpenResults: (runName: string, dataset: string, project: string) => void;
}) {
  const [section, setSection] = useState(() => {
    const requested = new URLSearchParams(window.location.search).get("step") || "scenario_configs";
    return normaliseSection(requested);
  });
  const [navigationTarget, setNavigationTarget] = useState<HTMLElement | null>(null);
  const editor = useScenarioConfigEditor(selection, section);
  useEffect(() => {
    setNavigationTarget(document.getElementById("config-section-tabs"));
  }, []);

  const chooseSection = (name: string) => {
    if (!confirmDiscardChanges()) return;
    setSection(name);
    const params = new URLSearchParams(window.location.search);
    params.set("view", "configure");
    params.set("step", name);
    params.delete("section");
    window.history.pushState(null, "", `?${params.toString()}`);
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const navigation =
    navigationTarget &&
    createPortal(
      <nav className={sidebarStyles["submenu-list"]} aria-label="Configuration and run pages">
        {Object.keys(SECTION_LABELS).map((name) => {
          const SectionIcon = icons[name as keyof typeof icons];
          return (
            <button
              key={name}
              className={`${sidebarStyles["submenu-item"]} ${section === name ? sidebarStyles["active"] : ""}`}
              onClick={() => chooseSection(name)}
              aria-current={section === name ? "page" : undefined}
            >
              <SectionIcon aria-hidden="true" />
              <b>{SECTION_LABELS[name]}</b>
            </button>
          );
        })}
      </nav>,
      navigationTarget,
    );

  if (section === "review_run")
    return (
      <>
        {navigation}
        <RunModel
          selection={selection}
          onEditConfiguration={() => chooseSection("scenario_configs")}
          runMonitor={runMonitor}
          onOpenResults={onOpenResults}
        />
      </>
    );

  return (
    <>
      {navigation}
      <PageHeader title={`Configure ${selection.scenario}`} />
      <div className={styles["config-layout"]}>
        <section
          className={[
            [editorPanelStyles["editor-panel"], styles["config-panel"]].join(" "),
            section === COMBINED_SECTION ? "" : "",
          ]
            .filter(Boolean)
            .join(" ")}
        >
          <header className={editorPanelStyles["editor-panel-head"]}>
            <div>
              <h2>{SECTION_LABELS[section]}</h2>
              {editor.config && <code>{editor.config.path}</code>}
            </div>
          </header>
          {editor.error && (
            <div className={[workspaceFeedbackStyles["notice"], workspaceFeedbackStyles["error"]].join(" ")}>
              {editor.error}
              <button onClick={() => void editor.reload()}>Reload</button>
            </div>
          )}
          {!editor.loading && editor.validationError && (
            <div className={[workspaceFeedbackStyles["notice"], workspaceFeedbackStyles["error"]].join(" ")}>
              {editor.validationError}
            </div>
          )}
          {editor.loading ? (
            <div className={editorPanelStyles["editor-loading"]}>
              <span className={workspaceFeedbackStyles["spinner"]} />
              Reading configuration…
            </div>
          ) : section === "scenario_configs" ? (
            <ScenarioSettings value={editor.draft} country={country} onChange={editor.setDraft} />
          ) : (
            <CombinedConstraintsEditor
              value={editor.draft}
              country={country}
              fuelRows={editor.fuelRows}
              fuelError={editor.fuelError}
              onFuelRowsChange={editor.setFuelRows}
              onChange={editor.setDraft}
            />
          )}
          <SaveDiscardActions
            hasChanges={editor.dirty}
            saving={editor.saving}
            saveDisabled={Boolean(editor.validationError)}
            status={editor.success}
            floating
            avoidSideControl={section === COMBINED_SECTION}
            onDiscard={editor.discard}
            onSave={() => void editor.save()}
          />
        </section>
      </div>
    </>
  );
}

function CombinedConstraintsEditor({
  value,
  country,
  fuelRows,
  fuelError,
  onFuelRowsChange,
  onChange,
}: {
  value: Record<string, unknown>;
  country: string;
  fuelRows: InputRow[];
  fuelError: string;
  onFuelRowsChange: (rows: InputRow[]) => void;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const co2 = (value.co2_management || {}) as Record<string, unknown>;
  const constraints = (value.custom_constraints || {}) as Record<string, unknown>;
  const visibleConstraints = country === "ALL" ? Object.values(constraints) : [constraints[country]];
  const constraintNames = [
    ...new Set([
      ...Object.keys(CONSTRAINT_DEFAULTS),
      ...visibleConstraints.flatMap((raw) => Object.keys(objectValue(raw))),
    ]),
  ];
  const tocItems = [
    { id: "config-co2-management", label: "CO₂ management" },
    ...constraintNames.map((name) => ({ id: constraintAnchor(name), label: prettyConfigLabel(name) })),
  ];
  return (
    <>
      <div className={styles["combined-config"]}>
        <section
          className={styles["combined-config-section"]}
          id="config-co2-management"
          aria-labelledby="co2-management-heading"
        >
          <header>
            <h3 id="co2-management-heading">CO₂ management</h3>
            <p>Country carbon cap or price by year.</p>
          </header>
          <Co2Editor value={co2} country={country} onChange={(next) => onChange({ ...value, co2_management: next })} />
        </section>
        <section className={styles["combined-config-section"]} aria-labelledby="custom-constraints-heading">
          <header>
            <h3 id="custom-constraints-heading">Custom constraints</h3>
            <p>Activate constraints and edit their parameters directly.</p>
          </header>
          <CustomConstraintsEditor
            value={constraints}
            country={country}
            countries={Object.keys(co2)}
            fuelRows={fuelRows}
            fuelError={fuelError}
            onFuelRowsChange={onFuelRowsChange}
            onChange={(next) => onChange({ ...value, custom_constraints: next })}
          />
        </section>
      </div>
      <ResultsToc
        id="config-section-list"
        heading="Sections"
        panelLabel="Configuration sections on this page"
        listName="section list"
        entries={tocItems.map((item) => ({ key: item.id, href: `#${item.id}`, label: item.label }))}
      />
    </>
  );
}

function CustomConstraintsEditor({
  value,
  country,
  countries,
  fuelRows,
  fuelError,
  onFuelRowsChange,
  onChange,
}: {
  value: Record<string, unknown>;
  country: string;
  countries: string[];
  fuelRows: InputRow[];
  fuelError: string;
  onFuelRowsChange: (rows: InputRow[]) => void;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const countryNames = country === "ALL" ? [...new Set([...countries, ...Object.keys(value)])] : [country];
  const constraintNames = [
    ...new Set([
      ...Object.keys(CONSTRAINT_DEFAULTS),
      ...countryNames.flatMap((name) => Object.keys(objectValue(value[name]))),
    ]),
  ];
  const updateConstraint = (countryName: string, constraintName: string, next: Record<string, unknown>) => {
    const countryConstraints = objectValue(value[countryName]);
    onChange({ ...value, [countryName]: { ...countryConstraints, [constraintName]: next } });
  };
  return (
    <div className={[scenarioConfigLayoutStyles["config-form"], styles["constraints-form"]].join(" ")}>
      {constraintNames.map((constraintName) => {
        const activeCountries = countryNames.filter(
          (countryName) => objectValue(objectValue(value[countryName])[constraintName]).activate === true,
        ).length;
        return (
          <section className={styles["constraint-group"]} id={constraintAnchor(constraintName)} key={constraintName}>
            <header>
              <h4>{prettyConfigLabel(constraintName)}</h4>
              <span>
                {activeCountries} of {countryNames.length} active
              </span>
            </header>
            <div className={styles["constraint-country-list"]}>
              {countryNames.map((countryName) => (
                <ConstraintCard
                  key={countryName}
                  name={constraintName}
                  country={countryName}
                  value={objectValue(objectValue(value[countryName])[constraintName])}
                  fuelRows={fuelRows}
                  fuelError={fuelError}
                  onFuelRowsChange={onFuelRowsChange}
                  onChange={(next) => updateConstraint(countryName, constraintName, next)}
                />
              ))}
            </div>
          </section>
        );
      })}
      {countryNames.length === 0 && (
        <div className={[editorPanelStyles["editor-empty"], editorPanelStyles["compact"]].join(" ")}>
          No countries are available for custom constraints.
        </div>
      )}
    </div>
  );
}
function ConstraintCard({
  name,
  country,
  value,
  fuelRows,
  fuelError,
  onFuelRowsChange,
  onChange,
}: {
  name: string;
  country: string;
  value: Record<string, unknown>;
  fuelRows: InputRow[];
  fuelError: string;
  onFuelRowsChange: (rows: InputRow[]) => void;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const editableValue = { ...(CONSTRAINT_DEFAULTS[name] || { activate: false }), ...value };
  const active = editableValue.activate === true;
  const fields = Object.keys(editableValue).filter((key) => key !== "activate");
  return (
    <section className={[styles["constraint-card"], active ? styles["active"] : ""].filter(Boolean).join(" ")}>
      <header>
        <h5>{country}</h5>
        <label className={styles["constraint-toggle"]}>
          <input
            type="checkbox"
            aria-label={`Activate ${prettyConfigLabel(name)} for ${country}`}
            checked={active}
            onChange={(event) => onChange({ ...editableValue, activate: event.target.checked })}
          />
          <i aria-hidden="true" />
          <span>{active ? "Active" : "Inactive"}</span>
        </label>
      </header>
      {active &&
        (name === "production_constraint_fuels" ? (
          <FuelProductionLimitsEditor
            country={country}
            value={editableValue}
            rows={fuelRows}
            error={fuelError}
            onRowsChange={onFuelRowsChange}
            onChange={onChange}
          />
        ) : fields.length ? (
          <ConstraintFields constraintName={name} value={editableValue} onChange={onChange} />
        ) : (
          <p className={styles["constraint-empty"]}>No parameters required.</p>
        ))}
    </section>
  );
}

function FuelProductionLimitsEditor({
  country,
  value,
  rows,
  error,
  onRowsChange,
  onChange,
}: {
  country: string;
  value: Record<string, unknown>;
  rows: InputRow[];
  error: string;
  onRowsChange: (rows: InputRow[]) => void;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const countryRows = rows.filter((row) => String(row.country) === country);
  const selected = new Set(Array.isArray(value.fuels) ? value.fuels.map(String) : []);
  const carriers = [...new Set([...countryRows.map((row) => String(row.carrier)), ...selected])].filter(Boolean).sort();
  const years = [...new Set(countryRows.map((row) => String(row.year)).filter(isYearKey))].sort(
    (left, right) => Number(left) - Number(right),
  );
  const rowsByCell = new Map(countryRows.map((row) => [`${String(row.carrier)}\u0000${String(row.year)}`, row]));
  const toggleCarrier = (carrier: string, checked: boolean) => {
    const next = new Set(selected);
    if (checked) next.add(carrier);
    else next.delete(carrier);
    onChange({ ...value, fuels: carriers.filter((name) => next.has(name)) });
  };
  const updateCell = (rowId: number, raw: string) =>
    onRowsChange(rows.map((row) => (row.__row_id === rowId ? { ...row, max_supply__mwh_year: raw } : row)));
  if (error)
    return (
      <div
        className={[
          workspaceFeedbackStyles["notice"],
          workspaceFeedbackStyles["error"],
          styles["fuel-limit-notice"],
        ].join(" ")}
      >
        {error}
      </div>
    );
  if (!countryRows.length && !selected.size)
    return <p className={styles["constraint-empty"]}>No fuel supply rows are available for {country}.</p>;
  return (
    <div className={[styles["constraint-fields"], styles["fuel-limit-fields"]].join(" ")}>
      <div className={[scenarioConfigControlsStyles["config-entry-block"], styles["constraint-wide"]].join(" ")}>
        <p className={scenarioConfigControlsStyles["field-label"]}>Maximum annual fuel supply (MWh/year)</p>
        <p className={editorPanelStyles["field-help"]}>
          Select the carriers to constrain, then enter their maximum supply for each year. Use <code>inf</code> for no
          numerical limit. Values are saved to <code>power/fuel_supplies.csv</code>.
        </p>
        <div className={scenarioConfigControlsStyles["config-table-wrap"]}>
          <table
            className={[
              dataTableStyles["table"],
              [
                scenarioConfigControlsStyles["config-entry-table"],
                styles["matrix-table"],
                styles["fuel-limit-table"],
              ].join(" "),
            ].join(" ")}
          >
            <thead>
              <tr>
                <th>Apply</th>
                <th>Fuel</th>
                {years.map((year) => (
                  <th key={year}>{year}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {carriers.map((carrier) => (
                <tr key={carrier}>
                  <td>
                    <input
                      className={styles["fuel-limit-check"]}
                      type="checkbox"
                      aria-label={`Apply fuel production limit to ${carrier} in ${country}`}
                      checked={selected.has(carrier)}
                      onChange={(event) => toggleCarrier(carrier, event.target.checked)}
                    />
                  </td>
                  <th scope="row">{carrier}</th>
                  {years.map((year) => {
                    const row = rowsByCell.get(`${carrier}\u0000${year}`);
                    return (
                      <td key={year}>
                        {row ? (
                          <input
                            aria-label={`${carrier} ${year} maximum supply in MWh per year`}
                            inputMode="decimal"
                            value={String(row.max_supply__mwh_year ?? "")}
                            disabled={!selected.has(carrier)}
                            onChange={(event) => updateCell(row.__row_id, event.target.value)}
                          />
                        ) : (
                          <span className={styles["fuel-limit-missing"]} title="No fuel supply row for this year">
                            —
                          </span>
                        )}
                      </td>
                    );
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}

function ConstraintFields({
  constraintName,
  value,
  onChange,
}: {
  constraintName: string;
  value: Record<string, unknown>;
  onChange: (value: Record<string, unknown>) => void;
}) {
  return (
    <div className={styles["constraint-fields"]}>
      {Object.entries(value)
        .filter(([key]) => key !== "activate")
        .map(([key, raw]) => (
          <ConstraintField
            key={key}
            constraintName={constraintName}
            name={key}
            value={raw}
            onChange={(next) => onChange({ ...value, [key]: next })}
          />
        ))}
    </div>
  );
}

function ConstraintField({
  constraintName,
  name,
  value,
  onChange,
}: {
  constraintName: string;
  name: string;
  value: unknown;
  onChange: (value: unknown) => void;
}) {
  if (Array.isArray(value)) {
    return (
      <Field className={styles["constraint-wide"]}>
        <span>{prettyConfigLabel(name)}</span>
        <input
          value={value.join(", ")}
          placeholder="Comma-separated values"
          onChange={(event) =>
            onChange(
              event.target.value
                .split(",")
                .map((item) => item.trim())
                .filter(Boolean),
            )
          }
        />
      </Field>
    );
  }
  if (value && typeof value === "object") {
    const mapping = objectValue(value);
    const matrix =
      (constraintName === "maximum_power_generation_constraint" && name === "value") || isYearMatrix(mapping);
    const years = YEAR_FIELDS.has(name) || (Object.keys(mapping).length > 0 && Object.keys(mapping).every(isYearKey));
    const label =
      name === "value" && constraintName === "maximum_power_generation_constraint"
        ? "Generation limits (TWh)"
        : name === "value" && constraintName === "capacity_factor_constraint"
          ? "Technology capacity factors"
          : prettyConfigLabel(name);
    if (matrix) return <YearMatrixEditor label={label} value={mapping} onChange={onChange} />;
    if (years || Object.values(mapping).every(isScalar))
      return <MappingTable label={years ? "Year" : label} value={mapping} yearKeys={years} onChange={onChange} />;
    return (
      <section className={styles["constraint-object"]}>
        <h6>{prettyConfigLabel(name)}</h6>
        <ConstraintFields constraintName={constraintName} value={mapping} onChange={onChange} />
      </section>
    );
  }
  if (typeof value === "boolean") {
    return (
      <ToggleField
        className={styles["constraint-boolean"]}
        checked={value}
        onChange={onChange}
        label={prettyConfigLabel(name)}
      />
    );
  }
  if (name === "method") {
    return (
      <Field>
        <span>Method</span>
        <select value={String(value || "static")} onChange={(event) => onChange(event.target.value)}>
          <option value="static">Static</option>
          <option value="dynamic">Dynamic</option>
        </select>
      </Field>
    );
  }
  if (name === "math_symbol") {
    return (
      <Field>
        <span>Comparison</span>
        <select value={String(value || "<=")} onChange={(event) => onChange(event.target.value)}>
          <option value="<=">At most (≤)</option>
          <option value=">=">At least (≥)</option>
          <option value="==">Exactly (=)</option>
        </select>
      </Field>
    );
  }
  if (typeof value === "string") {
    return (
      <Field>
        <span>{prettyConfigLabel(name)}</span>
        <input value={value} onChange={(event) => onChange(event.target.value)} />
      </Field>
    );
  }
  return (
    <Field>
      <span>{prettyConfigLabel(name)}</span>
      <input
        type="number"
        step="any"
        value={String(value ?? "")}
        onChange={(event) => onChange(inputNumber(event.target.value))}
      />
    </Field>
  );
}

function YearMatrixEditor({
  label,
  value,
  onChange,
}: {
  label: string;
  value: Record<string, unknown>;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const [newTechnology, setNewTechnology] = useState("");
  const [newYear, setNewYear] = useState("");
  const technologies = Object.keys(value).sort();
  const years = [
    ...new Set(
      Object.values(value)
        .flatMap((raw) => Object.keys(objectValue(raw)))
        .filter(isYearKey),
    ),
  ].sort((left, right) => Number(left) - Number(right));
  const setCell = (technology: string, year: string, raw: string) =>
    onChange({ ...value, [technology]: { ...objectValue(value[technology]), [year]: inputNumber(raw) } });
  const removeTechnology = (technology: string) => {
    const next = { ...value };
    delete next[technology];
    onChange(next);
  };
  const addTechnology = () => {
    const name = newTechnology.trim();
    if (!name || Object.hasOwn(value, name)) return;
    onChange({ ...value, [name]: Object.fromEntries(years.map((year) => [year, null])) });
    setNewTechnology("");
  };
  const addYear = () => {
    const year = newYear.trim();
    if (!isYearKey(year) || years.includes(year) || !technologies.length) return;
    onChange(
      Object.fromEntries(
        Object.entries(value).map(([technology, raw]) => [technology, { ...objectValue(raw), [year]: null }]),
      ),
    );
    setNewYear("");
  };
  return (
    <div className={[scenarioConfigControlsStyles["config-entry-block"], styles["matrix-entry-block"]].join(" ")}>
      <p className={scenarioConfigControlsStyles["field-label"]}>{label}</p>
      {technologies.length ? (
        <div className={scenarioConfigControlsStyles["config-table-wrap"]}>
          <table
            className={[
              dataTableStyles["table"],
              [scenarioConfigControlsStyles["config-entry-table"], styles["matrix-table"]].join(" "),
            ].join(" ")}
          >
            <thead>
              <tr>
                <th>Technology</th>
                {years.map((year) => (
                  <th key={year}>{year}</th>
                ))}
                <th>
                  <span className={workspaceUtilitiesStyles["sr-only"]}>Actions</span>
                </th>
              </tr>
            </thead>
            <tbody>
              {technologies.map((technology) => {
                const row = objectValue(value[technology]);
                return (
                  <tr key={technology}>
                    <th scope="row">{technology}</th>
                    {years.map((year) => (
                      <td key={year}>
                        <input
                          aria-label={`${technology} ${year} value`}
                          type="number"
                          step="any"
                          value={String(row[year] ?? "")}
                          onChange={(event) => setCell(technology, year, event.target.value)}
                        />
                      </td>
                    ))}
                    <td>
                      <IconButton
                        tone="danger"
                        aria-label={`Remove ${technology}`}
                        onClick={() => removeTechnology(technology)}
                      >
                        <Trash2 aria-hidden="true" />
                      </IconButton>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      ) : (
        <p className={styles["constraint-empty"]}>Add a technology to start this table.</p>
      )}
      <div className={styles["matrix-add-row"]}>
        <div className={scenarioConfigControlsStyles["config-table-add"]}>
          <Field>
            <span>New technology</span>
            <input
              value={newTechnology}
              placeholder="Technology"
              onChange={(event) => setNewTechnology(event.target.value)}
            />
          </Field>
          <Button
            disabled={!newTechnology.trim() || Object.hasOwn(value, newTechnology.trim())}
            onClick={addTechnology}
          >
            <Plus aria-hidden="true" />
            Add
          </Button>
        </div>
        <div className={scenarioConfigControlsStyles["config-table-add"]}>
          <Field>
            <span>New year</span>
            <input
              type="number"
              value={newYear}
              placeholder="2035"
              onChange={(event) => setNewYear(event.target.value)}
            />
          </Field>
          <Button
            disabled={!technologies.length || !isYearKey(newYear.trim()) || years.includes(newYear.trim())}
            onClick={addYear}
          >
            <Plus aria-hidden="true" />
            Add
          </Button>
        </div>
      </div>
    </div>
  );
}
