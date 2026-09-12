import scenarioConfigLayoutStyles from "./ScenarioConfigLayout.module.scss";
import editorPanelStyles from "../components/EditorPanel.module.scss";
import styles from "./ScenarioSettings.module.scss";
import IconButton from "../components/IconButton";
import { Field, SelectField, ToggleField } from "../components/FormControls";
import { useState } from "react";
import { Trash2 } from "lucide-react";

export default function ScenarioSettings({
  value,
  country,
  onChange,
}: {
  value: Record<string, unknown>;
  country: string;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const snapshots = (value.snapshots || {}) as Record<string, unknown>;
  const resolution = (value.resolution || {}) as Record<string, unknown>;
  const interest = (value.interest || {}) as Record<string, unknown>;
  const start = String(snapshots.start || "");
  const end = String(snapshots.end || "");
  const modelYear = Number(start.slice(0, 4)) || new Date().getFullYear();
  const [manual, setManual] = useState(
    !start.endsWith("-01-01") || end !== `${modelYear + 1}-01-01` || snapshots.inclusive !== "left",
  );
  const patch = (key: string, child: Record<string, unknown>) => onChange({ ...value, [key]: child });
  const toggleManual = (checked: boolean) => {
    setManual(checked);
    if (!checked) {
      patch("snapshots", { start: `${modelYear}-01-01`, end: `${modelYear + 1}-01-01`, inclusive: "left" });
    }
  };

  return (
    <div className={scenarioConfigLayoutStyles["config-form"]}>
      <div className={scenarioConfigLayoutStyles["form-grid"]}>
        <Field>
          <span>Model year</span>
          <input
            type="number"
            min="2010"
            max="3000"
            value={modelYear}
            onChange={(event) => {
              const year = Number(event.target.value);
              patch("snapshots", {
                ...snapshots,
                start: `${year}-01-01`,
                end: `${year + 1}-01-01`,
                inclusive: "left",
              });
            }}
          />
        </Field>
        <Field>
          <span>Remove assets below (MW)</span>
          <input
            type="number"
            min="0"
            step="0.1"
            value={String(value.remove_threshold ?? 0)}
            onChange={(event) => onChange({ ...value, remove_threshold: Number(event.target.value) })}
          />
        </Field>
      </div>
      <ToggleField label="Edit snapshot range manually" checked={manual} onChange={toggleManual} />
      {manual ? (
        <div className={[scenarioConfigLayoutStyles["form-grid"], scenarioConfigLayoutStyles["three"]].join(" ")}>
          <Field>
            <span>Snapshot start</span>
            <input
              type="date"
              value={start}
              onChange={(event) => patch("snapshots", { ...snapshots, start: event.target.value })}
            />
          </Field>
          <Field>
            <span>Snapshot end</span>
            <input
              type="date"
              value={end}
              onChange={(event) => patch("snapshots", { ...snapshots, end: event.target.value })}
            />
          </Field>
          <SelectField
            label="Inclusive"
            value={String(snapshots.inclusive || "left")}
            onChange={(inclusive) => patch("snapshots", { ...snapshots, inclusive })}
            options={["both", "neither", "left", "right"].map((option) => ({ value: option, label: option }))}
          />
        </div>
      ) : (
        <p className={editorPanelStyles["field-help"]}>
          Full-year hourly range: {modelYear}-01-01 to {modelYear + 1}-01-01, inclusive left.
        </p>
      )}
      <div className={scenarioConfigLayoutStyles["form-section"]}>
        <h3>Temporal resolution</h3>
        <div className={scenarioConfigLayoutStyles["form-grid"]}>
          <SelectField
            label="Method"
            value={String(resolution.method || "nth_hour")}
            onChange={(method) => patch("resolution", { ...resolution, method })}
            options={[
              { value: "nth_hour", label: "Every nth hour" },
              { value: "clustered", label: "Clustered representative days" },
            ]}
          />
          {resolution.method === "clustered" ? (
            <Field>
              <span>Number of days</span>
              <input
                type="number"
                min="1"
                value={String(resolution.number_of_days ?? 3)}
                onChange={(event) => patch("resolution", { ...resolution, number_of_days: Number(event.target.value) })}
              />
            </Field>
          ) : (
            <Field>
              <span>Step size</span>
              <input
                type="number"
                min="1"
                value={String(resolution.stepsize ?? 25)}
                onChange={(event) => patch("resolution", { ...resolution, stepsize: Number(event.target.value) })}
              />
            </Field>
          )}
        </div>
      </div>
      <div className={scenarioConfigLayoutStyles["form-section"]}>
        <h3>Interest rates</h3>
        <p className={editorPanelStyles["field-help"]}>Country-specific decimal rates; 0.05 means 5%.</p>
        <CountryValueEditor
          value={interest}
          country={country}
          onChange={(next) => onChange({ ...value, interest: next })}
        />
      </div>
    </div>
  );
}

function CountryValueEditor({
  value,
  country,
  onChange,
}: {
  value: Record<string, unknown>;
  country: string;
  onChange: (value: Record<string, unknown>) => void;
}) {
  if (country === "ALL") return <KeyValueEditor value={value} onChange={onChange} />;
  return (
    <div className={scenarioConfigLayoutStyles["key-value-grid"]}>
      <Field>
        <span>{country}</span>
        <input
          type="number"
          step="any"
          value={String(value[country] ?? "")}
          onChange={(event) => onChange({ ...value, [country]: Number(event.target.value) })}
        />
      </Field>
    </div>
  );
}

function KeyValueEditor({
  value,
  onChange,
}: {
  value: Record<string, unknown>;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const remove = (key: string) => {
    const next = { ...value };
    delete next[key];
    onChange(next);
  };
  return (
    <div>
      <div className={scenarioConfigLayoutStyles["key-value-grid"]}>
        {Object.entries(value).map(([key, raw]) => (
          <div className={styles["mapping-field"]} key={key}>
            <Field>
              <span>{key}</span>
              <input
                type="number"
                step="any"
                value={String(raw ?? "")}
                onChange={(event) => onChange({ ...value, [key]: Number(event.target.value) })}
              />
            </Field>
            <IconButton tone="danger" aria-label={`Remove ${key}`} onClick={() => remove(key)}>
              <Trash2 aria-hidden="true" />
            </IconButton>
          </div>
        ))}
      </div>
    </div>
  );
}
