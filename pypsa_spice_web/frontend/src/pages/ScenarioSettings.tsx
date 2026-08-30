import { useState } from "react";
import { Trash2 } from "lucide-react";
import { SelectField, ToggleField } from "../components/FormControls";

export default function ScenarioSettings({ value, country, onChange }: {
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

  return <div className="config-form">
    <div className="form-grid">
      <label className="field"><span>Model year</span><input
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
      /></label>
      <label className="field"><span>Remove assets below (MW)</span><input
        type="number"
        min="0"
        step="0.1"
        value={String(value.remove_threshold ?? 0)}
        onChange={(event) => onChange({ ...value, remove_threshold: Number(event.target.value) })}
      /></label>
    </div>
    <ToggleField label="Edit snapshot range manually" checked={manual} onChange={toggleManual} />
    {manual ? <div className="form-grid three">
      <label className="field"><span>Snapshot start</span><input type="date" value={start} onChange={(event) => patch("snapshots", { ...snapshots, start: event.target.value })} /></label>
      <label className="field"><span>Snapshot end</span><input type="date" value={end} onChange={(event) => patch("snapshots", { ...snapshots, end: event.target.value })} /></label>
      <SelectField
        label="Inclusive"
        value={String(snapshots.inclusive || "left")}
        onChange={(inclusive) => patch("snapshots", { ...snapshots, inclusive })}
        options={["both", "neither", "left", "right"].map((option) => ({ value: option, label: option }))}
      />
    </div> : <p className="field-help">Full-year hourly range: {modelYear}-01-01 to {modelYear + 1}-01-01, inclusive left.</p>}
    <div className="form-section">
      <h3>Temporal resolution</h3>
      <div className="form-grid">
        <SelectField
          label="Method"
          value={String(resolution.method || "nth_hour")}
          onChange={(method) => patch("resolution", { ...resolution, method })}
          options={[
            { value: "nth_hour", label: "Every nth hour" },
            { value: "clustered", label: "Clustered representative days" },
          ]}
        />
        {resolution.method === "clustered" ? <label className="field"><span>Number of days</span><input
          type="number"
          min="1"
          value={String(resolution.number_of_days ?? 3)}
          onChange={(event) => patch("resolution", { ...resolution, number_of_days: Number(event.target.value) })}
        /></label> : <label className="field"><span>Step size</span><input
          type="number"
          min="1"
          value={String(resolution.stepsize ?? 25)}
          onChange={(event) => patch("resolution", { ...resolution, stepsize: Number(event.target.value) })}
        /></label>}
      </div>
    </div>
    <div className="form-section">
      <h3>Interest rates</h3>
      <p className="field-help">Country-specific decimal rates; 0.05 means 5%.</p>
      <CountryValueEditor value={interest} country={country} onChange={(next) => onChange({ ...value, interest: next })} />
    </div>
  </div>;
}

function CountryValueEditor({ value, country, onChange }: {
  value: Record<string, unknown>;
  country: string;
  onChange: (value: Record<string, unknown>) => void;
}) {
  if (country === "ALL") return <KeyValueEditor value={value} onChange={onChange} />;
  return <div className="key-value-grid"><label className="field"><span>{country}</span><input
    type="number"
    step="any"
    value={String(value[country] ?? "")}
    onChange={(event) => onChange({ ...value, [country]: Number(event.target.value) })}
  /></label></div>;
}

function KeyValueEditor({ value, onChange }: {
  value: Record<string, unknown>;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const remove = (key: string) => {
    const next = { ...value };
    delete next[key];
    onChange(next);
  };
  return <div className="mapping-editor"><div className="key-value-grid">{Object.entries(value).map(([key, raw]) => <div className="mapping-field" key={key}>
    <label className="field"><span>{key}</span><input
      type="number"
      step="any"
      value={String(raw ?? "")}
      onChange={(event) => onChange({ ...value, [key]: Number(event.target.value) })}
    /></label>
    <button className="icon-button" aria-label={`Remove ${key}`} onClick={() => remove(key)}><Trash2 aria-hidden="true" /></button>
  </div>)}</div></div>;
}
