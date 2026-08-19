import { useEffect, useState } from "react";
import { Copy, X } from "lucide-react";
import "./NewScenarioDialog.css";
import { createScenario } from "../api";
import type { CreatedScenario, InputSelection } from "../types";

export default function NewScenarioDialog({ selection, scenarios, onClose, onCreated }: {
  selection: InputSelection;
  scenarios: string[];
  onClose: () => void;
  onCreated: (scenario: CreatedScenario) => void;
}) {
  const [name, setName] = useState("");
  const [source, setSource] = useState(selection.scenario);
  const [submitting, setSubmitting] = useState(false);
  const [error, setError] = useState("");
  const trimmedName = name.trim();
  const valid = Boolean(trimmedName && !trimmedName.startsWith(".") && !/[\\/\0]/.test(trimmedName) && trimmedName !== "global_input");

  useEffect(() => {
    const close = (event: KeyboardEvent) => { if (event.key === "Escape" && !submitting) onClose(); };
    window.addEventListener("keydown", close);
    return () => window.removeEventListener("keydown", close);
  }, [onClose, submitting]);

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!valid) return;
    setSubmitting(true); setError("");
    try {
      onCreated(await createScenario({
        dataset: selection.dataset,
        project: selection.project,
        source_scenario: source,
        new_scenario: trimmedName,
      }));
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not create the scenario.");
    } finally {
      setSubmitting(false);
    }
  };

  return <div className="dialog-backdrop" role="presentation" onMouseDown={(event) => { if (event.target === event.currentTarget && !submitting) onClose(); }}>
    <section className="dialog scenario-dialog" role="dialog" aria-modal="true" aria-labelledby="new-scenario-title">
      <header><h2 id="new-scenario-title">Create new scenario</h2><button className="icon-button" onClick={onClose} disabled={submitting} aria-label="Close"><X aria-hidden="true" /></button></header>
      <form onSubmit={submit}>
        <div className="scenario-dialog-body">
          <p>Duplicate a complete scenario as a new local starting point. Shared <code>global_input</code> files stay shared.</p>
          {error && <div className="notice error">{error}</div>}
          <label className="field"><span>New scenario name</span><input autoFocus value={name} onChange={(event) => setName(event.target.value)} maxLength={255} placeholder="for example: policy_high" aria-invalid={Boolean(name && !valid)} /></label>
          <label className="field"><span>Source scenario</span><select value={source} onChange={(event) => setSource(event.target.value)}>{scenarios.map((scenario) => <option key={scenario}>{scenario}</option>)}</select></label>
          <div className="scenario-path-preview"><span>Local destination</span><code>data/{selection.dataset}/{selection.project}/input/{trimmedName || "new_scenario"}</code></div>
        </div>
        <footer><button className="button secondary" type="button" onClick={onClose} disabled={submitting}>Cancel</button><button className="button primary" type="submit" disabled={!valid || submitting}><Copy aria-hidden="true" />{submitting ? "Creating…" : "Create scenario"}</button></footer>
      </form>
    </section>
  </div>;
}
