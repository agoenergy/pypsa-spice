// SPDX-FileCopyrightText: PyPSA-SPICE Developers
// SPDX-License-Identifier: GPL-2.0-or-later

import styles from "./RunModel.module.scss";
import workspaceFeedbackStyles from "./WorkspaceFeedback.module.scss";
import { useEffect, useId, useMemo, useState } from "react";

import { saveScenarioConfigSection } from "../api";
import { SelectField } from "./FormControls";
import SaveDiscardActions from "./SaveDiscardActions";
import { MappingTable } from "../pages/ScenarioConfigControls";
import { objectValue } from "../pages/ScenarioConfigEditorUtils";
import type { InputSelection, ScenarioConfigResponse } from "../types";
import { setEditorDirty } from "../utility";

export default function SolverSettings({
  selection,
  config,
  disabled,
  onSaved,
  onDirtyChange,
}: {
  selection: InputSelection;
  config: ScenarioConfigResponse;
  disabled: boolean;
  onSaved: (config: ScenarioConfigResponse) => void;
  onDirtyChange: (dirty: boolean) => void;
}) {
  const editorId = useId();
  const original = useMemo(() => structuredClone(config.sections.solving || {}), [config]);
  const [draft, setDraft] = useState<Record<string, unknown>>(original);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [status, setStatus] = useState("");
  useEffect(() => setDraft(original), [original]);

  const dirty = JSON.stringify(draft) !== JSON.stringify(original);
  useEffect(() => {
    onDirtyChange(dirty);
    setEditorDirty(editorId, dirty);
    return () => setEditorDirty(editorId, false);
  }, [editorId, dirty]);

  const solver = objectValue(draft.solver);
  const profiles = objectValue(draft.solver_options);
  const solverName = String(solver.name || "");
  const profileName = String(solver.options || "");
  const profile = objectValue(profiles[profileName]);

  // Both lists follow the "solving" section: solver.name plus the solvers named by
  // solver_options profiles ("gurobi-default" or an explicit "name"), and every profile.
  const profileSolver = (key: string) => String(objectValue(profiles[key]).name || key.split("-")[0]);
  const solverNames = [
    ...new Set([
      solverName,
      ...Object.keys(profiles)
        .filter((key) => key !== "default")
        .map(profileSolver),
    ]),
  ].filter(Boolean);
  const profileFits = (key: string, name: string) => key === "default" || profileSolver(key) === name;
  const profileNames = Object.keys(profiles);
  const validationError = Object.values(profile).some((value) => value === null)
    ? "Enter a value for every solver setting."
    : "";

  const updateSolver = (next: Record<string, unknown>) =>
    setDraft((current) => ({ ...current, solver: { ...objectValue(current.solver), ...next } }));
  const chooseSolver = (name: string) => {
    // Keep the profile if it suits the new solver, otherwise prefer "<solver>-default".
    const fitting = Object.keys(profiles).filter((key) => profileFits(key, name));
    const options = fitting.includes(profileName)
      ? profileName
      : fitting.find((key) => key === `${name}-default`) || fitting.find((key) => key !== "default") || "default";
    updateSolver({ name, options });
  };

  const save = async () => {
    setSaving(true);
    setError("");
    setStatus("");
    try {
      onSaved(await saveScenarioConfigSection(selection, "solving", config.revision, draft));
      setStatus("Solver settings saved.");
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not save the solver settings.");
    } finally {
      setSaving(false);
    }
  };

  return (
    <fieldset className={styles["solver-settings"]} disabled={disabled}>
      <legend>Solver settings</legend>
      <div className={styles["solver-choice"]}>
        <SelectField
          label="Solver"
          value={solverName}
          onChange={chooseSolver}
          options={solverNames.map((name) => ({ value: name, label: name }))}
        />
        <SelectField
          label="Solver options"
          value={profileName}
          onChange={(options) => updateSolver({ options })}
          options={profileNames.map((name) => ({ value: name, label: name }))}
        />
      </div>
      {Object.hasOwn(profiles, profileName) ? (
        <MappingTable
          label="Setting"
          value={profile}
          onChange={(next) =>
            setDraft((current) => ({
              ...current,
              solver_options: { ...objectValue(current.solver_options), [profileName]: next },
            }))
          }
        />
      ) : (
        <div className={[workspaceFeedbackStyles["notice"], workspaceFeedbackStyles["error"]].join(" ")}>
          The options profile "{profileName}" is not defined in solver_options.
        </div>
      )}
      {(error || (dirty && validationError)) && (
        <div className={[workspaceFeedbackStyles["notice"], workspaceFeedbackStyles["error"]].join(" ")}>
          {error || validationError}
        </div>
      )}
      <SaveDiscardActions
        hasChanges={dirty}
        saving={saving}
        onDiscard={() => setDraft(original)}
        onSave={() => void save()}
        saveDisabled={Boolean(validationError)}
        saveLabel="Save solver settings"
        status={status}
      />
    </fieldset>
  );
}
