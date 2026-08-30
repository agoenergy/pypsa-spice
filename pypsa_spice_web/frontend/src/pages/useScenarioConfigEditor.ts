import { useEffect, useId, useMemo, useState } from "react";
import {
  getFuelSupplyTable,
  getScenarioConfig,
  saveFuelSupplyLimits,
  saveScenarioConfigSection,
  saveScenarioConfigSections,
} from "../api";
import type { InputRow, InputSelection, InputTableResponse, ScenarioConfigResponse } from "../types";
import { setEditorDirty } from "../utility";
import {
  COMBINED_SECTION,
  draftForSection,
  fuelConstraintSnapshot,
  validateFuelLimits,
  validateSection,
} from "./ScenarioConfigEditorUtils";

export default function useScenarioConfigEditor(selection: InputSelection, section: string) {
  const editorId = useId();
  const [config, setConfig] = useState<ScenarioConfigResponse | null>(null);
  const [fuelTable, setFuelTable] = useState<InputTableResponse | null>(null);
  const [fuelRows, setFuelRows] = useState<InputRow[]>([]);
  const [fuelError, setFuelError] = useState("");
  const [draft, setDraft] = useState<Record<string, unknown>>({});
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [success, setSuccess] = useState("");

  const load = async (signal?: AbortSignal) => {
    setLoading(true);
    setError("");
    setFuelError("");
    setSuccess("");
    try {
      const [data, fuelResult] = await Promise.all([
        getScenarioConfig(selection, signal),
        getFuelSupplyTable(selection, signal)
          .then((value) => ({ value, error: "" }))
          .catch((reason) => ({
            value: null,
            error: reason instanceof Error ? reason.message : "Could not read fuel supply limits.",
          })),
      ]);
      setConfig(data);
      setDraft(draftForSection(data, section));
      setFuelTable(fuelResult.value);
      setFuelRows(structuredClone(fuelResult.value?.rows || []));
      setFuelError(fuelResult.error);
    } catch (reason) {
      if (!(reason instanceof DOMException && reason.name === "AbortError")) {
        setError(reason instanceof Error ? reason.message : "Could not load the scenario config.");
      }
    } finally {
      if (!signal?.aborted) setLoading(false);
    }
  };

  useEffect(() => {
    const controller = new AbortController();
    void load(controller.signal);
    return () => controller.abort();
  }, [selection.dataset, selection.project, selection.scenario]);

  useEffect(() => {
    if (!config) return;
    setDraft(draftForSection(config, section));
    setFuelRows(structuredClone(fuelTable?.rows || []));
    setSuccess("");
    setError("");
  }, [config, fuelTable, section]);

  const original = useMemo(() => config ? draftForSection(config, section) : {}, [config, section]);
  const configDirty = useMemo(() => JSON.stringify(draft) !== JSON.stringify(original), [draft, original]);
  const fuelChanges = useMemo(() => {
    if (section !== COMBINED_SECTION || !fuelTable) return [];
    const originals = new Map(fuelTable.rows.map((row) => [row.__row_id, row.max_supply__mwh_year]));
    return fuelRows.flatMap((row) => JSON.stringify(row.max_supply__mwh_year) === JSON.stringify(originals.get(row.__row_id))
      ? []
      : [{ row: row.__row_id, column: "max_supply__mwh_year", value: row.max_supply__mwh_year }]);
  }, [fuelRows, fuelTable, section]);
  const fuelDirty = fuelChanges.length > 0;
  const dirty = configDirty || fuelDirty;
  const fuelConstraintDirty = section === COMBINED_SECTION
    && JSON.stringify(fuelConstraintSnapshot(draft)) !== JSON.stringify(fuelConstraintSnapshot(original));
  const validationError = useMemo(
    () => validateSection(section, draft)
      || (fuelError ? fuelConstraintDirty ? fuelError : "" : validateFuelLimits(section, draft, fuelRows)),
    [section, draft, fuelRows, fuelError, fuelConstraintDirty],
  );

  useEffect(() => {
    const warn = (event: BeforeUnloadEvent) => { if (dirty) event.preventDefault(); };
    window.addEventListener("beforeunload", warn);
    return () => window.removeEventListener("beforeunload", warn);
  }, [dirty]);
  useEffect(() => {
    setEditorDirty(editorId, dirty);
    return () => setEditorDirty(editorId, false);
  }, [editorId, dirty]);

  const save = async () => {
    if (!config || !dirty) return;
    setSaving(true);
    setError("");
    setSuccess("");
    try {
      if (fuelDirty && fuelTable) {
        const savedFuel = await saveFuelSupplyLimits(selection, fuelTable.revision, fuelChanges);
        setFuelTable({ ...fuelTable, revision: savedFuel.revision, rows: structuredClone(fuelRows) });
      }
      const changedSections = section === COMBINED_SECTION ? Object.fromEntries(
        ["co2_management", "custom_constraints"]
          .filter((name) => JSON.stringify(draft[name]) !== JSON.stringify(original[name]))
          .map((name) => [name, draft[name] as Record<string, unknown>]),
      ) : {};
      const data = configDirty
        ? section === COMBINED_SECTION
          ? await saveScenarioConfigSections(selection, config.revision, changedSections)
          : await saveScenarioConfigSection(selection, section, config.revision, draft)
        : config;
      setConfig(data);
      setDraft(draftForSection(data, section));
      setSuccess("Changes saved.");
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Could not save this section.");
    } finally {
      setSaving(false);
    }
  };

  const discard = () => {
    setDraft(structuredClone(original));
    setFuelRows(structuredClone(fuelTable?.rows || []));
  };

  return {
    config,
    draft,
    setDraft,
    fuelRows,
    setFuelRows,
    fuelError,
    loading,
    saving,
    error,
    success,
    dirty,
    validationError,
    reload: load,
    save,
    discard,
  };
}
