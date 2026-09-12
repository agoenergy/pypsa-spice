import styles from "./ScenarioConfigControls.module.scss";
import dataTableStyles from "../components/DataTable.module.scss";
import workspaceUtilitiesStyles from "../components/WorkspaceUtilities.module.scss";
import IconButton from "../components/IconButton";
import { Field } from "../components/FormControls";
import { useState } from "react";
import { Plus, Trash2 } from "lucide-react";
import Button from "../components/Button";
import { inputNumber } from "./ScenarioConfigEditorUtils";

export function MappingTable({
  label,
  value,
  yearKeys = false,
  allowAdd = true,
  onChange,
}: {
  label: string;
  value: Record<string, unknown>;
  yearKeys?: boolean;
  allowAdd?: boolean;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const [newKey, setNewKey] = useState("");
  const entries = Object.entries(value).sort(([left], [right]) =>
    yearKeys ? Number(left) - Number(right) : left.localeCompare(right),
  );
  const remove = (key: string) => {
    const next = { ...value };
    delete next[key];
    onChange(next);
  };
  const add = () => {
    const key = newKey.trim();
    if (!key || Object.hasOwn(value, key)) return;
    onChange({ ...value, [key]: 0 });
    setNewKey("");
  };
  return (
    <div className={styles["config-entry-block"]}>
      <div className={styles["config-table-wrap"]}>
        <table className={[dataTableStyles["table"], styles["config-entry-table"]].join(" ")}>
          <thead>
            <tr>
              <th>{label}</th>
              <th>Value</th>
              <th>
                <span className={workspaceUtilitiesStyles["sr-only"]}>Actions</span>
              </th>
            </tr>
          </thead>
          <tbody>
            {entries.map(([key, raw]) => (
              <tr key={key}>
                <th scope="row">{key}</th>
                <td>
                  <input
                    aria-label={`${label} ${key} value`}
                    type={typeof raw === "string" ? "text" : "number"}
                    step="any"
                    value={String(raw ?? "")}
                    onChange={(event) =>
                      onChange({
                        ...value,
                        [key]: typeof raw === "string" ? event.target.value : inputNumber(event.target.value),
                      })
                    }
                  />
                </td>
                <td>
                  <IconButton tone="danger" aria-label={`Remove ${key}`} onClick={() => remove(key)}>
                    <Trash2 aria-hidden="true" />
                  </IconButton>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {allowAdd && (
        <div className={styles["config-table-add"]}>
          <Field>
            <span>New {label.toLowerCase()}</span>
            <input
              type={yearKeys ? "number" : "text"}
              value={newKey}
              placeholder={yearKeys ? "2035" : "Name"}
              onChange={(event) => setNewKey(event.target.value)}
              onKeyDown={(event) => {
                if (event.key === "Enter") {
                  event.preventDefault();
                  add();
                }
              }}
            />
          </Field>
          <Button disabled={!newKey.trim() || Object.hasOwn(value, newKey.trim())} onClick={add}>
            <Plus aria-hidden="true" />
            Add
          </Button>
        </div>
      )}
    </div>
  );
}
