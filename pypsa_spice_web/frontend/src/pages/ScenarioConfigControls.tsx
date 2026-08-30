import { useState } from "react";
import { Plus, Trash2 } from "lucide-react";
import { Button } from "../components/FormControls";
import { inputNumber } from "./ScenarioConfigEditorUtils";

export function MappingTable({ label, value, yearKeys = false, allowAdd = true, onChange }: {
  label: string;
  value: Record<string, unknown>;
  yearKeys?: boolean;
  allowAdd?: boolean;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const [newKey, setNewKey] = useState("");
  const entries = Object.entries(value).sort(([left], [right]) => yearKeys
    ? Number(left) - Number(right)
    : left.localeCompare(right));
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
  return <div className="config-entry-block">
    <div className="config-table-wrap"><table className="config-entry-table"><thead><tr>
      <th>{label}</th><th>Value</th><th><span className="sr-only">Actions</span></th>
    </tr></thead><tbody>{entries.map(([key, raw]) => <tr key={key}>
      <th scope="row">{key}</th>
      <td><input
        aria-label={`${label} ${key} value`}
        type={typeof raw === "string" ? "text" : "number"}
        step="any"
        value={String(raw ?? "")}
        onChange={(event) => onChange({
          ...value,
          [key]: typeof raw === "string" ? event.target.value : inputNumber(event.target.value),
        })}
      /></td>
      <td><button className="icon-button" aria-label={`Remove ${key}`} onClick={() => remove(key)}><Trash2 aria-hidden="true" /></button></td>
    </tr>)}</tbody></table></div>
    {allowAdd && <div className="config-table-add">
      <label className="field"><span>New {label.toLowerCase()}</span><input
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
      /></label>
      <Button disabled={!newKey.trim() || Object.hasOwn(value, newKey.trim())} onClick={add}><Plus aria-hidden="true" />Add</Button>
    </div>}
  </div>;
}
