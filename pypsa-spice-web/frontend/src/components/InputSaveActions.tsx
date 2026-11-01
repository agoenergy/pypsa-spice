import { createContext, useCallback, useContext, useRef, useState, type ReactNode } from "react";
import SaveDiscardActions from "./SaveDiscardActions";

interface InputEditorActions {
  changeCount: number;
  save: () => Promise<boolean>;
  discard: () => void;
}

type RegisterEditor = (id: string, actions: InputEditorActions) => () => void;
const InputSaveContext = createContext<RegisterEditor | null>(null);

export function useInputSaveActions() {
  const register = useContext(InputSaveContext);
  if (!register) throw new Error("Input tables must be inside InputSaveActions.");
  return register;
}

export default function InputSaveActions({ children }: { children: ReactNode }) {
  const [editors, setEditors] = useState(new Map<string, InputEditorActions>());
  const [saving, setSaving] = useState(false);
  const savingRef = useRef(false);
  const [status, setStatus] = useState("");
  const register = useCallback<RegisterEditor>((id, actions) => {
    setEditors((current) => new Map(current).set(id, actions));
    return () => {
      setEditors((current) => {
        const next = new Map(current);
        next.delete(id);
        return next;
      });
    };
  }, []);
  const pending = [...editors.values()].filter((editor) => editor.changeCount > 0);
  const changeCount = pending.reduce((total, editor) => total + editor.changeCount, 0);

  const save = async () => {
    if (savingRef.current || !pending.length) return;
    savingRef.current = true;
    setSaving(true);
    setStatus("");
    try {
      const results = await Promise.allSettled(pending.map((editor) => editor.save()));
      if (results.every((result) => result.status === "fulfilled" && result.value)) {
        setStatus(`Saved ${changeCount} ${changeCount === 1 ? "cell" : "cells"}.`);
      }
    } finally {
      savingRef.current = false;
      setSaving(false);
    }
  };

  return (
    <InputSaveContext.Provider value={register}>
      {children}
      <SaveDiscardActions
        floating
        hasChanges={changeCount > 0}
        saving={saving}
        saveLabel={`Save changes${changeCount ? ` (${changeCount})` : ""}`}
        status={changeCount ? "" : status}
        onSave={() => void save()}
        onDiscard={() => {
          pending.forEach((editor) => editor.discard());
          setStatus("");
        }}
      />
    </InputSaveContext.Provider>
  );
}
