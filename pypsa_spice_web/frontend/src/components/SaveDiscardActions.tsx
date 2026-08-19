import { Check, RotateCcw, Save } from "lucide-react";
import "./SaveDiscardActions.css";

interface SaveDiscardActionsProps {
  hasChanges: boolean;
  saving: boolean;
  onDiscard: () => void;
  onSave: () => void;
  saveDisabled?: boolean;
  saveLabel?: string;
  status?: string;
  floating?: boolean;
  avoidSideControl?: boolean;
}

export default function SaveDiscardActions({
  hasChanges,
  saving,
  onDiscard,
  onSave,
  saveDisabled = false,
  saveLabel = "Save changes",
  status = "",
  floating = false,
  avoidSideControl = false,
}: SaveDiscardActionsProps) {
  const className = [
    "save-discard-actions",
    floating ? "floating-save-discard-actions" : "",
    avoidSideControl ? "avoid-side-control" : "",
  ].filter(Boolean).join(" ");

  return <div className={className} role="group" aria-label="Save or discard changes">
    <button type="button" className="button secondary" disabled={!hasChanges || saving} onClick={onDiscard}>
      <RotateCcw aria-hidden="true" />Discard
    </button>
    <button type="button" className="button primary" disabled={!hasChanges || saving || saveDisabled} onClick={onSave}>
      <Save aria-hidden="true" />{saving ? "Saving…" : saveLabel}
    </button>
    {status && <span className="save-discard-status" role="status"><Check aria-hidden="true" />{status}</span>}
  </div>;
}
