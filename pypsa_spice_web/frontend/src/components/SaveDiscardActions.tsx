import styles from "./SaveDiscardActions.module.scss";
import Button from "./Button";
import { Check, RotateCcw, Save } from "lucide-react";

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
    styles["save-discard-actions"],
    floating ? styles["floating-save-discard-actions"] : "",
    avoidSideControl ? styles["avoid-side-control"] : "",
  ]
    .filter(Boolean)
    .join(" ");

  return (
    <div className={className} role="group" aria-label="Save or discard changes">
      <Button type="button" disabled={!hasChanges || saving} onClick={onDiscard}>
        <RotateCcw aria-hidden="true" />
        Discard
      </Button>
      <Button variant="primary" type="button" disabled={!hasChanges || saving || saveDisabled} onClick={onSave}>
        <Save aria-hidden="true" />
        {saving ? "Saving…" : saveLabel}
      </Button>
      {status && (
        <span className={styles["save-discard-status"]} role="status">
          <Check aria-hidden="true" />
          {status}
        </span>
      )}
    </div>
  );
}
