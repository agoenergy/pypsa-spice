import { useCallback, useRef, useState } from "react";
import { List, X } from "lucide-react";

import IconButton from "./IconButton";
import { useDismissOnEscapeOrOutside } from "../useDismiss";
import styles from "./ResultsToc.module.scss";

export interface TocEntry {
  key: string;
  href: string;
  label: string;
}

interface Props {
  /** DOM id of the panel, referenced by the trigger's aria-controls. */
  id: string;
  /** Panel heading, for example "Figures". */
  heading: string;
  /** Accessible name of the panel, for example "Figures on this page". */
  panelLabel: string;
  /** Noun used in the trigger and close labels, for example "figure list". */
  listName: string;
  entries: TocEntry[];
}

/** Floating jump-to list shared by the results, configuration, and comparison pages. */
export default function ResultsToc({ id, heading, panelLabel, listName, entries }: Props) {
  const [open, setOpen] = useState(false);
  const root = useRef<HTMLDivElement>(null);
  const close = useCallback(() => setOpen(false), []);
  useDismissOnEscapeOrOutside(open, root, close);

  if (!entries.length) return null;
  return (
    <div className={[styles["results-toc"], open ? styles["open"] : ""].filter(Boolean).join(" ")} ref={root}>
      {open && (
        <nav className={styles["results-toc-panel"]} id={id} aria-label={panelLabel}>
          <header>
            <h2>{heading}</h2>
            <IconButton alignEnd onClick={close} aria-label={`Close ${listName}`}>
              <X aria-hidden="true" />
            </IconButton>
          </header>
          <ol>
            {entries.map((entry, index) => (
              <li key={entry.key}>
                <a href={entry.href} onClick={close}>
                  <span>{String(index + 1).padStart(2, "0")}</span>
                  <b>{entry.label}</b>
                </a>
              </li>
            ))}
          </ol>
        </nav>
      )}
      <IconButton
        className={styles["results-toc-trigger"]}
        onClick={() => setOpen((current) => !current)}
        aria-label={`Open ${listName}`}
        aria-expanded={open}
        aria-controls={id}
      >
        <List aria-hidden="true" />
      </IconButton>
    </div>
  );
}
