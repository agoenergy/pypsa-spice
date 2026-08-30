import type { ReactNode } from "react";
import { ArrowRight } from "lucide-react";
import styles from "./WorkspaceActionCard.module.scss";

interface WorkspaceActionCardProps {
  icon: ReactNode;
  title: string;
  description: string;
  onClick?: () => void;
  disabled?: boolean;
}

export default function WorkspaceActionCard({
  icon,
  title,
  description,
  onClick,
  disabled = false,
}: WorkspaceActionCardProps) {
  return <li className={styles["action-card"]}>
    <button type="button" className={styles["action-button"]} onClick={onClick} disabled={disabled}>
      <span className={styles["action-icon"]}>{icon}</span>
      <span className={styles["action-copy"]}>
        <b>{title}{disabled && <em>Not implemented</em>}</b>
        <small>{description}</small>
      </span>
      {!disabled && <ArrowRight aria-hidden="true" />}
    </button>
  </li>;
}
