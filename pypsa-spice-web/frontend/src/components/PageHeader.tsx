import styles from "./PageHeader.module.scss";
import type { ReactNode } from "react";

export default function PageHeader({
  title,
  children,
  className = "",
}: {
  title: string;
  children?: ReactNode;
  className?: string;
}) {
  return (
    <section className={[styles["page-title"], className].filter(Boolean).join(" ")}>
      <h1>{title}</h1>
      {children}
    </section>
  );
}
