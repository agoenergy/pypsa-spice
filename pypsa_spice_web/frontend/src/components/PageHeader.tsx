import type { ReactNode } from "react";
import "./PageHeader.css";

export default function PageHeader({ title, children, className = "" }: { title: string; children?: ReactNode; className?: string }) {
  return <section className={`page-title ${className}`.trim()}>
    <h1>{title}</h1>
    {children}
  </section>;
}
