import type { ReactNode } from "react";
import {
  BarChart3,
  CarFront,
  CircleDollarSign,
  Cloud,
  Factory,
  FileInput,
  GitCompareArrows,
  House,
  LayoutDashboard,
  Moon,
  Settings2,
  Sun,
  Zap,
} from "lucide-react";
import type { Catalog, ViewMode } from "../types";
import styles from "./Sidebar.module.scss";

interface SidebarProps {
  open: boolean;
  view: ViewMode;
  sections: Catalog["sections"];
  activeSectionId?: string;
  darkMode: boolean;
  modelRunLabel?: string;
  onSelectView: (view: ViewMode) => void;
  onSelectSection: (sectionId: string) => void;
  onToggleDarkMode: () => void;
}

interface SidebarSectionProps {
  active: boolean;
  href: string;
  icon: ReactNode;
  label: string;
  onSelect: () => void;
  children?: ReactNode;
}

const sectionIcons = {
  power: Zap,
  industry: Factory,
  transport: CarFront,
  emissions: Cloud,
  costs: CircleDollarSign,
};

export default function Sidebar({
  open,
  view,
  sections,
  activeSectionId,
  darkMode,
  modelRunLabel,
  onSelectView,
  onSelectSection,
  onToggleDarkMode,
}: SidebarProps) {
  const navigation = [
    { view: "home" as const, href: "?view=home", icon: <House aria-hidden="true" />, label: "Home" },
    { view: "inputs" as const, href: "?view=inputs", icon: <FileInput aria-hidden="true" />, label: "Inputs" },
    { view: "compare" as const, href: "?view=compare", icon: <GitCompareArrows aria-hidden="true" />, label: "Scenario differences" },
    { view: "configure" as const, href: "?view=configure", icon: <Settings2 aria-hidden="true" />, label: "Configure & run" },
  ];

  return <aside className={`${styles["sidebar"]} ${open ? styles["open"] : ""}`}>
    <div className={styles["brand"]}>
      <a href="?view=home" onClick={(event) => { event.preventDefault(); onSelectView("home"); }} aria-label="Go to home">
        <img src="/ui/pypsa-logo.svg" alt="PyPSA-SPICE" />
      </a>
    </div>
    <div className={styles["workspace-nav"]} role="navigation" aria-label="Workflow">
      {navigation.map((item) => <SidebarSection
        key={item.view}
        active={view === item.view}
        href={item.href}
        icon={item.icon}
        label={item.label}
        onSelect={() => onSelectView(item.view)}
      >
        {item.view === "inputs" && view === "inputs" && <div className={styles["submenu-slot"]} id="input-table-menu" />}
        {item.view === "configure" && view === "configure" && <div className={styles["submenu-slot"]} id="config-section-tabs" />}
      </SidebarSection>)}
      <ResultsSidebarSection
        active={view === "outputs"}
        activeSectionId={activeSectionId}
        sections={sections}
        onSelect={() => onSelectView("outputs")}
        onSelectSection={onSelectSection}
      />
      <SidebarSection
        active={view === "dashboard"}
        href="?view=dashboard"
        icon={<LayoutDashboard aria-hidden="true" />}
        label="Dashboards"
        onSelect={() => onSelectView("dashboard")}
      />
    </div>
    <SidebarFooter darkMode={darkMode} modelRunLabel={modelRunLabel} onToggleDarkMode={onToggleDarkMode} />
  </aside>;
}

function SidebarSection({ active, href, icon, label, onSelect, children }: SidebarSectionProps) {
  return <div className={styles["section"]}>
    <a
      className={`${styles["primary"]} ${active ? styles["active"] : ""}`}
      href={href}
      onClick={(event) => { event.preventDefault(); onSelect(); }}
      aria-current={active ? "page" : undefined}
    >
      {icon}
      {label}
    </a>
    {children}
  </div>;
}

function ResultsSidebarSection({ active, activeSectionId, sections, onSelect, onSelectSection }: {
  active: boolean;
  activeSectionId?: string;
  sections: Catalog["sections"];
  onSelect: () => void;
  onSelectSection: (sectionId: string) => void;
}) {
  return <SidebarSection
    active={active}
    href={`?section=${activeSectionId || "power"}`}
    icon={<BarChart3 aria-hidden="true" />}
    label="Results"
    onSelect={onSelect}
  >
    {active && activeSectionId && <nav className={styles["submenu-list"]} aria-label="Result pages">
      {sections.map((section) => {
        const SectionIcon = sectionIcons[section.id as keyof typeof sectionIcons] || Zap;
        const sectionActive = section.id === activeSectionId;
        return <a
          className={`${styles["submenu-item"]} ${sectionActive ? styles["active"] : ""}`}
          href={`?section=${section.id}`}
          key={section.id}
          aria-current={sectionActive ? "page" : undefined}
          onClick={(event) => { event.preventDefault(); onSelectSection(section.id); }}
        >
          <SectionIcon aria-hidden="true" />
          <b>{section.label}</b>
          <small>{section.charts.length}</small>
        </a>;
      })}
    </nav>}
  </SidebarSection>;
}

function SidebarFooter({ darkMode, modelRunLabel, onToggleDarkMode }: {
  darkMode: boolean;
  modelRunLabel?: string;
  onToggleDarkMode: () => void;
}) {
  return <div className={styles["footer"]}>
    <div className={styles["statuses"]}>
      <span className={styles["status-indicator"]} data-label="Local files connected" aria-label="Local files connected" tabIndex={0}>
        <i className={styles["status-dot"]} aria-hidden="true" />
      </span>
      {modelRunLabel && <span
        className={`${styles["status-indicator"]} ${styles["model-run"]}`}
        data-label={modelRunLabel}
        aria-label={modelRunLabel}
        role="status"
        aria-live="polite"
        tabIndex={0}
      >
        <i className={styles["status-dot"]} aria-hidden="true" />
      </span>}
    </div>
    <a href="/docs" target="_blank">API</a>
    <button onClick={onToggleDarkMode} aria-label="Toggle dark mode">
      {darkMode ? <Sun aria-hidden="true" /> : <Moon aria-hidden="true" />}
    </button>
  </div>;
}
