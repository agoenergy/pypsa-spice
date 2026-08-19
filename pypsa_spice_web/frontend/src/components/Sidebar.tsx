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
import type { Catalog } from "../types";
import type { ViewMode } from "../utility";

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
  return <aside className={`sidebar ${open ? "open" : ""}`}>
    <div className="brand">
      <a
        href="?view=home"
        onClick={(event) => { event.preventDefault(); onSelectView("home"); }}
        aria-label="Go to home"
      >
        <img src="/brand/pypsa-logo.svg" alt="PyPSA-SPICE" />
      </a>
    </div>
    <div className="workspace-nav" role="navigation" aria-label="Workflow">
      <HomeSidebarSection active={view === "home"} onSelect={() => onSelectView("home")} />
      <InputsSidebarSection active={view === "inputs"} onSelect={() => onSelectView("inputs")} />
      <ScenarioDifferencesSidebarSection active={view === "compare"} onSelect={() => onSelectView("compare")} />
      <ConfigureSidebarSection active={view === "configure"} onSelect={() => onSelectView("configure")} />
      <ResultsSidebarSection
        active={view === "outputs"}
        activeSectionId={activeSectionId}
        sections={sections}
        onSelect={() => onSelectView("outputs")}
        onSelectSection={onSelectSection}
      />
      <DashboardsSidebarSection active={view === "dashboard"} onSelect={() => onSelectView("dashboard")} />
    </div>
    <SidebarFooter darkMode={darkMode} modelRunLabel={modelRunLabel} onToggleDarkMode={onToggleDarkMode} />
  </aside>;
}

function SidebarSection({ active, href, icon, label, onSelect, children }: SidebarSectionProps) {
  return <div className="sidebar-section">
    <a
      className={`sidebar-primary ${active ? "active" : ""}`}
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

function HomeSidebarSection({ active, onSelect }: SectionNavigationProps) {
  return <SidebarSection active={active} href="?view=home" icon={<House aria-hidden="true" />} label="Home" onSelect={onSelect} />;
}

function InputsSidebarSection({ active, onSelect }: SectionNavigationProps) {
  return <SidebarSection active={active} href="?view=inputs" icon={<FileInput aria-hidden="true" />} label="Inputs" onSelect={onSelect}>
    {active && <div className="sidebar-submenu-slot" id="input-table-menu" />}
  </SidebarSection>;
}

function ScenarioDifferencesSidebarSection({ active, onSelect }: SectionNavigationProps) {
  return <SidebarSection active={active} href="?view=compare" icon={<GitCompareArrows aria-hidden="true" />} label="Scenario differences" onSelect={onSelect} />;
}

function ConfigureSidebarSection({ active, onSelect }: SectionNavigationProps) {
  return <SidebarSection active={active} href="?view=configure" icon={<Settings2 aria-hidden="true" />} label="Configure & run" onSelect={onSelect}>
    {active && <div className="sidebar-submenu-slot" id="config-section-tabs" />}
  </SidebarSection>;
}

interface ResultsSidebarSectionProps extends SectionNavigationProps {
  activeSectionId?: string;
  sections: Catalog["sections"];
  onSelectSection: (sectionId: string) => void;
}

function ResultsSidebarSection({ active, activeSectionId, sections, onSelect, onSelectSection }: ResultsSidebarSectionProps) {
  return <SidebarSection
    active={active}
    href={`?section=${activeSectionId || "power"}`}
    icon={<BarChart3 aria-hidden="true" />}
    label="Results"
    onSelect={onSelect}
  >
    {active && activeSectionId && <nav className="sidebar-submenu-list" aria-label="Result pages">
      {sections.map((section) => {
        const SectionIcon = sectionIcons[section.id as keyof typeof sectionIcons] || Zap;
        const sectionActive = section.id === activeSectionId;
        return <a
          className={`sidebar-submenu-item ${sectionActive ? "active" : ""}`}
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

function DashboardsSidebarSection({ active, onSelect }: SectionNavigationProps) {
  return <SidebarSection active={active} href="?view=dashboard" icon={<LayoutDashboard aria-hidden="true" />} label="Dashboards" onSelect={onSelect} />;
}

interface SectionNavigationProps {
  active: boolean;
  onSelect: () => void;
}

function SidebarFooter({ darkMode, modelRunLabel, onToggleDarkMode }: {
  darkMode: boolean;
  modelRunLabel?: string;
  onToggleDarkMode: () => void;
}) {
  return <div className="sidebar-foot">
    <div className="sidebar-statuses">
      <span className="sidebar-status-indicator" data-label="Local files connected" aria-label="Local files connected" tabIndex={0}>
        <i className="status-dot" aria-hidden="true" />
      </span>
      {modelRunLabel && <span className="sidebar-status-indicator model-run" data-label={modelRunLabel} aria-label={modelRunLabel} role="status" aria-live="polite" tabIndex={0}>
        <i className="status-dot" aria-hidden="true" />
      </span>}
    </div>
    <a href="/docs" target="_blank">API</a>
    <button onClick={onToggleDarkMode} aria-label="Toggle dark mode">
      {darkMode ? <Sun aria-hidden="true" /> : <Moon aria-hidden="true" />}
    </button>
  </div>;
}
