import { useEffect, useMemo, useState } from "react";
import { createPortal } from "react-dom";
import { Cpu, Table2 } from "lucide-react";
import "./InputEditor.css";
import { SelectField } from "../components/FormControls";
import sidebarStyles from "../components/Sidebar.module.scss";
import type { InputCatalog, InputSelection } from "../types";
import { confirmDiscardChanges } from "../utility";
import { TableTitle, TableView } from "./TableEditor";
import TechnologyEditor, { TechnologyTitle } from "./TechnologyEditor";

type InputEditorView = "table" | "technology";

export default function InputEditor({ catalog, selection, onNavigate }: {
  catalog: InputCatalog;
  selection: InputSelection;
  onNavigate: () => void;
}) {
  const [view, setView] = useState<InputEditorView>("technology");
  const [sector, setSector] = useState("power");
  const project = catalog.datasets
    .find((item) => item.name === selection.dataset)
    ?.projects.find((item) => item.name === selection.project);
  const technologies = (project?.technologies || []).filter((item) => item.sector === sector);
  const [technologyId, setTechnologyId] = useState("");
  const [menuTarget, setMenuTarget] = useState<HTMLElement | null>(null);
  const [topbarTarget, setTopbarTarget] = useState<HTMLElement | null>(null);
  const definitions = useMemo(
    () => [...catalog.global_tables, ...(catalog.sector_tables[sector] || [])],
    [catalog, sector],
  );
  const defaultTable = catalog.global_tables.find((item) => item.id === "Technologies")?.id
    || definitions[0]?.id
    || "";
  const [tableId, setTableId] = useState(defaultTable);

  useEffect(() => {
    if (definitions.some((item) => item.id === tableId)) return;
    const preferred = ({
      power: "Power_generators",
      industry: "Heat_generators",
      transport: "Transport_loads",
    } as Record<string, string>)[sector];
    setTableId(definitions.find((item) => item.id === preferred)?.id || defaultTable);
  }, [defaultTable, definitions, sector, tableId]);

  useEffect(() => {
    if (!technologies.some((item) => item.id === technologyId)) setTechnologyId(technologies[0]?.id || "");
  }, [technologies, technologyId]);

  useEffect(() => {
    setMenuTarget(document.getElementById("input-table-menu"));
    setTopbarTarget(document.getElementById("input-topbar-controls"));
  }, []);

  const definition = definitions.find((item) => item.id === tableId);
  const technology = technologies.find((item) => item.id === technologyId) || technologies[0];
  const guarded = (action: () => void) => { if (confirmDiscardChanges()) action(); };
  const chooseView = (nextView: InputEditorView) => guarded(() => {
    setView(nextView);
    onNavigate();
  });

  return <>
    {menuTarget && createPortal(<nav className={sidebarStyles["submenu-list"]} aria-label="Input pages">
      <button className={`${sidebarStyles["submenu-item"]} ${view === "technology" ? sidebarStyles["active"] : ""}`} onClick={() => chooseView("technology")} aria-current={view === "technology" ? "page" : undefined}><Cpu aria-hidden="true" /><b>By technology</b></button>
      <button className={`${sidebarStyles["submenu-item"]} ${view === "table" ? sidebarStyles["active"] : ""}`} onClick={() => chooseView("table")} aria-current={view === "table" ? "page" : undefined}><Table2 aria-hidden="true" /><b>By table</b></button>
    </nav>, menuTarget)}
    {topbarTarget && createPortal(<>
      <SelectField
        variant="context"
        className="input-sector-control"
        label="Sector"
        value={sector}
        onChange={(value) => guarded(() => setSector(value))}
        options={["power", "industry", "transport"].map((value) => ({ value, label: value }))}
      />
      {view === "technology" ? <SelectField
        variant="context"
        className="input-technology-control"
        label="Technology"
        value={technologyId}
        onChange={(value) => guarded(() => setTechnologyId(value))}
        options={technologies.map((item) => ({ value: item.id, label: `${item.label} (${item.id})` }))}
      /> : <SelectField
        variant="context"
        className="input-table-control"
        label="Table"
        value={tableId}
        onChange={(value) => guarded(() => setTableId(value))}
        options={definitions.map((item) => ({ value: item.id, label: item.label }))}
      />}
    </>, topbarTarget)}
    {view === "technology"
      ? technology && <TechnologyTitle technology={technology} />
      : definition && <TableTitle definition={definition} />}
    {view === "table"
      ? definition
        ? <TableView
          key={`${selection.dataset}:${selection.project}:${selection.scenario}:${sector}:${definition.id}`}
          definition={definition}
          selection={selection}
        />
        : <div className="editor-empty">No configured table is available for this selection.</div>
      : technology
        ? <TechnologyEditor catalog={catalog} selection={selection} sector={sector} technology={technology} />
        : <div className="editor-empty">No mapped technology is available for this sector.</div>}
  </>;
}
