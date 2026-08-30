import PageHeader from "../components/PageHeader";
import type { InputCatalog, InputSelection, InputTechnology } from "../types";
import TableEditor from "./TableEditor";

export function TechnologyTitle({ technology }: { technology: InputTechnology }) {
  return <PageHeader title={technology.label} className="selection-title">
    <dl>
      <div><dt>PyPSA class</dt><dd>{technology.classes.join(", ") || "—"}</dd></div>
      <div><dt>Carrier</dt><dd>{technology.carriers.join(", ") || "—"}</dd></div>
    </dl>
  </PageHeader>;
}

export default function TechnologyEditor({ catalog, selection, sector, technology }: {
  catalog: InputCatalog;
  selection: InputSelection;
  sector: string;
  technology: InputTechnology;
}) {
  const globalDefinitions = catalog.global_tables.filter((item) => item.id !== "Demand_Profiles");
  const scenarioDefinitions = catalog.sector_tables[sector] || [];
  return <div className="technology-view">
    <section className="technology-group">
      <header><h2>Global input</h2><span>Changes here apply to every country and every scenario in this project.</span></header>
      <div className="technology-panels">{globalDefinitions.map((definition) => <TableEditor
        key={`${selection.dataset}:${selection.project}:global:${definition.id}:${technology.id}`}
        definition={definition}
        selection={selection}
        technology={technology}
        hideWhenEmpty
      />)}</div>
    </section>
    <section className="technology-group">
      <header><h2>Scenario input</h2><span>Assets and constraints for this scenario. Country filters appear only on tables with country-specific rows.</span></header>
      <div className="technology-panels">{scenarioDefinitions.map((definition) => <TableEditor
        key={`${selection.dataset}:${selection.project}:${selection.scenario}:${definition.id}:${technology.id}`}
        definition={definition}
        selection={selection}
        technology={technology}
        hideWhenEmpty
      />)}</div>
    </section>
  </div>;
}
