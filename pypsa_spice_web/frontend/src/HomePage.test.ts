import { describe, expect, it } from "vitest";
import { workspaceInventory } from "./HomePage";
import type { Catalog, InputCatalog } from "./types";

const catalog: Catalog = {
  datasets: [{
    name: "results",
    projects: [
      { name: "shared", scenarios: [{ name: "baseline", sectors: [] }, { name: "policy", sectors: [] }] },
      { name: "results-only", scenarios: [{ name: "latest", sectors: [] }] },
    ],
  }],
  sections: [
    { id: "power", label: "Power", title: "Power", charts: [{ id: "a", name: "A", summary: "sum", table_name: "a", leg_col: "carrier", type: "bar", hourly: false }] },
    { id: "costs", label: "Costs", title: "Costs", charts: [{ id: "b", name: "B", summary: "sum", table_name: "b", leg_col: "carrier", type: "bar", hourly: false }] },
  ],
  mappings: {},
};

const inputCatalog: InputCatalog = {
  table_query_version: 1,
  datasets: [{
    name: "inputs",
    projects: [
      { name: "shared", scenarios: ["baseline", "policy", "high-demand"], countries: [], technologies: [] },
      { name: "inputs-only", scenarios: ["baseline"], countries: [], technologies: [] },
    ],
  }],
  global_tables: [],
  sector_tables: {},
};

describe("home workspace inventory", () => {
  it("keeps workspaces distinct by dataset and project", () => {
    expect(workspaceInventory(catalog, inputCatalog)).toEqual([
      {
        key: "inputs::inputs-only",
        dataset: "inputs",
        project: "inputs-only",
        inputScenarios: ["baseline"],
        resultRuns: [],
        countries: [],
      },
      {
        key: "results::results-only",
        dataset: "results",
        project: "results-only",
        inputScenarios: [],
        resultRuns: ["latest"],
        countries: [],
      },
      {
        key: "inputs::shared",
        dataset: "inputs",
        project: "shared",
        inputScenarios: ["baseline", "policy", "high-demand"],
        resultRuns: [],
        countries: [],
      },
      {
        key: "results::shared",
        dataset: "results",
        project: "shared",
        inputScenarios: [],
        resultRuns: ["baseline", "policy"],
        countries: [],
      },
    ]);
  });

  it("merges inputs and results when they belong to the same workspace", () => {
    const matchingCatalog: Catalog = {
      ...catalog,
      datasets: [{
        name: "inputs",
        projects: [{ name: "shared", scenarios: [{ name: "baseline", sectors: [] }] }],
      }],
    };

    expect(workspaceInventory(matchingCatalog, inputCatalog).find((workspace) => workspace.key === "inputs::shared")).toEqual({
      key: "inputs::shared",
      dataset: "inputs",
      project: "shared",
      inputScenarios: ["baseline", "policy", "high-demand"],
      resultRuns: ["baseline"],
      countries: [],
    });
  });
});
