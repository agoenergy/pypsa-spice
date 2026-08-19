import { describe, expect, it } from "vitest";
import { filterComparisonSections } from "./ScenarioComparison";
import type { ScenarioDifferenceSection } from "../types";

const sections: ScenarioDifferenceSection[] = [
  {
    id: "config:co2",
    label: "CO₂ management",
    category: "configuration",
    kind: "config",
    changes: [
      {
        status: "changed",
        item: "DE",
        country: "DE",
        parameter: "Option",
        reference: "co2_cap",
        comparison: "co2_price",
        delta: null,
      },
    ],
  },
  {
    id: "input:power:generators",
    label: "Power generators",
    category: "power",
    kind: "input",
    changes: [
      {
        status: "added",
        item: "DE_solar",
        country: "DE",
        parameter: "P Nom",
        reference: null,
        comparison: 20,
        delta: null,
      },
      {
        status: "removed",
        item: "FR_coal",
        country: "FR",
        parameter: "P Nom",
        reference: 10,
        comparison: null,
        delta: null,
      },
    ],
  },
];

describe("filterComparisonSections", () => {
  it("removes groups without matching changed values", () => {
    const filtered = filterComparisonSections(sections, {
      query: "DE_solar",
    });

    expect(filtered).toHaveLength(1);
    expect(filtered[0].id).toBe("input:power:generators");
    expect(filtered[0].changes.map((change) => change.item)).toEqual(["DE_solar"]);
  });

  it("searches group labels, items, parameters, and values", () => {
    const filtered = filterComparisonSections(sections, {
      query: "co2_price",
    });

    expect(filtered.map((section) => section.id)).toEqual(["config:co2"]);
  });
});
