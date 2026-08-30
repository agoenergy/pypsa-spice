import type { InputRow, ScenarioConfigResponse } from "../types";

export const COMBINED_SECTION = "co2_constraints";

export const SECTION_LABELS: Record<string, string> = {
  scenario_configs: "Scenario settings",
  [COMBINED_SECTION]: "CO₂ & constraints",
  review_run: "Review & run",
};

export const YEAR_FIELDS = new Set(["ei_fraction", "res_generation_share"]);

export const CONFIG_LABELS: Record<string, string> = {
  energy_independence: "Energy independence",
  production_constraint_fuels: "Fuel production limits",
  reserve_margin: "Reserve margin",
  res_generation: "Renewable generation share",
  thermal_must_run: "Thermal must-run",
  capacity_factor_constraint: "Capacity factor",
  maximum_power_generation_constraint: "Maximum power generation",
  pe_conv_fraction: "Primary energy conversion fractions",
  ei_fraction: "Energy independence fraction",
  fuels: "Restricted fuels",
  epsilon_load: "Load reserve fraction",
  epsilon_vre: "VRE reserve contribution",
  contingency: "Contingency (MW)",
  res_generation_share: "Renewable generation share",
  min_must_run_ratio: "Minimum must-run ratio",
};

export const CONSTRAINT_DEFAULTS: Record<string, Record<string, unknown>> = {
  energy_independence: {
    activate: false,
    pe_conv_fraction: { Solar: null, Wind: null, Geothermal: null, Water: null },
    ei_fraction: {},
  },
  production_constraint_fuels: { activate: false, fuels: [] },
  reserve_margin: { activate: false, epsilon_load: null, epsilon_vre: null, contingency: null, method: "static" },
  res_generation: { activate: false, math_symbol: "<=", res_generation_share: {} },
  thermal_must_run: { activate: false, min_must_run_ratio: null },
  capacity_factor_constraint: { activate: false, value: {} },
  maximum_power_generation_constraint: { activate: false, value: {} },
};

export function normaliseSection(requested: string): string {
  if (requested === "co2_management" || requested === "custom_constraints") return COMBINED_SECTION;
  return Object.hasOwn(SECTION_LABELS, requested) ? requested : "scenario_configs";
}

export function draftForSection(config: ScenarioConfigResponse, section: string): Record<string, unknown> {
  if (section === COMBINED_SECTION) {
    return {
      co2_management: structuredClone(config.sections.co2_management || {}),
      custom_constraints: structuredClone(config.sections.custom_constraints || {}),
    };
  }
  return structuredClone(config.sections[section] || {});
}

export function fuelConstraintSnapshot(value: Record<string, unknown>): Record<string, unknown> {
  return Object.fromEntries(Object.entries(objectValue(value.custom_constraints)).map(([country, raw]) => [
    country,
    objectValue(objectValue(raw).production_constraint_fuels),
  ]));
}

export function objectValue(value: unknown): Record<string, unknown> {
  return value && typeof value === "object" && !Array.isArray(value) ? value as Record<string, unknown> : {};
}

export function prettyConfigLabel(value: string): string {
  return CONFIG_LABELS[value] || value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

export function constraintAnchor(name: string): string {
  return `config-constraint-${name.replace(/[^a-zA-Z0-9_-]/g, "-")}`;
}

export function isYearKey(value: string): boolean {
  return /^\d{4}$/.test(value);
}

export function isScalar(value: unknown): boolean {
  return value === null || ["string", "number", "boolean"].includes(typeof value);
}

export function isYearMatrix(value: Record<string, unknown>): boolean {
  const entries = Object.values(value);
  return entries.length > 0 && entries.every((raw) => {
    if (!raw || typeof raw !== "object" || Array.isArray(raw)) return false;
    return Object.keys(objectValue(raw)).every(isYearKey);
  });
}

export function inputNumber(value: string): number | null {
  return value === "" ? null : Number(value);
}

export function validateFuelLimits(section: string, value: Record<string, unknown>, rows: InputRow[]): string {
  if (section !== COMBINED_SECTION) return "";
  const constraints = objectValue(value.custom_constraints);
  for (const [country, rawConstraints] of Object.entries(constraints)) {
    const production = objectValue(objectValue(rawConstraints).production_constraint_fuels);
    if (production.activate !== true) continue;
    const fuels = Array.isArray(production.fuels) ? production.fuels.map(String) : [];
    if (!fuels.length) return `Select at least one fuel production limit for ${country}.`;
    for (const fuel of fuels) {
      const fuelRows = rows.filter((row) => String(row.country) === country && String(row.carrier) === fuel);
      if (!fuelRows.length) return `No fuel supply data is available for ${fuel} in ${country}.`;
      for (const row of fuelRows) {
        const raw = row.max_supply__mwh_year;
        if (typeof raw === "string" && raw.trim().toLowerCase() === "inf") continue;
        if (raw === "" || raw === null || !Number.isFinite(Number(raw)) || Number(raw) < 0) {
          return `${fuel} limits for ${country} must be zero or greater, or inf.`;
        }
      }
    }
  }
  return "";
}

export function validateSection(section: string, value: Record<string, unknown>): string {
  if (section === "scenario_configs") {
    const snapshots = objectValue(value.snapshots);
    if (!snapshots.start || !snapshots.end) return "Snapshot start and end dates are required.";
    if (String(snapshots.end) < String(snapshots.start)) return "Snapshot end must be on or after snapshot start.";
    const resolution = objectValue(value.resolution);
    if (!["nth_hour", "clustered"].includes(String(resolution.method))) {
      return "Choose a supported temporal resolution method.";
    }
    if (resolution.method === "clustered" && Number(resolution.number_of_days) < 1) {
      return "Number of days must be at least 1.";
    }
    if (resolution.method !== "clustered" && Number(resolution.stepsize) < 1) return "Step size must be at least 1.";
    if (!Number.isFinite(Number(value.remove_threshold)) || Number(value.remove_threshold) < 0) {
      return "Remove threshold must be zero or greater.";
    }
  }
  if (section === COMBINED_SECTION || section === "co2_management") {
    const co2 = section === COMBINED_SECTION ? objectValue(value.co2_management) : value;
    for (const raw of Object.values(co2)) {
      const config = objectValue(raw);
      if (!["co2_cap", "co2_price"].includes(String(config.option))) {
        return "Choose a supported CO₂ instrument for every country.";
      }
      if (Object.values(objectValue(config.value)).some((item) => !Number.isFinite(Number(item)))) {
        return "CO₂ year values must be numbers.";
      }
    }
  }
  return "";
}
