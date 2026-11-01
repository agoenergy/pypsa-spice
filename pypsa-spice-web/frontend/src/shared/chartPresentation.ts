import type { Catalog } from "../types";

const fallbackColors = [
  "#e6007e",
  "#005ca9",
  "#60a917",
  "#ec6608",
  "#7553a6",
  "#009e8e",
  "#c33c54",
  "#79848d",
  "#d5a400",
  "#3f7c85",
  "#9b4b96",
  "#86a6c2",
];

export function formatLegendLabel(value: string, mappings: Catalog["mappings"]): string {
  return mappings[value]?.label || value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

export function getLegendColour(value: string, index: number, mappings: Catalog["mappings"]): string {
  return mappings[value]?.color || fallbackColors[Math.max(0, index) % fallbackColors.length];
}

// Plotly draws its own text, so it cannot read the CSS type scale in global.scss.
// These mirror --text-xs and --text-sm so chart type matches the surrounding interface.
export const chartFont = { body: 12, hover: 13 };
