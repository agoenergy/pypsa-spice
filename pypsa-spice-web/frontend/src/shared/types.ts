import type { Data, PlotData } from "plotly.js";

export type AxisRange = [number, number];
export type TimeAxisRange = [string, string];

export interface SharedYAxisRanges {
  primary?: AxisRange;
  secondary?: AxisRange;
}

// @types/plotly.js 3.0.13 omits per-point scatter text positions supported by Plotly.
// https://plotly.com/javascript/reference/scatter/#scatter-textposition
export type ScatterTextTrace = Omit<Partial<PlotData>, "type" | "textposition"> & {
  type: "scatter";
  textposition?: PlotData["textposition"] | PlotData["textposition"][];
};

export type PlotTrace = Data | ScatterTextTrace;
