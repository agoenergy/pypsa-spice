import type { Data, PlotData } from "plotly.js";

// @types/plotly.js 3.0.13 omits per-point scatter text positions supported by Plotly.
// https://plotly.com/javascript/reference/scatter/#scatter-textposition
export type ScatterTextTrace = Omit<Partial<PlotData>, "type" | "textposition"> & {
  type: "scatter";
  textposition?: PlotData["textposition"] | PlotData["textposition"][];
};

export type PlotTrace = Data | ScatterTextTrace;
