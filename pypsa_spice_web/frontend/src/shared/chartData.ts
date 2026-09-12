import type { ChartDefinition, ChartResponse, ResultRow } from "../types";
import type { AxisRange, SharedYAxisRanges, TimeAxisRange } from "./types";

export function aggregate(rows: ResultRow[], chart: ChartDefinition): Map<string, { x: string | number; y: number }[]> {
  const xKey = chart.hourly ? "snapshot" : "year";
  const values = new Map<string, number>();
  for (const row of rows) {
    const x = row[xKey];
    const legend = String(row[chart.leg_col] ?? "Series");
    if (x === null || x === undefined) continue;
    const key = `${x}\u0000${legend}`;
    values.set(key, (values.get(key) || 0) + Number(row.value || 0));
  }
  const groups = new Map<string, { x: string | number; y: number }[]>();
  for (const [key, y] of values) {
    const [rawX, legend] = key.split("\u0000");
    if (!groups.has(legend)) groups.set(legend, []);
    groups.get(legend)!.push({ x: chart.hourly ? rawX : Number(rawX), y });
  }
  for (const points of groups.values()) points.sort((a, b) => String(a.x).localeCompare(String(b.x)));
  return groups;
}

export function getLegendValues(chart: ChartDefinition, ...responses: (ChartResponse | null)[]): string[] {
  const values = new Set<string>();
  for (const response of responses) {
    if (!response) continue;
    for (const value of aggregate(response.rows, chart).keys()) values.add(value);
  }
  return [...values];
}

function paddedRange(minimum: number, maximum: number, includeZero: boolean): AxisRange {
  const minimumWithZero = includeZero ? Math.min(0, minimum) : minimum;
  const maximumWithZero = includeZero ? Math.max(0, maximum) : maximum;
  const span = maximumWithZero - minimumWithZero;
  const padding = span > 0 ? span * 0.05 : Math.max(Math.abs(minimumWithZero) * 0.05, 1);

  return [
    includeZero && minimumWithZero === 0 ? 0 : minimumWithZero - padding,
    includeZero && maximumWithZero === 0 ? 0 : maximumWithZero + padding,
  ];
}

function responseAxisExtents(
  response: ChartResponse,
  chart: ChartDefinition,
  hiddenLegendValues: ReadonlySet<string>,
  secondary: boolean,
): [number, number] | null {
  const exactExtent = secondary ? response.meta.axis_extents?.secondary : response.meta.axis_extents?.primary;
  if (exactExtent) return exactExtent;

  const aggregates = aggregate(response.rows, chart);
  const shouldStack = (chart.type.includes("bar") && chart.type !== "grouped_bar") || chart.type === "area_share";
  let minimum = Number.POSITIVE_INFINITY;
  let maximum = Number.NEGATIVE_INFINITY;

  if (shouldStack) {
    const totals = new Map<string, { positive: number; negative: number }>();
    for (const [legend, points] of aggregates) {
      const isSecondary = Boolean(chart.secondary_y_lab?.includes(legend));
      if (hiddenLegendValues.has(legend) || isSecondary !== secondary) continue;
      for (const point of points) {
        if (!Number.isFinite(point.y)) continue;
        const key = String(point.x);
        const total = totals.get(key) || { positive: 0, negative: 0 };
        if (point.y >= 0) total.positive += point.y;
        else total.negative += point.y;
        totals.set(key, total);
      }
    }
    for (const total of totals.values()) {
      minimum = Math.min(minimum, total.negative);
      maximum = Math.max(maximum, total.positive);
    }
  } else {
    for (const [legend, points] of aggregates) {
      const isSecondary = Boolean(chart.secondary_y_lab?.includes(legend));
      if (hiddenLegendValues.has(legend) || isSecondary !== secondary) continue;
      for (const point of points) {
        if (!Number.isFinite(point.y)) continue;
        minimum = Math.min(minimum, point.y);
        maximum = Math.max(maximum, point.y);
      }
    }
  }

  return Number.isFinite(minimum) && Number.isFinite(maximum) ? [minimum, maximum] : null;
}

export function getSharedYAxisRanges(
  responses: ChartResponse[],
  chart: ChartDefinition,
  hiddenLegendValues: ReadonlySet<string>,
): SharedYAxisRanges {
  const rangeForAxis = (secondary: boolean): AxisRange | undefined => {
    const extents = responses
      .map((response) => responseAxisExtents(response, chart, hiddenLegendValues, secondary))
      .filter((extent): extent is [number, number] => extent !== null);
    if (extents.length === 0) return undefined;
    const minimum = Math.min(...extents.map(([value]) => value));
    const maximum = Math.max(...extents.map(([, value]) => value));
    return paddedRange(minimum, maximum, !secondary);
  };

  return {
    primary: rangeForAxis(false),
    secondary: chart.secondary_y_lab?.length ? rangeForAxis(true) : undefined,
  };
}

function timestampValue(value: string): number {
  return new Date(value.replace(" ", "T")).getTime();
}

function earlierTimestamp(left: string, right: string): string {
  const leftValue = timestampValue(left);
  const rightValue = timestampValue(right);
  if (Number.isFinite(leftValue) && Number.isFinite(rightValue)) return leftValue <= rightValue ? left : right;
  return left.localeCompare(right) <= 0 ? left : right;
}

function laterTimestamp(left: string, right: string): string {
  return earlierTimestamp(left, right) === left ? right : left;
}

function responseTimeRange(response: ChartResponse, selectedStart: string, selectedEnd: string): TimeAxisRange | null {
  let availableStart = response.meta.available_start;
  let availableEnd = response.meta.available_end;
  if (!availableStart || !availableEnd) {
    for (const row of response.rows) {
      if (row.snapshot === null || row.snapshot === undefined) continue;
      const snapshot = String(row.snapshot);
      availableStart = availableStart ? earlierTimestamp(availableStart, snapshot) : snapshot;
      availableEnd = availableEnd ? laterTimestamp(availableEnd, snapshot) : snapshot;
    }
  }
  if (!availableStart || !availableEnd || response.meta.returned_rows === 0) return null;

  const start = selectedStart ? laterTimestamp(availableStart, selectedStart) : availableStart;
  const end = selectedEnd ? earlierTimestamp(availableEnd, selectedEnd) : availableEnd;
  return earlierTimestamp(start, end) === start ? [start, end] : null;
}

export function getSharedXAxisRange(
  responses: ChartResponse[],
  chart: ChartDefinition,
  selectedStart = "",
  selectedEnd = "",
): TimeAxisRange | undefined {
  if (!chart.hourly) return undefined;
  const ranges = responses
    .map((response) => responseTimeRange(response, selectedStart, selectedEnd))
    .filter((range): range is TimeAxisRange => range !== null);
  if (ranges.length === 0) return undefined;
  const start = ranges.map(([value]) => value).reduce(earlierTimestamp);
  const end = ranges.map(([, value]) => value).reduce(laterTimestamp);
  if (timestampValue(start) !== timestampValue(end)) return [start, end];

  const centre = timestampValue(start);
  if (!Number.isFinite(centre)) return undefined;
  const halfHour = 30 * 60 * 1000;
  return [new Date(centre - halfHour).toISOString(), new Date(centre + halfHour).toISOString()];
}

export function differenceAggregates(
  primary: ChartResponse,
  comparison: ChartResponse,
  chart: ChartDefinition,
): Map<string, { x: string | number; y: number }[]> {
  const first = aggregate(primary.rows, chart);
  const second = aggregate(comparison.rows, chart);
  const legends = [...new Set([...first.keys(), ...second.keys()])];
  const differences = new Map<string, { x: string | number; y: number }[]>();
  for (const legend of legends) {
    const firstPoints = new Map((first.get(legend) || []).map((point) => [String(point.x), point.y]));
    const secondPoints = new Map((second.get(legend) || []).map((point) => [String(point.x), point.y]));
    const xValues = [...new Set([...firstPoints.keys(), ...secondPoints.keys()])].sort();
    differences.set(
      legend,
      xValues.map((value) => ({
        x: chart.hourly ? value : Number(value),
        y: (secondPoints.get(value) || 0) - (firstPoints.get(value) || 0),
      })),
    );
  }
  return differences;
}

export function buildDifferenceRows(
  primary: ChartResponse,
  comparison: ChartResponse,
  chart: ChartDefinition,
  primaryName: string,
  comparisonName: string,
): ResultRow[] {
  const xKey = chart.hourly ? "snapshot" : "year";
  const unitValues = new Set(
    [...primary.rows, ...comparison.rows]
      .map((row) => row.unit)
      .filter((unit): unit is string | number => unit !== null && unit !== undefined),
  );
  const unit = unitValues.size === 1 ? [...unitValues][0] : undefined;
  const rows: ResultRow[] = [];
  for (const [legend, points] of differenceAggregates(primary, comparison, chart)) {
    if (points.every((point) => point.y === 0)) continue;
    for (const point of points) {
      rows.push({
        scenario: `${comparisonName} − ${primaryName}`,
        [chart.leg_col]: legend,
        ...(unit !== undefined ? { unit } : {}),
        [xKey]: point.x,
        value: point.y,
      });
    }
  }
  return rows;
}
