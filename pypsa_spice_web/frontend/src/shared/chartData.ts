import type { ChartDefinition, ChartResponse, ResultRow } from "../types";

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
