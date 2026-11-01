import styles from "./ChartLegend.module.scss";
import type { Catalog } from "../types";
import { formatLegendLabel, getLegendColour } from "../shared/chartPresentation";

export function ChartLegend({
  values,
  mappings,
  hiddenValues,
  onToggle,
}: {
  values: string[];
  mappings: Catalog["mappings"];
  hiddenValues: ReadonlySet<string>;
  onToggle: (value: string) => void;
}) {
  return (
    <div className={styles["html-legend"]} aria-label="Chart legend">
      {values.map((value, index) => (
        <button
          type="button"
          className={[styles["html-legend-item"], hiddenValues.has(value) ? styles["is-hidden"] : ""]
            .filter(Boolean)
            .join(" ")}
          key={value}
          aria-pressed={!hiddenValues.has(value)}
          title={`${hiddenValues.has(value) ? "Show" : "Hide"} ${formatLegendLabel(value, mappings)}`}
          onClick={() => onToggle(value)}
        >
          <i style={{ backgroundColor: getLegendColour(value, index, mappings) }} aria-hidden="true" />
          <span>{formatLegendLabel(value, mappings)}</span>
        </button>
      ))}
    </div>
  );
}
