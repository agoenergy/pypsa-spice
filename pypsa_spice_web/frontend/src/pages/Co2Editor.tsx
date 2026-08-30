import { SelectField } from "../components/FormControls";
import { MappingTable } from "./ScenarioConfigControls";

export default function Co2Editor({ value, country, onChange }: {
  value: Record<string, unknown>;
  country: string;
  onChange: (value: Record<string, unknown>) => void;
}) {
  const updateCountry = (countryName: string, next: Record<string, unknown>) => {
    onChange({ ...value, [countryName]: next });
  };
  const entries = Object.entries(value).filter(([name]) => country === "ALL" || name === country);
  return <div className="config-form">
    {entries.map(([name, raw]) => {
      const config = (raw || {}) as Record<string, unknown>;
      return <article className="country-config" key={name}>
        <header>
          <h3>{name}</h3>
          <SelectField
            label="Instrument"
            value={String(config.option || "co2_cap")}
            onChange={(option) => updateCountry(name, { ...config, option })}
            compact
            options={[
              { value: "co2_cap", label: "CO₂ cap" },
              { value: "co2_price", label: "CO₂ price" },
            ]}
          />
        </header>
        <MappingTable
          label="Year"
          value={(config.value || {}) as Record<string, unknown>}
          yearKeys
          allowAdd={false}
          onChange={(next) => updateCountry(name, { ...config, value: next })}
        />
      </article>;
    })}
    {entries.length === 0 && <div className="editor-empty compact">No CO₂ configuration is defined for {country}.</div>}
  </div>;
}
