import type { ChangeEvent, LabelHTMLAttributes } from "react";
import { Search } from "lucide-react";
import styles from "./FormControls.module.scss";

export interface SelectOption {
  value: string;
  label: string;
}

interface SelectFieldProps {
  label: string;
  value: string;
  options: SelectOption[];
  onChange: (value: string) => void;
  className?: string;
  compact?: boolean;
  variant?: "context" | "field";
}

interface FieldProps extends LabelHTMLAttributes<HTMLLabelElement> {
  compact?: boolean;
  variant?: "context" | "field";
}

export function Field({ compact = false, variant = "field", className = "", ...props }: FieldProps) {
  return (
    <label
      {...props}
      data-control={variant === "context" ? "context" : "field"}
      data-compact={compact || undefined}
      className={[
        styles[variant === "context" ? "context-control" : "field"],
        compact ? styles["compact"] : "",
        className,
      ]
        .filter(Boolean)
        .join(" ")}
    />
  );
}

export function SelectField({
  label,
  value,
  options,
  onChange,
  className = "",
  compact = false,
  variant = "field",
}: SelectFieldProps) {
  return (
    <Field variant={variant} compact={compact} className={className}>
      <span>{label}</span>
      <select value={value} onChange={(event) => onChange(event.target.value)}>
        {options.map((option) => (
          <option value={option.value} key={option.value || "empty"}>
            {option.label}
          </option>
        ))}
      </select>
    </Field>
  );
}

interface SearchFieldProps {
  value: string;
  onChange: (value: string) => void;
  placeholder: string;
  label?: string;
  className?: string;
}

export function SearchField({ value, onChange, placeholder, label = placeholder, className = "" }: SearchFieldProps) {
  return (
    <label data-control="search" className={`${styles["search"]} ${className}`.trim()}>
      <Search aria-hidden="true" />
      <input
        value={value}
        onChange={(event: ChangeEvent<HTMLInputElement>) => onChange(event.target.value)}
        aria-label={label}
        type="search"
        placeholder={placeholder}
      />
    </label>
  );
}

interface ToggleFieldProps {
  label: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
  className?: string;
}

export function ToggleField({ label, checked, onChange, className = "" }: ToggleFieldProps) {
  return (
    <label className={`${styles["toggle-row"]} ${className}`.trim()}>
      <input type="checkbox" checked={checked} onChange={(event) => onChange(event.target.checked)} />
      <span>{label}</span>
    </label>
  );
}
