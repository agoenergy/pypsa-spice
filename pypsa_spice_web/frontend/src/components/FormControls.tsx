import type { ButtonHTMLAttributes, ChangeEvent, ReactNode } from "react";
import { Search } from "lucide-react";

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

export function SelectField({
  label,
  value,
  options,
  onChange,
  className = "",
  compact = false,
  variant = "field",
}: SelectFieldProps) {
  const baseClass = variant === "context" ? "context-control" : `field${compact ? " compact" : ""}`;
  return <label className={`${baseClass} ${className}`.trim()}>
    <span>{label}</span>
    <select value={value} onChange={(event) => onChange(event.target.value)}>
      {options.map((option) => <option value={option.value} key={option.value || "empty"}>{option.label}</option>)}
    </select>
  </label>;
}

interface SearchFieldProps {
  value: string;
  onChange: (value: string) => void;
  placeholder: string;
  className?: string;
}

export function SearchField({ value, onChange, placeholder, className = "" }: SearchFieldProps) {
  return <label className={`search ${className}`.trim()}>
    <Search aria-hidden="true" />
    <input
      value={value}
      onChange={(event: ChangeEvent<HTMLInputElement>) => onChange(event.target.value)}
      type="search"
      placeholder={placeholder}
    />
  </label>;
}

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: "primary" | "secondary";
  children: ReactNode;
}

export function Button({ variant = "secondary", className = "", children, type = "button", ...props }: ButtonProps) {
  return <button type={type} className={`button ${variant} ${className}`.trim()} {...props}>{children}</button>;
}

interface ToggleFieldProps {
  label: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
  className?: string;
}

export function ToggleField({ label, checked, onChange, className = "" }: ToggleFieldProps) {
  return <label className={`toggle-row ${className}`.trim()}>
    <input type="checkbox" checked={checked} onChange={(event) => onChange(event.target.checked)} />
    <span>{label}</span>
  </label>;
}
