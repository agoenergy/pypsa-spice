import type { ButtonHTMLAttributes } from "react";
import styles from "./Button.module.scss";

interface ButtonProps extends ButtonHTMLAttributes<HTMLButtonElement> {
  variant?: "primary" | "secondary" | "danger";
}

export default function Button({ variant = "secondary", className = "", type = "button", ...props }: ButtonProps) {
  return (
    <button
      {...props}
      type={type}
      data-control="button"
      className={`${styles["button"]} ${styles[variant]} ${className}`.trim()}
    />
  );
}
