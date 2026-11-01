import type { AnchorHTMLAttributes, ButtonHTMLAttributes } from "react";
import styles from "./IconButton.module.scss";

interface IconControlProps {
  "aria-label": string;
  variant?: "plain" | "surface" | "toolbar";
  tone?: "default" | "danger";
  alignEnd?: boolean;
}

type IconButtonProps = Omit<ButtonHTMLAttributes<HTMLButtonElement>, "aria-label"> & IconControlProps;
type IconLinkProps = Omit<AnchorHTMLAttributes<HTMLAnchorElement>, "aria-label"> & IconControlProps & { href: string };

function controlClasses(variant: string, tone: string, alignEnd: boolean, className: string) {
  return [
    styles["icon-control"],
    styles[variant],
    tone === "danger" ? styles["danger"] : "",
    alignEnd ? styles["align-end"] : "",
    className,
  ]
    .filter(Boolean)
    .join(" ");
}

export default function IconButton({
  "aria-label": label,
  title = label,
  variant = "plain",
  tone = "default",
  alignEnd = false,
  className = "",
  type = "button",
  ...props
}: IconButtonProps) {
  return (
    <button
      {...props}
      type={type}
      aria-label={label}
      title={title}
      data-control="icon"
      className={controlClasses(variant, tone, alignEnd, className)}
    />
  );
}

export function IconLink({
  "aria-label": label,
  title = label,
  variant = "plain",
  tone = "default",
  alignEnd = false,
  className = "",
  ...props
}: IconLinkProps) {
  return (
    <a
      {...props}
      aria-label={label}
      title={title}
      data-control="icon"
      className={controlClasses(variant, tone, alignEnd, className)}
    />
  );
}
