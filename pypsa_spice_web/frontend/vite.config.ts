import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import { relative, resolve } from "node:path";
import { createHash } from "node:crypto";

export default defineConfig({
  root: resolve(__dirname),
  base: "/ui/",
  plugins: [react()],
  css: {
    modules: {
      // Shared @value imports and direct SCSS imports must resolve to the same
      // class, independently of Sass formatting or the importing stylesheet.
      generateScopedName: (name, filename) => {
        const path = relative(__dirname, filename).replaceAll("\\", "/");
        const hash = createHash("sha256").update(path).digest("hex").slice(0, 8);
        return `_${name}_${hash}`;
      },
    },
  },
  build: {
    outDir: resolve(__dirname, "dist"),
    emptyOutDir: true,
  },
  server: {
    port: 5173,
    proxy: {
      "/api": "http://127.0.0.1:8000",
      "/vendor": "http://127.0.0.1:8000",
    },
  },
});
