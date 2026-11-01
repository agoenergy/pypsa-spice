import { readFileSync, readdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import ts from "typescript";
import { preprocessCSS, resolveConfig } from "vite";
import { describe, expect, it } from "vitest";

const sourceRoot = resolve(__dirname);
const files = readdirSync(sourceRoot, { recursive: true }) as string[];

describe("stylesheet boundaries", () => {
  it("keeps component styles in SCSS modules and global styles free of classes", () => {
    const stylesheets = files.filter((file) => /\.(css|scss)$/.test(file));
    expect(stylesheets.filter((file) => file !== "global.scss" && !file.endsWith(".module.scss"))).toEqual([]);
    const global = readFileSync(resolve(sourceRoot, "global.scss"), "utf8");
    expect(global).not.toMatch(/\.[a-zA-Z_-][\w-]*\s*[{,:]/);
    for (const file of stylesheets) {
      expect(readFileSync(resolve(sourceRoot, file), "utf8"), file).not.toContain(":global");
    }
  });

  it("resolves JSX classes and shared selector imports to the same scoped exports", async () => {
    const config = await resolveConfig({}, "build");
    const modules = new Map<string, Record<string, string>>();
    for (const file of files.filter((file) => file.endsWith(".module.scss"))) {
      const filename = resolve(sourceRoot, file);
      const compiled = await preprocessCSS(readFileSync(filename, "utf8"), filename, config);
      modules.set(filename, compiled.modules ?? {});
    }
    for (const [filename, exports] of modules) {
      const source = readFileSync(filename, "utf8");
      for (const match of source.matchAll(/@value\s+([^;]+?)\s+from\s+"([^"]+)";/g)) {
        const shared = modules.get(resolve(dirname(filename), match[2]));
        expect(shared, `${filename}: ${match[2]}`).toBeDefined();
        for (const name of match[1].split(",").map((name) => name.trim())) {
          expect(shared?.[name], `${filename}: shared ${name}`).toBeTruthy();
          expect(exports[name], `${filename}: shared ${name}`).toBe(shared?.[name]);
        }
      }
    }
    for (const file of files.filter((file) => file.endsWith(".tsx"))) {
      const filename = resolve(sourceRoot, file);
      const ast = ts.createSourceFile(
        filename,
        readFileSync(filename, "utf8"),
        ts.ScriptTarget.Latest,
        true,
        ts.ScriptKind.TSX,
      );
      const bindings = new Map<string, Record<string, string>>();
      for (const node of ast.statements) {
        if (ts.isImportDeclaration(node) && node.importClause?.name && ts.isStringLiteral(node.moduleSpecifier)) {
          const exports = modules.get(resolve(dirname(filename), node.moduleSpecifier.text));
          if (exports) bindings.set(node.importClause.name.text, exports);
        }
      }
      const visit = (node: ts.Node) => {
        if (
          ts.isElementAccessExpression(node) &&
          ts.isIdentifier(node.expression) &&
          ts.isStringLiteral(node.argumentExpression)
        ) {
          const exports = bindings.get(node.expression.text);
          if (exports) expect(exports[node.argumentExpression.text], `${file}: ${node.getText(ast)}`).toBeTruthy();
        }
        ts.forEachChild(node, visit);
      };
      visit(ast);
    }
  }, 15000);
});
