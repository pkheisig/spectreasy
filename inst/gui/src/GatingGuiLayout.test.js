import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";

const styles = readFileSync(new URL("./GatingGui.css", import.meta.url), "utf8");

test("control sidebar cards retain their two-row four-corner layout", () => {
  assert.match(styles, /\.app-shell \.file-row\s*\{[^}]*display:\s*grid;/s);
  assert.match(styles, /grid-template-areas:\s*"fluorophore channel"\s*"marker control-type";/s);
  assert.match(styles, /\.file-row-fluorophore\s*\{[^}]*grid-area:\s*fluorophore;/s);
  assert.match(styles, /\.file-row-channel\s*\{[^}]*grid-area:\s*channel;/s);
  assert.match(styles, /\.file-row-marker\s*\{[^}]*grid-area:\s*marker;/s);
  assert.match(styles, /\.file-row-control-type\s*\{[^}]*grid-area:\s*control-type;/s);
});
