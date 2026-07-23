/// <reference types="node" />
import test from "node:test";
import assert from "node:assert/strict";
import {
  createRefreshSequence,
  mergeProjectRefresh,
  normalizeProjectPickerPayload,
} from "./projectRefresh.ts";
import type { ProjectState } from "./types.ts";

function project(overrides: Partial<ProjectState> = {}): ProjectState {
  return {
    projectName: "project",
    projectPath: "/project",
    cytometer: "Aurora",
    method: "AutoSpectral",
    lastAction: "",
    lastActionAt: "",
    artifacts: [],
    mapping: [],
    mappingDirty: false,
    gatesDirty: false,
    missingInputDirs: [],
    controlInputDir: "scc",
    sampleInputDir: "samples",
    dataRevision: "empty",
    matrixFiles: null,
    sampleFiles: null,
    gatingFiles: null,
    gatingMetadata: {},
    scan: { controls: 0, samples: 0, matrices: 0, reports: 0, qcMetrics: 0, spectralVariants: 0, gates: 0 },
    ...overrides,
  };
}

test("backend refresh preserves an unsaved mapping draft and stale artifacts", () => {
  const draft = project({
    mappingDirty: true,
    mapping: [{ id: "one", file: "one.fcs", fluorophore: "FITC", marker: "CD3", channel: "B1-A", controlType: "cell", universalNegative: "" }],
    artifacts: [{ id: "reference", name: "matrix", type: "matrix", group: "controls", status: "stale", detail: "", path: "", updated: "" }],
  });
  const refreshed = project({
    projectName: "updated",
    mapping: [],
    artifacts: [{ id: "reference", name: "matrix", type: "matrix", group: "controls", status: "current", detail: "", path: "", updated: "" }],
  });

  const merged = mergeProjectRefresh(draft, refreshed);
  assert.equal(merged.projectName, "updated");
  assert.deepEqual(merged.mapping, draft.mapping);
  assert.equal(merged.mappingDirty, true);
  assert.equal(merged.artifacts[0]?.status, "stale");
});

test("backend refresh replaces mapping after the draft has been saved", () => {
  const refreshed = project({ projectName: "updated" });
  assert.equal(mergeProjectRefresh(project(), refreshed), refreshed);
});

test("a project switch never carries a draft into a different project", () => {
  const draft = project({ mappingDirty: true });
  const other = project({ projectPath: "/other", projectName: "other" });
  assert.equal(mergeProjectRefresh(draft, other), other);
});

test("refresh sequence accepts only the most recently started request", () => {
  const sequence = createRefreshSequence();
  const older = sequence.begin();
  const newer = sequence.begin();
  assert.equal(sequence.isCurrent(older), false);
  assert.equal(sequence.isCurrent(newer), true);
});

test("project picker payloads accept direct and legacy R scalar wrappers", () => {
  assert.deepEqual(
    normalizeProjectPickerPayload({ success: true, cancelled: false, project_path: "/projects/new" }),
    { success: true, cancelled: false, message: "", projectPath: "/projects/new" },
  );
  assert.deepEqual(
    normalizeProjectPickerPayload({ success: [true], cancelled: [false], project: { project_path: ["/projects/existing"] } }),
    { success: true, cancelled: false, message: "", projectPath: "/projects/existing" },
  );
});

test("project picker payloads reject cancellation and missing paths", () => {
  assert.equal(normalizeProjectPickerPayload({ success: [false], cancelled: [true] }).cancelled, true);
  assert.deepEqual(
    normalizeProjectPickerPayload({ success: [true], cancelled: [false] }),
    { success: false, cancelled: false, message: "The folder picker returned no project path." },
  );
});
