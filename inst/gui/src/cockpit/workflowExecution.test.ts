/// <reference types="node" />
import assert from "node:assert/strict";
import test from "node:test";
import { createWorkflowPayload, workflowHasExistingResults } from "./workflowExecution.ts";
import { defaultWorkflowSettings } from "./types.ts";

test("control workflow requests a versioned complete stage", () => {
  const settings = defaultWorkflowSettings("/project");
  const payload = createWorkflowPayload("control", settings);
  assert.equal(payload.output_collision, "version");
});

test("control rerun detection does not depend on a QC report", () => {
  assert.equal(workflowHasExistingResults("control", {
    artifacts: [],
    matrixFiles: ["spectreasy_outputs/unmix_controls/scc_reference_matrix.csv"],
  }), true);
  assert.equal(workflowHasExistingResults("control", {
    artifacts: [],
    matrixFiles: ["spectreasy_outputs/unmix_controls_10/scc_unmixing_matrix.csv"],
  }), true);
  assert.equal(workflowHasExistingResults("control", {
    artifacts: [],
    matrixFiles: [],
  }), false);
});
