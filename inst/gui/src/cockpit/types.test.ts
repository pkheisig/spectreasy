import assert from "node:assert/strict";
import test from "node:test";

import {
  defaultWorkflowSettings,
  normalizeWorkflowUnmixingMethod,
  WORKFLOW_UNMIXING_METHODS,
} from "./types.ts";

test("workflow method registry keeps AutoSpectral as the default", () => {
  const settings = defaultWorkflowSettings("");

  assert.equal(settings.control.method, "AutoSpectral");
  assert.equal(settings.sample.method, "AutoSpectral");
  assert.equal(WORKFLOW_UNMIXING_METHODS[0], "AutoSpectral");
  assert.deepEqual(
    [...WORKFLOW_UNMIXING_METHODS],
    ["AutoSpectral", "OLS", "WLS", "RWLS", "NNLS"],
  );
});

test("unsupported persisted workflow methods migrate to AutoSpectral", () => {
  assert.equal(normalizeWorkflowUnmixingMethod("unsupported"), "AutoSpectral");
  assert.equal(normalizeWorkflowUnmixingMethod("WLS"), "WLS");
});
