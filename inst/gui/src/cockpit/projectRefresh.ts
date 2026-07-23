import type { ProjectState } from "./types";

export type ProjectPickerResult = {
  success: boolean;
  cancelled: boolean;
  message: string;
  projectPath?: string;
};

function pickerScalar(value: unknown, fallback = ""): string {
  if (Array.isArray(value)) return pickerScalar(value[0], fallback);
  if (value == null) return fallback;
  return String(value);
}

/** Normalize scalar and single-element-array payloads emitted by R/Plumber. */
export function normalizeProjectPickerPayload(value: unknown): ProjectPickerResult {
  const payload = value && typeof value === "object"
    ? value as Record<string, unknown>
    : {};
  const project = payload.project && typeof payload.project === "object"
    ? payload.project as Record<string, unknown>
    : {};
  const cancelled = pickerScalar(payload.cancelled, "false") === "true";
  if (cancelled) {
    return { success: false, cancelled: true, message: "Project selection cancelled." };
  }
  const error = pickerScalar(payload.error, "");
  const success = pickerScalar(payload.success, "false") === "true";
  const projectPath = pickerScalar(payload.project_path ?? project.project_path, "").trim();
  if (!success || error) {
    return {
      success: false,
      cancelled: false,
      message: error || "The project folder could not be opened.",
    };
  }
  if (!projectPath) {
    return {
      success: false,
      cancelled: false,
      message: "The folder picker returned no project path.",
    };
  }
  return { success: true, cancelled: false, message: "", projectPath };
}

/**
 * Merge a backend refresh without discarding an in-progress mapping draft.
 * Backend-derived fields still refresh while the locally edited mapping and
 * its stale artifact markers remain authoritative until the user saves.
 */
export function mergeProjectRefresh(
  current: ProjectState,
  incoming: ProjectState,
): ProjectState {
  if (!current.mappingDirty || current.projectPath !== incoming.projectPath) {
    return incoming;
  }

  const staleArtifactIds = new Set(
    current.artifacts
      .filter((artifact) => artifact.status === "stale")
      .map((artifact) => artifact.id),
  );

  return {
    ...incoming,
    mapping: current.mapping,
    mappingDirty: true,
    artifacts: incoming.artifacts.map((artifact) =>
      staleArtifactIds.has(artifact.id)
        ? { ...artifact, status: "stale" }
        : artifact,
    ),
  };
}

/** A monotonic sequence prevents an older async refresh from winning a race. */
export function createRefreshSequence() {
  let latest = 0;
  return {
    begin(): number {
      latest += 1;
      return latest;
    },
    isCurrent(sequence: number): boolean {
      return sequence === latest;
    },
  };
}
