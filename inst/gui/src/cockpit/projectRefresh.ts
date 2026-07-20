import type { ProjectState } from "./types";

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
