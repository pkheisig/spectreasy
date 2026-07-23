import type { ProjectState } from './types'

export const emptyProject: ProjectState = {
  projectName: 'Choose a project folder',
  projectPath: '',
  cytometer: 'Auto',
  method: 'AutoSpectral',
  lastAction: '',
  lastActionAt: '',
  artifacts: [],
  mapping: [],
  mappingDirty: false,
  gatesDirty: false,
  missingInputDirs: [],
  controlInputDir: 'scc',
  sampleInputDir: 'samples',
  dataRevision: 'empty',
  matrixFiles: null,
  sampleFiles: null,
  gatingFiles: null,
  gatingMetadata: {},
  scan: { controls: 0, samples: 0, matrices: 0, reports: 0, qcMetrics: 0, spectralVariants: 0, gates: 0 },
}
