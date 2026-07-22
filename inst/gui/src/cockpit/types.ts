export type SectionId =
  | 'controls'
  | 'samples'
  | 'matrix'
  | 'panel'
  | 'af'
  | 'settings'

export type CockpitAppletId =
  | 'control-gating'
  | 'panel-builder'
  | 'matrix-adjustment'
  | 'control-qc-report'
  | 'sample-qc-report'

export type ArtifactStatus = 'current' | 'stale' | 'user' | 'missing' | 'unknown'
export type JobState = 'idle' | 'running' | 'complete' | 'failed'
export type StepState = 'complete' | 'ready' | 'warning' | 'blocked' | 'stale' | 'idle'

export const WORKFLOW_UNMIXING_METHODS = ['AutoSpectral', 'OLS', 'WLS', 'RWLS', 'NNLS'] as const

export function normalizeWorkflowUnmixingMethod(value: unknown): string {
  const method = String(value ?? '')
  return WORKFLOW_UNMIXING_METHODS.includes(method as typeof WORKFLOW_UNMIXING_METHODS[number])
    ? method
    : 'AutoSpectral'
}

export type Artifact = {
  id: string
  name: string
  type: string
  group: string
  status: ArtifactStatus
  detail: string
  path: string
  relativePath?: string
  updated: string
  updatedEpoch?: number
  run?: string
}

export type MappingRow = {
  id: string
  file: string
  fluorophore: string
  marker: string
  channel: string
  controlType: 'cell' | 'bead'
  universalNegative: string
  isViability?: boolean
  warning?: string
  ignored?: boolean
  ignoredReason?: string
}

export type ProjectState = {
  projectName: string
  projectPath: string
  cytometer: string
  method: string
  lastAction: string
  lastActionAt: string
  artifacts: Artifact[]
  mapping: MappingRow[]
  mappingDirty: boolean
  gatesDirty: boolean
  missingInputDirs: string[]
  controlInputDir: string
  sampleInputDir: string
  dataRevision: string
  matrixFiles: string[] | null
  sampleFiles: string[] | null
  gatingFiles: Array<Record<string, unknown>> | null
  gatingMetadata: Record<string, unknown>
  scan: {
    controls: number
    samples: number
    matrices: number
    reports: number
    qcMetrics: number
    spectralVariants: number
    gates: number
  }
}

export type BackendStatus = {
  connected: boolean
  version: string
  mode: string
  message: string
  apiPort: string
  packageReady: boolean
}

export type Job = {
  label: string
  state: JobState
  progress: number
  subtask: string
  startedAt?: string
  finishedAt?: string
  output?: string
}

export type ExecutionLogEntry = {
  id: string
  kind: 'command' | 'info' | 'success' | 'warning' | 'error'
  text: string
}

export type PanelPayload = {
  cytometer: string
  configuration: string
  fluorophores: Array<{ fluorophore: string; peak_detector: string; peak_laser: string }>
  selected: string[]
  complexity_index: number | null
  peak_detectors: string[]
}

export type Report = {
  id: string
  title: string
  type: 'Control QC' | 'Sample QC' | 'Panel overview'
  format: 'HTML' | 'PDF'
  run: string
  created: string
  createdEpoch?: number
  status: 'current' | 'stale'
  matrix: string
  path?: string
  promptPath?: string
  promptBytes?: number
}

export type ControlSettings = {
  sccDir: string
  controlFile: string
  outputDir: string
  method: string
  cytometer: string
  autoCreateMapping: boolean
  autoUnknownFluorPolicy: 'by_channel' | 'empty' | 'filename'
  manualGateFile: string
  afNBands: number
  afMaxCells: number
  defaultSampleType: 'beads' | 'cells'
  histogramPctBeads: number
  histogramDirectionBeads: 'right' | 'both' | 'left'
  histogramPctCells: number
  histogramDirectionCells: 'right' | 'both' | 'left'
  outlierPercentile: number
  debrisPercentile: number
  beadGateScale: number
  maxClusters: number
  minClusterProportion: number
  gateContourBeads: number
  gateContourCells: number
  subsampleN: number
  unmixScatterPanelSizeMm: number
  rwlsMaxIter: number
  unmixThreads: number
  seed: number
  saveQcPlots: boolean
  saveReport: boolean
  outputFormat: 'html' | 'pdf'
  sccBackgroundMethod: 'scatter_knn' | 'none'
  sccBackgroundK: number
  spectralVariantSomNodes: number
  spectralVariantTopK: number
  spectralVariantCosineThreshold: number
  spectralVariantMaxVariants: number
  spectralVariantMinEvents: number
  autospectralNCandidates: number
  autospectralNSpectral: number
  autospectralMinEvents: number
  refine: boolean
}

export type SampleSettings = {
  sampleDir: string
  matrixFile: string
  detectorNoiseFile: string
  outputDir: string
  method: string
  rwlsMaxIter: number
  nThreads: number
  spectralVariantTopK: number
  spectralVariantMinAbundance: number
  spectralVariantPositiveFraction: number
  spectralVariantMinImprovement: number
  spectralVariantLibraryFile: string
  estimateAf: boolean
  writeFcs: boolean
  saveReport: boolean
  outputFormat: 'html' | 'pdf'
  reportPerSample: boolean
  saveQcPlots: boolean
  plotNEvents: number
  chunkSize: number
  seed: number
  returnType: 'list' | 'flowSet' | 'SingleCellExperiment'
}

export type AfSettings = {
  fcsFile: string
  saveName: string
  saveOverwrite: boolean
  afNBands: number
  afMaxCells: number
  seed: number
}

export type AppearanceSettings = {
  theme: 'light' | 'dark'
  density: 'compact' | 'comfortable' | 'spacious'
  fontScale: number
  fontFamily: 'avenir' | 'futura' | 'atkinson' | 'charter' | 'palatino' | 'monaco' | 'system'
  sidebarWidth: number
  cornerRadius: number
  shadows: 'none' | 'subtle' | 'raised'
  backgroundTexture: boolean
}

export const INTERFACE_SCALE_MIN = 55
export const INTERFACE_SCALE_MAX = 95
export const INTERFACE_SCALE_STEP = 5
export const INTERFACE_SCALE_DEFAULT = 75

export function normalizeInterfaceScale(value: unknown): number {
  const numeric = typeof value === 'number' && Number.isFinite(value)
    ? value
    : INTERFACE_SCALE_DEFAULT
  const rounded = Math.round(numeric / INTERFACE_SCALE_STEP) * INTERFACE_SCALE_STEP
  return Math.min(INTERFACE_SCALE_MAX, Math.max(INTERFACE_SCALE_MIN, rounded))
}

export function interfaceScaleLevel(value: unknown): number {
  return Math.round((normalizeInterfaceScale(value) - INTERFACE_SCALE_MIN) / INTERFACE_SCALE_STEP) + 1
}

export type WorkflowSettings = {
  projectPath: string
  control: ControlSettings
  sample: SampleSettings
  af: AfSettings
  appearance: AppearanceSettings
}

export function defaultWorkflowSettings(projectPath: string): WorkflowSettings {
  return {
    projectPath,
    control: {
      sccDir: 'scc',
      controlFile: 'fcs_mapping.csv',
      outputDir: 'spectreasy_outputs',
      method: 'AutoSpectral',
      cytometer: 'auto',
      autoCreateMapping: true,
      autoUnknownFluorPolicy: 'by_channel',
      manualGateFile: '',
      afNBands: 100,
      afMaxCells: 50000,
      defaultSampleType: 'beads',
      histogramPctBeads: 0.98,
      histogramDirectionBeads: 'right',
      histogramPctCells: 0.35,
      histogramDirectionCells: 'right',
      outlierPercentile: 0.02,
      debrisPercentile: 0.08,
      beadGateScale: 1.3,
      maxClusters: 10,
      minClusterProportion: 0.03,
      gateContourBeads: 0.95,
      gateContourCells: 0.9,
      subsampleN: 5000,
      unmixScatterPanelSizeMm: 30,
      rwlsMaxIter: 1,
      unmixThreads: 1,
      seed: 1,
      saveQcPlots: false,
      saveReport: true,
      outputFormat: 'html',
      sccBackgroundMethod: 'scatter_knn',
      sccBackgroundK: 2,
      spectralVariantSomNodes: 16,
      spectralVariantTopK: 3,
      spectralVariantCosineThreshold: 0.98,
      spectralVariantMaxVariants: 8,
      spectralVariantMinEvents: 50,
      autospectralNCandidates: 1000,
      autospectralNSpectral: 200,
      autospectralMinEvents: 10,
      refine: false,
    },
    sample: {
      sampleDir: 'samples',
      matrixFile: 'spectreasy_outputs/unmix_controls/scc_reference_matrix.csv',
      detectorNoiseFile: 'spectreasy_outputs/unmix_controls/scc_detector_noise.csv',
      outputDir: 'spectreasy_outputs',
      method: 'AutoSpectral',
      rwlsMaxIter: 1,
      nThreads: 1,
      spectralVariantTopK: 3,
      spectralVariantMinAbundance: 1,
      spectralVariantPositiveFraction: 0.02,
      spectralVariantMinImprovement: 0.01,
      spectralVariantLibraryFile: '',
      estimateAf: false,
      writeFcs: true,
      saveReport: true,
      outputFormat: 'html',
      reportPerSample: false,
      saveQcPlots: false,
      plotNEvents: 10000,
      chunkSize: 50000,
      seed: 1,
      returnType: 'list',
    },
    af: {
      fcsFile: 'scc/unstained_cells.fcs',
      saveName: '',
      saveOverwrite: false,
      afNBands: 100,
      afMaxCells: 50000,
      seed: 1,
    },
    appearance: {
      theme: 'light',
      density: 'comfortable',
      fontScale: INTERFACE_SCALE_DEFAULT,
      fontFamily: 'avenir',
      sidebarWidth: 242,
      cornerRadius: 8,
      shadows: 'subtle',
      backgroundTexture: false,
    },
  }
}
