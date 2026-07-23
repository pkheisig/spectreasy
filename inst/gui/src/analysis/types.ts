export type AnalysisRole = 'positive' | 'negative' | 'root' | 'terminal' | null
export type GateType = 'root' | 'rectangle' | 'ellipse' | 'polygon' | 'range'
export type AxisTransform = 'linear' | 'asinh' | 'biexponential'
export type PlotColor = 'density' | 'marker'
export type AnalysisPlotType = 'scatter' | 'histogram' | 'contour' | 'hexbin'

export type GatePoint = { x: number; y: number }

export type GateGeometry = {
  x_min?: number
  x_max?: number
  y_min?: number
  y_max?: number
  min?: number
  max?: number
  center_x?: number
  center_y?: number
  radius_x?: number
  radius_y?: number
  points?: GatePoint[]
}

export type AnalysisPopulation = {
  id: string
  name: string
  parent_id: string | null
  type: GateType
  role: AnalysisRole
  source_file?: string | null
  x?: string
  y?: string
  geometry?: GateGeometry | null
}

export type AnalysisPlotState = {
  id: string
  type: AnalysisPlotType
  population_id: string
  x: string
  y: string
  color_by: PlotColor
  color_marker?: string
  color_palette?: string
  overlay_population_id?: string
  x_transform: AxisTransform
  y_transform: AxisTransform
  point_size?: number
  opacity?: number
  x_min?: number | null
  x_max?: number | null
  y_min?: number | null
  y_max?: number | null
}

export type AnalysisWorkspaceState = {
  schema_version: number
  updated_at: string | null
  source_path: string
  selected_file: string
  active_population_id: string
  seed: number
  populations: AnalysisPopulation[]
  plots: AnalysisPlotState[]
  annotations: Array<{ population_id: string; label: string; rule?: string }>
  root_event_id: number | null
  root_population_id: string | null
  root_source_file: string | null
}

export type GateSetImportSummary = {
  current_gate_count: number
  imported_gate_count: number
  plot_references_redirected: number
  overlays_cleared: number
  annotations_removed: number
  trajectory_root_cleared: boolean
}

export type PreparedGateSetImport = {
  workspace: AnalysisWorkspaceState
  warnings: string[]
  source_name: string
  summary: GateSetImportSummary
}

export type AnalysisFile = {
  id: string
  name: string
  path: string
  size: number
  modified: string
  event_count: number | null
  parameter_count: number | null
  channels: string[]
  descriptions: string[]
  written_by: string
  error?: string
}

export type AnalysisSource = {
  id: string
  label: string
  role: 'samples' | 'controls' | 'unmixed'
  path: string
  files: AnalysisFile[]
}

export type ChannelLabel = { channel: string; marker: string }
export type AnalysisEvent = { event_id: number; x: number; y: number; color?: number }

export type AnalysisEventPayload = {
  success: boolean
  error?: string
  file: string
  population_id: string
  x: string
  y: string
  color: string
  events: AnalysisEvent[]
  population_count: number
  parent_count: number
  total_count: number
  displayed_count: number
  channels: ChannelLabel[]
}

export type AnalysisMethod = {
  id: string
  name: string
  family: 'reduction' | 'clustering' | 'trajectory'
  package: string
  runtime?: 'R' | 'python'
  available: boolean
  availability_state: 'ready' | 'setup-required' | 'unavailable'
  availability_label: string
  installed: boolean
  adapter_verified: boolean
  research_only: boolean
  version: string
  citation: string
  doi: string
  requirements: string | string[]
  blocker: string
  next_action: string
  prerequisites: string[]
  automatic_prerequisites: string[]
  user_prerequisites: string[]
  outputs: string[]
  pipeline: string[]
  supports_3d: boolean
  visible?: boolean
  parameters: AnalysisMethodParameter[]
}

export type AnalysisMethodParameterChoice = {
  value: string
  label: string
}

export type AnalysisMethodParameter = {
  id: string
  label: string
  type: 'number' | 'integer' | 'boolean' | 'select'
  default: number | string | boolean
  minimum?: number | null
  maximum?: number | null
  step?: number | null
  choices?: AnalysisMethodParameterChoice[] | null
  description?: string
}

export type AnalysisParameterValue = number | string | boolean
export type AnalysisAdvancedSettings = Record<string, Record<string, AnalysisParameterValue>>

export type PopulationStatistics = {
  population_id: string
  population_name: string
  count: number
  parent_count: number
  total_count: number
  percent_parent: number
  percent_total: number
  markers: Array<{ marker: string; median: number; mean: number; robust_sd: number }>
}

export type AnalysisRunEvent = {
  event_id: number
  source_file?: string
  source_event_id?: number
  sample_id?: number
  dimension_1: number
  dimension_2: number
  dimension_3?: number
  cluster_id?: number
  trajectory_branch?: number | string
  pseudotime?: number
  predicted_identity?: string
  identity_score?: number
  identity_margin?: number
  [dimension: `dimension_${number}`]: number | undefined
  [marker: `marker_${number}`]: number | undefined
}

export type AnalysisIdentitySignature = {
  id: string
  name: string
  color: string
  positive_markers: string[]
  negative_markers: string[]
}

export type AnalysisIdentityAnnotation = {
  identity_id: string
  analysis_id: string
  method: 'robust-signed-marker-score'
  signatures: Array<Omit<AnalysisIdentitySignature, 'id'>>
  thresholds: { min_score: number; min_margin: number; evidence_sensitivity?: number }
  counts: Record<string, number>
  event_count: number
  assigned_count: number
  unassigned_count: number
  model: AnalysisArtifactFile
  scores: AnalysisArtifactFile
  metadata: AnalysisArtifactFile
  created_at: string
}

export type AnalysisArtifactFile = {
  path: string
  sha256: string
}

export type AnalysisRunResult = {
  metadata: {
    analysis_id: string
    display_name?: string
    method: AnalysisMethod
    cluster_method?: AnalysisMethod | null
    reduction_method?: AnalysisMethod | null
    source_file: string
    source_files?: string[]
    population_id: string
    population_path: string
    markers: string[]
    marker_channels?: string[]
    marker_columns?: Array<{ marker: string; channel?: string; column: `marker_${number}` }>
    coordinate_count: number
    coordinate_labels: string[]
    seed: number
    event_count: number
    root_event_id: number | null
    runtime_seconds: number
    events_file: string
    pipeline?: Array<{ id: string; name: string; inputs: string[]; outputs: string[]; artifact_id?: string; reused?: boolean }>
    artifacts?: {
      clustering?: AnalysisArtifact | null
      embedding?: AnalysisArtifact | null
      trajectory?: AnalysisArtifact | null
    }
    identity_annotation?: AnalysisIdentityAnnotation | null
    intermediate_files?: Array<{ id: string; path: string; sha256: string }>
    plot_files?: Array<{ format: string; path: string; sha256: string }>
  }
  events: AnalysisRunEvent[]
  metadata_file: string
}

export type AnalysisArtifact = {
  id: string
  type: 'clustering' | 'embedding' | 'trajectory'
  method: AnalysisMethod
  reused: boolean
  parameters: Record<string, AnalysisParameterValue>
  object: { path: string; sha256: string }
  values: { path: string; sha256: string }
}

export type GateDraft = {
  type: 'rectangle' | 'ellipse' | 'polygon' | 'range'
  x: string
  y?: string
  geometry: GateGeometry
}
