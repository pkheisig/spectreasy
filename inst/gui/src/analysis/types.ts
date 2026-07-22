export type AnalysisRole = 'positive' | 'negative' | 'root' | 'terminal' | null
export type GateType = 'root' | 'rectangle' | 'polygon' | 'range'
export type AxisTransform = 'linear' | 'asinh'
export type PlotColor = 'density' | 'marker'

export type GatePoint = { x: number; y: number }

export type GateGeometry = {
  x_min?: number
  x_max?: number
  y_min?: number
  y_max?: number
  min?: number
  max?: number
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
  type: 'scatter' | 'density'
  population_id: string
  x: string
  y: string
  color_by: PlotColor
  color_marker?: string
  overlay_population_id?: string
  x_transform: AxisTransform
  y_transform: AxisTransform
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
  available: boolean
  installed: boolean
  adapter_verified: boolean
  research_only: boolean
  version: string
  citation: string
  doi: string
  requirements: string | string[]
}

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
  dimension_1: number
  dimension_2: number
  cluster_id?: number
  pseudotime?: number
}

export type AnalysisRunResult = {
  metadata: {
    analysis_id: string
    method: AnalysisMethod
    source_file: string
    population_id: string
    population_path: string
    markers: string[]
    seed: number
    event_count: number
    root_event_id: number | null
    runtime_seconds: number
    events_file: string
    plot_files?: Array<{ format: string; path: string; sha256: string }>
  }
  events: AnalysisRunEvent[]
  metadata_file: string
}

export type GateDraft = {
  type: 'rectangle' | 'polygon'
  x: string
  y: string
  geometry: GateGeometry
}
