import { useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  Boxes,
  Check,
  ChevronRight,
  ChevronDown,
  CircleAlert,
  CircleDashed,
  CircleHelp,
  Crosshair,
  Download,
  FileDown,
  FileUp,
  FlaskConical,
  FolderOpen,
  MousePointer2,
  Pentagon,
  Plus,
  Save,
  Search,
  SlidersHorizontal,
  SquareDashed,
  Trash2,
  Undo2,
  Redo2,
  X,
} from 'lucide-react'
import { analysisRequest } from './api'
import { AnalysisAdvancedSettingsPanel } from './AnalysisAdvancedSettings'
import { advancedSettingsIssues, resolvedAdvancedSettings } from './advancedSettings'
import { AnalysisGuideDialog } from './AnalysisGuideDialog'
import { AnalysisResultPlots } from './AnalysisResultPlots'
import { CellIdentityPanel } from './CellIdentityPanel'
import { CONTINUOUS_PALETTES, paletteGradient } from './resultColor'
import { PlotCard } from './PlotCard'
import type { PlotTool } from './AnalysisPlot'
import type {
  AnalysisAdvancedSettings,
  AnalysisMethod,
  AnalysisParameterValue,
  AnalysisPlotState,
  AnalysisPopulation,
  AnalysisRunResult,
  AnalysisSource,
  AnalysisWorkspaceState,
  ChannelLabel,
  GateDraft,
  PopulationStatistics,
} from './types'
import './analysis.css'

type AnalysisWorkspaceProps = {
  projectPath?: string
  cockpitTheme?: 'light' | 'dark'
  onRequestExit?: () => void
}

type SaveState = 'saved' | 'saving' | 'error'
type InspectorTab = 'plot' | 'population' | 'export'
type AnalysisMode = 'explore' | 'trajectory'
type AnalysisResultTab = 'plots' | 'identities'

function asArray<T>(value: T | T[] | null | undefined): T[] {
  if (value == null) return []
  return Array.isArray(value) ? value : [value]
}

function scalarString(value: unknown, fallback = ''): string {
  if (Array.isArray(value)) return scalarString(value[0], fallback)
  return typeof value === 'string' || typeof value === 'number' ? String(value) : fallback
}

function scalarNumber(value: unknown, fallback = 0): number {
  const parsed = Number(Array.isArray(value) ? value[0] : value)
  return Number.isFinite(parsed) ? parsed : fallback
}

function scalarBoolean(value: unknown): boolean {
  const unboxed = Array.isArray(value) ? value[0] : value
  return unboxed === true || unboxed === 1 || unboxed === 'true'
}

function normalizeSource(source: AnalysisSource): AnalysisSource {
  return {
    ...source,
    id: scalarString(source.id),
    label: scalarString(source.label),
    role: scalarString(source.role, 'samples') as AnalysisSource['role'],
    path: scalarString(source.path),
    files: asArray(source.files).map((file) => ({
      ...file,
      id: scalarString(file.id),
      name: scalarString(file.name),
      path: scalarString(file.path),
      size: scalarNumber(file.size),
      modified: scalarString(file.modified),
      event_count: scalarNumber(file.event_count),
      parameter_count: scalarNumber(file.parameter_count),
      channels: asArray(file.channels).map((value) => scalarString(value)),
      descriptions: asArray(file.descriptions).map((value) => scalarString(value)),
      written_by: scalarString(file.written_by),
      error: scalarString(file.error),
    })),
  }
}

function normalizeWorkspace(value: AnalysisWorkspaceState): AnalysisWorkspaceState {
  const role = (candidate: unknown): AnalysisPopulation['role'] => {
    const text = scalarString(candidate)
    return text === 'positive' || text === 'negative' || text === 'root' || text === 'terminal' ? text : null
  }
  const rootEvent = scalarNumber(value.root_event_id, Number.NaN)
  return {
    ...value,
    schema_version: scalarNumber(value.schema_version, 2),
    source_path: scalarString(value.source_path),
    selected_file: scalarString(value.selected_file),
    active_population_id: scalarString(value.active_population_id, 'root'),
    seed: scalarNumber(value.seed, 20260723),
    root_event_id: Number.isInteger(rootEvent) && rootEvent > 0 ? rootEvent : null,
    root_population_id: scalarString(value.root_population_id) || null,
    root_source_file: scalarString(value.root_source_file) || null,
    populations: asArray(value.populations).map((population) => ({
      ...population,
      id: scalarString(population.id),
      name: scalarString(population.name),
      parent_id: scalarString(population.parent_id) || null,
      type: scalarString(population.type, 'polygon') as AnalysisPopulation['type'],
      role: role(population.role),
      source_file: scalarString(population.source_file) || null,
      x: scalarString(population.x),
      y: scalarString(population.y),
    })),
    plots: asArray(value.plots).map((plot) => ({
      ...plot,
      id: scalarString(plot.id),
      type: (scalarString(plot.type, 'scatter') === 'density' ? 'scatter' : scalarString(plot.type, 'scatter')) as AnalysisPlotState['type'],
      population_id: scalarString(plot.population_id, 'root'),
      x: scalarString(plot.x),
      y: scalarString(plot.y),
      color_by: scalarString(plot.color_by, 'density') as AnalysisPlotState['color_by'],
      color_marker: scalarString(plot.color_marker),
      color_palette: scalarString(plot.color_palette, 'viridis'),
      overlay_population_id: scalarString(plot.overlay_population_id),
      x_transform: scalarString(plot.x_transform, 'linear') as AnalysisPlotState['x_transform'],
      y_transform: scalarString(plot.y_transform, 'linear') as AnalysisPlotState['y_transform'],
      point_size: scalarNumber(plot.point_size, 2.4),
      opacity: scalarNumber(plot.opacity, 0.82),
      x_min: plot.x_min == null ? null : scalarNumber(plot.x_min),
      x_max: plot.x_max == null ? null : scalarNumber(plot.x_max),
      y_min: plot.y_min == null ? null : scalarNumber(plot.y_min),
      y_max: plot.y_max == null ? null : scalarNumber(plot.y_max),
    })),
    annotations: asArray(value.annotations),
  }
}

function cloneWorkspace(value: AnalysisWorkspaceState): AnalysisWorkspaceState {
  return typeof structuredClone === 'function'
    ? structuredClone(value)
    : JSON.parse(JSON.stringify(value)) as AnalysisWorkspaceState
}

function normalizeMethod(method: AnalysisMethod): AnalysisMethod {
  return {
    ...method,
    id: scalarString(method.id),
    name: scalarString(method.name),
    family: scalarString(method.family) as AnalysisMethod['family'],
    package: scalarString(method.package),
    available: scalarBoolean(method.available),
    availability_state: scalarString(method.availability_state, method.available ? 'ready' : 'unavailable') as AnalysisMethod['availability_state'],
    availability_label: scalarString(method.availability_label, method.available ? 'Ready' : 'Unavailable'),
    installed: scalarBoolean(method.installed),
    adapter_verified: scalarBoolean(method.adapter_verified),
    research_only: scalarBoolean(method.research_only),
    runtime: scalarString(method.runtime, 'R') as AnalysisMethod['runtime'],
    visible: method.visible == null ? true : scalarBoolean(method.visible),
    version: scalarString(method.version),
    citation: scalarString(method.citation),
    doi: scalarString(method.doi),
    requirements: asArray(method.requirements).map((value) => scalarString(value)),
    blocker: scalarString(method.blocker, methodRequirement(method)),
    next_action: scalarString(method.next_action),
    prerequisites: asArray(method.prerequisites).map((value) => scalarString(value)),
    automatic_prerequisites: asArray(method.automatic_prerequisites).map((value) => scalarString(value)),
    user_prerequisites: asArray(method.user_prerequisites).map((value) => scalarString(value)),
    outputs: asArray(method.outputs).map((value) => scalarString(value)),
    pipeline: asArray(method.pipeline).map((value) => scalarString(value)),
    supports_3d: scalarBoolean(method.supports_3d),
    parameters: asArray(method.parameters).map((parameter) => ({
      ...parameter,
      id: scalarString(parameter.id),
      label: scalarString(parameter.label),
      type: scalarString(parameter.type, 'number') as typeof parameter.type,
      default: parameter.default as AnalysisParameterValue,
      minimum: parameter.minimum == null ? null : scalarNumber(parameter.minimum),
      maximum: parameter.maximum == null ? null : scalarNumber(parameter.maximum),
      step: parameter.step == null ? null : scalarNumber(parameter.step),
      choices: asArray(parameter.choices).map((choice) => ({
        value: scalarString(choice.value),
        label: scalarString(choice.label),
      })),
      description: scalarString(parameter.description),
    })),
  }
}

function populationChildren(populations: AnalysisPopulation[], parentId: string) {
  return populations.filter((population) => population.parent_id === parentId)
}

function descendants(populations: AnalysisPopulation[], populationId: string): Set<string> {
  const result = new Set<string>([populationId])
  let changed = true
  while (changed) {
    changed = false
    for (const population of populations) {
      if (population.parent_id && result.has(population.parent_id) && !result.has(population.id)) {
        result.add(population.id)
        changed = true
      }
    }
  }
  return result
}

function populationContainsRoot(populations: AnalysisPopulation[], populationId: string, rootPopulationId: string | null) {
  let current = rootPopulationId
  const visited = new Set<string>()
  while (current && !visited.has(current)) {
    if (current === populationId) return true
    visited.add(current)
    current = populations.find((population) => population.id === current)?.parent_id ?? null
  }
  return false
}

function channelLabels(sourceFile: AnalysisSource['files'][number] | undefined): ChannelLabel[] {
  if (!sourceFile) return []
  return sourceFile.channels.map((channel, index) => ({ channel, marker: sourceFile.descriptions[index] || '' }))
}

function channelLabel(channel: ChannelLabel) {
  return channel.marker && channel.marker !== channel.channel
    ? `${channel.marker} · ${channel.channel}`
    : channel.channel
}

function preferredAxes(channels: string[]): [string, string] {
  const pick = (patterns: RegExp[], fallback: string) => patterns.map((pattern) => channels.find((channel) => pattern.test(channel))).find(Boolean) || fallback
  const x = pick([/^FSC-A$/i, /^FSC.*-A$/i], channels[0] || '')
  const y = pick([/^(?:[A-Z])?SSC-A$/i, /SSC.*-A$/i], channels.find((channel) => channel !== x) || channels[1] || x)
  return [x, y]
}

function preferredMarkers(file: AnalysisSource['files'][number] | undefined) {
  if (!file) return []
  const markerChannels = file.channels.filter((channel, index) => file.descriptions[index] && !/^(time|fsc|ssc)/i.test(channel))
  return (markerChannels.length ? markerChannels : file.channels).slice(0, 8)
}

function methodRequirement(method: AnalysisMethod) {
  const requirements = asArray(method.requirements).filter(Boolean)
  return requirements.join('; ')
}

function PopulationTree({
  populations,
  parentId,
  activeId,
  depth = 0,
  onSelect,
  collapsed,
  onToggle,
  visibleIds,
}: {
  populations: AnalysisPopulation[]
  parentId: string
  activeId: string
  depth?: number
  onSelect: (id: string) => void
  collapsed: Set<string>
  onToggle: (id: string) => void
  visibleIds: Set<string> | null
}) {
  return populationChildren(populations, parentId).filter((population) => !visibleIds || visibleIds.has(population.id)).map((population) => {
    const children = populationChildren(populations, population.id)
    const isCollapsed = collapsed.has(population.id)
    return (
      <div key={population.id}>
        <button
          type="button"
          className={`analysis-population-row ${activeId === population.id ? 'is-active' : ''}`}
          style={{ paddingLeft: 12 + depth * 16 }}
          onClick={() => onSelect(population.id)}
          onDoubleClick={() => onSelect(population.id)}
        >
          <span
            className={`analysis-tree-toggle ${children.length ? '' : 'is-leaf'}`}
            role={children.length ? 'button' : undefined}
            tabIndex={children.length ? 0 : -1}
            aria-label={children.length ? `${isCollapsed ? 'Expand' : 'Collapse'} ${population.name}` : undefined}
            onClick={(event) => { if (children.length) { event.stopPropagation(); onToggle(population.id) } }}
            onKeyDown={(event) => {
              if (children.length && (event.key === 'Enter' || event.key === ' ')) {
                event.preventDefault()
                event.stopPropagation()
                onToggle(population.id)
              }
            }}
          >
            {isCollapsed ? <ChevronRight size={12} /> : <ChevronDown size={12} />}
          </span>
          <span className={`population-swatch role-${population.role || 'none'}`} />
          <span>{population.name}</span>
          {population.role ? <small>{population.role === 'positive' ? 'POS' : population.role === 'negative' ? 'NEG' : population.role}</small> : null}
        </button>
        {children.length && !isCollapsed ? (
          <PopulationTree populations={populations} parentId={population.id} activeId={activeId} depth={depth + 1} onSelect={onSelect} collapsed={collapsed} onToggle={onToggle} visibleIds={visibleIds} />
        ) : null}
      </div>
    )
  })
}

export default function AnalysisWorkspace({ projectPath = '', cockpitTheme = 'light', onRequestExit }: AnalysisWorkspaceProps) {
  const [sources, setSources] = useState<AnalysisSource[]>([])
  const [methods, setMethods] = useState<AnalysisMethod[]>([])
  const [workspace, setWorkspace] = useState<AnalysisWorkspaceState | null>(null)
  const [ready, setReady] = useState(false)
  const [status, setStatus] = useState('Opening analysis workspace')
  const [saveState, setSaveState] = useState<SaveState>('saved')
  const [tool, setTool] = useState<PlotTool>('select')
  const [tab, setTab] = useState<InspectorTab>('plot')
  const [selectedPlotId, setSelectedPlotId] = useState('plot-1')
  const [analysisDialogOpen, setAnalysisDialogOpen] = useState(false)
  const [analysisGuideOpen, setAnalysisGuideOpen] = useState(false)
  const [gateDraft, setGateDraft] = useState<GateDraft | null>(null)
  const [gateName, setGateName] = useState('Population')
  const [gateRole, setGateRole] = useState('')
  const [fileSpecificGate, setFileSpecificGate] = useState(false)
  const [statistics, setStatistics] = useState<PopulationStatistics | null>(null)
  const [statisticsError, setStatisticsError] = useState('')
  const [stainingMarker, setStainingMarker] = useState('')
  const [stainingResult, setStainingResult] = useState<Record<string, number | string> | null>(null)
  const [actionMessage, setActionMessage] = useState('')
  const [exportFormat, setExportFormat] = useState('fcs')
  const [exportCount, setExportCount] = useState(0)
  const [exportFolder, setExportFolder] = useState('spectreasy_outputs/analysis/exports')
  const [exportAllFiles, setExportAllFiles] = useState(false)
  const [analysisMode, setAnalysisMode] = useState<AnalysisMode>('explore')
  const [selectedClusterMethod, setSelectedClusterMethod] = useState('')
  const [selectedReductionMethod, setSelectedReductionMethod] = useState('')
  const [selectedTrajectoryMethod, setSelectedTrajectoryMethod] = useState('')
  const [selectedMarkers, setSelectedMarkers] = useState<string[]>([])
  const [methodRunning, setMethodRunning] = useState(false)
  const [methodJobId, setMethodJobId] = useState('')
  const [methodProgress, setMethodProgress] = useState('')
  const [methodResult, setMethodResult] = useState<AnalysisRunResult | null>(null)
  const [analysisResultTab, setAnalysisResultTab] = useState<AnalysisResultTab>('plots')
  const [advancedSettings, setAdvancedSettings] = useState<AnalysisAdvancedSettings>({})
  const [analysisMaxEvents, setAnalysisMaxEvents] = useState(20000)
  const [analysisCofactor, setAnalysisCofactor] = useState(150)
  const [analysisFiles, setAnalysisFiles] = useState<string[]>([])
  const [populationQuery, setPopulationQuery] = useState('')
  const [collapsedPopulationIds, setCollapsedPopulationIds] = useState<Set<string>>(() => new Set())
  const [inspectorOpen, setInspectorOpen] = useState(
    () => typeof window === 'undefined' || !window.matchMedia('(max-width: 850px)').matches,
  )
  const undoStack = useRef<AnalysisWorkspaceState[]>([])
  const redoStack = useRef<AnalysisWorkspaceState[]>([])
  const importWorkspaceInput = useRef<HTMLInputElement>(null)
  const [, setHistoryRevision] = useState(0)

  const bumpHistory = useCallback(() => setHistoryRevision((revision) => revision + 1), [])

  useEffect(() => {
    const root = document.documentElement
    const body = document.body
    const host = document.getElementById('root')
    const previous = {
      rootOverflow: root.style.overflow,
      rootWidth: root.style.width,
      bodyOverflow: body.style.overflow,
      bodyWidth: body.style.width,
      bodyMaxWidth: body.style.maxWidth,
      bodyMinWidth: body.style.minWidth,
      hostCssText: host?.style.cssText ?? '',
    }
    root.style.overflow = 'hidden'
    root.style.width = '100%'
    body.style.overflow = 'hidden'
    body.style.width = '100vw'
    body.style.maxWidth = '100vw'
    body.style.minWidth = '0'
    if (host) {
      host.style.position = 'fixed'
      host.style.inset = '0'
      host.style.width = '100vw'
      host.style.maxWidth = '100vw'
      host.style.overflow = 'hidden'
    }
    return () => {
      root.style.overflow = previous.rootOverflow
      root.style.width = previous.rootWidth
      body.style.overflow = previous.bodyOverflow
      body.style.width = previous.bodyWidth
      body.style.maxWidth = previous.bodyMaxWidth
      body.style.minWidth = previous.bodyMinWidth
      if (host) host.style.cssText = previous.hostCssText
    }
  }, [])

  const updateWorkspace = useCallback((update: (current: AnalysisWorkspaceState) => AnalysisWorkspaceState) => {
    setWorkspace((current) => {
      if (!current) return current
      const next = update(current)
      if (next === current) return current
      undoStack.current.push(cloneWorkspace(current))
      if (undoStack.current.length > 100) undoStack.current.shift()
      redoStack.current = []
      bumpHistory()
      setSaveState('saving')
      return next
    })
  }, [bumpHistory])

  const replaceWorkspace = useCallback((next: AnalysisWorkspaceState, recordHistory = true) => {
    setWorkspace((current) => {
      if (recordHistory && current) {
        undoStack.current.push(cloneWorkspace(current))
        if (undoStack.current.length > 100) undoStack.current.shift()
      } else {
        undoStack.current = []
      }
      redoStack.current = []
      bumpHistory()
      setSaveState('saving')
      return cloneWorkspace(next)
    })
  }, [bumpHistory])

  const undoWorkspace = useCallback(() => {
    setWorkspace((current) => {
      const previous = undoStack.current.pop()
      if (!current || !previous) return current
      redoStack.current.push(cloneWorkspace(current))
      bumpHistory()
      setSaveState('saving')
      return previous
    })
  }, [bumpHistory])

  const redoWorkspace = useCallback(() => {
    setWorkspace((current) => {
      const next = redoStack.current.pop()
      if (!current || !next) return current
      undoStack.current.push(cloneWorkspace(current))
      bumpHistory()
      setSaveState('saving')
      return next
    })
  }, [bumpHistory])

  useEffect(() => {
    let cancelled = false
    void Promise.all([
      analysisRequest<{ success: boolean; sources: AnalysisSource[] }>('/analysis/sources', projectPath),
      analysisRequest<{ success: boolean; workspace: AnalysisWorkspaceState }>('/analysis/workspace', projectPath),
      analysisRequest<{ success: boolean; methods: AnalysisMethod[] }>('/analysis/methods', projectPath),
    ]).then(([sourceResponse, workspaceResponse, methodResponse]) => {
      if (cancelled) return
      const loadedSources = asArray(sourceResponse.sources).map(normalizeSource)
      const loadedMethods = asArray(methodResponse.methods).map(normalizeMethod)
      let loaded = normalizeWorkspace(workspaceResponse.workspace)
      const preferredSource = loadedSources.find((source) => source.path === loaded.source_path)
        ?? loadedSources.find((source) => source.role === 'unmixed')
        ?? loadedSources.find((source) => source.role === 'samples')
        ?? loadedSources[0]
      const selected = preferredSource?.files.find((file) => file.path === loaded.selected_file) ?? preferredSource?.files[0]
      const [x, y] = preferredAxes(selected?.channels ?? [])
      loaded = {
        ...loaded,
        source_path: preferredSource?.path ?? '',
        selected_file: selected?.path ?? '',
        plots: loaded.plots.map((plot) => ({
          ...plot,
          x: selected?.channels.includes(plot.x) ? plot.x : x,
          y: selected?.channels.includes(plot.y) ? plot.y : y,
        })),
      }
      setSources(loadedSources)
      setMethods(loadedMethods)
      setSelectedClusterMethod('')
      setSelectedReductionMethod('')
      setSelectedTrajectoryMethod('')
      setSelectedMarkers(preferredMarkers(selected))
      setStainingMarker(preferredMarkers(selected)[0] ?? selected?.channels[0] ?? '')
      setAnalysisFiles(selected?.path ? [selected.path] : [])
      setWorkspace(loaded)
      undoStack.current = []
      redoStack.current = []
      bumpHistory()
      setSelectedPlotId(loaded.plots[0]?.id ?? 'plot-1')
      setReady(true)
      setStatus(loadedSources.length ? 'Ready' : 'No FCS sources found')
    }).catch((reason: unknown) => {
      if (!cancelled) setStatus(reason instanceof Error ? reason.message : String(reason))
    })
    return () => { cancelled = true }
  }, [bumpHistory, projectPath])

  useEffect(() => {
    const onKeyDown = (event: KeyboardEvent) => {
      const target = event.target as HTMLElement | null
      if (target?.matches('input, textarea, select, [contenteditable="true"]')) return
      const modifier = event.metaKey || event.ctrlKey
      if (!modifier) return
      if (event.key.toLocaleLowerCase() === 'z') {
        event.preventDefault()
        if (event.shiftKey) redoWorkspace()
        else undoWorkspace()
      } else if (event.ctrlKey && event.key.toLocaleLowerCase() === 'y') {
        event.preventDefault()
        redoWorkspace()
      }
    }
    window.addEventListener('keydown', onKeyDown)
    return () => window.removeEventListener('keydown', onKeyDown)
  }, [redoWorkspace, undoWorkspace])

  useEffect(() => {
    if (!ready || !workspace || saveState !== 'saving') return
    const timer = window.setTimeout(() => {
      void analysisRequest<{ success: boolean }>('/analysis/workspace', projectPath, {
        method: 'POST', body: JSON.stringify({ workspace }),
      }).then(() => setSaveState('saved')).catch(() => setSaveState('error'))
    }, 450)
    return () => window.clearTimeout(timer)
  }, [projectPath, ready, saveState, workspace])

  const selectedSource = useMemo(
    () => sources.find((source) => source.path === workspace?.source_path) ?? sources[0],
    [sources, workspace?.source_path],
  )
  const selectedFile = useMemo(
    () => selectedSource?.files.find((file) => file.path === workspace?.selected_file) ?? selectedSource?.files[0],
    [selectedSource, workspace?.selected_file],
  )
  const channels = useMemo(() => channelLabels(selectedFile), [selectedFile])
  const activePopulation = useMemo(
    () => workspace?.populations.find((population) => population.id === workspace.active_population_id) ?? workspace?.populations[0],
    [workspace],
  )
  const selectedPlot = useMemo(
    () => workspace?.plots.find((plot) => plot.id === selectedPlotId) ?? workspace?.plots[0],
    [selectedPlotId, workspace],
  )
  const selectedCluster = useMemo(
    () => methods.find((method) => method.id === selectedClusterMethod),
    [methods, selectedClusterMethod],
  )
  const selectedReduction = useMemo(
    () => methods.find((method) => method.id === selectedReductionMethod),
    [methods, selectedReductionMethod],
  )
  const selectedTrajectory = useMemo(
    () => methods.find((method) => method.id === selectedTrajectoryMethod),
    [methods, selectedTrajectoryMethod],
  )
  const selectedAnalysisMethods = useMemo(() => {
    if (analysisMode === 'explore') {
      return [selectedClusterMethod === 'none' ? undefined : selectedCluster, selectedReduction]
        .filter((method): method is AnalysisMethod => Boolean(method))
    }
    return [
      selectedTrajectory,
      selectedTrajectory?.prerequisites.includes('clustering') ? selectedCluster : undefined,
      selectedTrajectory?.prerequisites.includes('reduction') ? selectedReduction : undefined,
    ].filter((method): method is AnalysisMethod => Boolean(method))
  }, [analysisMode, selectedCluster, selectedClusterMethod, selectedReduction, selectedTrajectory])
  const visiblePopulationIds = useMemo(() => {
    const query = populationQuery.trim().toLocaleLowerCase()
    if (!query || !workspace) return null
    const result = new Set<string>()
    for (const population of workspace.populations) {
      if (!population.name.toLocaleLowerCase().includes(query)) continue
      let current: AnalysisPopulation | undefined = population
      while (current) {
        result.add(current.id)
        current = current.parent_id
          ? workspace.populations.find((candidate) => candidate.id === current?.parent_id)
          : undefined
      }
    }
    return result
  }, [populationQuery, workspace])
  const methodIssues = (() => {
    const issues: string[] = []
    const eventCountKnown = statistics?.population_id === activePopulation?.id
    const eventCountReady = !eventCountKnown || (statistics?.count ?? 0) >= 20
    const rootReady = Boolean(
      workspace?.root_event_id
      && workspace.root_source_file === selectedFile?.path
      && activePopulation
      && populationContainsRoot(workspace.populations, activePopulation.id, workspace.root_population_id),
    )
    if (analysisMode === 'explore') {
      if (!selectedClusterMethod) issues.push('Choose a clustering method or explicitly skip clustering.')
      if (selectedClusterMethod !== 'none' && !selectedCluster?.available) issues.push('Choose an enabled clustering method.')
      if (!selectedReductionMethod) issues.push('Choose a dimensional-reduction method.')
      if (selectedReductionMethod && !selectedReduction?.available) issues.push('Choose an enabled dimensional-reduction method.')
    } else {
      if (!selectedTrajectoryMethod) issues.push('Choose a trajectory method.')
      if (selectedTrajectoryMethod && !selectedTrajectory?.available) issues.push('Choose an enabled trajectory method.')
      if (selectedTrajectory?.prerequisites.includes('clustering')) {
        if (!selectedClusterMethod || selectedClusterMethod === 'none') issues.push(`${selectedTrajectory.name} requires clustering.`)
        else if (!selectedCluster?.available) issues.push('Choose an enabled clustering method.')
      }
      if (selectedTrajectory?.prerequisites.includes('reduction')) {
        if (!selectedReductionMethod) issues.push(`${selectedTrajectory.name} requires a map.`)
        else if (!selectedReduction?.available) issues.push('Choose an enabled dimensional-reduction method.')
      }
      const rootRequired = selectedTrajectory?.user_prerequisites.includes('trajectory-root') ?? false
      if (rootRequired && !rootReady) {
        issues.push(workspace?.root_event_id
          ? 'The trajectory root belongs to another file or population. Set a new root here.'
          : 'Use Set root in the gating toolbar, then click an event in this population.')
      }
    }
    if (selectedMarkers.length < 2) issues.push('Select at least two markers.')
    if (!analysisFiles.length) issues.push('Select at least one sample file.')
    if (!eventCountReady) issues.push('The selected population needs at least 20 events.')
    issues.push(...advancedSettingsIssues(selectedAnalysisMethods, advancedSettings))
    return issues
  })()
  const methodCanRun = methodIssues.length === 0

  useEffect(() => {
    if (!workspace || !selectedFile || !activePopulation) return
    const controller = new AbortController()
    const markers = Array.from(new Set(workspace.plots.flatMap((plot) => [plot.x, plot.y]))).filter(Boolean)
    void analysisRequest<{ success: boolean; result: PopulationStatistics }>('/analysis/statistics', projectPath, {
      method: 'POST', signal: controller.signal,
      body: JSON.stringify({ file: selectedFile.path, populationId: activePopulation.id, markers }),
    }).then((response) => {
      setStatistics(response.result)
      setStatisticsError('')
    }).catch((reason: unknown) => {
      if (!controller.signal.aborted) setStatisticsError(reason instanceof Error ? reason.message : String(reason))
    })
    return () => controller.abort()
  }, [activePopulation, projectPath, selectedFile, workspace])

  function selectPopulation(id: string, assignToPlot = true) {
    updateWorkspace((current) => ({
      ...current,
      active_population_id: id,
      plots: assignToPlot
        ? current.plots.map((plot) => plot.id === selectedPlotId ? { ...plot, population_id: id } : plot)
        : current.plots,
    }))
    setStainingResult(null)
  }

  function updateSelectedPlot(patch: Partial<AnalysisPlotState>) {
    if (!selectedPlot) return
    updateWorkspace((current) => ({
      ...current,
      plots: current.plots.map((plot) => plot.id === selectedPlot.id ? { ...plot, ...patch } : plot),
    }))
  }

  function addPlot() {
    if (!workspace) return
    const [x, y] = preferredAxes(channels.map((channel) => channel.channel))
    const id = `plot-${crypto.randomUUID()}`
    const plot: AnalysisPlotState = {
      id,
      type: 'scatter',
      population_id: workspace.active_population_id,
      x,
      y,
      color_by: 'density',
      color_palette: 'viridis',
      x_transform: 'linear',
      y_transform: 'linear',
      point_size: 2.4,
      opacity: 0.82,
    }
    updateWorkspace((current) => ({ ...current, plots: [...current.plots, plot] }))
    setSelectedPlotId(id)
    setTab('plot')
  }

  function duplicatePlot(plot: AnalysisPlotState) {
    const id = `plot-${crypto.randomUUID()}`
    updateWorkspace((current) => {
      const index = current.plots.findIndex((candidate) => candidate.id === plot.id)
      const plots = [...current.plots]
      plots.splice(index + 1, 0, { ...plot, id })
      return { ...current, plots }
    })
    setSelectedPlotId(id)
    setTab('plot')
  }

  function movePlot(plotId: string, direction: -1 | 1) {
    updateWorkspace((current) => {
      const index = current.plots.findIndex((plot) => plot.id === plotId)
      const target = index + direction
      if (index < 0 || target < 0 || target >= current.plots.length) return current
      const plots = [...current.plots]
      const [moved] = plots.splice(index, 1)
      plots.splice(target, 0, moved)
      return { ...current, plots }
    })
  }

  function changeSource(path: string) {
    const source = sources.find((candidate) => candidate.path === path)
    const file = source?.files[0]
    const [x, y] = preferredAxes(file?.channels ?? [])
    updateWorkspace((current) => ({
      ...current,
      source_path: path,
      selected_file: file?.path ?? '',
      active_population_id: 'root',
      root_event_id: null,
      root_population_id: null,
      root_source_file: null,
      plots: current.plots.map((plot) => ({ ...plot, population_id: 'root', x, y, color_marker: x })),
    }))
    setSelectedMarkers(preferredMarkers(file))
    setStainingMarker(preferredMarkers(file)[0] ?? file?.channels[0] ?? '')
    setMethodResult(null)
    setAnalysisFiles(file?.path ? [file.path] : [])
  }

  function changeFile(path: string) {
    const file = selectedSource?.files.find((candidate) => candidate.path === path)
    updateWorkspace((current) => ({
      ...current,
      selected_file: path,
      root_event_id: null,
      root_population_id: null,
      root_source_file: null,
    }))
    setSelectedMarkers(preferredMarkers(file))
    setStainingMarker(preferredMarkers(file)[0] ?? file?.channels[0] ?? '')
    setMethodResult(null)
    setAnalysisFiles(path ? [path] : [])
    setActionMessage('')
  }

  function updatePopulation(id: string, patch: Partial<AnalysisPopulation>) {
    updateWorkspace((current) => ({
      ...current,
      populations: current.populations.map((population) => population.id === id ? { ...population, ...patch } : population),
    }))
  }

  function updateActiveGateGeometry(patch: NonNullable<AnalysisPopulation['geometry']>) {
    if (!activePopulation || activePopulation.id === 'root') return
    updatePopulation(activePopulation.id, {
      geometry: { ...(activePopulation.geometry ?? {}), ...patch },
    })
  }

  function createGate() {
    if (!workspace || !gateDraft || !selectedFile) return
    const id = `population-${crypto.randomUUID()}`
    const node: AnalysisPopulation = {
      id,
      name: gateName.trim() || 'Population',
      parent_id: workspace.active_population_id,
      type: gateDraft.type,
      role: gateRole === 'positive' || gateRole === 'negative' || gateRole === 'root' || gateRole === 'terminal' ? gateRole : null,
      source_file: fileSpecificGate ? selectedFile.path : null,
      x: gateDraft.x,
      y: gateDraft.y,
      geometry: gateDraft.geometry,
    }
    updateWorkspace((current) => ({ ...current, populations: [...current.populations, node] }))
    setGateDraft(null)
    setGateName('Population')
    setGateRole('')
    setFileSpecificGate(false)
    setTool('select')
  }

  function deleteActivePopulation() {
    if (!workspace || !activePopulation || activePopulation.id === 'root') return
    const removing = descendants(workspace.populations, activePopulation.id)
    if (!window.confirm(`Delete ${activePopulation.name} and ${removing.size - 1} descendant population${removing.size === 2 ? '' : 's'}?`)) return
    updateWorkspace((current) => ({
      ...current,
      active_population_id: activePopulation.parent_id || 'root',
      root_event_id: current.root_population_id && removing.has(current.root_population_id) ? null : current.root_event_id,
      root_population_id: current.root_population_id && removing.has(current.root_population_id) ? null : current.root_population_id,
      root_source_file: current.root_population_id && removing.has(current.root_population_id) ? null : current.root_source_file,
      populations: current.populations.filter((population) => !removing.has(population.id)),
      plots: current.plots.map((plot) => ({
        ...plot,
        population_id: removing.has(plot.population_id) ? (activePopulation.parent_id || 'root') : plot.population_id,
        overlay_population_id: plot.overlay_population_id && removing.has(plot.overlay_population_id) ? '' : plot.overlay_population_id,
      })),
    }))
  }

  async function calculateStainingIndex() {
    if (!selectedFile || !stainingMarker) return
    setActionMessage('Calculating staining index…')
    try {
      const response = await analysisRequest<{ success: boolean; result: Record<string, number | string> }>('/analysis/statistics', projectPath, {
        method: 'POST', body: JSON.stringify({ mode: 'staining_index', file: selectedFile.path, marker: stainingMarker }),
      })
      setStainingResult(response.result)
      setActionMessage('')
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    }
  }

  async function exportPopulation() {
    if (!workspace || !selectedFile || !selectedSource) return
    setActionMessage('Exporting selected population…')
    try {
      const files = exportAllFiles ? selectedSource.files.map((file) => file.path) : [selectedFile.path]
      const response = await analysisRequest<{ success: boolean; result: { files: Array<{ path: string; events: number }> } }>('/analysis/export', projectPath, {
        method: 'POST',
        body: JSON.stringify({
          files,
          populationId: workspace.active_population_id,
          format: exportFormat,
          maxEvents: exportCount,
          seed: workspace.seed,
          outputFolder: exportFolder,
        }),
      })
      const exported = asArray(response.result.files)
      setActionMessage(`Exported ${exported.length} file${exported.length === 1 ? '' : 's'} · ${exported.map((file) => file.path).join(', ')}`)
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    }
  }

  async function exportStatistics() {
    if (!workspace || !selectedFile || !selectedSource) return
    setActionMessage('Exporting population statistics…')
    try {
      const files = exportAllFiles ? selectedSource.files.map((file) => file.path) : [selectedFile.path]
      const response = await analysisRequest<{ success: boolean; result: { path: string; rows: number } }>('/analysis/export', projectPath, {
        method: 'POST',
        body: JSON.stringify({
          mode: 'statistics', files,
          populationIds: workspace.populations.map((population) => population.id),
          markers: channels.map((channel) => channel.channel),
          outputFolder: exportFolder,
        }),
      })
      setActionMessage(`Exported ${response.result.rows.toLocaleString()} statistic rows · ${response.result.path}`)
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    }
  }

  async function chooseExportFolder() {
    setActionMessage('Opening folder picker…')
    try {
      const response = await analysisRequest<{ success: boolean; cancelled: boolean; path?: string }>('/analysis/export-folder', projectPath, {
        method: 'POST',
        body: JSON.stringify({ current: exportFolder }),
      })
      if (response.cancelled) {
        setActionMessage('Folder selection cancelled.')
        return
      }
      if (response.path) setExportFolder(response.path)
      setActionMessage(response.path ? `Export folder: ${response.path}` : '')
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    }
  }

  function downloadWorkspace() {
    if (!workspace) return
    const envelope = {
      format: 'spectreasy-analysis-workspace',
      schema_version: 1,
      exported_at: new Date().toISOString(),
      workspace,
    }
    const blob = new Blob([`${JSON.stringify(envelope, null, 2)}\n`], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const anchor = document.createElement('a')
    anchor.href = url
    anchor.download = `spectreasy-workspace-${new Date().toISOString().slice(0, 10)}.json`
    anchor.click()
    URL.revokeObjectURL(url)
    setActionMessage('Workspace JSON exported.')
  }

  async function importWorkspaceFile(file: File) {
    try {
      const parsed = JSON.parse(await file.text()) as {
        format?: string
        workspace?: AnalysisWorkspaceState
      } | AnalysisWorkspaceState
      const candidate = 'workspace' in parsed && parsed.workspace ? parsed.workspace : parsed as AnalysisWorkspaceState
      const response = await analysisRequest<{ success: boolean; workspace: AnalysisWorkspaceState }>('/analysis/workspace', projectPath, {
        method: 'POST',
        body: JSON.stringify({ workspace: candidate }),
      })
      const imported = normalizeWorkspace(response.workspace)
      replaceWorkspace(imported)
      setSelectedPlotId(imported.plots[0]?.id ?? '')
      setActionMessage(`Imported ${imported.populations.length} populations and ${imported.plots.length} plots. Undo is available.`)
    } catch (reason) {
      setActionMessage(`Workspace import failed: ${reason instanceof Error ? reason.message : String(reason)}`)
    } finally {
      if (importWorkspaceInput.current) importWorkspaceInput.current.value = ''
    }
  }

  async function runMethod() {
    if (!workspace || !selectedFile || !methodCanRun) return
    setMethodRunning(true)
    setMethodProgress('Queueing analysis')
    setActionMessage('')
    try {
      const requestBody = {
        file: analysisFiles[0] ?? selectedFile.path,
        files: analysisFiles,
        populationId: workspace.active_population_id,
        ...(analysisMode === 'explore'
          ? {
              clusterMethod: selectedClusterMethod,
              reductionMethod: selectedReductionMethod,
            }
          : { method: selectedTrajectoryMethod }),
        ...(analysisMode === 'trajectory'
          ? {
              clusterMethod: selectedClusterMethod,
              reductionMethod: selectedReductionMethod,
            }
          : {}),
        markers: selectedMarkers,
        maxEvents: analysisMaxEvents,
        seed: workspace.seed,
        rootEventId: workspace.root_event_id,
        rootSourceFile: workspace.root_source_file,
        cofactor: analysisCofactor,
        advancedSettings: resolvedAdvancedSettings(selectedAnalysisMethods, advancedSettings),
      }
      const response = await analysisRequest<{ success: boolean; job: { job_id: string; state: string; message: string } }>('/analysis/jobs', projectPath, {
        method: 'POST',
        body: JSON.stringify(requestBody),
      })
      const jobId = response.job.job_id
      setMethodJobId(jobId)
      setMethodProgress(response.job.message || 'Queued')
      for (;;) {
        await new Promise((resolve) => window.setTimeout(resolve, 450))
        const poll = await analysisRequest<{
          success: boolean
          job: { state: string; message: string; error?: string; result?: AnalysisRunResult }
        }>(`/analysis/jobs?job_id=${encodeURIComponent(jobId)}`, projectPath)
        setMethodProgress(poll.job.message || poll.job.state)
        if (poll.job.state === 'completed' && poll.job.result) {
          setMethodResult(poll.job.result)
          setAnalysisResultTab('plots')
          break
        }
        if (poll.job.state === 'failed' || poll.job.state === 'cancelled') {
          throw new Error(poll.job.error || poll.job.message || 'Analysis stopped')
        }
      }
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    } finally {
      setMethodRunning(false)
      setMethodJobId('')
      setMethodProgress('')
    }
  }

  async function cancelMethod() {
    if (!methodJobId) return
    setMethodProgress('Cancelling…')
    try {
      await analysisRequest(`/analysis/jobs?job_id=${encodeURIComponent(methodJobId)}`, projectPath, { method: 'DELETE' })
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    }
  }

  if (!workspace) {
    return <div className={`analysis-shell theme-${cockpitTheme}`}><div className="analysis-boot"><FlaskConical size={28} /><strong>{status}</strong></div></div>
  }

  return (
    <div className={`analysis-shell theme-${cockpitTheme}`}>
      <header className="analysis-topbar">
        <div className="analysis-brand"><FlaskConical size={18} /><strong>Population analysis</strong><span>v2 workspace</span></div>
        <label><span>Source</span><select value={selectedSource?.path ?? ''} onChange={(event) => changeSource(event.target.value)}>{sources.map((source) => <option key={source.id} value={source.path}>{source.label} · {source.files.length}</option>)}</select></label>
        <label className="analysis-file-select"><span>File</span><select value={selectedFile?.path ?? ''} onChange={(event) => changeFile(event.target.value)}>{selectedSource?.files.map((file) => <option key={file.id} value={file.path}>{file.name} · {(file.event_count ?? 0).toLocaleString()}</option>)}</select></label>
        <div className={`analysis-save-state is-${saveState}`}><Save size={13} />{saveState === 'saved' ? 'Saved' : saveState === 'saving' ? 'Saving' : 'Save failed'}</div>
        <button type="button" className="analysis-close" onClick={onRequestExit} aria-label="Close population analysis"><X size={17} /></button>
      </header>

      <div className="analysis-body">
        <aside className="analysis-populations">
          <div className="analysis-panel-title"><Boxes size={15} /><strong>Populations</strong></div>
          <label className="analysis-population-search">
            <Search size={13} />
            <input aria-label="Filter populations" value={populationQuery} onChange={(event) => setPopulationQuery(event.target.value)} placeholder="Filter populations" />
            {populationQuery ? <button type="button" aria-label="Clear population filter" onClick={() => setPopulationQuery('')}><X size={12} /></button> : null}
          </label>
          <button type="button" className={`analysis-population-row root-row ${workspace.active_population_id === 'root' ? 'is-active' : ''}`} onClick={() => selectPopulation('root')}>
            <span className="population-swatch role-none" /><span>All events</span>
          </button>
          <PopulationTree
            populations={workspace.populations}
            parentId="root"
            activeId={workspace.active_population_id}
            onSelect={selectPopulation}
            collapsed={populationQuery ? new Set() : collapsedPopulationIds}
            onToggle={(id) => setCollapsedPopulationIds((current) => {
              const next = new Set(current)
              if (next.has(id)) next.delete(id)
              else next.add(id)
              return next
            })}
            visibleIds={visiblePopulationIds}
          />
          <div className="analysis-population-footer">
            <small>{selectedSource?.role === 'unmixed' ? 'Unmixed sample data' : selectedSource?.role === 'controls' ? 'SCC/control data' : 'Raw sample data'}</small>
            <span>{selectedFile?.parameter_count ?? 0} parameters</span>
          </div>
        </aside>

        <main className="analysis-main">
          <div className="analysis-toolstrip" role="toolbar" aria-label="Gating tools">
            <button type="button" disabled={undoStack.current.length === 0} onClick={undoWorkspace} title="Undo workspace change (Ctrl/Cmd+Z)"><Undo2 size={15} /> Undo</button>
            <button type="button" disabled={redoStack.current.length === 0} onClick={redoWorkspace} title="Redo workspace change (Ctrl/Cmd+Shift+Z)"><Redo2 size={15} /> Redo</button>
            <span className="analysis-tool-divider" />
            <button type="button" className={tool === 'select' ? 'is-active' : ''} onClick={() => setTool('select')} title="Select"><MousePointer2 size={15} /> Select</button>
            <button type="button" className={tool === 'rectangle' ? 'is-active' : ''} onClick={() => setTool('rectangle')} title="Rectangle gate"><SquareDashed size={15} /> Rectangle</button>
            <button type="button" className={tool === 'ellipse' ? 'is-active' : ''} onClick={() => setTool('ellipse')} title="Ellipse gate"><CircleDashed size={15} /> Ellipse</button>
            <button type="button" className={tool === 'polygon' ? 'is-active' : ''} onClick={() => setTool('polygon')} title="Polygon gate"><Pentagon size={15} /> Polygon</button>
            <button type="button" className={tool === 'range' ? 'is-active' : ''} disabled={selectedPlot?.type !== 'histogram'} onClick={() => setTool('range')} title="Histogram range gate"><SquareDashed size={15} /> Range</button>
            <button type="button" className={tool === 'root' ? 'is-active' : ''} onClick={() => setTool('root')} title="Select trajectory root"><Crosshair size={15} /> Set root</button>
            <span className="analysis-tool-divider" />
            <button type="button" onClick={addPlot} title="Add a fixed-size plot"><Plus size={15} /> Add plot</button>
            <button type="button" onClick={() => setAnalysisDialogOpen(true)} title={`Analyze ${activePopulation?.name}`}><FlaskConical size={15} /> Analyze population</button>
            <button type="button" onClick={() => setAnalysisGuideOpen(true)} title="Open population analysis guide"><CircleHelp size={15} /> Guide</button>
            <button type="button" className="analysis-compact-inspector-toggle" onClick={() => setInspectorOpen((open) => !open)} aria-label="Toggle plot settings"><SlidersHorizontal size={15} /> Settings</button>
            <span className="analysis-tool-context">Parent: <strong>{activePopulation?.name}</strong></span>
          </div>

          <section className="analysis-plot-grid">
            {workspace.plots.map((plot, plotIndex) => (
              <PlotCard
                key={plot.id}
                projectPath={projectPath}
                file={selectedFile?.path ?? ''}
                plot={plot}
                channels={channels}
                populations={workspace.populations}
                seed={workspace.seed}
                tool={tool}
                rootEventId={workspace.root_event_id}
                selected={selectedPlot?.id === plot.id}
                onSelect={() => { setSelectedPlotId(plot.id); setTab('plot') }}
                onChange={(patch) => updateWorkspace((current) => ({
                  ...current,
                  plots: current.plots.map((candidate) => candidate.id === plot.id ? { ...candidate, ...patch } : candidate),
                }))}
                onRemove={() => {
                  const remaining = workspace.plots.filter((candidate) => candidate.id !== plot.id)
                  updateWorkspace((current) => ({ ...current, plots: remaining }))
                  if (selectedPlot?.id === plot.id) setSelectedPlotId(remaining[0]?.id ?? '')
                }}
                onDuplicate={() => duplicatePlot(plot)}
                onMove={(direction) => movePlot(plot.id, direction)}
                canMoveLeft={plotIndex > 0}
                canMoveRight={plotIndex < workspace.plots.length - 1}
                onGateDraft={(draft) => { setGateDraft(draft); setGateName(`Population ${workspace.populations.length}`) }}
                onRootEvent={(eventId) => {
                  updateWorkspace((current) => ({
                    ...current,
                    root_event_id: eventId,
                    root_population_id: plot.population_id,
                    root_source_file: selectedFile?.path ?? null,
                  }))
                  setTool('select')
                }}
                onSelectGate={(populationId) => {
                  selectPopulation(populationId, false)
                  setTab('population')
                }}
                onUpdateGate={(populationId, geometry) => updatePopulation(populationId, { geometry })}
              />
            ))}
            {workspace.plots.length === 0 ? (
              <div className="analysis-empty-plots">
                <strong>No plots</strong>
                <span>Add a plot when you are ready to inspect or gate this population.</span>
                <button type="button" className="analysis-primary" onClick={addPlot}><Plus size={14} /> Add plot</button>
              </div>
            ) : null}
          </section>
        </main>

        <aside className={`analysis-inspector ${inspectorOpen ? 'is-open' : ''}`}>
          <nav>
            <button type="button" className={tab === 'plot' ? 'is-active' : ''} onClick={() => setTab('plot')}>Plot settings</button>
            <button type="button" className={tab === 'population' ? 'is-active' : ''} onClick={() => setTab('population')}>Population</button>
            <button type="button" className={tab === 'export' ? 'is-active' : ''} onClick={() => setTab('export')}>Export</button>
            <button type="button" className="analysis-inspector-close" onClick={() => setInspectorOpen(false)} aria-label="Close plot settings"><X size={14} /></button>
          </nav>

          {tab === 'plot' && selectedPlot ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section analysis-plot-settings-title">
                <div><SlidersHorizontal size={15} /><h3>Plot settings</h3></div>
                <p>Settings apply to the highlighted plot. Plot cards keep a fixed square size as more plots are added.</p>
              </section>
              <section className="analysis-inspector-section analysis-setting-grid">
                <label>Population
                  <select value={selectedPlot.population_id} onChange={(event) => {
                    updateSelectedPlot({ population_id: event.target.value })
                    updateWorkspace((current) => ({ ...current, active_population_id: event.target.value }))
                  }}>
                    {workspace.populations.map((population) => <option key={population.id} value={population.id}>{population.name}</option>)}
                  </select>
                </label>
                <label>Plot type
                  <select aria-label="Plot type" value={selectedPlot.type} onChange={(event) => updateSelectedPlot({ type: event.target.value as AnalysisPlotState['type'] })}>
                    <option value="scatter">Scatter</option>
                    <option value="contour">Contour</option>
                    <option value="histogram">Histogram</option>
                    <option value="hexbin">Hexbin</option>
                  </select>
                </label>
              </section>
              {selectedPlot.type === 'scatter' ? <section className="analysis-inspector-section">
                <h3>Display</h3>
                <label>Point color<select aria-label="Plot color" value={selectedPlot.color_by} onChange={(event) => updateSelectedPlot({ color_by: event.target.value as AnalysisPlotState['color_by'] })}><option value="density">Local density</option><option value="marker">Continuous marker value</option></select></label>
                {selectedPlot.color_by === 'marker' ? <label>Color marker<select aria-label="Color marker" value={selectedPlot.color_marker || selectedPlot.x} onChange={(event) => updateSelectedPlot({ color_marker: event.target.value })}>{channels.map((channel) => <option key={channel.channel} value={channel.channel}>{channelLabel(channel)}</option>)}</select></label> : null}
                {selectedPlot.color_by === 'marker' ? <label>Color gradient<select aria-label="Marker color palette" value={selectedPlot.color_palette || 'viridis'} onChange={(event) => updateSelectedPlot({ color_palette: event.target.value })}>{CONTINUOUS_PALETTES.map((palette) => <option key={palette.id} value={palette.id}>{palette.name}</option>)}</select><span className="analysis-palette-swatch analysis-gating-palette" style={{ backgroundImage: paletteGradient((selectedPlot.color_palette || 'viridis') as Parameters<typeof paletteGradient>[0]) }} /></label> : null}
              </section> : null}
              <section className="analysis-inspector-section">
                <h3>Transforms</h3>
                <div className="analysis-transform-grid">
                  <label>X transform<select aria-label="X transform" value={selectedPlot.x_transform} onChange={(event) => updateSelectedPlot({ x_transform: event.target.value as AnalysisPlotState['x_transform'] })}><option value="linear">Linear</option><option value="asinh">Asinh</option><option value="biexponential">Biexponential</option></select></label>
                  {selectedPlot.type !== 'histogram' ? <label>Y transform<select aria-label="Y transform" value={selectedPlot.y_transform} onChange={(event) => updateSelectedPlot({ y_transform: event.target.value as AnalysisPlotState['y_transform'] })}><option value="linear">Linear</option><option value="asinh">Asinh</option><option value="biexponential">Biexponential</option></select></label> : null}
                </div>
              </section>
              {selectedPlot.type === 'scatter' ? <section className="analysis-inspector-section">
                <h3>Points</h3>
                <label>Point size<input aria-label="Point size" type="range" min="1" max="6" step=".2" value={selectedPlot.point_size ?? 2.4} onChange={(event) => updateSelectedPlot({ point_size: Number(event.target.value) })} /></label>
                <label>Opacity<input aria-label="Point opacity" type="range" min=".1" max="1" step=".05" value={selectedPlot.opacity ?? .82} onChange={(event) => updateSelectedPlot({ opacity: Number(event.target.value) })} /></label>
              </section> : null}
              <section className="analysis-inspector-section">
                <h3>Axis limits</h3>
                <p>Leave blank for robust automatic limits.</p>
                <div className="analysis-axis-limit-grid">
                  <label>X min<input type="number" value={selectedPlot.x_min ?? ''} onChange={(event) => updateSelectedPlot({ x_min: event.target.value === '' ? null : Number(event.target.value) })} /></label>
                  <label>X max<input type="number" value={selectedPlot.x_max ?? ''} onChange={(event) => updateSelectedPlot({ x_max: event.target.value === '' ? null : Number(event.target.value) })} /></label>
                  {selectedPlot.type !== 'histogram' ? <>
                    <label>Y min<input type="number" value={selectedPlot.y_min ?? ''} onChange={(event) => updateSelectedPlot({ y_min: event.target.value === '' ? null : Number(event.target.value) })} /></label>
                    <label>Y max<input type="number" value={selectedPlot.y_max ?? ''} onChange={(event) => updateSelectedPlot({ y_max: event.target.value === '' ? null : Number(event.target.value) })} /></label>
                  </> : null}
                </div>
              </section>
              <section className="analysis-inspector-section">
                <h3>Backgating</h3>
                <p>Overlay another population on this plot without changing its active filter.</p>
                <label>Backgate population<select aria-label="Backgate population" value={selectedPlot.overlay_population_id ?? ''} onChange={(event) => updateSelectedPlot({ overlay_population_id: event.target.value })}><option value="">No backgate</option>{workspace.populations.filter((population) => population.id !== selectedPlot.population_id).map((population) => <option key={population.id} value={population.id}>{population.name}</option>)}</select></label>
              </section>
              <section className="analysis-inspector-section">
                <button type="button" className="analysis-primary" onClick={addPlot}><Plus size={14} /> Add another plot</button>
              </section>
            </div>
          ) : null}

          {tab === 'population' && activePopulation ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section">
                <label>Population name<input value={activePopulation.name} disabled={activePopulation.id === 'root'} onChange={(event) => updatePopulation(activePopulation.id, { name: event.target.value })} /></label>
                <label>Semantic role<select value={activePopulation.role ?? ''} disabled={activePopulation.id === 'root'} onChange={(event) => updatePopulation(activePopulation.id, { role: (event.target.value || null) as AnalysisPopulation['role'] })}><option value="">None</option><option value="positive">POS</option><option value="negative">NEG</option><option value="root">Trajectory root population</option><option value="terminal">Terminal population</option></select></label>
                {activePopulation.id !== 'root' ? <button type="button" className="analysis-danger" onClick={deleteActivePopulation}><Trash2 size={14} /> Delete population</button> : null}
              </section>
              {activePopulation.id !== 'root' && activePopulation.geometry ? (
                <section className="analysis-inspector-section">
                  <h3>Gate coordinates</h3>
                  <p>Edit exact raw values. The gate updates immediately in every matching plot.</p>
                  <div className="analysis-gate-coordinate-grid">
                    {activePopulation.type === 'range' ? <>
                      <label>Minimum<input type="number" value={activePopulation.geometry.min ?? ''} onChange={(event) => updateActiveGateGeometry({ min: Number(event.target.value) })} /></label>
                      <label>Maximum<input type="number" value={activePopulation.geometry.max ?? ''} onChange={(event) => updateActiveGateGeometry({ max: Number(event.target.value) })} /></label>
                    </> : null}
                    {activePopulation.type === 'rectangle' ? <>
                      <label>X minimum<input type="number" value={activePopulation.geometry.x_min ?? ''} onChange={(event) => updateActiveGateGeometry({ x_min: Number(event.target.value) })} /></label>
                      <label>X maximum<input type="number" value={activePopulation.geometry.x_max ?? ''} onChange={(event) => updateActiveGateGeometry({ x_max: Number(event.target.value) })} /></label>
                      <label>Y minimum<input type="number" value={activePopulation.geometry.y_min ?? ''} onChange={(event) => updateActiveGateGeometry({ y_min: Number(event.target.value) })} /></label>
                      <label>Y maximum<input type="number" value={activePopulation.geometry.y_max ?? ''} onChange={(event) => updateActiveGateGeometry({ y_max: Number(event.target.value) })} /></label>
                    </> : null}
                    {activePopulation.type === 'ellipse' ? <>
                      <label>Center X<input type="number" value={activePopulation.geometry.center_x ?? ''} onChange={(event) => updateActiveGateGeometry({ center_x: Number(event.target.value) })} /></label>
                      <label>Center Y<input type="number" value={activePopulation.geometry.center_y ?? ''} onChange={(event) => updateActiveGateGeometry({ center_y: Number(event.target.value) })} /></label>
                      <label>Radius X<input type="number" min="0" value={activePopulation.geometry.radius_x ?? ''} onChange={(event) => updateActiveGateGeometry({ radius_x: Math.max(0, Number(event.target.value)) })} /></label>
                      <label>Radius Y<input type="number" min="0" value={activePopulation.geometry.radius_y ?? ''} onChange={(event) => updateActiveGateGeometry({ radius_y: Math.max(0, Number(event.target.value)) })} /></label>
                    </> : null}
                  </div>
                  {activePopulation.type === 'polygon' ? (
                    <div className="analysis-polygon-coordinate-list">
                      {(activePopulation.geometry.points ?? []).map((point, index) => (
                        <div key={index}>
                          <span>{index + 1}</span>
                          <label>X<input aria-label={`Polygon point ${index + 1} X`} type="number" value={point.x} onChange={(event) => {
                            const points = [...(activePopulation.geometry?.points ?? [])]
                            points[index] = { ...point, x: Number(event.target.value) }
                            updateActiveGateGeometry({ points })
                          }} /></label>
                          <label>Y<input aria-label={`Polygon point ${index + 1} Y`} type="number" value={point.y} onChange={(event) => {
                            const points = [...(activePopulation.geometry?.points ?? [])]
                            points[index] = { ...point, y: Number(event.target.value) }
                            updateActiveGateGeometry({ points })
                          }} /></label>
                        </div>
                      ))}
                    </div>
                  ) : null}
                </section>
              ) : null}
              <section className="analysis-inspector-section">
                <h3>Statistics</h3>
                {statistics ? <div className="analysis-stat-grid"><div><span>Events</span><strong>{statistics.count.toLocaleString()}</strong></div><div><span>% parent</span><strong>{statistics.percent_parent.toFixed(2)}</strong></div><div><span>% total</span><strong>{statistics.percent_total.toFixed(2)}</strong></div></div> : <p>{statisticsError || 'Calculating…'}</p>}
                {statistics?.markers.map((marker) => <div className="analysis-marker-stat" key={marker.marker}><strong>{marker.marker}</strong><span>median {Number(marker.median).toFixed(2)}</span><span>robust SD {Number(marker.robust_sd).toFixed(2)}</span></div>)}
              </section>
              <section className="analysis-inspector-section">
                <h3>Staining index</h3>
                <p>Uses the POS and NEG gates with a shared parent: (median POS − median NEG) / (2 × 1.4826 MAD NEG).</p>
                <select value={stainingMarker} onChange={(event) => setStainingMarker(event.target.value)}>{channels.map((channel) => <option key={channel.channel} value={channel.channel}>{channel.marker ? `${channel.marker} · ${channel.channel}` : channel.channel}</option>)}</select>
                <button type="button" className="analysis-primary" onClick={() => void calculateStainingIndex()}>Calculate</button>
                {stainingResult ? <div className="analysis-staining-result"><span>Staining index</span><strong>{Number(stainingResult.staining_index).toFixed(3)}</strong><small>POS n={Number(stainingResult.positive_count).toLocaleString()} · NEG n={Number(stainingResult.negative_count).toLocaleString()}</small></div> : null}
              </section>
              <section className="analysis-inspector-section"><h3>Trajectory root</h3><p>{workspace.root_event_id ? `Event ${workspace.root_event_id.toLocaleString()} is selected.` : 'Choose Set root, then click an event in a plot.'}</p></section>
              <section className="analysis-inspector-section analysis-population-analysis-callout">
                <h3>Population analysis</h3>
                <p>Dimensional reduction, clustering and trajectory methods run in a separate workspace for {activePopulation.name}.</p>
                <button type="button" className="analysis-primary" onClick={() => setAnalysisDialogOpen(true)}><FlaskConical size={14} /> Analyze {activePopulation.name}</button>
              </section>
            </div>
          ) : null}

          {tab === 'export' ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section">
                <h3>Export {activePopulation?.name}</h3>
                <label>Format<select value={exportFormat} onChange={(event) => setExportFormat(event.target.value)}><option value="fcs">FCS</option><option value="csv">CSV</option><option value="both">FCS + CSV</option></select></label>
                <label>Maximum events<input type="number" min="0" value={exportCount} onChange={(event) => setExportCount(Math.max(0, Number(event.target.value)))} /><small>0 exports every event; otherwise seeded random sampling without replacement.</small></label>
                <label>Master seed<input type="number" min="1" value={workspace.seed} onChange={(event) => updateWorkspace((current) => ({ ...current, seed: Math.max(1, Number(event.target.value)) }))} /></label>
                <label>Output folder<span className="analysis-folder-field"><input value={exportFolder} onChange={(event) => setExportFolder(event.target.value)} /><button type="button" onClick={() => void chooseExportFolder()} aria-label="Choose export folder"><FolderOpen size={14} /></button></span><small>The folder must be inside the active project.</small></label>
                <label className="analysis-check"><input type="checkbox" checked={exportAllFiles} onChange={(event) => setExportAllFiles(event.target.checked)} /><span>Export this population from every file in the selected source</span></label>
                <button type="button" className="analysis-primary" onClick={() => void exportPopulation()}><Download size={14} /> Export population</button>
                <button type="button" className="analysis-secondary" onClick={() => void exportStatistics()}><Download size={14} /> Export hierarchy statistics CSV</button>
              </section>
              <section className="analysis-inspector-section">
                <h3>Workspace</h3>
                <p>Move the complete hierarchy, plots, gates, roles, transforms, annotations and trajectory root between projects.</p>
                <button type="button" className="analysis-secondary" onClick={downloadWorkspace}><FileDown size={14} /> Export workspace JSON</button>
                <button type="button" className="analysis-secondary" onClick={() => importWorkspaceInput.current?.click()}><FileUp size={14} /> Import workspace JSON</button>
                <input
                  ref={importWorkspaceInput}
                  className="analysis-hidden-input"
                  type="file"
                  accept=".json,application/json"
                  aria-label="Import analysis workspace"
                  onChange={(event) => {
                    const file = event.target.files?.[0]
                    if (file) void importWorkspaceFile(file)
                  }}
                />
              </section>
              <section className="analysis-inspector-section analysis-provenance"><h3>Recorded provenance</h3><p>Source file, source event row, population path, seed, event count, software writer and SHA-256 checksum.</p></section>
            </div>
          ) : null}

          {actionMessage ? <div className="analysis-action-message" role="status">{actionMessage}</div> : null}
        </aside>
      </div>

      {gateDraft ? (
        <div className="analysis-modal-backdrop" role="presentation">
          <form className="analysis-gate-dialog" onSubmit={(event) => { event.preventDefault(); createGate() }}>
            <header><div><strong>Create subpopulation</strong><span>{gateDraft.type} gate under {activePopulation?.name}</span></div><button type="button" onClick={() => setGateDraft(null)}><X size={16} /></button></header>
            <label>Name<input autoFocus value={gateName} onChange={(event) => setGateName(event.target.value)} /></label>
            <label>Role<select value={gateRole} onChange={(event) => setGateRole(event.target.value)}><option value="">None</option><option value="positive">POS for staining index</option><option value="negative">NEG for staining index</option><option value="root">Trajectory root population</option><option value="terminal">Terminal population</option></select></label>
            <label className="analysis-check"><input type="checkbox" checked={fileSpecificGate} onChange={(event) => setFileSpecificGate(event.target.checked)} /><span>Apply only to {selectedFile?.name}</span></label>
            <footer><button type="button" onClick={() => setGateDraft(null)}>Cancel</button><button type="submit" className="analysis-primary">Create population</button></footer>
          </form>
        </div>
      ) : null}

      {analysisDialogOpen && activePopulation ? (
        <div className="analysis-modal-backdrop analysis-method-backdrop" role="presentation">
          <section className="analysis-method-dialog" role="dialog" aria-modal="true" aria-label={`Analyze ${activePopulation.name}`}>
            <header className="analysis-method-dialog-header">
              <div>
                <FlaskConical size={18} />
                <span>
                  <strong>Population analysis</strong>
                  <small>{activePopulation.name} · {analysisFiles.length || 1} file{(analysisFiles.length || 1) === 1 ? '' : 's'}</small>
                </span>
              </div>
              <nav className="analysis-result-tabs" aria-label="Population analysis results">
                <button
                  type="button"
                  className={analysisResultTab === 'plots' ? 'is-active' : ''}
                  disabled={!methodResult}
                  onClick={() => setAnalysisResultTab('plots')}
                >
                  Plots
                </button>
                <button
                  type="button"
                  className={analysisResultTab === 'identities' ? 'is-active' : ''}
                  disabled={!methodResult}
                  onClick={() => setAnalysisResultTab('identities')}
                  title={methodResult ? 'Configure marker-score cell identities' : 'Run a map or trajectory first'}
                >
                  Cell identities
                </button>
              </nav>
              <div className="analysis-method-dialog-actions">
                <button type="button" className="analysis-secondary" onClick={() => setAnalysisGuideOpen(true)}><CircleHelp size={14} /> Guide</button>
                <button type="button" className="analysis-icon-button" aria-label="Close population analysis" onClick={() => setAnalysisDialogOpen(false)}><X size={16} /></button>
              </div>
            </header>

            <div className="analysis-method-dialog-body">
              <aside className="analysis-method-markers">
                <div className="analysis-method-column-title">
                  <strong>Files</strong>
                  <small>{analysisFiles.length} pooled</small>
                </div>
                <div className="analysis-analysis-file-picker">
                  {selectedSource?.files.map((file) => (
                    <label key={file.id}>
                      <input
                        type="checkbox"
                        checked={analysisFiles.includes(file.path)}
                        onChange={(event) => setAnalysisFiles((current) => event.target.checked
                          ? [...new Set([...current, file.path])]
                          : current.filter((path) => path !== file.path))}
                      />
                      <span>{file.name}</span>
                      <small>{(file.event_count ?? 0).toLocaleString()}</small>
                    </label>
                  ))}
                </div>
                <div className="analysis-method-column-title">
                  <strong>Markers</strong>
                  <small>{selectedMarkers.length} selected</small>
                </div>
                <div className="analysis-marker-picker">
                  {channels.map((channel) => (
                    <label key={channel.channel}>
                      <input
                        type="checkbox"
                        checked={selectedMarkers.includes(channel.channel)}
                        onChange={(event) => setSelectedMarkers((current) => event.target.checked
                          ? [...new Set([...current, channel.channel])]
                          : current.filter((marker) => marker !== channel.channel))}
                      />
                      <span>{channel.marker || channel.channel}</span>
                      <small>{channel.marker ? channel.channel : ''}</small>
                    </label>
                  ))}
                </div>
              </aside>

              <main className="analysis-method-catalog">
                <nav className="analysis-mode-switch" aria-label="Analysis workflow">
                  <button type="button" className={analysisMode === 'explore' ? 'is-active' : ''} onClick={() => setAnalysisMode('explore')}>Cluster + map</button>
                  <button type="button" className={analysisMode === 'trajectory' ? 'is-active' : ''} onClick={() => setAnalysisMode('trajectory')}>Trajectory</button>
                </nav>

                {analysisMode === 'explore' ? (
                  <div className="analysis-pipeline-flow">
                    <section className="analysis-pipeline-stage">
                      <header>
                        <span>1</span>
                        <div><strong>Cluster</strong><small>Saved and reused across maps</small></div>
                        <label className="analysis-skip-clustering">
                          <input
                            type="checkbox"
                            checked={selectedClusterMethod === 'none'}
                            onChange={(event) => setSelectedClusterMethod(event.target.checked
                              ? 'none'
                              : methods.find((method) => method.family === 'clustering' && method.available)?.id ?? '')}
                          />
                          Skip
                        </label>
                      </header>
                      <label className="analysis-method-select">
                        Algorithm
                        <select
                          aria-label="Clustering method"
                          value={selectedClusterMethod === 'none' ? '' : selectedClusterMethod}
                          disabled={selectedClusterMethod === 'none'}
                          onChange={(event) => setSelectedClusterMethod(event.target.value)}
                        >
                          <option value="">Choose clustering…</option>
                          {methods.filter((method) => method.family === 'clustering' && method.visible !== false).map((method) => (
                            <option key={method.id} value={method.id} disabled={!method.available}>
                              {method.name}{method.available ? '' : ` — ${method.availability_label}`}
                            </option>
                          ))}
                        </select>
                      </label>
                    </section>

                    <ChevronRight className="analysis-pipeline-arrow" size={18} aria-hidden="true" />

                    <section className="analysis-pipeline-stage">
                      <header><span>2</span><div><strong>Map</strong><small>Choose independently</small></div></header>
                      <label className="analysis-method-select">
                        Coordinates
                        <select aria-label="Dimensional-reduction method" value={selectedReductionMethod} onChange={(event) => setSelectedReductionMethod(event.target.value)}>
                          <option value="">Choose map…</option>
                          {methods.filter((method) => method.family === 'reduction' && method.visible !== false).map((method) => (
                            <option key={method.id} value={method.id} disabled={!method.available}>
                              {method.name}{method.available ? '' : ` — ${method.availability_label}`}
                            </option>
                          ))}
                        </select>
                      </label>
                    </section>
                  </div>
                ) : (
                  <section className="analysis-trajectory-stage">
                    <header><strong>Trajectory method</strong><small>Prerequisites are checked before execution</small></header>
                    <div className="analysis-trajectory-selectors">
                      <label className="analysis-method-select">
                        Method
                        <select aria-label="Trajectory method" value={selectedTrajectoryMethod} onChange={(event) => setSelectedTrajectoryMethod(event.target.value)}>
                          <option value="">Choose trajectory…</option>
                          {methods.filter((method) => method.family === 'trajectory' && method.visible !== false).map((method) => (
                            <option key={method.id} value={method.id} disabled={!method.available}>
                              {method.name}{method.available ? '' : ` — ${method.availability_label}`}
                            </option>
                          ))}
                        </select>
                      </label>
                      {selectedTrajectory?.prerequisites.includes('clustering') ? (
                        <label className="analysis-method-select">
                          Cluster object
                          <select aria-label="Trajectory clustering method" value={selectedClusterMethod === 'none' ? '' : selectedClusterMethod} onChange={(event) => setSelectedClusterMethod(event.target.value)}>
                            <option value="">Choose clustering…</option>
                            {methods.filter((method) => method.family === 'clustering' && method.visible !== false).map((method) => (
                              <option key={method.id} value={method.id} disabled={!method.available}>{method.name}</option>
                            ))}
                          </select>
                        </label>
                      ) : null}
                      {selectedTrajectory?.prerequisites.includes('reduction') ? (
                        <label className="analysis-method-select">
                          Map object
                          <select aria-label="Trajectory dimensional-reduction method" value={selectedReductionMethod} onChange={(event) => setSelectedReductionMethod(event.target.value)}>
                            <option value="">Choose map…</option>
                            {methods.filter((method) => method.family === 'reduction' && method.visible !== false).map((method) => (
                              <option key={method.id} value={method.id} disabled={!method.available}>{method.name}</option>
                            ))}
                          </select>
                        </label>
                      ) : null}
                    </div>
                    {selectedTrajectory ? (
                      <p className="analysis-prerequisite-note">
                        Root event required
                        {selectedTrajectory.prerequisites.length ? ` · Automatic input: ${selectedTrajectory.prerequisites.join(' + ')}` : ''}
                      </p>
                    ) : null}
                  </section>
                )}

                <AnalysisAdvancedSettingsPanel
                  methods={selectedAnalysisMethods}
                  values={advancedSettings}
                  onChange={(methodId, parameterId, value) => setAdvancedSettings((current) => ({
                    ...current,
                    [methodId]: { ...current[methodId], [parameterId]: value },
                  }))}
                  onReset={(methodId) => setAdvancedSettings((current) => {
                    const next = { ...current }
                    delete next[methodId]
                    return next
                  })}
                />

                <footer className="analysis-pipeline-run">
                  <details className="analysis-method-settings">
                    <summary>Run settings</summary>
                    <div className="analysis-method-parameters">
                      <label>Maximum events<input type="number" min="100" max="200000" value={analysisMaxEvents} onChange={(event) => setAnalysisMaxEvents(Math.max(100, Number(event.target.value)))} /></label>
                      <label>Asinh cofactor<input type="number" min="1" value={analysisCofactor} onChange={(event) => setAnalysisCofactor(Math.max(1, Number(event.target.value)))} /></label>
                      <label>Master seed<input type="number" min="1" value={workspace.seed} onChange={(event) => updateWorkspace((current) => ({ ...current, seed: Math.max(1, Number(event.target.value)) }))} /></label>
                    </div>
                  </details>
                  <div className={`analysis-run-status ${methodCanRun ? 'is-ready' : 'is-blocked'}`} role="status">
                    {methodRunning || methodCanRun ? <Check size={14} /> : <CircleAlert size={14} />}
                    <span>{methodRunning ? methodProgress : methodCanRun ? 'Ready to run' : methodIssues[0]}</span>
                  </div>
                  {methodRunning ? (
                    <button type="button" className="analysis-danger" disabled={!methodJobId} onClick={() => void cancelMethod()}><X size={14} /> Cancel</button>
                  ) : (
                    <button type="button" className="analysis-primary" disabled={!methodCanRun} onClick={() => void runMethod()}>
                      <FlaskConical size={14} /> {analysisMode === 'explore' ? 'Run pipeline' : 'Run trajectory'}
                    </button>
                  )}
                </footer>
                {actionMessage ? <div className="analysis-action-message" role="status">{actionMessage}</div> : null}
              </main>
            </div>

            {methodResult ? (
              <div className="analysis-method-result-area">
                <div className="analysis-result-workbench">
                  <section className="analysis-result-card analysis-method-result">
                    <header>
                      <div>
                        <strong>{methodResult.metadata.display_name || methodResult.metadata.method.name}</strong>
                        <span>{methodResult.metadata.event_count.toLocaleString()} events · {methodResult.metadata.runtime_seconds.toFixed(2)} s</span>
                        {methodResult.metadata.artifacts?.clustering?.reused ? <em className="analysis-cache-badge">Cluster reused</em> : null}
                        {methodResult.metadata.artifacts?.embedding?.reused ? <em className="analysis-cache-badge">Map reused</em> : null}
                        {methodResult.metadata.artifacts?.trajectory?.reused ? <em className="analysis-cache-badge">Trajectory reused</em> : null}
                      </div>
                      <button type="button" onClick={() => { setMethodResult(null); setAnalysisResultTab('plots') }}><X size={13} /> Clear result</button>
                    </header>
                    {analysisResultTab === 'plots' ? (
                      <AnalysisResultPlots key={`${methodResult.metadata.analysis_id}-${methodResult.metadata.identity_annotation?.identity_id ?? 'unannotated'}`} result={methodResult} />
                    ) : (
                      <CellIdentityPanel
                        key={`${methodResult.metadata.analysis_id}-identities`}
                        projectPath={projectPath}
                        result={methodResult}
                        onAnnotated={(annotation, annotatedEvents) => {
                          const byEvent = new Map(annotatedEvents.map((event) => [Number(event.event_id), event]))
                          setMethodResult((current) => current ? {
                            ...current,
                            metadata: { ...current.metadata, identity_annotation: annotation },
                            events: current.events.map((event) => ({ ...event, ...byEvent.get(Number(event.event_id)) })),
                          } : current)
                        }}
                      />
                    )}
                  </section>
                </div>
              </div>
            ) : null}
          </section>
        </div>
      ) : null}
      {analysisGuideOpen ? <AnalysisGuideDialog methods={methods} onClose={() => setAnalysisGuideOpen(false)} /> : null}
    </div>
  )
}
