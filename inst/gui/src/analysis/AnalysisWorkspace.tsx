import { useCallback, useEffect, useMemo, useState } from 'react'
import {
  Boxes,
  ChevronRight,
  Crosshair,
  Download,
  FlaskConical,
  MousePointer2,
  Pentagon,
  Plus,
  Save,
  SquareDashed,
  Trash2,
  X,
} from 'lucide-react'
import { analysisRequest } from './api'
import { AnalysisResultPlot } from './AnalysisResultPlot'
import { PlotCard } from './PlotCard'
import type { PlotTool } from './AnalysisPlot'
import type {
  AnalysisMethod,
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
type InspectorTab = 'population' | 'export' | 'methods'

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
    schema_version: scalarNumber(value.schema_version, 1),
    source_path: scalarString(value.source_path),
    selected_file: scalarString(value.selected_file),
    active_population_id: scalarString(value.active_population_id, 'root'),
    seed: scalarNumber(value.seed, 20260723),
    root_event_id: Number.isInteger(rootEvent) && rootEvent > 0 ? rootEvent : null,
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
      type: scalarString(plot.type, 'scatter') as AnalysisPlotState['type'],
      population_id: scalarString(plot.population_id, 'root'),
      x: scalarString(plot.x),
      y: scalarString(plot.y),
      color_by: scalarString(plot.color_by, 'density') as AnalysisPlotState['color_by'],
      color_marker: scalarString(plot.color_marker),
      overlay_population_id: scalarString(plot.overlay_population_id),
      x_transform: scalarString(plot.x_transform, 'linear') as AnalysisPlotState['x_transform'],
      y_transform: scalarString(plot.y_transform, 'linear') as AnalysisPlotState['y_transform'],
    })),
    annotations: asArray(value.annotations),
  }
}

function normalizeMethod(method: AnalysisMethod): AnalysisMethod {
  return {
    ...method,
    id: scalarString(method.id),
    name: scalarString(method.name),
    family: scalarString(method.family) as AnalysisMethod['family'],
    package: scalarString(method.package),
    available: scalarBoolean(method.available),
    installed: scalarBoolean(method.installed),
    adapter_verified: scalarBoolean(method.adapter_verified),
    research_only: scalarBoolean(method.research_only),
    version: scalarString(method.version),
    citation: scalarString(method.citation),
    doi: scalarString(method.doi),
    requirements: asArray(method.requirements).map((value) => scalarString(value)),
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

function channelLabels(sourceFile: AnalysisSource['files'][number] | undefined): ChannelLabel[] {
  if (!sourceFile) return []
  return sourceFile.channels.map((channel, index) => ({ channel, marker: sourceFile.descriptions[index] || '' }))
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
}: {
  populations: AnalysisPopulation[]
  parentId: string
  activeId: string
  depth?: number
  onSelect: (id: string) => void
}) {
  return populationChildren(populations, parentId).map((population) => {
    const children = populationChildren(populations, population.id)
    return (
      <div key={population.id}>
        <button
          type="button"
          className={`analysis-population-row ${activeId === population.id ? 'is-active' : ''}`}
          style={{ paddingLeft: 12 + depth * 16 }}
          onClick={() => onSelect(population.id)}
          onDoubleClick={() => onSelect(population.id)}
        >
          <ChevronRight size={12} className={children.length ? '' : 'is-leaf'} />
          <span className={`population-swatch role-${population.role || 'none'}`} />
          <span>{population.name}</span>
          {population.role ? <small>{population.role === 'positive' ? 'POS' : population.role === 'negative' ? 'NEG' : population.role}</small> : null}
        </button>
        {children.length ? (
          <PopulationTree populations={populations} parentId={population.id} activeId={activeId} depth={depth + 1} onSelect={onSelect} />
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
  const [tab, setTab] = useState<InspectorTab>('population')
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
  const [selectedMethod, setSelectedMethod] = useState('pca')
  const [selectedMarkers, setSelectedMarkers] = useState<string[]>([])
  const [methodRunning, setMethodRunning] = useState(false)
  const [methodResult, setMethodResult] = useState<AnalysisRunResult | null>(null)
  const [analysisMaxEvents, setAnalysisMaxEvents] = useState(20000)
  const [analysisCofactor, setAnalysisCofactor] = useState(150)
  const [analysisNeighbors, setAnalysisNeighbors] = useState(15)
  const [analysisClusters, setAnalysisClusters] = useState(10)
  const [analysisPerplexity, setAnalysisPerplexity] = useState(30)

  const updateWorkspace = useCallback((update: (current: AnalysisWorkspaceState) => AnalysisWorkspaceState) => {
    setWorkspace((current) => current ? update(current) : current)
    setSaveState('saving')
  }, [])

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
        plots: loaded.plots.length ? loaded.plots.map((plot) => ({
          ...plot,
          x: selected?.channels.includes(plot.x) ? plot.x : x,
          y: selected?.channels.includes(plot.y) ? plot.y : y,
        })) : [{ id: 'plot-1', type: 'scatter', population_id: 'root', x, y, color_by: 'density', x_transform: 'linear', y_transform: 'linear' }],
      }
      const availableDefault = loadedMethods.find((method) => method.available && method.family === 'reduction')
      setSources(loadedSources)
      setMethods(loadedMethods)
      setSelectedMethod(availableDefault?.id ?? 'pca')
      setSelectedMarkers(preferredMarkers(selected))
      setStainingMarker(preferredMarkers(selected)[0] ?? selected?.channels[0] ?? '')
      setWorkspace(loaded)
      setReady(true)
      setStatus(loadedSources.length ? 'Ready' : 'No FCS sources found')
    }).catch((reason: unknown) => {
      if (!cancelled) setStatus(reason instanceof Error ? reason.message : String(reason))
    })
    return () => { cancelled = true }
  }, [projectPath])

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

  function selectPopulation(id: string) {
    updateWorkspace((current) => ({
      ...current,
      active_population_id: id,
      plots: current.plots.map((plot) => ({ ...plot, population_id: id })),
    }))
    setStainingResult(null)
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
      plots: current.plots.map((plot) => ({ ...plot, population_id: 'root', x, y, color_marker: x })),
    }))
    setSelectedMarkers(preferredMarkers(file))
    setStainingMarker(preferredMarkers(file)[0] ?? file?.channels[0] ?? '')
    setMethodResult(null)
  }

  function changeFile(path: string) {
    const file = selectedSource?.files.find((candidate) => candidate.path === path)
    updateWorkspace((current) => ({ ...current, selected_file: path }))
    setSelectedMarkers(preferredMarkers(file))
    setStainingMarker(preferredMarkers(file)[0] ?? file?.channels[0] ?? '')
    setMethodResult(null)
    setActionMessage('')
  }

  function updatePopulation(id: string, patch: Partial<AnalysisPopulation>) {
    updateWorkspace((current) => ({
      ...current,
      populations: current.populations.map((population) => population.id === id ? { ...population, ...patch } : population),
    }))
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

  async function runMethod() {
    if (!workspace || !selectedFile || !selectedMethod) return
    setMethodRunning(true)
    setActionMessage('')
    try {
      const response = await analysisRequest<{ success: boolean; result: AnalysisRunResult }>('/analysis/run', projectPath, {
        method: 'POST',
        body: JSON.stringify({
          file: selectedFile.path,
          populationId: workspace.active_population_id,
          method: selectedMethod,
          markers: selectedMarkers,
          maxEvents: analysisMaxEvents,
          seed: workspace.seed,
          rootEventId: workspace.root_event_id,
          cofactor: analysisCofactor,
          neighbors: analysisNeighbors,
          clusters: analysisClusters,
          perplexity: analysisPerplexity,
        }),
      })
      setMethodResult(response.result)
    } catch (reason) {
      setActionMessage(reason instanceof Error ? reason.message : String(reason))
    } finally {
      setMethodRunning(false)
    }
  }

  if (!workspace) {
    return <div className={`analysis-shell theme-${cockpitTheme}`}><div className="analysis-boot"><FlaskConical size={28} /><strong>{status}</strong></div></div>
  }

  return (
    <div className={`analysis-shell theme-${cockpitTheme}`}>
      <header className="analysis-topbar">
        <div className="analysis-brand"><FlaskConical size={18} /><strong>Sample analysis</strong><span>v2 workspace</span></div>
        <label><span>Source</span><select value={selectedSource?.path ?? ''} onChange={(event) => changeSource(event.target.value)}>{sources.map((source) => <option key={source.id} value={source.path}>{source.label} · {source.files.length}</option>)}</select></label>
        <label className="analysis-file-select"><span>File</span><select value={selectedFile?.path ?? ''} onChange={(event) => changeFile(event.target.value)}>{selectedSource?.files.map((file) => <option key={file.id} value={file.path}>{file.name} · {(file.event_count ?? 0).toLocaleString()}</option>)}</select></label>
        <div className={`analysis-save-state is-${saveState}`}><Save size={13} />{saveState === 'saved' ? 'Saved' : saveState === 'saving' ? 'Saving' : 'Save failed'}</div>
        <button type="button" className="analysis-close" onClick={onRequestExit} aria-label="Close sample analysis"><X size={17} /></button>
      </header>

      <div className="analysis-body">
        <aside className="analysis-populations">
          <div className="analysis-panel-title"><Boxes size={15} /><strong>Populations</strong></div>
          <button type="button" className={`analysis-population-row root-row ${workspace.active_population_id === 'root' ? 'is-active' : ''}`} onClick={() => selectPopulation('root')}>
            <span className="population-swatch role-none" /><span>All events</span>
          </button>
          <PopulationTree populations={workspace.populations} parentId="root" activeId={workspace.active_population_id} onSelect={selectPopulation} />
          <div className="analysis-population-footer">
            <small>{selectedSource?.role === 'unmixed' ? 'Unmixed sample data' : selectedSource?.role === 'controls' ? 'SCC/control data' : 'Raw sample data'}</small>
            <span>{selectedFile?.parameter_count ?? 0} parameters</span>
          </div>
        </aside>

        <main className="analysis-main">
          <div className="analysis-toolstrip" role="toolbar" aria-label="Gating tools">
            <button type="button" className={tool === 'select' ? 'is-active' : ''} onClick={() => setTool('select')} title="Select"><MousePointer2 size={15} /> Select</button>
            <button type="button" className={tool === 'rectangle' ? 'is-active' : ''} onClick={() => setTool('rectangle')} title="Rectangle gate"><SquareDashed size={15} /> Rectangle</button>
            <button type="button" className={tool === 'polygon' ? 'is-active' : ''} onClick={() => setTool('polygon')} title="Polygon gate"><Pentagon size={15} /> Polygon</button>
            <button type="button" className={tool === 'root' ? 'is-active' : ''} onClick={() => setTool('root')} title="Select trajectory root"><Crosshair size={15} /> Set root</button>
            <span className="analysis-tool-context">Parent: <strong>{activePopulation?.name}</strong></span>
          </div>

          {methodResult ? (
            <section className="analysis-result-card">
              <header><div><strong>{methodResult.metadata.method.name}</strong><span>{methodResult.metadata.event_count.toLocaleString()} events · {methodResult.metadata.runtime_seconds.toFixed(2)} s</span></div><button type="button" onClick={() => setMethodResult(null)}><X size={14} /> Close result</button></header>
              <AnalysisResultPlot result={methodResult} />
              {methodResult.metadata.plot_files?.length ? <footer className="analysis-result-exports"><strong>Saved plot formats</strong>{methodResult.metadata.plot_files.map((file) => <span key={file.format}>{file.format.toUpperCase()} · {file.path}</span>)}</footer> : null}
            </section>
          ) : null}

          <section className={`analysis-plot-grid plots-${Math.min(workspace.plots.length, 4)}`}>
            {workspace.plots.map((plot) => (
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
                onChange={(patch) => updateWorkspace((current) => ({ ...current, plots: current.plots.map((candidate) => candidate.id === plot.id ? { ...candidate, ...patch } : candidate) }))}
                onRemove={() => updateWorkspace((current) => ({ ...current, plots: current.plots.filter((candidate) => candidate.id !== plot.id) }))}
                onGateDraft={(draft) => { setGateDraft(draft); setGateName(`Population ${workspace.populations.length}`) }}
                onRootEvent={(eventId) => { updateWorkspace((current) => ({ ...current, root_event_id: eventId })); setTool('select') }}
              />
            ))}
            <button
              type="button"
              className="analysis-add-plot"
              onClick={() => {
                const [x, y] = preferredAxes(channels.map((channel) => channel.channel))
                const plot: AnalysisPlotState = { id: `plot-${crypto.randomUUID()}`, type: 'scatter', population_id: workspace.active_population_id, x, y, color_by: 'density', x_transform: 'linear', y_transform: 'linear' }
                updateWorkspace((current) => ({ ...current, plots: [...current.plots, plot] }))
              }}
            >
              <Plus size={22} /><strong>Add plot</strong><span>Choose any parameters after adding</span>
            </button>
          </section>
        </main>

        <aside className="analysis-inspector">
          <nav>
            <button type="button" className={tab === 'population' ? 'is-active' : ''} onClick={() => setTab('population')}>Population</button>
            <button type="button" className={tab === 'export' ? 'is-active' : ''} onClick={() => setTab('export')}>Export</button>
            <button type="button" className={tab === 'methods' ? 'is-active' : ''} onClick={() => setTab('methods')}>Discover</button>
          </nav>

          {tab === 'population' && activePopulation ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section">
                <label>Population name<input value={activePopulation.name} disabled={activePopulation.id === 'root'} onChange={(event) => updatePopulation(activePopulation.id, { name: event.target.value })} /></label>
                <label>Semantic role<select value={activePopulation.role ?? ''} disabled={activePopulation.id === 'root'} onChange={(event) => updatePopulation(activePopulation.id, { role: (event.target.value || null) as AnalysisPopulation['role'] })}><option value="">None</option><option value="positive">POS</option><option value="negative">NEG</option><option value="root">Trajectory root population</option><option value="terminal">Terminal population</option></select></label>
                {activePopulation.id !== 'root' ? <button type="button" className="analysis-danger" onClick={deleteActivePopulation}><Trash2 size={14} /> Delete population</button> : null}
              </section>
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
            </div>
          ) : null}

          {tab === 'export' ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section">
                <h3>Export {activePopulation?.name}</h3>
                <label>Format<select value={exportFormat} onChange={(event) => setExportFormat(event.target.value)}><option value="fcs">FCS</option><option value="csv">CSV</option><option value="both">FCS + CSV</option></select></label>
                <label>Maximum events<input type="number" min="0" value={exportCount} onChange={(event) => setExportCount(Math.max(0, Number(event.target.value)))} /><small>0 exports every event; otherwise seeded random sampling without replacement.</small></label>
                <label>Master seed<input type="number" min="1" value={workspace.seed} onChange={(event) => updateWorkspace((current) => ({ ...current, seed: Math.max(1, Number(event.target.value)) }))} /></label>
                <label>Output folder<input value={exportFolder} onChange={(event) => setExportFolder(event.target.value)} /></label>
                <label className="analysis-check"><input type="checkbox" checked={exportAllFiles} onChange={(event) => setExportAllFiles(event.target.checked)} /><span>Export this population from every file in the selected source</span></label>
                <button type="button" className="analysis-primary" onClick={() => void exportPopulation()}><Download size={14} /> Export population</button>
                <button type="button" className="analysis-secondary" onClick={() => void exportStatistics()}><Download size={14} /> Export hierarchy statistics CSV</button>
              </section>
              <section className="analysis-inspector-section analysis-provenance"><h3>Recorded provenance</h3><p>Source file, source event row, population path, seed, event count, software writer and SHA-256 checksum.</p></section>
            </div>
          ) : null}

          {tab === 'methods' ? (
            <div className="analysis-inspector-scroll">
              <section className="analysis-inspector-section">
                <h3>Markers</h3>
                <div className="analysis-marker-picker">{channels.map((channel) => <label key={channel.channel}><input type="checkbox" checked={selectedMarkers.includes(channel.channel)} onChange={(event) => setSelectedMarkers((current) => event.target.checked ? [...new Set([...current, channel.channel])] : current.filter((marker) => marker !== channel.channel))} /><span>{channel.marker || channel.channel}</span><small>{channel.marker ? channel.channel : ''}</small></label>)}</div>
              </section>
              {(['reduction', 'clustering', 'trajectory'] as const).map((family) => (
                <section className="analysis-inspector-section" key={family}>
                  <h3>{family === 'reduction' ? 'Dimensional reduction' : family === 'clustering' ? 'Clustering' : 'Pseudotime & trajectory'}</h3>
                  <div className="analysis-method-list">{methods.filter((method) => method.family === family).map((method) => <button type="button" key={method.id} disabled={!method.available} className={`${selectedMethod === method.id ? 'is-selected' : ''} ${method.available ? '' : 'is-unavailable'}`} onClick={() => setSelectedMethod(method.id)}><span><strong>{method.name}</strong>{method.available ? <em>ready</em> : <em>{method.research_only || !method.adapter_verified ? 'validation required' : 'not installed'}</em>}</span><small>{method.citation}{method.version ? ` · ${method.package} ${method.version}` : ''}</small>{!method.available ? <small>{methodRequirement(method)}</small> : null}</button>)}</div>
                </section>
              ))}
              <section className="analysis-inspector-section analysis-run-method">
                <div className="analysis-method-parameters">
                  <label>Max events<input type="number" min="100" max="200000" value={analysisMaxEvents} onChange={(event) => setAnalysisMaxEvents(Math.max(100, Number(event.target.value)))} /></label>
                  <label>Asinh cofactor<input type="number" min="1" value={analysisCofactor} onChange={(event) => setAnalysisCofactor(Math.max(1, Number(event.target.value)))} /></label>
                  {['umap', 'diffusion-map', 'dpt', 'phenograph'].includes(selectedMethod) ? <label>Neighbors (k)<input type="number" min="2" max="200" value={analysisNeighbors} onChange={(event) => setAnalysisNeighbors(Math.max(2, Number(event.target.value)))} /></label> : null}
                  {selectedMethod === 'flowsom' ? <label>Metaclusters<input type="number" min="2" max="100" value={analysisClusters} onChange={(event) => setAnalysisClusters(Math.max(2, Number(event.target.value)))} /></label> : null}
                  {selectedMethod === 'tsne' ? <label>Perplexity<input type="number" min="2" max="200" value={analysisPerplexity} onChange={(event) => setAnalysisPerplexity(Math.max(2, Number(event.target.value)))} /></label> : null}
                </div>
                <button type="button" className="analysis-primary" disabled={methodRunning || selectedMarkers.length < 2 || (selectedMethod === 'dpt' && !workspace.root_event_id)} onClick={() => void runMethod()}><FlaskConical size={14} /> {methodRunning ? 'Running…' : `Run ${methods.find((method) => method.id === selectedMethod)?.name ?? 'method'}`}</button>
                {selectedMethod === 'dpt' && !workspace.root_event_id ? <small>Select a trajectory root event first.</small> : null}
              </section>
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
    </div>
  )
}
