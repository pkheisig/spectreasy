import { useEffect, useMemo, useState } from 'react'
import { ArrowLeft, ArrowRight, Copy, X } from 'lucide-react'
import { analysisRequest } from './api'
import { AnalysisPlot } from './AnalysisPlot'
import type { PlotTool } from './AnalysisPlot'
import type {
  AnalysisEventPayload,
  AnalysisPlotState,
  AnalysisPopulation,
  ChannelLabel,
  GateDraft,
} from './types'

type PlotCardProps = {
  projectPath: string
  file: string
  plot: AnalysisPlotState
  channels: ChannelLabel[]
  populations: AnalysisPopulation[]
  seed: number
  tool: PlotTool
  rootEventId: number | null
  selected: boolean
  onSelect: () => void
  onRemove: () => void
  onDuplicate: () => void
  onMove: (direction: -1 | 1) => void
  canMoveLeft: boolean
  canMoveRight: boolean
  onGateDraft: (draft: GateDraft) => void
  onRootEvent: (eventId: number) => void
  onSelectGate: (populationId: string) => void
  onUpdateGate: (populationId: string, geometry: NonNullable<AnalysisPopulation['geometry']>) => void
}

function channelLabel(channel: ChannelLabel) {
  return channel.marker && channel.marker !== channel.channel
    ? `${channel.marker} · ${channel.channel}`
    : channel.channel
}

function eventUrl(file: string, plot: AnalysisPlotState, seed: number, populationId = plot.population_id) {
  const color = plot.color_by === 'marker' ? plot.color_marker ?? '' : ''
  const params = new URLSearchParams({
    file,
    population_id: populationId,
    x: plot.x,
    y: plot.y,
    color,
    max_points: '8000',
    seed: String(seed),
  })
  return `/analysis/events?${params}`
}

export function PlotCard({
  projectPath,
  file,
  plot,
  channels,
  populations,
  seed,
  tool,
  rootEventId,
  selected,
  onSelect,
  onRemove,
  onDuplicate,
  onMove,
  canMoveLeft,
  canMoveRight,
  onGateDraft,
  onRootEvent,
  onSelectGate,
  onUpdateGate,
}: PlotCardProps) {
  const [payload, setPayload] = useState<AnalysisEventPayload | null>(null)
  const [overlay, setOverlay] = useState<AnalysisEventPayload | null>(null)
  const [error, setError] = useState('')

  useEffect(() => {
    if (!file) return
    const controller = new AbortController()
    void analysisRequest<AnalysisEventPayload>(eventUrl(file, plot, seed), projectPath, { signal: controller.signal })
      .then((result) => { setPayload(result); setError('') })
      .catch((reason: unknown) => {
        if (!controller.signal.aborted) setError(reason instanceof Error ? reason.message : String(reason))
      })
    return () => controller.abort()
  }, [file, plot, projectPath, seed])

  useEffect(() => {
    const overlayId = plot.overlay_population_id
    if (!file || !overlayId) return
    const controller = new AbortController()
    void analysisRequest<AnalysisEventPayload>(eventUrl(file, plot, seed, overlayId), projectPath, { signal: controller.signal })
      .then(setOverlay)
      .catch(() => { if (!controller.signal.aborted) setOverlay(null) })
    return () => controller.abort()
  }, [file, plot, plot.overlay_population_id, projectPath, seed])

  const children = useMemo(
    () => populations.filter((population) => population.parent_id === plot.population_id),
    [plot.population_id, populations],
  )
  const population = populations.find((candidate) => candidate.id === plot.population_id)
  const overlayPopulation = populations.find((candidate) => candidate.id === plot.overlay_population_id)

  return (
    <article
      className={`analysis-plot-card ${selected ? 'is-selected' : ''}`}
      aria-label={`${population?.name ?? 'Population'} plot`}
      onPointerDownCapture={onSelect}
    >
      <header className="analysis-plot-header">
        <div>
          <strong>{population?.name ?? 'All events'}</strong>
          <span>{payload ? `${payload.population_count.toLocaleString()} events` : 'Loading events'}</span>
        </div>
        <div className="analysis-plot-header-actions">
          <button type="button" className="analysis-icon-button" aria-label="Move plot left" disabled={!canMoveLeft} onClick={() => onMove(-1)}><ArrowLeft size={13} /></button>
          <button type="button" className="analysis-icon-button" aria-label="Move plot right" disabled={!canMoveRight} onClick={() => onMove(1)}><ArrowRight size={13} /></button>
          <button type="button" className="analysis-icon-button" aria-label="Duplicate plot" onClick={onDuplicate}><Copy size={13} /></button>
          <button type="button" className="analysis-icon-button" aria-label="Remove plot" onClick={onRemove}><X size={14} /></button>
        </div>
      </header>
      <div className="analysis-plot-meta">
        {plot.type !== 'histogram' ? <span>{channelLabel(channels.find((channel) => channel.channel === plot.y) ?? { channel: plot.y, marker: '' })}</span> : <span>Histogram</span>}
        <span className="analysis-plot-versus">{plot.type === 'histogram' ? '·' : '×'}</span>
        <span>{channelLabel(channels.find((channel) => channel.channel === plot.x) ?? { channel: plot.x, marker: '' })}</span>
        {overlayPopulation ? <em>Backgate · {overlayPopulation.name}</em> : null}
      </div>
      {error ? <div className="analysis-plot-error" role="alert">{error}</div> : (
        <AnalysisPlot
          plot={plot}
          payload={payload}
          overlayPayload={plot.overlay_population_id ? overlay : null}
          childGates={children}
          tool={tool}
          rootEventId={rootEventId}
          onGateDraft={onGateDraft}
          onRootEvent={onRootEvent}
          onSelectGate={onSelectGate}
          onUpdateGate={onUpdateGate}
        />
      )}
      <footer className="analysis-plot-footer">
        <span>{plot.type === 'histogram' ? 'Histogram' : plot.type === 'contour' ? 'Density contours' : plot.color_by === 'marker' ? `Continuous · ${plot.color_marker || plot.x}` : 'Density scatter'}</span>
        <span>{plot.type === 'histogram' ? plot.x_transform : `${plot.x_transform} / ${plot.y_transform}`}</span>
      </footer>
    </article>
  )
}
