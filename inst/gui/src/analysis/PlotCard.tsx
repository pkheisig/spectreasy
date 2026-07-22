import { useEffect, useMemo, useState } from 'react'
import { X } from 'lucide-react'
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
  onChange: (patch: Partial<AnalysisPlotState>) => void
  onRemove: () => void
  onGateDraft: (draft: GateDraft) => void
  onRootEvent: (eventId: number) => void
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
  onChange,
  onRemove,
  onGateDraft,
  onRootEvent,
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

  return (
    <article className="analysis-plot-card">
      <header className="analysis-plot-header">
        <select aria-label="X axis" value={plot.x} onChange={(event) => onChange({ x: event.target.value })}>
          {channels.map((channel) => <option key={channel.channel} value={channel.channel}>{channelLabel(channel)}</option>)}
        </select>
        <span className="analysis-plot-versus">×</span>
        <select aria-label="Y axis" value={plot.y} onChange={(event) => onChange({ y: event.target.value })}>
          {channels.map((channel) => <option key={channel.channel} value={channel.channel}>{channelLabel(channel)}</option>)}
        </select>
        <select aria-label="Plot color" value={plot.color_by} onChange={(event) => onChange({ color_by: event.target.value as AnalysisPlotState['color_by'] })}>
          <option value="density">Density</option>
          <option value="marker">Marker</option>
        </select>
        {plot.color_by === 'marker' ? (
          <select aria-label="Color marker" value={plot.color_marker || plot.x} onChange={(event) => onChange({ color_marker: event.target.value })}>
            {channels.map((channel) => <option key={channel.channel} value={channel.channel}>{channelLabel(channel)}</option>)}
          </select>
        ) : null}
        <select aria-label="Backgate population" value={plot.overlay_population_id ?? ''} onChange={(event) => onChange({ overlay_population_id: event.target.value })}>
          <option value="">No backgate</option>
          {populations.filter((population) => population.id !== plot.population_id).map((population) => (
            <option key={population.id} value={population.id}>{population.name}</option>
          ))}
        </select>
        <button type="button" className="analysis-icon-button" aria-label="Remove plot" onClick={onRemove}><X size={14} /></button>
      </header>
      <div className="analysis-plot-meta">
        <span>{payload ? `${payload.population_count.toLocaleString()} events` : 'Loading'}</span>
        <label> X <select value={plot.x_transform} onChange={(event) => onChange({ x_transform: event.target.value as AnalysisPlotState['x_transform'] })}><option value="linear">linear</option><option value="asinh">asinh</option></select></label>
        <label> Y <select value={plot.y_transform} onChange={(event) => onChange({ y_transform: event.target.value as AnalysisPlotState['y_transform'] })}><option value="linear">linear</option><option value="asinh">asinh</option></select></label>
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
        />
      )}
    </article>
  )
}
