import { useEffect, useMemo, useRef, useState, type RefObject } from 'react'
import { Box, ChevronDown, Download, Square, Trash2 } from 'lucide-react'
import { computeDensityBuckets } from '../gating/GatingData.js'
import type { AnalysisRunResult } from './types'
import { linearScale, PLOT_MARGIN, robustExtent } from './geometry'
import {
  CATEGORICAL_COLORS,
  CONTINUOUS_PALETTES,
  identityColor,
  markerValue,
  paletteColors,
  paletteGradient,
  resultColorModes,
  type ContinuousPaletteId,
  type ResultColorMode,
} from './resultColor'
import { coordinateLabel, eventValue, resultCoordinates, type CoordinateKey } from './resultPlot'
import {
  downloadDataUrl,
  downloadText,
  resultPlotHtml,
  resultPlotSvg,
} from './resultExport'

function continuousColor(values: number[], value: number, palette: string[]) {
  const domain = robustExtent(values)
  const position = (value - domain[0]) / (domain[1] - domain[0] || 1)
  const index = Math.max(0, Math.min(palette.length - 1, Math.floor(position * (palette.length - 1))))
  return palette[index]
}

function Result2D({
  result,
  xKey,
  yKey,
  colors,
  canvasRef,
}: {
  result: AnalysisRunResult
  xKey: CoordinateKey
  yKey: CoordinateKey
  colors: string[]
  canvasRef: RefObject<HTMLCanvasElement | null>
}) {
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const width = 372
    const height = 372
    const ratio = window.devicePixelRatio || 1
    canvas.width = width * ratio
    canvas.height = height * ratio
    canvas.style.width = '100%'
    canvas.style.height = '100%'
    const context = canvas.getContext('2d')
    if (!context) return
    context.setTransform(ratio, 0, 0, ratio, 0, 0)
    context.fillStyle = '#f8f7f3'
    context.fillRect(0, 0, width, height)
    const xDomain = robustExtent(result.events.map((event) => eventValue(event, xKey)))
    const yDomain = robustExtent(result.events.map((event) => eventValue(event, yKey)))
    const x = linearScale(xDomain, [PLOT_MARGIN.left, width - PLOT_MARGIN.right])
    const y = linearScale(yDomain, [height - PLOT_MARGIN.bottom, PLOT_MARGIN.top])
    context.globalAlpha = 0.82
    for (const [index, event] of result.events.entries()) {
      context.fillStyle = colors[index] ?? '#247f9e'
      context.fillRect(x.toScreen(eventValue(event, xKey)) - 1.2, y.toScreen(eventValue(event, yKey)) - 1.2, 2.4, 2.4)
    }
    context.globalAlpha = 1
    context.strokeStyle = '#c7c3ba'
    context.strokeRect(PLOT_MARGIN.left, PLOT_MARGIN.top, width - PLOT_MARGIN.left - PLOT_MARGIN.right, height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)
    context.fillStyle = '#746f68'
    context.font = '11px "Avenir Next", sans-serif'
    context.fillText(coordinateLabel(result, xKey), width / 2 - 28, height - 15)
    context.save()
    context.translate(16, height / 2 + 25)
    context.rotate(-Math.PI / 2)
    context.fillText(coordinateLabel(result, yKey), 0, 0)
    context.restore()
  }, [canvasRef, colors, result, xKey, yKey])

  return <canvas ref={canvasRef} aria-label={`${result.metadata.method.name} 2D result plot`} />
}

function Result3D({
  result,
  xKey,
  yKey,
  zKey,
  colorMode,
  paletteId,
  colorLabel,
  plotRef,
}: {
  result: AnalysisRunResult
  xKey: CoordinateKey
  yKey: CoordinateKey
  zKey: CoordinateKey
  colorMode: ResultColorMode
  paletteId: ContinuousPaletteId
  colorLabel: string
  plotRef: RefObject<HTMLDivElement | null>
}) {
  const [error, setError] = useState('')

  useEffect(() => {
    const element = plotRef.current
    if (!element) return
    let disposed = false
    let purge: (() => void) | undefined
    setError('')
    void import('plotly.js-gl3d-dist-min').then((module) => {
      if (disposed) return
      const Plotly = module.default
      const palette = paletteColors(paletteId)
      const continuous = colorMode === 'pseudotime' || colorMode.startsWith('marker:')
      let colors: string | number[] | string[] = '#247f9e'
      if (colorMode === 'density') {
        const xDomain = robustExtent(result.events.map((event) => eventValue(event, xKey)))
        const yDomain = robustExtent(result.events.map((event) => eventValue(event, yKey)))
        const buckets = computeDensityBuckets(
          result.events.map((event) => ({ x: eventValue(event, xKey), y: eventValue(event, yKey) })),
          'x', 'y', xDomain, yDomain,
        )
        const byEvent = Array<string>(result.events.length).fill(palette[0])
        buckets.forEach((bucket, colorIndex) => bucket.forEach((eventIndex) => { byEvent[eventIndex] = palette[colorIndex] ?? palette[0] }))
        colors = byEvent
      } else if (colorMode === 'cluster') {
        colors = result.events.map((event) => CATEGORICAL_COLORS[((event.cluster_id ?? 1) - 1) % CATEGORICAL_COLORS.length])
      } else if (colorMode === 'identity') {
        colors = result.events.map((event) => identityColor(result, event.predicted_identity))
      } else {
        colors = result.events.map((event) => markerValue(event, colorMode))
      }
      const marker = {
        size: 2.6,
        opacity: 0.82,
        color: colors,
        ...(continuous ? {
          colorscale: palette.map((color, index) => [index / Math.max(1, palette.length - 1), color]),
          colorbar: { title: { text: colorLabel, font: { size: 9 } }, tickfont: { size: 8 }, thickness: 8, len: 0.55, x: 0.93 },
        } : {}),
      }
      return Plotly.newPlot(element, [{
        type: 'scatter3d',
        mode: 'markers',
        x: result.events.map((event) => eventValue(event, xKey)),
        y: result.events.map((event) => eventValue(event, yKey)),
        z: result.events.map((event) => eventValue(event, zKey)),
        customdata: result.events.map((event) => [
          event.event_id,
          colorMode === 'identity' ? (event.predicted_identity ?? 'Unassigned') : '',
        ]),
        hovertemplate: colorMode === 'identity'
          ? 'Event %{customdata[0]}<br>%{customdata[1]}<br>%{x:.3f}<br>%{y:.3f}<br>%{z:.3f}<extra></extra>'
          : 'Event %{customdata[0]}<br>%{x:.3f}<br>%{y:.3f}<br>%{z:.3f}<extra></extra>',
        marker,
      }], {
        autosize: true,
        margin: { l: 8, r: 8, t: 8, b: 8 },
        paper_bgcolor: '#f8f7f3',
        plot_bgcolor: '#f8f7f3',
        scene: {
          bgcolor: '#f8f7f3',
          aspectmode: 'cube',
          xaxis: { title: { text: coordinateLabel(result, xKey), font: { size: 10 } }, tickfont: { size: 8 }, gridcolor: '#ddd9d1', zerolinecolor: '#c7c3ba' },
          yaxis: { title: { text: coordinateLabel(result, yKey), font: { size: 10 } }, tickfont: { size: 8 }, gridcolor: '#ddd9d1', zerolinecolor: '#c7c3ba' },
          zaxis: { title: { text: coordinateLabel(result, zKey), font: { size: 10 } }, tickfont: { size: 8 }, gridcolor: '#ddd9d1', zerolinecolor: '#c7c3ba' },
        },
        showlegend: false,
      }, {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ['toImage', 'sendDataToCloud'],
      }).then(() => {
        purge = () => Plotly.purge(element)
      })
    }).catch((reason: unknown) => {
      if (!disposed) setError(reason instanceof Error ? reason.message : String(reason))
    })
    return () => {
      disposed = true
      purge?.()
    }
  }, [colorLabel, colorMode, paletteId, plotRef, result, xKey, yKey, zKey])

  if (error) return <div className="analysis-result-plot-error" role="alert">3D rendering could not start: {error}</div>
  return <div ref={plotRef} className="analysis-result-plotly" aria-label={`${result.metadata.method.name} interactive 3D result plot`} />
}

function CoordinatePicker({
  result,
  coordinates,
  mode,
  setMode,
  xKey,
  setXKey,
  yKey,
  setYKey,
  zKey,
  setZKey,
}: {
  result: AnalysisRunResult
  coordinates: CoordinateKey[]
  mode: '2d' | '3d'
  setMode: (mode: '2d' | '3d') => void
  xKey: CoordinateKey
  setXKey: (key: CoordinateKey) => void
  yKey: CoordinateKey
  setYKey: (key: CoordinateKey) => void
  zKey: CoordinateKey
  setZKey: (key: CoordinateKey) => void
}) {
  const canRender3D = coordinates.length >= 3
  const summary = mode === '3d'
    ? `X ${coordinateLabel(result, xKey)} · Y ${coordinateLabel(result, yKey)} · Z ${coordinateLabel(result, zKey)}`
    : `X ${coordinateLabel(result, xKey)} · Y ${coordinateLabel(result, yKey)}`
  return (
    <details className="analysis-coordinate-picker">
      <summary aria-label="Choose result coordinates"><span>{summary}</span><ChevronDown size={12} /></summary>
      <div className="analysis-coordinate-popover">
        {canRender3D ? (
          <div className="analysis-result-mode" role="group" aria-label="Result plot dimensionality">
            <button type="button" className={mode === '2d' ? 'is-active' : ''} onClick={() => setMode('2d')}><Square size={12} /> 2D</button>
            <button type="button" className={mode === '3d' ? 'is-active' : ''} onClick={() => setMode('3d')}><Box size={12} /> 3D</button>
          </div>
        ) : null}
        <table>
          <tbody>
            <tr>
              <th scope="row">X</th>
              <td><select aria-label="Result X coordinate" value={xKey} onChange={(event) => setXKey(event.target.value as CoordinateKey)}>{coordinates.map((key) => <option key={key} value={key} disabled={key === yKey || (mode === '3d' && key === zKey)}>{coordinateLabel(result, key)}</option>)}</select></td>
            </tr>
            <tr>
              <th scope="row">Y</th>
              <td><select aria-label="Result Y coordinate" value={yKey} onChange={(event) => setYKey(event.target.value as CoordinateKey)}>{coordinates.map((key) => <option key={key} value={key} disabled={key === xKey || (mode === '3d' && key === zKey)}>{coordinateLabel(result, key)}</option>)}</select></td>
            </tr>
            {mode === '3d' ? (
              <tr>
                <th scope="row">Z</th>
                <td><select aria-label="Result Z coordinate" value={zKey} onChange={(event) => setZKey(event.target.value as CoordinateKey)}>{coordinates.map((key) => <option key={key} value={key} disabled={key === xKey || key === yKey}>{coordinateLabel(result, key)}</option>)}</select></td>
              </tr>
            ) : null}
          </tbody>
        </table>
      </div>
    </details>
  )
}

function pointColors(
  result: AnalysisRunResult,
  xKey: CoordinateKey,
  yKey: CoordinateKey,
  colorMode: ResultColorMode,
  paletteId: ContinuousPaletteId,
): string[] {
  const palette = paletteColors(paletteId)
  if (colorMode === 'density') {
    const xDomain = robustExtent(result.events.map((event) => eventValue(event, xKey)))
    const yDomain = robustExtent(result.events.map((event) => eventValue(event, yKey)))
    const buckets = computeDensityBuckets(
      result.events.map((event) => ({ x: eventValue(event, xKey), y: eventValue(event, yKey) })),
      'x', 'y', xDomain, yDomain,
    )
    const colors = Array<string>(result.events.length).fill(palette[0])
    buckets.forEach((bucket, colorIndex) => bucket.forEach((index) => {
      colors[index] = palette[colorIndex] ?? palette[0]
    }))
    return colors
  }
  if (colorMode === 'cluster') {
    return result.events.map((event) => CATEGORICAL_COLORS[((event.cluster_id ?? 1) - 1) % CATEGORICAL_COLORS.length])
  }
  if (colorMode === 'identity') {
    return result.events.map((event) => identityColor(result, event.predicted_identity))
  }
  const values = result.events.map((event) => markerValue(event, colorMode))
  return result.events.map((event) => continuousColor(values, markerValue(event, colorMode), palette))
}

export function AnalysisResultPlot({
  result,
  canRemove = false,
  onRemove,
}: {
  result: AnalysisRunResult
  canRemove?: boolean
  onRemove?: () => void
}) {
  const coordinates = useMemo(() => resultCoordinates(result), [result])
  const colorModes = useMemo(() => resultColorModes(result), [result])
  const defaultColorMode: ResultColorMode = colorModes.some((item) => item.id === 'identity')
    ? 'identity'
    : colorModes.some((item) => item.id === 'pseudotime') ? 'pseudotime'
    : colorModes.some((mode) => mode.id === 'cluster') ? 'cluster' : 'density'
  const [mode, setMode] = useState<'2d' | '3d'>('2d')
  const [xKey, setXKey] = useState<CoordinateKey>(coordinates[0] ?? 'dimension_1')
  const [yKey, setYKey] = useState<CoordinateKey>(coordinates[1] ?? 'dimension_2')
  const [zKey, setZKey] = useState<CoordinateKey>(coordinates[2] ?? 'dimension_3')
  const [colorMode, setColorMode] = useState<ResultColorMode>(defaultColorMode)
  const [paletteId, setPaletteId] = useState<ContinuousPaletteId>(defaultColorMode === 'density' ? 'control-density' : 'viridis')
  const [exportOpen, setExportOpen] = useState(false)
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const plotRef = useRef<HTMLDivElement>(null)
  const canRender3D = coordinates.length >= 3
  const effectiveColorMode = colorModes.some((candidate) => candidate.id === colorMode) ? colorMode : defaultColorMode
  const selectedColor = colorModes.find((candidate) => candidate.id === effectiveColorMode) ?? colorModes[0]
  const colors = useMemo(
    () => pointColors(result, xKey, yKey, effectiveColorMode, paletteId),
    [effectiveColorMode, paletteId, result, xKey, yKey],
  )

  function chooseColor(next: ResultColorMode) {
    setColorMode(next)
    if (next === 'density') setPaletteId('control-density')
    else if (next === 'pseudotime' || next.startsWith('marker:')) setPaletteId((current) => current === 'control-density' ? 'viridis' : current)
  }

  async function exportPlot(format: 'png' | 'svg' | 'html') {
    const base = `${result.metadata.analysis_id}-${coordinateLabel(result, xKey)}-${coordinateLabel(result, yKey)}`.replaceAll(/[^a-zA-Z0-9._-]+/g, '-')
    if (format === 'html') {
      downloadText(`${base}.html`, resultPlotHtml({ result, xKey, yKey, zKey, colors, mode }), 'text/html;charset=utf-8')
    } else if (format === 'svg' && mode === '2d') {
      downloadText(`${base}.svg`, resultPlotSvg({ result, xKey, yKey, colors }), 'image/svg+xml;charset=utf-8')
    } else if (mode === '2d') {
      const url = canvasRef.current?.toDataURL('image/png')
      if (url) downloadDataUrl(`${base}.png`, url)
    } else if (plotRef.current) {
      const module = await import('plotly.js-gl3d-dist-min')
      const url = await module.default.toImage(plotRef.current, {
        format,
        width: 900,
        height: 900,
      })
      downloadDataUrl(`${base}.${format}`, url)
    }
    setExportOpen(false)
  }

  return (
    <div className="analysis-result-visual">
      <div className="analysis-result-controls">
        <CoordinatePicker
          result={result}
          coordinates={coordinates}
          mode={mode}
          setMode={setMode}
          xKey={xKey}
          setXKey={setXKey}
          yKey={yKey}
          setYKey={setYKey}
          zKey={zKey}
          setZKey={setZKey}
        />
        <label className="analysis-result-color-control">Color<select aria-label="Result color" value={effectiveColorMode} onChange={(event) => chooseColor(event.target.value as ResultColorMode)}>{colorModes.map((item) => <option key={item.id} value={item.id}>{item.label}</option>)}</select></label>
        {selectedColor?.continuous ? (
          <label className="analysis-result-palette-control">
            Palette
            <select aria-label="Result color palette" value={paletteId} onChange={(event) => setPaletteId(event.target.value as ContinuousPaletteId)}>
              {CONTINUOUS_PALETTES.map((palette) => <option key={palette.id} value={palette.id}>{palette.name}</option>)}
            </select>
            <span className="analysis-palette-swatch" style={{ backgroundImage: paletteGradient(paletteId) }} aria-hidden="true" />
          </label>
        ) : null}
        <div className="analysis-result-export">
          <button type="button" aria-expanded={exportOpen} onClick={() => setExportOpen((current) => !current)}><Download size={12} /> Export</button>
          {exportOpen ? (
            <div role="menu">
              <button type="button" onClick={() => void exportPlot('png')}>PNG image</button>
              <button type="button" onClick={() => void exportPlot('svg')}>SVG vector</button>
              <button type="button" onClick={() => void exportPlot('html')}>Interactive HTML</button>
            </div>
          ) : null}
        </div>
        {canRemove ? <button type="button" className="analysis-result-remove" aria-label="Remove result plot" onClick={onRemove}><Trash2 size={12} /></button> : null}
      </div>
      <div className="analysis-result-plot">
        {mode === '3d' && canRender3D
          ? <Result3D result={result} xKey={xKey} yKey={yKey} zKey={zKey} colorMode={effectiveColorMode} paletteId={paletteId} colorLabel={selectedColor?.label ?? 'Value'} plotRef={plotRef} />
          : <Result2D result={result} xKey={xKey} yKey={yKey} colors={colors} canvasRef={canvasRef} />}
      </div>
    </div>
  )
}
