import type { AnalysisRunResult } from './types.ts'
import { coordinateLabel, eventValue, resultCoordinates, type CoordinateKey } from './resultPlot.ts'

function csvCell(value: unknown): string {
  const text = value == null ? '' : String(value)
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text
}

export function resultTableCsv(result: AnalysisRunResult): string {
  const coordinates = resultCoordinates(result)
  const markers = result.metadata.marker_columns ?? []
  const optional = [
    ['source_file', result.events.some((event) => Boolean(event.source_file))],
    ['source_event_id', result.events.some((event) => event.source_event_id != null)],
    ['sample_id', result.events.some((event) => event.sample_id != null)],
    ['cluster_id', result.events.some((event) => event.cluster_id != null)],
    ['trajectory_branch', result.events.some((event) => event.trajectory_branch != null)],
    ['pseudotime', result.events.some((event) => event.pseudotime != null)],
    ['predicted_identity', result.events.some((event) => event.predicted_identity != null)],
    ['identity_score', result.events.some((event) => event.identity_score != null)],
    ['identity_margin', result.events.some((event) => event.identity_margin != null)],
  ] as const
  const includedOptional = optional.filter(([, included]) => included).map(([key]) => key)
  const headers = [
    'event_id',
    ...coordinates.map((key) => coordinateLabel(result, key)),
    ...markers.map((marker) => marker.marker),
    ...includedOptional,
  ]
  const rows = result.events.map((event) => [
    event.event_id,
    ...coordinates.map((key) => eventValue(event, key)),
    ...markers.map((marker) => event[marker.column] ?? ''),
    ...includedOptional.map((key) => event[key] ?? ''),
  ])
  return [headers, ...rows].map((row) => row.map(csvCell).join(',')).join('\r\n')
}

export function resultPlotSvg({
  result,
  xKey,
  yKey,
  colors,
}: {
  result: AnalysisRunResult
  xKey: CoordinateKey
  yKey: CoordinateKey
  colors: string[]
}): string {
  const width = 720
  const height = 720
  const margin = { left: 72, right: 24, top: 24, bottom: 64 }
  const xValues = result.events.map((event) => eventValue(event, xKey))
  const yValues = result.events.map((event) => eventValue(event, yKey))
  const extent = (values: number[]) => {
    const minimum = Math.min(...values)
    const maximum = Math.max(...values)
    const pad = (maximum - minimum || 1) * 0.03
    return [minimum - pad, maximum + pad]
  }
  const [xMin, xMax] = extent(xValues)
  const [yMin, yMax] = extent(yValues)
  const sx = (value: number) => margin.left + ((value - xMin) / (xMax - xMin || 1)) * (width - margin.left - margin.right)
  const sy = (value: number) => height - margin.bottom - ((value - yMin) / (yMax - yMin || 1)) * (height - margin.top - margin.bottom)
  const points = result.events.map((event, index) => (
    `<circle cx="${sx(eventValue(event, xKey)).toFixed(2)}" cy="${sy(eventValue(event, yKey)).toFixed(2)}" r="1.7" fill="${colors[index] ?? '#247f9e'}" fill-opacity=".82"/>`
  )).join('')
  return [
    '<?xml version="1.0" encoding="UTF-8"?>',
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`,
    '<rect width="100%" height="100%" fill="#f8f7f3"/>',
    `<rect x="${margin.left}" y="${margin.top}" width="${width - margin.left - margin.right}" height="${height - margin.top - margin.bottom}" fill="none" stroke="#c7c3ba"/>`,
    points,
    `<text x="${width / 2}" y="${height - 20}" text-anchor="middle" font-family="Arial,sans-serif" font-size="16" fill="#746f68">${coordinateLabel(result, xKey)}</text>`,
    `<text x="20" y="${height / 2}" text-anchor="middle" transform="rotate(-90 20 ${height / 2})" font-family="Arial,sans-serif" font-size="16" fill="#746f68">${coordinateLabel(result, yKey)}</text>`,
    '</svg>',
  ].join('')
}

function safeJson(value: unknown): string {
  return JSON.stringify(value).replaceAll('</', '<\\/')
}

export function resultPlotHtml({
  result,
  xKey,
  yKey,
  zKey,
  colors,
  mode,
}: {
  result: AnalysisRunResult
  xKey: CoordinateKey
  yKey: CoordinateKey
  zKey: CoordinateKey
  colors: string[]
  mode: '2d' | '3d'
}): string {
  const trace = {
    type: mode === '3d' ? 'scatter3d' : 'scattergl',
    mode: 'markers',
    x: result.events.map((event) => eventValue(event, xKey)),
    y: result.events.map((event) => eventValue(event, yKey)),
    ...(mode === '3d' ? { z: result.events.map((event) => eventValue(event, zKey)) } : {}),
    text: result.events.map((event) => `Event ${event.event_id}`),
    marker: { size: mode === '3d' ? 3 : 4, opacity: 0.82, color: colors },
    hovertemplate: '%{text}<extra></extra>',
  }
  const layout = {
    paper_bgcolor: '#f8f7f3',
    plot_bgcolor: '#f8f7f3',
    margin: { l: 55, r: 20, t: 20, b: 50 },
    ...(mode === '3d'
      ? { scene: { aspectmode: 'cube', xaxis: { title: coordinateLabel(result, xKey), showgrid: false, zeroline: false }, yaxis: { title: coordinateLabel(result, yKey), showgrid: false, zeroline: false }, zaxis: { title: coordinateLabel(result, zKey), showgrid: false, zeroline: false } } }
      : { xaxis: { title: coordinateLabel(result, xKey), showgrid: false, zeroline: false }, yaxis: { title: coordinateLabel(result, yKey), showgrid: false, zeroline: false } }),
  }
  return `<!doctype html><html><head><meta charset="utf-8"><title>${result.metadata.display_name ?? result.metadata.method.name}</title><script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script></head><body style="margin:0;background:#f8f7f3"><div id="plot" style="width:100vw;height:100vh"></div><script>Plotly.newPlot("plot",[${safeJson(trace)}],${safeJson(layout)},{responsive:true,displaylogo:false});</script></body></html>`
}

export function downloadText(filename: string, content: string, type: string) {
  const url = URL.createObjectURL(new Blob([content], { type }))
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  anchor.click()
  window.setTimeout(() => URL.revokeObjectURL(url), 1000)
}

export function downloadDataUrl(filename: string, url: string) {
  const anchor = document.createElement('a')
  anchor.href = url
  anchor.download = filename
  anchor.click()
}
