import { useEffect, useMemo, useRef, useState } from 'react'
import type {
  AnalysisEventPayload,
  AnalysisPlotState,
  AnalysisPopulation,
  GateDraft,
  GatePoint,
} from './types'
import {
  densityBuckets,
  formatAxisValue,
  inverseTransformValue,
  linearScale,
  nearestEvent,
  PLOT_MARGIN,
  robustExtent,
  transformValue,
} from './geometry'

type PlotTool = 'select' | 'rectangle' | 'polygon' | 'root'

type ScaleContext = {
  x: ReturnType<typeof linearScale>
  y: ReturnType<typeof linearScale>
  width: number
  height: number
}

type AnalysisPlotProps = {
  plot: AnalysisPlotState
  payload: AnalysisEventPayload | null
  overlayPayload?: AnalysisEventPayload | null
  childGates: AnalysisPopulation[]
  tool: PlotTool
  rootEventId: number | null
  onGateDraft: (draft: GateDraft) => void
  onRootEvent: (eventId: number) => void
}

const densityPalette = ['#315a86', '#247f9e', '#19a899', '#76bf71', '#dfc34f', '#f1843d', '#d94b43']

function colorForMarker(value: number, low: number, high: number) {
  const normalized = Math.max(0, Math.min(1, (value - low) / (high - low || 1)))
  const red = Math.round(45 + normalized * 218)
  const green = Math.round(104 + (1 - Math.abs(normalized - 0.45) * 2) * 90)
  const blue = Math.round(179 - normalized * 126)
  return `rgb(${red}, ${Math.max(58, green)}, ${Math.max(50, blue)})`
}

function gateScreenPoints(node: AnalysisPopulation, scale: ScaleContext, plot: AnalysisPlotState): GatePoint[] {
  const geometry = node.geometry
  if (!geometry) return []
  const raw = node.type === 'rectangle'
    ? [
        { x: geometry.x_min ?? 0, y: geometry.y_min ?? 0 },
        { x: geometry.x_max ?? 0, y: geometry.y_min ?? 0 },
        { x: geometry.x_max ?? 0, y: geometry.y_max ?? 0 },
        { x: geometry.x_min ?? 0, y: geometry.y_max ?? 0 },
      ]
    : geometry.points ?? []
  return raw.map((point) => ({
    x: scale.x.toScreen(transformValue(point.x, plot.x_transform)),
    y: scale.y.toScreen(transformValue(point.y, plot.y_transform)),
  }))
}

export function AnalysisPlot({
  plot,
  payload,
  overlayPayload,
  childGates,
  tool,
  rootEventId,
  onGateDraft,
  onRootEvent,
}: AnalysisPlotProps) {
  const hostRef = useRef<HTMLDivElement>(null)
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const scaleRef = useRef<ScaleContext | null>(null)
  const [size, setSize] = useState({ width: 640, height: 430 })
  const [dragStart, setDragStart] = useState<GatePoint | null>(null)
  const [dragCurrent, setDragCurrent] = useState<GatePoint | null>(null)
  const [polygon, setPolygon] = useState<GatePoint[]>([])

  useEffect(() => {
    const host = hostRef.current
    if (!host) return
    const update = () => setSize({ width: Math.max(360, host.clientWidth), height: Math.max(300, host.clientHeight) })
    update()
    const observer = new ResizeObserver(update)
    observer.observe(host)
    return () => observer.disconnect()
  }, [])

  const transformed = useMemo(() => (payload?.events ?? []).map((event) => ({
    ...event,
    x: transformValue(event.x, plot.x_transform),
    y: transformValue(event.y, plot.y_transform),
  })), [payload?.events, plot.x_transform, plot.y_transform])

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const ratio = window.devicePixelRatio || 1
    canvas.width = Math.round(size.width * ratio)
    canvas.height = Math.round(size.height * ratio)
    canvas.style.width = `${size.width}px`
    canvas.style.height = `${size.height}px`
    const context = canvas.getContext('2d')
    if (!context) return
    context.setTransform(ratio, 0, 0, ratio, 0, 0)
    context.clearRect(0, 0, size.width, size.height)
    context.fillStyle = '#07131d'
    context.fillRect(0, 0, size.width, size.height)

    const xDomain = robustExtent(transformed.map((event) => event.x))
    const yDomain = robustExtent(transformed.map((event) => event.y))
    const xScale = linearScale(xDomain, [PLOT_MARGIN.left, size.width - PLOT_MARGIN.right])
    const yScale = linearScale(yDomain, [size.height - PLOT_MARGIN.bottom, PLOT_MARGIN.top])
    const scale = { x: xScale, y: yScale, width: size.width, height: size.height }
    scaleRef.current = scale

    context.strokeStyle = 'rgba(148, 170, 184, 0.16)'
    context.fillStyle = '#9bb2c0'
    context.font = '10px "Avenir Next", sans-serif'
    context.lineWidth = 1
    for (let index = 0; index < 5; index += 1) {
      const fraction = index / 4
      const x = PLOT_MARGIN.left + fraction * (size.width - PLOT_MARGIN.left - PLOT_MARGIN.right)
      const y = size.height - PLOT_MARGIN.bottom - fraction * (size.height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)
      context.beginPath()
      context.moveTo(x, PLOT_MARGIN.top)
      context.lineTo(x, size.height - PLOT_MARGIN.bottom)
      context.stroke()
      context.beginPath()
      context.moveTo(PLOT_MARGIN.left, y)
      context.lineTo(size.width - PLOT_MARGIN.right, y)
      context.stroke()
      const rawX = inverseTransformValue(xDomain[0] + fraction * (xDomain[1] - xDomain[0]), plot.x_transform)
      const rawY = inverseTransformValue(yDomain[0] + fraction * (yDomain[1] - yDomain[0]), plot.y_transform)
      context.fillText(formatAxisValue(rawX), x - 12, size.height - 18)
      context.fillText(formatAxisValue(rawY), 6, y + 3)
    }

    if (transformed.length) {
      const density = densityBuckets(transformed, xDomain, yDomain)
      const markerValues = transformed.map((event) => event.color ?? 0).filter(Number.isFinite).sort((a, b) => a - b)
      const markerLow = markerValues[Math.floor(markerValues.length * 0.02)] ?? 0
      const markerHigh = markerValues[Math.floor(markerValues.length * 0.98)] ?? 1
      const order = transformed.map((_, index) => index).sort((a, b) => density[a] - density[b])
      for (const index of order) {
        const event = transformed[index]
        const bucket = Math.min(densityPalette.length - 1, Math.floor(density[index] * densityPalette.length))
        context.fillStyle = plot.color_by === 'marker' && Number.isFinite(event.color)
          ? colorForMarker(event.color ?? 0, markerLow, markerHigh)
          : densityPalette[bucket]
        context.globalAlpha = 0.76
        context.beginPath()
        context.arc(xScale.toScreen(event.x), yScale.toScreen(event.y), 1.45, 0, Math.PI * 2)
        context.fill()
      }
      context.globalAlpha = 1
    }

    if (overlayPayload?.events?.length) {
      context.fillStyle = '#ff4f9a'
      context.globalAlpha = 0.75
      for (const event of overlayPayload.events) {
        context.beginPath()
        context.arc(
          xScale.toScreen(transformValue(event.x, plot.x_transform)),
          yScale.toScreen(transformValue(event.y, plot.y_transform)),
          2.1,
          0,
          Math.PI * 2,
        )
        context.fill()
      }
      context.globalAlpha = 1
    }

    for (const gate of childGates) {
      if (gate.x !== plot.x || gate.y !== plot.y) continue
      const points = gateScreenPoints(gate, scale, plot)
      if (points.length < 3) continue
      context.strokeStyle = gate.role === 'positive' ? '#58cf82' : gate.role === 'negative' ? '#ed6a67' : '#f4cf61'
      context.fillStyle = 'rgba(244, 207, 97, 0.08)'
      context.lineWidth = 1.6
      context.beginPath()
      context.moveTo(points[0].x, points[0].y)
      points.slice(1).forEach((point) => context.lineTo(point.x, point.y))
      context.closePath()
      context.fill()
      context.stroke()
      context.fillStyle = '#dce8ee'
      context.font = '600 10px "Avenir Next", sans-serif'
      context.fillText(gate.name, points[0].x + 4, points[0].y - 5)
    }

    const root = payload?.events.find((event) => event.event_id === rootEventId)
    if (root) {
      const x = xScale.toScreen(transformValue(root.x, plot.x_transform))
      const y = yScale.toScreen(transformValue(root.y, plot.y_transform))
      context.strokeStyle = '#f6f1df'
      context.lineWidth = 2
      context.beginPath()
      context.arc(x, y, 6, 0, Math.PI * 2)
      context.stroke()
      context.beginPath()
      context.moveTo(x - 9, y)
      context.lineTo(x + 9, y)
      context.moveTo(x, y - 9)
      context.lineTo(x, y + 9)
      context.stroke()
    }
  }, [childGates, overlayPayload, payload?.events, plot, rootEventId, size, transformed])

  function localPoint(event: React.PointerEvent<HTMLCanvasElement> | React.MouseEvent<HTMLCanvasElement>): GatePoint {
    const bounds = event.currentTarget.getBoundingClientRect()
    return { x: event.clientX - bounds.left, y: event.clientY - bounds.top }
  }

  function finishRectangle(end: GatePoint) {
    const scale = scaleRef.current
    if (!scale || !dragStart) return
    const x0 = inverseTransformValue(scale.x.fromScreen(dragStart.x), plot.x_transform)
    const x1 = inverseTransformValue(scale.x.fromScreen(end.x), plot.x_transform)
    const y0 = inverseTransformValue(scale.y.fromScreen(dragStart.y), plot.y_transform)
    const y1 = inverseTransformValue(scale.y.fromScreen(end.y), plot.y_transform)
    onGateDraft({
      type: 'rectangle', x: plot.x, y: plot.y,
      geometry: { x_min: Math.min(x0, x1), x_max: Math.max(x0, x1), y_min: Math.min(y0, y1), y_max: Math.max(y0, y1) },
    })
    setDragStart(null)
    setDragCurrent(null)
  }

  function finishPolygon(points: GatePoint[]) {
    const scale = scaleRef.current
    if (!scale || points.length < 3) return
    onGateDraft({
      type: 'polygon', x: plot.x, y: plot.y,
      geometry: {
        points: points.map((point) => ({
          x: inverseTransformValue(scale.x.fromScreen(point.x), plot.x_transform),
          y: inverseTransformValue(scale.y.fromScreen(point.y), plot.y_transform),
        })),
      },
    })
    setPolygon([])
  }

  return (
    <div ref={hostRef} className={`analysis-plot-surface tool-${tool}`}>
      <canvas
        ref={canvasRef}
        aria-label={`${plot.x} by ${plot.y} event plot`}
        onPointerDown={(event) => {
          if (tool !== 'rectangle') return
          const point = localPoint(event)
          setDragStart(point)
          setDragCurrent(point)
          event.currentTarget.setPointerCapture(event.pointerId)
        }}
        onPointerMove={(event) => {
          if (tool === 'rectangle' && dragStart) setDragCurrent(localPoint(event))
        }}
        onPointerUp={(event) => {
          if (tool === 'rectangle' && dragStart) finishRectangle(localPoint(event))
        }}
        onClick={(event) => {
          const point = localPoint(event)
          if (tool === 'polygon') setPolygon((current) => [...current, point])
          if (tool === 'root' && payload && scaleRef.current) {
            const found = nearestEvent(payload.events, point, scaleRef.current.x, scaleRef.current.y, plot.x_transform, plot.y_transform)
            if (found) onRootEvent(found.event_id)
          }
        }}
        onDoubleClick={(event) => {
          if (tool !== 'polygon') return
          event.preventDefault()
          const point = localPoint(event)
          finishPolygon([...polygon, point])
        }}
      />
      <svg className="analysis-drawing-layer" width={size.width} height={size.height} aria-hidden="true">
        {dragStart && dragCurrent ? (
          <rect
            x={Math.min(dragStart.x, dragCurrent.x)}
            y={Math.min(dragStart.y, dragCurrent.y)}
            width={Math.abs(dragCurrent.x - dragStart.x)}
            height={Math.abs(dragCurrent.y - dragStart.y)}
          />
        ) : null}
        {polygon.length ? <polyline points={polygon.map((point) => `${point.x},${point.y}`).join(' ')} /> : null}
      </svg>
      {!payload ? <div className="analysis-plot-loading">Loading events…</div> : null}
      {payload?.events.length === 0 ? <div className="analysis-plot-loading">No events in this population</div> : null}
      <span className="analysis-axis-label analysis-axis-x">{plot.x}</span>
      <span className="analysis-axis-label analysis-axis-y">{plot.y}</span>
    </div>
  )
}

export type { PlotTool }
