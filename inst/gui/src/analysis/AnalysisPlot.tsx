import { useEffect, useMemo, useRef, useState } from 'react'
import type {
  AnalysisEventPayload,
  AnalysisPlotState,
  AnalysisPopulation,
  GateDraft,
  GatePoint,
} from './types'
import { computeDensityBuckets } from '../gating/GatingData.js'
import { DENSITY_PALETTE } from '../gating/GatingPalette.js'
import { paletteColors } from './resultColor'
import {
  densityContourSegments,
  formatAxisValue,
  histogramCounts,
  inverseTransformValue,
  linearScale,
  nearestEvent,
  PLOT_MARGIN,
  robustExtent,
  transformValue,
} from './geometry'

type PlotTool = 'select' | 'rectangle' | 'ellipse' | 'polygon' | 'range' | 'root'

type ScaleContext = {
  x: ReturnType<typeof linearScale>
  y: ReturnType<typeof linearScale>
  width: number
  height: number
}

type TwoDimensionalGateEdit = {
  gate: AnalysisPopulation
  mode: 'move' | 'x-min' | 'x-max' | 'y-min' | 'y-max' | 'radius-x' | 'radius-y' | 'vertex'
  vertex?: number
  startRaw: GatePoint
  geometry: NonNullable<AnalysisPopulation['geometry']>
}

function containsScreenPoint(point: GatePoint, polygon: GatePoint[]) {
  let inside = false
  for (let current = 0, previous = polygon.length - 1; current < polygon.length; previous = current, current += 1) {
    const a = polygon[current]
    const b = polygon[previous]
    const crosses = ((a.y > point.y) !== (b.y > point.y))
      && point.x < ((b.x - a.x) * (point.y - a.y)) / (b.y - a.y || 1) + a.x
    if (crosses) inside = !inside
  }
  return inside
}

function screenDistance(a: GatePoint, b: GatePoint) {
  return Math.hypot(a.x - b.x, a.y - b.y)
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
  onSelectGate: (populationId: string) => void
  onUpdateGate: (populationId: string, geometry: NonNullable<AnalysisPopulation['geometry']>) => void
}

function colorForMarker(value: number, low: number, high: number, colors: string[]) {
  const normalized = Math.max(0, Math.min(1, (value - low) / (high - low || 1)))
  return colors[Math.min(colors.length - 1, Math.floor(normalized * colors.length))] ?? '#197783'
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
    : node.type === 'ellipse'
      ? Array.from({ length: 64 }, (_, index) => {
          const angle = (index / 64) * Math.PI * 2 - Math.PI / 2
          return {
            x: (geometry.center_x ?? 0) + Math.cos(angle) * (geometry.radius_x ?? 0),
            y: (geometry.center_y ?? 0) + Math.sin(angle) * (geometry.radius_y ?? 0),
          }
        })
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
  onSelectGate,
  onUpdateGate,
}: AnalysisPlotProps) {
  const hostRef = useRef<HTMLDivElement>(null)
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const scaleRef = useRef<ScaleContext | null>(null)
  const [size, setSize] = useState({ width: 640, height: 430 })
  const [dragStart, setDragStart] = useState<GatePoint | null>(null)
  const [dragCurrent, setDragCurrent] = useState<GatePoint | null>(null)
  const [polygon, setPolygon] = useState<GatePoint[]>([])
  const [rangeEdit, setRangeEdit] = useState<{
    gateId: string
    mode: 'minimum' | 'maximum' | 'move'
    startX: number
    min: number
    max: number
  } | null>(null)
  const [twoDimensionalEdit, setTwoDimensionalEdit] = useState<TwoDimensionalGateEdit | null>(null)

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
    context.fillStyle = '#f8f7f3'
    context.fillRect(0, 0, size.width, size.height)

    const automaticXDomain = robustExtent(transformed.map((event) => event.x))
    const xDomain: [number, number] = [
      plot.x_min == null ? automaticXDomain[0] : transformValue(plot.x_min, plot.x_transform),
      plot.x_max == null ? automaticXDomain[1] : transformValue(plot.x_max, plot.x_transform),
    ]
    const histogram = histogramCounts(transformed.map((event) => event.x), xDomain)
    const automaticYDomain: [number, number] = plot.type === 'histogram'
      ? [0, Math.max(...histogram, 1) * 1.08]
      : robustExtent(transformed.map((event) => event.y))
    const yDomain: [number, number] = plot.type === 'histogram'
      ? automaticYDomain
      : [
          plot.y_min == null ? automaticYDomain[0] : transformValue(plot.y_min, plot.y_transform),
          plot.y_max == null ? automaticYDomain[1] : transformValue(plot.y_max, plot.y_transform),
        ]
    const xScale = linearScale(xDomain, [PLOT_MARGIN.left, size.width - PLOT_MARGIN.right])
    const yScale = linearScale(yDomain, [size.height - PLOT_MARGIN.bottom, PLOT_MARGIN.top])
    const scale = { x: xScale, y: yScale, width: size.width, height: size.height }
    scaleRef.current = scale

    context.strokeStyle = 'rgba(135, 130, 120, 0.16)'
    context.fillStyle = '#746f68'
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
      const rawY = plot.type === 'histogram'
        ? yDomain[0] + fraction * (yDomain[1] - yDomain[0])
        : inverseTransformValue(yDomain[0] + fraction * (yDomain[1] - yDomain[0]), plot.y_transform)
      context.fillText(formatAxisValue(rawX), x - 12, size.height - 18)
      context.fillText(formatAxisValue(rawY), 6, y + 3)
    }

    if (transformed.length && plot.type === 'scatter') {
      const densityBuckets = computeDensityBuckets(transformed, 'x', 'y', xDomain, yDomain)
      const markerValues = transformed.map((event) => event.color ?? 0).filter(Number.isFinite).sort((a, b) => a - b)
      const markerLow = markerValues[Math.floor(markerValues.length * 0.02)] ?? 0
      const markerHigh = markerValues[Math.floor(markerValues.length * 0.98)] ?? 1
      const markerPalette = paletteColors((plot.color_palette || 'viridis') as Parameters<typeof paletteColors>[0], 128)
      densityBuckets.forEach((bucket, colorIndex) => {
        if (!bucket.length) return
        for (const index of bucket) {
          const event = transformed[index]
          if (!event) continue
          context.fillStyle = plot.color_by === 'marker' && Number.isFinite(event.color)
            ? colorForMarker(event.color ?? 0, markerLow, markerHigh, markerPalette)
            : DENSITY_PALETTE[colorIndex] ?? 'rgba(0, 0, 255, 0.7)'
          context.globalAlpha = plot.opacity ?? 0.82
          const pointSize = plot.point_size ?? 2.4
          context.fillRect(xScale.toScreen(event.x) - pointSize / 2, yScale.toScreen(event.y) - pointSize / 2, pointSize, pointSize)
        }
      })
      context.globalAlpha = 1
    }

    if (transformed.length && plot.type === 'contour') {
      const contours = densityContourSegments(transformed, xDomain, yDomain)
      context.lineWidth = 1.35
      context.lineCap = 'round'
      context.lineJoin = 'round'
      for (const segment of contours) {
        const colorIndex = Math.min(DENSITY_PALETTE.length - 1, Math.round((segment.level / 0.64) * (DENSITY_PALETTE.length - 1)))
        context.strokeStyle = DENSITY_PALETTE[colorIndex] ?? '#167f95'
        context.beginPath()
        context.moveTo(xScale.toScreen(segment.x1), yScale.toScreen(segment.y1))
        context.lineTo(xScale.toScreen(segment.x2), yScale.toScreen(segment.y2))
        context.stroke()
      }
    }

    if (transformed.length && plot.type === 'hexbin') {
      const columns = 34
      const rows = 34
      const bins = new Map<string, { x: number; y: number; count: number }>()
      for (const event of transformed) {
        const column = Math.max(0, Math.min(columns - 1, Math.floor(((event.x - xDomain[0]) / (xDomain[1] - xDomain[0] || 1)) * columns)))
        const row = Math.max(0, Math.min(rows - 1, Math.floor(((event.y - yDomain[0]) / (yDomain[1] - yDomain[0] || 1)) * rows)))
        const key = `${column}:${row}`
        const bin = bins.get(key) ?? { x: column, y: row, count: 0 }
        bin.count += 1
        bins.set(key, bin)
      }
      const maximum = Math.max(1, ...Array.from(bins.values(), (bin) => bin.count))
      const cellWidth = (size.width - PLOT_MARGIN.left - PLOT_MARGIN.right) / columns
      const cellHeight = (size.height - PLOT_MARGIN.top - PLOT_MARGIN.bottom) / rows
      for (const bin of bins.values()) {
        const centerX = PLOT_MARGIN.left + (bin.x + 0.5) * cellWidth
        const centerY = size.height - PLOT_MARGIN.bottom - (bin.y + 0.5) * cellHeight
        const radius = Math.max(2, Math.min(cellWidth, cellHeight) * 0.62)
        const colorIndex = Math.min(DENSITY_PALETTE.length - 1, Math.floor((bin.count / maximum) * DENSITY_PALETTE.length))
        context.fillStyle = DENSITY_PALETTE[colorIndex] ?? '#197783'
        context.beginPath()
        for (let vertex = 0; vertex < 6; vertex += 1) {
          const angle = (Math.PI / 3) * vertex
          const x = centerX + Math.cos(angle) * radius
          const y = centerY + Math.sin(angle) * radius
          if (vertex === 0) context.moveTo(x, y)
          else context.lineTo(x, y)
        }
        context.closePath()
        context.fill()
      }
    }

    if (transformed.length && plot.type === 'histogram') {
      const binWidth = (size.width - PLOT_MARGIN.left - PLOT_MARGIN.right) / histogram.length
      context.fillStyle = 'rgba(25, 119, 131, 0.24)'
      context.strokeStyle = '#197783'
      context.lineWidth = 1.2
      context.beginPath()
      context.moveTo(PLOT_MARGIN.left, size.height - PLOT_MARGIN.bottom)
      histogram.forEach((count, index) => {
        const x = PLOT_MARGIN.left + index * binWidth
        const y = yScale.toScreen(count)
        context.lineTo(x, y)
        context.lineTo(x + binWidth, y)
      })
      context.lineTo(size.width - PLOT_MARGIN.right, size.height - PLOT_MARGIN.bottom)
      context.closePath()
      context.fill()
      context.stroke()
    }

    if (overlayPayload?.events?.length && plot.type !== 'histogram') {
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

    if (plot.type === 'histogram') {
      for (const gate of childGates) {
        if (gate.type !== 'range' || gate.x !== plot.x || !gate.geometry) continue
        const left = xScale.toScreen(transformValue(gate.geometry.min ?? 0, plot.x_transform))
        const right = xScale.toScreen(transformValue(gate.geometry.max ?? 0, plot.x_transform))
        context.fillStyle = gate.role === 'positive' ? 'rgba(88, 207, 130, .12)' : gate.role === 'negative' ? 'rgba(237, 106, 103, .12)' : 'rgba(244, 207, 97, .12)'
        context.strokeStyle = gate.role === 'positive' ? '#58cf82' : gate.role === 'negative' ? '#ed6a67' : '#d3a91d'
        context.lineWidth = 1.6
        context.fillRect(Math.min(left, right), PLOT_MARGIN.top, Math.abs(right - left), size.height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)
        context.strokeRect(Math.min(left, right), PLOT_MARGIN.top, Math.abs(right - left), size.height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)
        context.fillStyle = '#393a37'
        context.font = '600 10px "Avenir Next", sans-serif'
        context.fillText(gate.name, Math.min(left, right) + 4, PLOT_MARGIN.top + 12)
      }
    }

    for (const gate of plot.type === 'histogram' ? [] : childGates) {
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
      context.fillStyle = '#393a37'
      context.font = '600 10px "Avenir Next", sans-serif'
      context.fillText(gate.name, points[0].x + 4, points[0].y - 5)
    }

    const root = plot.type === 'histogram' ? undefined : payload?.events.find((event) => event.event_id === rootEventId)
    if (root) {
      const x = xScale.toScreen(transformValue(root.x, plot.x_transform))
      const y = yScale.toScreen(transformValue(root.y, plot.y_transform))
      context.strokeStyle = '#1f3737'
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

  function finishEllipse(end: GatePoint) {
    const scale = scaleRef.current
    if (!scale || !dragStart) return
    const x0 = inverseTransformValue(scale.x.fromScreen(dragStart.x), plot.x_transform)
    const x1 = inverseTransformValue(scale.x.fromScreen(end.x), plot.x_transform)
    const y0 = inverseTransformValue(scale.y.fromScreen(dragStart.y), plot.y_transform)
    const y1 = inverseTransformValue(scale.y.fromScreen(end.y), plot.y_transform)
    onGateDraft({
      type: 'ellipse', x: plot.x, y: plot.y,
      geometry: {
        center_x: (x0 + x1) / 2,
        center_y: (y0 + y1) / 2,
        radius_x: Math.abs(x1 - x0) / 2,
        radius_y: Math.abs(y1 - y0) / 2,
      },
    })
    setDragStart(null)
    setDragCurrent(null)
  }

  function finishRange(end: GatePoint) {
    const scale = scaleRef.current
    if (!scale || !dragStart) return
    const first = inverseTransformValue(scale.x.fromScreen(dragStart.x), plot.x_transform)
    const second = inverseTransformValue(scale.x.fromScreen(end.x), plot.x_transform)
    onGateDraft({
      type: 'range',
      x: plot.x,
      geometry: { min: Math.min(first, second), max: Math.max(first, second) },
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
          if (tool === 'select' && plot.type === 'histogram' && scaleRef.current) {
            const point = localPoint(event)
            const matching = childGates.filter((gate) => gate.type === 'range' && gate.x === plot.x && gate.geometry)
            for (const gate of matching.reverse()) {
              const minimum = gate.geometry?.min ?? 0
              const maximum = gate.geometry?.max ?? 0
              const left = scaleRef.current.x.toScreen(transformValue(minimum, plot.x_transform))
              const right = scaleRef.current.x.toScreen(transformValue(maximum, plot.x_transform))
              const low = Math.min(left, right)
              const high = Math.max(left, right)
              const mode = Math.abs(point.x - low) <= 7
                ? 'minimum'
                : Math.abs(point.x - high) <= 7
                  ? 'maximum'
                  : point.x >= low && point.x <= high
                    ? 'move'
                    : null
              if (mode) {
                onSelectGate(gate.id)
                setRangeEdit({ gateId: gate.id, mode, startX: point.x, min: minimum, max: maximum })
                event.currentTarget.setPointerCapture(event.pointerId)
                return
              }
            }
          }
          if (tool === 'select' && plot.type !== 'histogram' && scaleRef.current) {
            const point = localPoint(event)
            const rawPoint = {
              x: inverseTransformValue(scaleRef.current.x.fromScreen(point.x), plot.x_transform),
              y: inverseTransformValue(scaleRef.current.y.fromScreen(point.y), plot.y_transform),
            }
            const matching = childGates.filter((gate) => gate.x === plot.x && gate.y === plot.y && gate.geometry)
            for (const gate of matching.reverse()) {
              const screenPoints = gateScreenPoints(gate, scaleRef.current, plot)
              if (screenPoints.length < 3) continue
              let mode: TwoDimensionalGateEdit['mode'] | null = null
              let vertex: number | undefined
              if (gate.type === 'polygon') {
                const nearest = screenPoints.reduce((best, candidate, index) => {
                  const distance = screenDistance(point, candidate)
                  return distance < best.distance ? { index, distance } : best
                }, { index: -1, distance: Number.POSITIVE_INFINITY })
                if (nearest.distance <= 8) {
                  mode = 'vertex'
                  vertex = nearest.index
                }
              } else if (gate.type === 'rectangle') {
                const left = Math.min(...screenPoints.map((candidate) => candidate.x))
                const right = Math.max(...screenPoints.map((candidate) => candidate.x))
                const top = Math.min(...screenPoints.map((candidate) => candidate.y))
                const bottom = Math.max(...screenPoints.map((candidate) => candidate.y))
                if (Math.abs(point.x - left) <= 7) mode = 'x-min'
                else if (Math.abs(point.x - right) <= 7) mode = 'x-max'
                else if (Math.abs(point.y - bottom) <= 7) mode = 'y-min'
                else if (Math.abs(point.y - top) <= 7) mode = 'y-max'
              } else if (gate.type === 'ellipse') {
                const geometry = gate.geometry ?? {}
                const center = {
                  x: scaleRef.current.x.toScreen(transformValue(geometry.center_x ?? 0, plot.x_transform)),
                  y: scaleRef.current.y.toScreen(transformValue(geometry.center_y ?? 0, plot.y_transform)),
                }
                const xHandle = {
                  x: scaleRef.current.x.toScreen(transformValue((geometry.center_x ?? 0) + (geometry.radius_x ?? 0), plot.x_transform)),
                  y: center.y,
                }
                const yHandle = {
                  x: center.x,
                  y: scaleRef.current.y.toScreen(transformValue((geometry.center_y ?? 0) + (geometry.radius_y ?? 0), plot.y_transform)),
                }
                if (screenDistance(point, xHandle) <= 9) mode = 'radius-x'
                else if (screenDistance(point, yHandle) <= 9) mode = 'radius-y'
              }
              if (!mode && containsScreenPoint(point, screenPoints)) mode = 'move'
              if (mode) {
                onSelectGate(gate.id)
                setTwoDimensionalEdit({
                  gate,
                  mode,
                  vertex,
                  startRaw: rawPoint,
                  geometry: { ...(gate.geometry ?? {}) },
                })
                event.currentTarget.setPointerCapture(event.pointerId)
                return
              }
            }
          }
          const canDrawTwoDimensional = plot.type !== 'histogram' && (tool === 'rectangle' || tool === 'ellipse')
          const canDrawRange = plot.type === 'histogram' && tool === 'range'
          if (!canDrawTwoDimensional && !canDrawRange) return
          const point = localPoint(event)
          setDragStart(point)
          setDragCurrent(point)
          event.currentTarget.setPointerCapture(event.pointerId)
        }}
        onPointerMove={(event) => {
          if (rangeEdit && scaleRef.current) {
            const raw = inverseTransformValue(scaleRef.current.x.fromScreen(localPoint(event).x), plot.x_transform)
            const startRaw = inverseTransformValue(scaleRef.current.x.fromScreen(rangeEdit.startX), plot.x_transform)
            if (rangeEdit.mode === 'minimum') {
              onUpdateGate(rangeEdit.gateId, { min: Math.min(raw, rangeEdit.max), max: rangeEdit.max })
            } else if (rangeEdit.mode === 'maximum') {
              onUpdateGate(rangeEdit.gateId, { min: rangeEdit.min, max: Math.max(raw, rangeEdit.min) })
            } else {
              const delta = raw - startRaw
              onUpdateGate(rangeEdit.gateId, { min: rangeEdit.min + delta, max: rangeEdit.max + delta })
            }
            return
          }
          if (twoDimensionalEdit && scaleRef.current) {
            const screen = localPoint(event)
            const raw = {
              x: inverseTransformValue(scaleRef.current.x.fromScreen(screen.x), plot.x_transform),
              y: inverseTransformValue(scaleRef.current.y.fromScreen(screen.y), plot.y_transform),
            }
            const dx = raw.x - twoDimensionalEdit.startRaw.x
            const dy = raw.y - twoDimensionalEdit.startRaw.y
            const geometry = { ...twoDimensionalEdit.geometry }
            if (twoDimensionalEdit.mode === 'move') {
              if (twoDimensionalEdit.gate.type === 'rectangle') {
                geometry.x_min = (geometry.x_min ?? 0) + dx
                geometry.x_max = (geometry.x_max ?? 0) + dx
                geometry.y_min = (geometry.y_min ?? 0) + dy
                geometry.y_max = (geometry.y_max ?? 0) + dy
              } else if (twoDimensionalEdit.gate.type === 'ellipse') {
                geometry.center_x = (geometry.center_x ?? 0) + dx
                geometry.center_y = (geometry.center_y ?? 0) + dy
              } else if (twoDimensionalEdit.gate.type === 'polygon') {
                geometry.points = (geometry.points ?? []).map((point) => ({ x: point.x + dx, y: point.y + dy }))
              }
            } else if (twoDimensionalEdit.mode === 'x-min') geometry.x_min = Math.min(raw.x, geometry.x_max ?? raw.x)
            else if (twoDimensionalEdit.mode === 'x-max') geometry.x_max = Math.max(raw.x, geometry.x_min ?? raw.x)
            else if (twoDimensionalEdit.mode === 'y-min') geometry.y_min = Math.min(raw.y, geometry.y_max ?? raw.y)
            else if (twoDimensionalEdit.mode === 'y-max') geometry.y_max = Math.max(raw.y, geometry.y_min ?? raw.y)
            else if (twoDimensionalEdit.mode === 'radius-x') geometry.radius_x = Math.abs(raw.x - (geometry.center_x ?? 0))
            else if (twoDimensionalEdit.mode === 'radius-y') geometry.radius_y = Math.abs(raw.y - (geometry.center_y ?? 0))
            else if (twoDimensionalEdit.mode === 'vertex' && twoDimensionalEdit.vertex != null) {
              const points = [...(geometry.points ?? [])]
              points[twoDimensionalEdit.vertex] = raw
              geometry.points = points
            }
            onUpdateGate(twoDimensionalEdit.gate.id, geometry)
            return
          }
          if ((tool === 'rectangle' || tool === 'ellipse' || tool === 'range') && dragStart) setDragCurrent(localPoint(event))
        }}
        onPointerUp={(event) => {
          if (rangeEdit) {
            setRangeEdit(null)
            return
          }
          if (twoDimensionalEdit) {
            setTwoDimensionalEdit(null)
            return
          }
          if (tool === 'rectangle' && dragStart) finishRectangle(localPoint(event))
          if (tool === 'ellipse' && dragStart) finishEllipse(localPoint(event))
          if (tool === 'range' && dragStart) finishRange(localPoint(event))
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
        {dragStart && dragCurrent && tool === 'rectangle' ? (
          <rect
            x={Math.min(dragStart.x, dragCurrent.x)}
            y={Math.min(dragStart.y, dragCurrent.y)}
            width={Math.abs(dragCurrent.x - dragStart.x)}
            height={Math.abs(dragCurrent.y - dragStart.y)}
          />
        ) : null}
        {dragStart && dragCurrent && tool === 'ellipse' ? (
          <ellipse
            cx={(dragStart.x + dragCurrent.x) / 2}
            cy={(dragStart.y + dragCurrent.y) / 2}
            rx={Math.abs(dragCurrent.x - dragStart.x) / 2}
            ry={Math.abs(dragCurrent.y - dragStart.y) / 2}
          />
        ) : null}
        {dragStart && dragCurrent && tool === 'range' ? (
          <rect
            x={Math.min(dragStart.x, dragCurrent.x)}
            y={PLOT_MARGIN.top}
            width={Math.abs(dragCurrent.x - dragStart.x)}
            height={Math.max(0, size.height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)}
          />
        ) : null}
        {polygon.length ? <polyline points={polygon.map((point) => `${point.x},${point.y}`).join(' ')} /> : null}
      </svg>
      {!payload ? <div className="analysis-plot-loading">Loading events…</div> : null}
      {payload?.events.length === 0 ? <div className="analysis-plot-loading">No events in this population</div> : null}
      <span className="analysis-axis-label analysis-axis-x">{plot.x}</span>
      <span className="analysis-axis-label analysis-axis-y">{plot.type === 'histogram' ? 'Count' : plot.y}</span>
    </div>
  )
}

export type { PlotTool }
