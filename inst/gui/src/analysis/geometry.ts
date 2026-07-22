import type { AnalysisEvent, AxisTransform, GatePoint } from './types'

export const PLOT_MARGIN = { top: 14, right: 16, bottom: 42, left: 54 }

export function transformValue(value: number, transform: AxisTransform): number {
  if (transform === 'asinh') return Math.asinh(value / 150)
  return value
}

export function inverseTransformValue(value: number, transform: AxisTransform): number {
  if (transform === 'asinh') return Math.sinh(value) * 150
  return value
}

export function robustExtent(values: number[]): [number, number] {
  const finite = values.filter(Number.isFinite).sort((a, b) => a - b)
  if (!finite.length) return [0, 1]
  const low = finite[Math.floor((finite.length - 1) * 0.002)] ?? finite[0]
  const high = finite[Math.floor((finite.length - 1) * 0.998)] ?? finite.at(-1) ?? 1
  if (low === high) return [low - 0.5, high + 0.5]
  const padding = (high - low) * 0.045
  return [low - padding, high + padding]
}

export function linearScale(domain: [number, number], range: [number, number]) {
  const span = domain[1] - domain[0] || 1
  const rangeSpan = range[1] - range[0]
  return {
    toScreen: (value: number) => range[0] + ((value - domain[0]) / span) * rangeSpan,
    fromScreen: (value: number) => domain[0] + ((value - range[0]) / rangeSpan) * span,
  }
}

export function densityBuckets(events: AnalysisEvent[], xDomain: [number, number], yDomain: [number, number], bins = 52): number[] {
  const grid = new Uint32Array(bins * bins)
  const xSpan = xDomain[1] - xDomain[0] || 1
  const ySpan = yDomain[1] - yDomain[0] || 1
  const indexes = events.map((event) => {
    const x = Math.max(0, Math.min(bins - 1, Math.floor(((event.x - xDomain[0]) / xSpan) * bins)))
    const y = Math.max(0, Math.min(bins - 1, Math.floor(((event.y - yDomain[0]) / ySpan) * bins)))
    const index = y * bins + x
    grid[index] += 1
    return index
  })
  let max = 1
  for (const count of grid) if (count > max) max = count
  return indexes.map((index) => grid[index] / max)
}

export function nearestEvent(
  events: AnalysisEvent[],
  point: GatePoint,
  xScale: ReturnType<typeof linearScale>,
  yScale: ReturnType<typeof linearScale>,
  xTransform: AxisTransform,
  yTransform: AxisTransform,
  radius = 12,
): AnalysisEvent | null {
  let best: AnalysisEvent | null = null
  let bestDistance = radius * radius
  for (const event of events) {
    const dx = xScale.toScreen(transformValue(event.x, xTransform)) - point.x
    const dy = yScale.toScreen(transformValue(event.y, yTransform)) - point.y
    const distance = dx * dx + dy * dy
    if (distance <= bestDistance) {
      best = event
      bestDistance = distance
    }
  }
  return best
}

export function formatAxisValue(value: number): string {
  const absolute = Math.abs(value)
  if (absolute >= 1_000_000) return `${(value / 1_000_000).toFixed(1)}M`
  if (absolute >= 1_000) return `${(value / 1_000).toFixed(1).replace('.0', '')}K`
  if (absolute > 0 && absolute < 0.01) return value.toExponential(1)
  return value.toFixed(absolute < 10 ? 1 : 0)
}
