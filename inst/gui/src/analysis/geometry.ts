import type { AnalysisEvent, AxisTransform, GatePoint } from './types'

export const PLOT_MARGIN = { top: 14, right: 16, bottom: 42, left: 54 }

export function transformValue(value: number, transform: AxisTransform): number {
  if (transform === 'asinh') return Math.asinh(value / 150)
  if (transform === 'biexponential') return Math.sign(value) * Math.log10(1 + Math.abs(value) / 50)
  return value
}

export function inverseTransformValue(value: number, transform: AxisTransform): number {
  if (transform === 'asinh') return Math.sinh(value) * 150
  if (transform === 'biexponential') return Math.sign(value) * 50 * (10 ** Math.abs(value) - 1)
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

export function histogramCounts(values: number[], domain: [number, number], bins = 64): number[] {
  const counts = Array.from({ length: bins }, () => 0)
  const span = domain[1] - domain[0] || 1
  for (const value of values) {
    if (!Number.isFinite(value) || value < domain[0] || value > domain[1]) continue
    const index = Math.max(0, Math.min(bins - 1, Math.floor(((value - domain[0]) / span) * bins)))
    counts[index] += 1
  }
  return counts
}

export type ContourSegment = {
  x1: number
  y1: number
  x2: number
  y2: number
  level: number
}

export function densityContourSegments(
  events: Array<{ x: number; y: number }>,
  xDomain: [number, number],
  yDomain: [number, number],
  bins = 96,
): ContourSegment[] {
  if (!events.length || bins < 4) return []
  const side = bins + 1
  const field = new Float32Array(side * side)
  const xSpan = xDomain[1] - xDomain[0] || 1
  const ySpan = yDomain[1] - yDomain[0] || 1
  for (const event of events) {
    if (!Number.isFinite(event.x) || !Number.isFinite(event.y)) continue
    const x = Math.max(0, Math.min(bins, Math.round(((event.x - xDomain[0]) / xSpan) * bins)))
    const y = Math.max(0, Math.min(bins, Math.round(((event.y - yDomain[0]) / ySpan) * bins)))
    field[y * side + x] += 1
  }
  let source = field
  for (let pass = 0; pass < 6; pass += 1) {
    const target = new Float32Array(source.length)
    for (let y = 0; y < side; y += 1) {
      for (let x = 0; x < side; x += 1) {
        let sum = 0
        let weight = 0
        for (let dy = -1; dy <= 1; dy += 1) {
          for (let dx = -1; dx <= 1; dx += 1) {
            const nx = Math.max(0, Math.min(bins, x + dx))
            const ny = Math.max(0, Math.min(bins, y + dy))
            const localWeight = dx === 0 && dy === 0 ? 4 : dx === 0 || dy === 0 ? 2 : 1
            sum += source[ny * side + nx] * localWeight
            weight += localWeight
          }
        }
        target[y * side + x] = sum / weight
      }
    }
    source = target
  }
  let maximum = 0
  for (const value of source) maximum = Math.max(maximum, value)
  if (maximum <= 0) return []
  const cases: Record<number, Array<[number, number]>> = {
    1: [[3, 0]], 2: [[0, 1]], 3: [[3, 1]], 4: [[1, 2]],
    5: [[3, 0], [1, 2]], 6: [[0, 2]], 7: [[3, 2]], 8: [[2, 3]],
    9: [[0, 2]], 10: [[0, 1], [2, 3]], 11: [[1, 2]], 12: [[3, 1]],
    13: [[0, 1]], 14: [[3, 0]],
  }
  const edgePoint = (edge: number, x: number, y: number, threshold: number): [number, number] => {
    const topLeft = source[y * side + x]
    const topRight = source[y * side + x + 1]
    const bottomRight = source[(y + 1) * side + x + 1]
    const bottomLeft = source[(y + 1) * side + x]
    const interpolate = (first: number, second: number) => Math.max(0, Math.min(1, (threshold - first) / (second - first || 1)))
    const coordinates: Array<[number, number]> = [
      [x + interpolate(topLeft, topRight), y],
      [x + 1, y + interpolate(topRight, bottomRight)],
      [x + 1 - interpolate(bottomRight, bottomLeft), y + 1],
      [x, y + 1 - interpolate(bottomLeft, topLeft)],
    ]
    const point = coordinates[edge] ?? coordinates[0]
    return [
      xDomain[0] + (point[0] / bins) * xSpan,
      yDomain[0] + (point[1] / bins) * ySpan,
    ]
  }
  const segments: ContourSegment[] = []
  for (const level of [0.08, 0.16, 0.28, 0.44, 0.64]) {
    const threshold = maximum * level
    for (let y = 0; y < bins; y += 1) {
      for (let x = 0; x < bins; x += 1) {
        const topLeft = source[y * side + x] >= threshold ? 1 : 0
        const topRight = source[y * side + x + 1] >= threshold ? 2 : 0
        const bottomRight = source[(y + 1) * side + x + 1] >= threshold ? 4 : 0
        const bottomLeft = source[(y + 1) * side + x] >= threshold ? 8 : 0
        for (const [firstEdge, secondEdge] of cases[topLeft | topRight | bottomRight | bottomLeft] ?? []) {
          const first = edgePoint(firstEdge, x, y, threshold)
          const second = edgePoint(secondEdge, x, y, threshold)
          segments.push({ x1: first[0], y1: first[1], x2: second[0], y2: second[1], level })
        }
      }
    }
  }
  return segments
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
