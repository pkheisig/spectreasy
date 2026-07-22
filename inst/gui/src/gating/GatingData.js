import { resolveApiBase, resolveApiToken } from '../apiBase'
import { withGatingApiToken } from '../gatingApi.js'
import {
  DEFAULT_EVENT_COUNT,
  DEFAULT_HISTOGRAM_BINS,
  EVENT_COUNT_STEPS,
  clamp,
  normalizeEventCount,
} from './GatingGeometry.jsx'

const API_BASE = resolveApiBase()
const API_TOKEN = resolveApiToken()

function gatingApiRequest(path, options = {}) {
  const projectPath = window.sessionStorage.getItem('spectreasy-project-path') || ''
  const requestOptions = { ...options }
  let requestPath = path
  if (projectPath && String(requestOptions.method || 'GET').toUpperCase() === 'GET') {
    requestPath += `${requestPath.includes('?') ? '&' : '?'}project_path=${encodeURIComponent(projectPath)}`
  } else if (projectPath && requestOptions.body != null) {
    try {
      requestOptions.body = JSON.stringify({ ...JSON.parse(String(requestOptions.body)), projectPath })
    } catch {
      // Non-JSON uploads do not use project-scoped gating endpoints.
    }
  }
  return fetch(`${API_BASE}${requestPath}`, withGatingApiToken(requestOptions, API_TOKEN)).then((res) => {
    if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
    return res.json()
  })
}

function eventStepIndex(value) {
  const idx = EVENT_COUNT_STEPS.indexOf(normalizeEventCount(value))
  return idx >= 0 ? idx : EVENT_COUNT_STEPS.indexOf(DEFAULT_EVENT_COUNT)
}

function eventStepLabel(value, totalEvents) {
  const numeric = normalizeEventCount(value)
  if (Number.isFinite(totalEvents) && totalEvents > 0 && totalEvents < numeric) {
    return `${totalEvents.toLocaleString()} events`
  }
  return `${numeric.toLocaleString()} events`
}

function getDensityCurve(values, domain, steps = DEFAULT_HISTOGRAM_BINS) {
  const [min, max] = domain
  if (min === max) return []
  const step = (max - min) / steps
  const points = []
  const n = values.length
  if (n === 0) return []
  const bins = new Array(steps + 1).fill(0)
  values.forEach((value) => {
    if (!Number.isFinite(value)) return
    const idx = clamp(Math.floor(((value - min) / (max - min)) * steps), 0, steps)
    bins[idx] += 1
  })
  const smoothRadius = Math.max(1, Math.min(3, Math.floor(steps / 40)))
  for (let i = 0; i <= steps; i++) {
    const x = min + i * step
    let sum = 0
    let weightSum = 0
    for (let j = -smoothRadius; j <= smoothRadius; j++) {
      const idx = i + j
      if (idx < 0 || idx > steps) continue
      const weight = Math.exp(-(j * j) / 4)
      sum += bins[idx] * weight
      weightSum += weight
    }
    points.push({ x: x, y: sum / Math.max(weightSum, 1) })
  }
  return points
}

function densityPalette(size = 64) {
  const jetColors = (v) => {
    const r = clamp(Math.min(4 * v - 1.5, -4 * v + 4.5), 0, 1)
    const g = clamp(Math.min(4 * v - 0.5, -4 * v + 3.5), 0, 1)
    const b = clamp(Math.min(4 * v + 0.5, -4 * v + 2.5), 0, 1)
    return `rgba(${Math.floor(r * 255)}, ${Math.floor(g * 255)}, ${Math.floor(b * 255)}, 0.7)`
  }
  return Array.from({ length: size }, (_, i) => jetColors(i / Math.max(size - 1, 1)))
}

const DENSITY_PALETTE = densityPalette()

function computeDensityBuckets(points, xField, yField, xDomain, yDomain) {
  const n = points.length
  const buckets = Array.from({ length: DENSITY_PALETTE.length }, () => [])
  if (n === 0) return buckets

  const xs = points.map((p) => p[xField])
  const ys = points.map((p) => p[yField])
  const [minX, maxX] = xDomain || [Math.min(...xs), Math.max(...xs)]
  const [minY, maxY] = yDomain || [Math.min(...ys), Math.max(...ys)]

  const rx = maxX - minX || 1
  const ry = maxY - minY || 1

  // Build a smooth density field, then sample it continuously at every event.
  // Assigning one density to every event in a hard grid cell makes the colour
  // field look tiled when a long scatter tail compresses the main population
  // into a relatively small part of the plot.
  const numBins = 256
  const gridSide = numBins + 1
  const grid = new Float32Array(gridSide * gridSide)

  points.forEach((p) => {
    if (p[xField] < minX || p[xField] > maxX || p[yField] < minY || p[yField] > maxY) return
    const bx = clamp(Math.floor(((p[xField] - minX) / rx) * numBins), 0, numBins)
    const by = clamp(Math.floor(((p[yField] - minY) / ry) * numBins), 0, numBins)
    grid[by * gridSide + bx] += 1
  })

  // A separable Gaussian blur gives every grid node a local density estimate.
  // Bilinear sampling below removes the remaining cell boundaries.
  const sigma = 2.0
  const radius = Math.ceil(sigma * 3)
  const kernel = new Float32Array(radius * 2 + 1)
  let kernelSum = 0
  for (let offset = -radius; offset <= radius; offset++) {
    const weight = Math.exp(-(offset * offset) / (2 * sigma * sigma))
    kernel[offset + radius] = weight
    kernelSum += weight
  }
  for (let i = 0; i < kernel.length; i++) kernel[i] /= kernelSum

  const horizontal = new Float32Array(grid.length)
  const smoothed = new Float32Array(grid.length)
  for (let y = 0; y < gridSide; y++) {
    const row = y * gridSide
    for (let x = 0; x < gridSide; x++) {
      let sum = 0
      for (let offset = -radius; offset <= radius; offset++) {
        const nx = clamp(x + offset, 0, numBins)
        sum += grid[row + nx] * kernel[offset + radius]
      }
      horizontal[row + x] = sum
    }
  }
  for (let y = 0; y < gridSide; y++) {
    for (let x = 0; x < gridSide; x++) {
      let sum = 0
      for (let offset = -radius; offset <= radius; offset++) {
        const ny = clamp(y + offset, 0, numBins)
        sum += horizontal[ny * gridSide + x] * kernel[offset + radius]
      }
      smoothed[y * gridSide + x] = sum
    }
  }

  const densities = new Float32Array(n)
  let maxD = 1
  for (let pointIndex = 0; pointIndex < n; pointIndex++) {
    const p = points[pointIndex]
    if (p[xField] < minX || p[xField] > maxX || p[yField] < minY || p[yField] > maxY) continue
    const gx = clamp(((p[xField] - minX) / rx) * numBins, 0, numBins)
    const gy = clamp(((p[yField] - minY) / ry) * numBins, 0, numBins)
    const x0 = Math.floor(gx)
    const y0 = Math.floor(gy)
    const x1 = Math.min(x0 + 1, numBins)
    const y1 = Math.min(y0 + 1, numBins)
    const tx = gx - x0
    const ty = gy - y0
    const top = smoothed[y0 * gridSide + x0] * (1 - tx) + smoothed[y0 * gridSide + x1] * tx
    const bottom = smoothed[y1 * gridSide + x0] * (1 - tx) + smoothed[y1 * gridSide + x1] * tx
    const density = top * (1 - ty) + bottom * ty
    densities[pointIndex] = density
    if (density > maxD) maxD = density
  }

  for (let i = 0; i < n; i++) {
    const bucket = clamp(Math.floor((densities[i] / maxD) * (DENSITY_PALETTE.length - 1)), 0, DENSITY_PALETTE.length - 1)
    buckets[bucket].push(i)
  }
  return buckets
}

export {
  DENSITY_PALETTE,
  computeDensityBuckets,
  eventStepIndex,
  eventStepLabel,
  getDensityCurve,
  normalizeEventCount,
  gatingApiRequest,
}
