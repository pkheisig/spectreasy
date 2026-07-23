function finiteExtent(values) {
  let minimum = Infinity
  let maximum = -Infinity
  for (const item of values || []) {
    const value = Number(item)
    if (!Number.isFinite(value)) continue
    minimum = Math.min(minimum, value)
    maximum = Math.max(maximum, value)
  }
  return Number.isFinite(minimum) ? [minimum, maximum] : null
}

function finiteSorted(values) {
  return (values || [])
    .map(Number)
    .filter(Number.isFinite)
    .sort((left, right) => left - right)
}

function quantile(sorted, probability) {
  if (!sorted.length) return 0
  const position = Math.max(0, Math.min(1, probability)) * (sorted.length - 1)
  const lower = Math.floor(position)
  const upper = Math.ceil(position)
  if (lower === upper) return sorted[lower]
  const weight = position - lower
  return sorted[lower] * (1 - weight) + sorted[upper] * weight
}

export function adaptiveHistogramCofactor(values) {
  const sorted = finiteSorted(values)
  if (sorted.length < 2) return 1

  const q001 = quantile(sorted, 0.001)
  const q05 = quantile(sorted, 0.05)
  const q10 = quantile(sorted, 0.10)
  const q999 = quantile(sorted, 0.999)
  const robustSpan = Math.max(q999 - q001, Number.EPSILON)
  const lowerPopulationScale = Math.max(
    Math.abs(q05),
    Math.abs(q10),
    Math.abs(q10 - q05),
    robustSpan / 1000,
  )
  const cofactor = lowerPopulationScale / Math.sinh(2.5)
  return Math.max(robustSpan / 1e6, Math.min(robustSpan / 2, cofactor))
}

export function adaptiveHistogramTransform(values) {
  const cofactor = adaptiveHistogramCofactor(values)
  return {
    cofactor,
    forward: (value) => Math.asinh(Number(value) / cofactor),
    inverse: (value) => Math.sinh(Number(value)) * cofactor,
  }
}

export function dragHistogramGateInDisplaySpace(vertices, displayDelta, transform, vertexIndex = null) {
  const delta = Number(displayDelta)
  if (!Array.isArray(vertices) || !Number.isFinite(delta) || !transform?.forward || !transform?.inverse) {
    return Array.isArray(vertices) ? vertices.map((vertex) => ({ ...vertex })) : []
  }

  return vertices.map((vertex, index) => {
    if (vertexIndex !== null && index !== vertexIndex) return { ...vertex }
    const displayX = Number(transform.forward(vertex?.x))
    const rawX = Number(transform.inverse(displayX + delta))
    return Number.isFinite(rawX) ? { ...vertex, x: rawX, y: 0 } : { ...vertex }
  })
}

export function histogramDomainIncludingGates(rawDomain, eventValues, gates) {
  const domain = finiteExtent(rawDomain || []) || finiteExtent(eventValues || []) || [0, 1]
  const eventDomain = finiteExtent(eventValues || []) || domain
  const gateValues = (gates || [])
    .flatMap((gate) => gate?.vertices || [])
    .map((vertex) => Number(vertex?.x))
    .filter(Number.isFinite)
    .map((value) => Math.max(eventDomain[0], Math.min(eventDomain[1], value)))

  const expanded = [
    Math.min(domain[0], eventDomain[0], ...gateValues),
    Math.max(domain[1], eventDomain[1], ...gateValues),
  ]
  return expanded[0] === expanded[1] ? [expanded[0], expanded[0] + 1] : expanded
}

export function padHistogramDomain(domain, fraction = 0.04) {
  const finite = finiteExtent(domain || []) || [0, 1]
  const span = finite[1] - finite[0]
  const padding = (span > 0 ? span : Math.max(Math.abs(finite[0]), 1)) * Math.max(0, Number(fraction) || 0)
  return [finite[0] - padding, finite[1] + padding]
}

export function panHistogramDomain(rawDomain, pixelDelta, plotWidth, rawLimits) {
  const domain = finiteExtent(rawDomain || [])
  const limits = finiteExtent(rawLimits || [])
  const width = Number(plotWidth)
  const delta = Number(pixelDelta)
  if (!domain || !limits || !Number.isFinite(width) || width <= 0 || !Number.isFinite(delta)) {
    return domain || [0, 1]
  }

  const span = domain[1] - domain[0]
  const limitSpan = limits[1] - limits[0]
  if (span <= 0 || limitSpan < span) return domain

  const shift = -(delta / width) * span
  let nextMin = domain[0] + shift
  let nextMax = domain[1] + shift
  if (nextMin < limits[0]) {
    nextMin = limits[0]
    nextMax = nextMin + span
  }
  if (nextMax > limits[1]) {
    nextMax = limits[1]
    nextMin = nextMax - span
  }
  return [nextMin, nextMax]
}

export function clearHistogramGatesForFile(gates, filename) {
  const next = { ...gates }
  delete next[`positive:${filename}`]
  delete next[`negative:${filename}`]
  return next
}
