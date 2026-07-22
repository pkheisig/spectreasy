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
