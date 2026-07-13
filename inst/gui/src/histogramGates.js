function finiteExtent(values) {
  const finite = values.map(Number).filter(Number.isFinite)
  if (!finite.length) return null
  return [Math.min(...finite), Math.max(...finite)]
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
    Math.min(domain[0], ...gateValues),
    Math.max(domain[1], ...gateValues),
  ]
  return expanded[0] === expanded[1] ? [expanded[0], expanded[0] + 1] : expanded
}

export function clearHistogramGatesForFile(gates, filename) {
  const next = { ...gates }
  delete next[`positive:${filename}`]
  delete next[`negative:${filename}`]
  return next
}
