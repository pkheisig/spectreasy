export function fileUsesHistogramGates(file) {
  const explicit = file?.uses_histogram_gates
  if (typeof explicit === 'boolean') return explicit
  if (typeof explicit === 'string' && explicit.trim()) {
    return ['true', 't', '1', 'yes'].includes(explicit.trim().toLowerCase())
  }
  return !file?.is_af
}

export function fileUsesNegativeHistogramGate(file) {
  if (!fileUsesHistogramGates(file)) return false
  const explicit = file?.uses_negative_histogram_gate
  if (typeof explicit === 'boolean') return explicit
  if (typeof explicit === 'string' && explicit.trim()) {
    return ['true', 't', '1', 'yes'].includes(explicit.trim().toLowerCase())
  }
  return true
}

export function pruneUnavailableNegativeGates(gates, files) {
  const next = { ...(gates || {}) }
  ;(files || []).forEach((file) => {
    if (!fileUsesNegativeHistogramGate(file)) {
      delete next[`negative:${file?.filename || ''}`]
    }
  })
  return next
}
