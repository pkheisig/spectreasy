export function reconcileGateCsvRows(rows, files) {
  const knownFiles = new Set(
    (files || []).map((file) => String(file?.filename || '').trim()).filter(Boolean),
  )
  const ignoredFiles = new Set()
  const negativeDisabledFiles = new Set(
    (files || [])
      .filter((file) => file?.uses_negative_histogram_gate === false || String(file?.uses_negative_histogram_gate || '').toLowerCase() === 'false')
      .map((file) => String(file?.filename || '').trim())
      .filter(Boolean),
  )
  let ignoredRowCount = 0

  const compatibleRows = (rows || []).filter((row) => {
    const gateType = String(row?.gate_type || '').trim()
    const filename = String(row?.filename || '').trim()
    if (gateType === 'negative' && negativeDisabledFiles.has(filename)) return false
    if (gateType === 'setting' || !filename || knownFiles.has(filename)) return true

    ignoredFiles.add(filename)
    ignoredRowCount += 1
    return false
  })

  return {
    rows: compatibleRows,
    ignoredFiles: Array.from(ignoredFiles),
    ignoredRowCount,
  }
}
