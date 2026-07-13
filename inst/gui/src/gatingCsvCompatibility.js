export function reconcileGateCsvRows(rows, files) {
  const knownFiles = new Set(
    (files || []).map((file) => String(file?.filename || '').trim()).filter(Boolean),
  )
  const ignoredFiles = new Set()
  let ignoredRowCount = 0

  const compatibleRows = (rows || []).filter((row) => {
    const gateType = String(row?.gate_type || '').trim()
    const filename = String(row?.filename || '').trim()
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
