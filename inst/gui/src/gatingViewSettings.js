const SCATTER_VIEW_PLOTS = ['cell', 'singlet']
const CONTROL_SCOPES = ['cells', 'beads']

export function normalizeViewDomain(value) {
  if (!Array.isArray(value) || value.length < 2) return null
  const lo = Number(value[0])
  const hi = Number(value[1])
  return Number.isFinite(lo) && Number.isFinite(hi) && hi > lo ? [lo, hi] : null
}

export function normalizePlotView(value, includeY = true) {
  const x = normalizeViewDomain(value?.x)
  const y = includeY ? normalizeViewDomain(value?.y) : null
  if (!x && !y) return null
  return { x, y }
}

export function buildViewSettingRows(viewSettings = {}, files = []) {
  const rows = []
  const addDomain = (plot, scope, filename, axis, domain) => {
    const clean = normalizeViewDomain(domain)
    if (!clean) return
    rows.push({
      gate_type: 'setting',
      scope,
      filename,
      x_channel: `view_${plot}`,
      y_channel: `${axis}_domain`,
      plot_mode: 'setting',
      vertex_index: 0,
      x: clean[0],
      y: clean[1],
    })
  }

  SCATTER_VIEW_PLOTS.forEach((plot) => {
    CONTROL_SCOPES.forEach((scope) => {
      const view = normalizePlotView(viewSettings?.[plot]?.[scope])
      if (!view) return
      addDomain(plot, scope, '', 'x', view.x)
      addDomain(plot, scope, '', 'y', view.y)
    })
  })

  const knownFiles = new Set(files.map((file) => String(file?.filename || '')).filter(Boolean))
  Object.entries(viewSettings?.histogram || {}).forEach(([filename, rawView]) => {
    if (!knownFiles.has(filename)) return
    const view = normalizePlotView(rawView, false)
    if (view?.x) addDomain('histogram', 'file', filename, 'x', view.x)
  })
  return rows
}

export function parseViewSettings(rows = []) {
  const views = { cell: {}, singlet: {}, histogram: {} }
  rows.forEach((row) => {
    if (String(row?.gate_type || '') !== 'setting') return
    const match = /^view_(cell|singlet|histogram)$/.exec(String(row?.x_channel || ''))
    const axisMatch = /^(x|y)_domain$/.exec(String(row?.y_channel || ''))
    if (!match || !axisMatch) return
    const plot = match[1]
    const axis = axisMatch[1]
    if (plot === 'histogram' && axis !== 'x') return
    const target = plot === 'histogram' ? String(row?.filename || '') : String(row?.scope || '')
    if (!target || (plot !== 'histogram' && !CONTROL_SCOPES.includes(target))) return
    const domain = normalizeViewDomain([row.x, row.y])
    if (!domain) return
    if (!views[plot][target]) views[plot][target] = { x: null, y: null }
    views[plot][target][axis] = domain
  })
  return views
}
