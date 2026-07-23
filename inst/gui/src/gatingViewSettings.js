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

export function buildViewSettingRows(viewSettings = {}) {
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
      const view = normalizePlotView(viewSettings?.[plot]?.[scope] || viewSettings?.[plot]?.global)
      if (!view) return
      addDomain(plot, scope, '', 'x', view.x)
      addDomain(plot, scope, '', 'y', view.y)
    })
  })

  return rows
}

export function parseViewSettings(rows = []) {
  const views = { cell: {}, singlet: {}, histogram: {} }
  const legacyViews = { cell: null, singlet: null }

  const parseRow = (row) => {
    if (String(row?.gate_type || '') !== 'setting') return
    const match = /^view_(cell|singlet)$/.exec(String(row?.x_channel || ''))
    const axisMatch = /^(x|y)_domain$/.exec(String(row?.y_channel || ''))
    if (!match || !axisMatch) return
    const plot = match[1]
    const axis = axisMatch[1]
    const scope = String(row?.scope || '')
    if (scope !== 'global' && !CONTROL_SCOPES.includes(scope)) return
    const domain = normalizeViewDomain([row.x, row.y])
    if (!domain) return
    return { plot, axis, scope, domain }
  }

  const assignDomain = (target, axis, domain) => {
    const previous = target[axis]
    target[axis] = previous
      ? [Math.min(previous[0], domain[0]), Math.max(previous[1], domain[1])]
      : domain
  }

  rows.forEach((row) => {
    const parsed = parseRow(row)
    if (!parsed) return
    const { plot, axis, scope, domain } = parsed
    if (scope === 'global') {
      if (!legacyViews[plot]) legacyViews[plot] = { x: null, y: null }
      assignDomain(legacyViews[plot], axis, domain)
      return
    }
    if (!views[plot][scope]) views[plot][scope] = { x: null, y: null }
    assignDomain(views[plot][scope], axis, domain)
  })

  SCATTER_VIEW_PLOTS.forEach((plot) => {
    const legacy = legacyViews[plot]
    if (!legacy) return
    CONTROL_SCOPES.forEach((scope) => {
      if (!views[plot][scope]) views[plot][scope] = { x: null, y: null }
      if (!views[plot][scope].x && legacy.x) views[plot][scope].x = [...legacy.x]
      if (!views[plot][scope].y && legacy.y) views[plot][scope].y = [...legacy.y]
    })
  })
  return views
}
