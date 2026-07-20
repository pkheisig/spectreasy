import { buildViewSettingRows, parseViewSettings } from '../gatingViewSettings.js'
import {
  fileUsesHistogramGates,
  fileUsesNegativeHistogramGate,
} from '../gatingEligibility.js'
import {
  MIN_CONFIRM_EVENTS,
  REQUIRED_GATE_CSV_COLUMNS,
  fileControlType,
  normalizeEventCount,
  normalizeHistogramBins,
  normalizeHistogramTransform,
  pointInPolygonValues,
} from './GatingGeometry.jsx'

function cleanVertices(vertices) {
  if (!Array.isArray(vertices)) return []
  return vertices
    .map((vertex) => ({ x: Number(vertex?.x), y: Number(vertex?.y) }))
    .filter((vertex) => Number.isFinite(vertex.x) && Number.isFinite(vertex.y))
}

function csvEscape(value) {
  const clean = value === null || value === undefined ? '' : String(value)
  return `"${clean.replace(/"/g, '""')}"`
}

function gateRowsToCsv(rows) {
  const columns = ['gate_type', 'scope', 'filename', 'x_channel', 'y_channel', 'plot_mode', 'vertex_index', 'x', 'y']
  const lines = [
    columns.map(csvEscape).join(','),
    ...rows.map((row) => columns.map((column) => csvEscape(row?.[column] ?? '')).join(',')),
  ]
  return `${lines.join('\n')}\n`
}

function parseCsvDocument(text) {
  const rows = []
  let row = []
  let field = ''
  let inQuotes = false

  for (let i = 0; i < text.length; i += 1) {
    const char = text[i]
    const next = text[i + 1]

    if (inQuotes) {
      if (char === '"' && next === '"') {
        field += '"'
        i += 1
      } else if (char === '"') {
        inQuotes = false
      } else {
        field += char
      }
      continue
    }

    if (char === '"') {
      inQuotes = true
    } else if (char === ',') {
      row.push(field)
      field = ''
    } else if (char === '\n') {
      row.push(field)
      rows.push(row)
      row = []
      field = ''
    } else if (char !== '\r') {
      field += char
    }
  }

  if (field.length || row.length) {
    row.push(field)
    rows.push(row)
  }
  if (rows.length === 0) return { headers: [], rows: [] }

  const headers = rows[0].map((header) => header.trim())
  const dataRows = rows.slice(1)
    .filter((values) => values.some((value) => String(value || '').trim().length > 0))
    .map((values) => {
      const out = {}
      headers.forEach((header, index) => {
        if (header) out[header] = values[index] ?? ''
      })
      return out
    })
  return { headers, rows: dataRows }
}

function validateGateCsvRows(rows, headers = []) {
  const missingColumns = REQUIRED_GATE_CSV_COLUMNS.filter((column) => !headers.includes(column))
  if (missingColumns.length) {
    throw new Error(`CSV is missing required columns: ${missingColumns.join(', ')}`)
  }

  const allowedGateTypes = new Set(['cell', 'singlet', 'positive', 'negative', 'setting'])
  const allowedPlotModes = new Set(['scatter', 'histogram', 'separator', 'positive_1d', 'negative_1d', 'missing', 'blocked'])

  rows.forEach((row, index) => {
    const rowNumber = index + 2
    const gateType = String(row.gate_type || '').trim()
    const plotMode = String(row.plot_mode || '').trim()
    const scope = String(row.scope || '').trim()
    const filename = String(row.filename || '').trim()

    if (!gateType) throw new Error(`CSV row ${rowNumber} has no gate_type`)
    if (!allowedGateTypes.has(gateType)) throw new Error(`CSV row ${rowNumber} has unknown gate_type "${gateType}"`)
    if (gateType === 'setting') return
    if (!plotMode) throw new Error(`CSV row ${rowNumber} has no plot_mode`)
    if (!allowedPlotModes.has(plotMode)) throw new Error(`CSV row ${rowNumber} has unknown plot_mode "${plotMode}"`)
    if (scope === 'file' && !filename) throw new Error(`CSV row ${rowNumber} is file-specific but has no filename`)
    if ((gateType === 'positive' || gateType === 'negative') && plotMode !== 'missing' && !filename) {
      throw new Error(`CSV row ${rowNumber} is a ${gateType} gate but has no filename`)
    }
    if (plotMode === 'missing' || plotMode === 'blocked') return

    const vertexIndex = Number(row.vertex_index)
    const x = Number(row.x)
    const y = Number(row.y)
    if (!Number.isInteger(vertexIndex) || vertexIndex < 1) {
      throw new Error(`CSV row ${rowNumber} has an invalid vertex_index`)
    }
    if (!Number.isFinite(x) || !Number.isFinite(y)) {
      throw new Error(`CSV row ${rowNumber} has invalid gate coordinates`)
    }
  })

}

async function saveCsvWithSystemPicker(filename, csvText) {
  const savePicker = window.showSaveFilePicker
  if (typeof savePicker === 'function') {
    const handle = await savePicker({
      suggestedName: filename,
      types: [{
        description: 'CSV files',
        accept: { 'text/csv': ['.csv'] },
      }],
    })
    const writable = await handle.createWritable()
    await writable.write(new Blob([csvText], { type: 'text/csv;charset=utf-8' }))
    await writable.close()
    return handle.name || filename
  }

  const blob = new Blob([csvText], { type: 'text/csv;charset=utf-8' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  link.remove()
  URL.revokeObjectURL(url)
  return filename
}

function gateHasValidShape(gate) {
  if (gate?.mode === 'blocked') return true
  const vertices = cleanVertices(gate?.vertices)
  if (vertices.length === 0) return false
  if (gate.mode === 'separator') return vertices.length >= 1
  if (gate.mode === 'positive_1d' || gate.mode === 'negative_1d') return vertices.length >= 2
  return vertices.length >= 3
}

function gateIsFinalized(gate) {
  return Boolean(gate) && gate.mode !== 'blocked' && gateHasValidShape(gate)
}

function resolveGateForFile(gates, type, file) {
  const filename = file?.filename || ''
  const controlType = fileControlType(file)
  const fileGate = gates[`${type}:${filename}`]
  if (fileGate?.mode === 'blocked') return null
  if (fileGate) return fileGate
  if (type === 'positive' || type === 'negative') return null
  return gates[`${type}:${controlType}`] || null
}

function spectrumGateState(gates, file) {
  const cell = resolveGateForFile(gates, 'cell', file)
  const singlet = resolveGateForFile(gates, 'singlet', file)
  const positive = resolveGateForFile(gates, 'positive', file)
  const usesHistogram = fileUsesHistogramGates(file)
  return {
    cell,
    singlet,
    positive,
    usesHistogram,
    eligible: gateIsFinalized(cell) &&
      gateIsFinalized(singlet) &&
      (!usesHistogram || gateIsFinalized(positive)),
  }
}

function spectrumCacheKeyForFile(gates, file) {
  const state = spectrumGateState(gates, file)
  if (!state.eligible) return ''
  const relevantGates = {
    cell: state.cell,
    singlet: state.singlet,
    positive: state.usesHistogram ? state.positive : null,
  }
  return `${file.filename}:${JSON.stringify(relevantGates)}`
}

function payloadEventSource(payload) {
  if (Array.isArray(payload?.events)) {
    return {
      rows: payload.events.length,
      value: (rowIndex, field) => payload.events[rowIndex]?.[field],
    }
  }
  const compact = payload?.events_compact
  const rows = Number(compact?.rows) || 0
  const fields = compact?.fields || []
  const values = compact?.values
  if (!rows || !fields.length || !(values instanceof Float32Array)) return null
  const fieldIndices = new Map(fields.map((field, index) => [field, index]))
  return {
    rows,
    value(rowIndex, field) {
      const fieldIndex = fieldIndices.get(field)
      return fieldIndex === undefined ? NaN : values[fieldIndex * rows + rowIndex]
    },
  }
}

function filterPayloadPolygon(source, inputIndices, gate, xField, yField) {
  if (!source || !gateIsFinalized(gate)) return []
  const output = []
  const count = inputIndices ? inputIndices.length : source.rows
  for (let i = 0; i < count; i++) {
    const rowIndex = inputIndices ? inputIndices[i] : i
    const x = source.value(rowIndex, xField)
    const y = source.value(rowIndex, yField)
    if (Number.isFinite(x) && Number.isFinite(y) && pointInPolygonValues(x, y, gate.vertices)) {
      output.push(rowIndex)
    }
  }
  return output
}

function countPayloadHistogram(source, inputIndices, gate, field = 'peak') {
  if (!source || !gateIsFinalized(gate)) return 0
  const vertices = gate.vertices || []
  const limits = vertices.map((point) => Number(point.x)).filter(Number.isFinite)
  if (!limits.length) return 0
  const lo = Math.min(...limits)
  const hi = Math.max(...limits)
  let count = 0
  for (const rowIndex of inputIndices || []) {
    const value = source.value(rowIndex, field)
    if (gate.mode === 'separator' ? value >= lo : value >= lo && value <= hi) count += 1
  }
  return count
}

function validateFileForConfirm(file, payload, gates) {
  const filename = file?.filename || ''
  const issues = []
  const source = payloadEventSource(payload)
  if (!filename) return issues
  if (!source?.rows) {
    issues.push({ filename, message: 'events are still loading' })
    return issues
  }

  const cellGate = resolveGateForFile(gates, 'cell', file)
  const singletGate = resolveGateForFile(gates, 'singlet', file)
  if (!gateIsFinalized(cellGate)) issues.push({ filename, message: 'cell gate unfinished' })
  if (!gateIsFinalized(singletGate)) issues.push({ filename, message: 'singlet gate unfinished' })

  const cells = filterPayloadPolygon(source, null, cellGate, cellGate?.xChannel, cellGate?.yChannel)
  const singlets = filterPayloadPolygon(source, cells, singletGate, singletGate?.xChannel, singletGate?.yChannel)

  if (!fileUsesHistogramGates(file)) {
    if (gateIsFinalized(singletGate) && singlets.length < MIN_CONFIRM_EVENTS) {
      issues.push({ filename, message: `singlet gate has only ${singlets.length.toLocaleString()} events` })
    }
    return issues
  }

  const positiveGate = resolveGateForFile(gates, 'positive', file)
  if (!gateIsFinalized(positiveGate)) {
    issues.push({ filename, message: 'positive gate unfinished' })
    return issues
  }

  const positiveCount = countPayloadHistogram(source, singlets, positiveGate)
  if (positiveCount < MIN_CONFIRM_EVENTS) {
    issues.push({ filename, message: `positive gate has only ${positiveCount.toLocaleString()} events` })
  }
  return issues
}

function formatConfirmIssues(issues) {
  const grouped = new Map()
  issues.forEach((issue) => {
    if (!grouped.has(issue.filename)) grouped.set(issue.filename, [])
    grouped.get(issue.filename).push(issue.message)
  })
  const lines = Array.from(grouped.entries()).map(([filename, messages]) => `${filename}: ${messages.join(', ')}`)
  const shown = lines.slice(0, 8)
  if (lines.length > shown.length) shown.push(`+${lines.length - shown.length} more files`)
  return shown
}

function buildRows(gates, files, settings = {}) {
  const rows = [
    {
      gate_type: 'setting',
      scope: 'global',
      filename: '',
      x_channel: 'point_size',
      y_channel: '',
      plot_mode: 'setting',
      vertex_index: 0,
      x: Number(settings.pointSize || 1.5),
      y: '',
    },
    {
      gate_type: 'setting',
      scope: 'global',
      filename: '',
      x_channel: 'max_points',
      y_channel: '',
      plot_mode: 'setting',
      vertex_index: 0,
      x: normalizeEventCount(settings.maxPoints),
      y: '',
    },
    {
      gate_type: 'setting',
      scope: 'global',
      filename: '',
      x_channel: 'histogram_bins',
      y_channel: '',
      plot_mode: 'setting',
      vertex_index: 0,
      x: normalizeHistogramBins(settings.histogramBins),
      y: '',
    },
    {
      gate_type: 'setting',
      scope: 'global',
      filename: '',
      x_channel: 'histogram_transform',
      y_channel: '',
      plot_mode: 'setting',
      vertex_index: 0,
      x: normalizeHistogramTransform(settings.histogramTransform),
      y: '',
    },
  ]
  rows.push(...buildViewSettingRows(settings.viewSettings, files))
  const negativeDisabledFiles = new Set(
    files.filter((file) => !fileUsesNegativeHistogramGate(file)).map((file) => file.filename),
  )
  Object.values(gates).forEach((gate) => {
    if (gate?.type === 'negative' && negativeDisabledFiles.has(gate.filename)) return
    if (gate?.mode === 'blocked') {
      rows.push({
        gate_type: gate.type,
        scope: 'file',
        filename: gate.filename || '',
        x_channel: gate.xChannel || '',
        y_channel: gate.yChannel || '',
        plot_mode: 'blocked',
        vertex_index: 0,
        x: '',
        y: '',
      })
      return
    }
    if (!gateHasValidShape(gate)) return
    cleanVertices(gate.vertices).forEach((vertex, index) => {
      rows.push({
        gate_type: gate.type,
        scope: gate.scope,
        filename: gate.scope === 'file' ? gate.filename : '',
        x_channel: gate.xChannel || '',
        y_channel: gate.yChannel || '',
        plot_mode: gate.mode,
        vertex_index: index + 1,
        x: vertex.x,
        y: vertex.y,
      })
    })
  })
  files.forEach((file) => {
    if (!gates[`positive:${file.filename}`] && fileUsesHistogramGates(file)) {
      rows.push({
        gate_type: 'positive',
        scope: 'file',
        filename: file.filename,
        x_channel: file.channel,
        y_channel: '',
        plot_mode: 'missing',
        vertex_index: 0,
        x: '',
        y: '',
      })
    }
  })
  return rows
}

function parseConfigRows(rows) {
  const gates = {}
  const viewSettings = parseViewSettings(rows)
  let pointSize = null
  let maxPoints = null
  let histogramBins = null
  let histogramTransform = null
  rows.forEach((row) => {
    if (row.gate_type === 'setting') {
      if (row.x_channel === 'point_size') {
        pointSize = Number(row.x) || 1.5
      } else if (row.x_channel === 'max_points') {
        maxPoints = normalizeEventCount(row.x)
      } else if (row.x_channel === 'histogram_bins') {
        histogramBins = normalizeHistogramBins(row.x)
      } else if (row.x_channel === 'histogram_transform') {
        histogramTransform = normalizeHistogramTransform(row.x)
      }
      return
    }
    if (!row.gate_type || row.plot_mode === 'missing') return
    const scope = row.scope || (row.gate_type === 'positive' ? 'file' : 'cells')
    const filename = row.filename || ''
    const key = scope === 'file' ? `${row.gate_type}:${filename}` : `${row.gate_type}:${scope}`
    if (!gates[key]) {
      gates[key] = {
        type: row.gate_type,
        scope,
        filename: filename,
        xChannel: row.x_channel,
        yChannel: row.y_channel,
        mode: row.plot_mode || 'scatter',
        vertices: [],
      }
    }
    if (row.plot_mode === 'blocked') return
    const vertexIndex = Number(row.vertex_index)
    if (vertexIndex < 1) return
    const vertex = { x: Number(row.x), y: Number(row.y), vertexIndex }
    if (Number.isFinite(vertex.x) && Number.isFinite(vertex.y)) {
      gates[key].vertices.push(vertex)
    }
  })
  Object.keys(gates).forEach((key) => {
    gates[key].vertices = gates[key].vertices
      .sort((a, b) => a.vertexIndex - b.vertexIndex)
      .map(({ x, y }) => ({ x, y }))
    if (!gateHasValidShape(gates[key])) delete gates[key]
  })
  return {
    gates,
    pointSize,
    maxPoints: maxPoints === null ? null : normalizeEventCount(maxPoints),
    histogramBins,
    histogramTransform,
    viewSettings,
  }
}

function hasViewSettings(value) {
  return ['cell', 'singlet', 'histogram'].some((plot) => Object.keys(value?.[plot] || {}).length > 0)
}

export {
  buildRows,
  formatConfirmIssues,
  gateHasValidShape,
  gateIsFinalized,
  gateRowsToCsv,
  hasViewSettings,
  parseConfigRows,
  parseCsvDocument,
  resolveGateForFile,
  saveCsvWithSystemPicker,
  spectrumCacheKeyForFile,
  spectrumGateState,
  validateFileForConfirm,
  validateGateCsvRows,
}
