import React, { useEffect, useMemo, useRef, useState } from 'react'
import {
  CheckCircle2,
  Check,
  Eraser,
  FileDown,
  FolderOpen,
  Save,
  ScatterChart,
  Upload,
  Hexagon,
  ChevronLeft,
  ChevronRight,
  Sun,
  Moon,
  Settings,
  WandSparkles,
} from 'lucide-react'
import './GatingGui.css'
import { reconcileGateCsvRows } from './gatingCsvCompatibility.js'
import {
  fileUsesHistogramGates,
  fileUsesNegativeHistogramGate,
  pruneUnavailableNegativeGates,
} from './gatingEligibility.js'
import { clearHistogramGatesForFile, histogramDomainIncludingGates } from './histogramGates.js'

const API_BASE = (() => {
  const envBase = import.meta.env.VITE_API_BASE?.trim()
  if (envBase) return envBase.replace(/\/$/, '')
  if (typeof window !== 'undefined') {
    if (window.location.port === '5174') return 'http://localhost:8000'
    return window.location.origin.replace(/\/$/, '')
  }
  return 'http://localhost:8000'
})()
const CONFIG_NAME = 'ssc_gate_config.csv'
const GUI_MODULE = 'control_gating'
const DEFAULT_EVENT_COUNT = 50000
const EVENT_COUNT_VERSION = 2
const AXIS_SETTINGS_VERSION = 2
const EVENT_COUNT_STEPS = [1000, 2000, 3000, 5000, 10000, 20000, 50000, 100000]
const DEFAULT_HISTOGRAM_BINS = 100
const HISTOGRAM_BIN_MIN = 5
const HISTOGRAM_BIN_MAX = 500
const DEFAULT_HISTOGRAM_TRANSFORM = 'asinh'
const HISTOGRAM_TRANSFORMS = [
  { value: 'asinh', label: 'Asinh' },
  { value: 'linear', label: 'Linear' },
  { value: 'log10', label: 'Log10' },
  { value: 'biexponential', label: 'Biexponential' },
]
const PRELOAD_POINTS = 100000
const MIN_CONFIRM_EVENTS = 200
const REQUIRED_GATE_CSV_COLUMNS = ['gate_type', 'scope', 'filename', 'x_channel', 'y_channel', 'plot_mode', 'vertex_index', 'x', 'y']
const PLOT_WIDTH = 520
const PLOT_HEIGHT = 420
const PAD = { left: 54, right: 18, top: 18, bottom: 46 }

function axisLabel(labels, channel) {
  const desc = labels?.[channel]
  return desc && desc !== channel ? `${channel} (${desc})` : channel
}

function channelTitle(channel) {
  return channel || ''
}

function clamp(value, lo, hi) {
  return Math.max(lo, Math.min(hi, value))
}

function normalizeHistogramBins(value) {
  const numeric = Math.round(Number(value))
  if (!Number.isFinite(numeric)) return DEFAULT_HISTOGRAM_BINS
  return clamp(numeric, HISTOGRAM_BIN_MIN, HISTOGRAM_BIN_MAX)
}

function normalizeHistogramTransform(value) {
  const clean = String(value || '').toLowerCase()
  return HISTOGRAM_TRANSFORMS.some((item) => item.value === clean) ? clean : DEFAULT_HISTOGRAM_TRANSFORM
}

function histogramTransformFns(name) {
  const transform = normalizeHistogramTransform(name)
  if (transform === 'log10') {
    return {
      forward: (value) => Math.log10(Math.max(Number(value), 1)),
      inverse: (value) => Math.pow(10, Number(value)),
    }
  }
  if (transform === 'asinh') {
    const cofactor = 150
    return {
      forward: (value) => Math.asinh(Number(value) / cofactor),
      inverse: (value) => Math.sinh(Number(value)) * cofactor,
    }
  }
  if (transform === 'biexponential') {
    const cofactor = 50
    return {
      forward: (value) => Math.sign(Number(value)) * Math.log10(1 + Math.abs(Number(value)) / cofactor),
      inverse: (value) => Math.sign(Number(value)) * cofactor * (Math.pow(10, Math.abs(Number(value))) - 1),
    }
  }
  return {
    forward: (value) => Number(value),
    inverse: (value) => Number(value),
  }
}

function histogramTransformLabel(value) {
  const selected = HISTOGRAM_TRANSFORMS.find((item) => item.value === normalizeHistogramTransform(value))
  return selected?.label || 'Linear'
}

function TransformDropdown({ value, onChange }) {
  const [open, setOpen] = useState(false)
  const rootRef = useRef(null)
  const normalizedValue = normalizeHistogramTransform(value)

  useEffect(() => {
    const closeOnOutsideClick = (event) => {
      if (rootRef.current && !rootRef.current.contains(event.target)) setOpen(false)
    }
    window.addEventListener('pointerdown', closeOnOutsideClick)
    return () => window.removeEventListener('pointerdown', closeOnOutsideClick)
  }, [])

  return (
    <div className="transform-dropdown-control">
      <span>Transform</span>
      <div className="transform-dropdown" ref={rootRef}>
        <button
          type="button"
          className="transform-select"
          onClick={(event) => {
            event.stopPropagation()
            setOpen((current) => !current)
          }}
        >
          {histogramTransformLabel(normalizedValue)}
        </button>
        {open && (
          <div className="transform-menu" role="listbox">
            {HISTOGRAM_TRANSFORMS.map((item) => {
              const selected = item.value === normalizedValue
              return (
                <button
                  key={item.value}
                  type="button"
                  className={selected ? 'selected' : ''}
                  onClick={(event) => {
                    event.stopPropagation()
                    onChange?.(item.value)
                    setOpen(false)
                  }}
                >
                  <span className="transform-check">{selected ? <Check size={14} /> : null}</span>
                  {item.label}
                </button>
              )
            })}
          </div>
        )}
      </div>
    </div>
  )
}

function extent(values) {
  const finite = values.filter((v) => Number.isFinite(v))
  if (!finite.length) return [0, 1]
  const sorted = [...finite].sort((a, b) => a - b)
  let max = sorted[Math.floor(sorted.length * 0.998)]
  if (max <= 0) {
    max = Math.max(...finite)
  }
  if (max <= 0) {
    max = 1
  }
  const pad = max * 0.04
  return [0, max + pad]
}

function getTicks(domain, count = 5) {
  const [min, max] = domain
  if (min === max) return [min]
  const step = (max - min) / (count - 1)
  const ticks = []
  for (let i = 0; i < count; i++) {
    ticks.push(min + i * step)
  }
  return ticks
}

function formatTickValue(val) {
  if (Math.abs(val) >= 1000000) {
    return (val / 1000000).toFixed(1).replace(/\.0$/, '') + 'M'
  }
  if (Math.abs(val) >= 1000) {
    return (val / 1000).toFixed(0) + 'K'
  }
  return val.toFixed(0)
}

function makeScale(domain, range) {
  const [d0, d1] = domain
  const [r0, r1] = range
  return {
    toScreen(value) {
      return r0 + ((value - d0) / (d1 - d0)) * (r1 - r0)
    },
    fromScreen(value) {
      return d0 + ((value - r0) / (r1 - r0)) * (d1 - d0)
    },
  }
}

function pointInPolygon(point, polygon) {
  if (!polygon || polygon.length < 3) return false
  let inside = false
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i].x
    const yi = polygon[i].y
    const xj = polygon[j].x
    const yj = polygon[j].y
    const intersect = yi > point.y !== yj > point.y && point.x < ((xj - xi) * (point.y - yi)) / (yj - yi) + xi
    if (intersect) inside = !inside
  }
  return inside
}

function gateKey(type, filename, controlType, forceFile = false) {
  if (type === 'positive' || forceFile) return `${type}:${filename}`
  return `${type}:${controlType || 'cells'}`
}

function histogramGateKey(type, filename) {
  return type === 'negative' ? `negative:${filename}` : `positive:${filename}`
}

function fileControlType(file) {
  const value = String(file?.['control.type'] || 'cells').toLowerCase()
  return value.startsWith('bead') ? 'beads' : 'cells'
}

function displayEventsForLimit(events, maxPoints) {
  if (!Array.isArray(events) || events.length === 0) return []
  const limit = Number(maxPoints)
  if (!Number.isFinite(limit) || limit <= 0 || events.length <= limit) return events
  const out = []
  const step = events.length / limit
  for (let i = 0; i < limit; i++) {
    out.push(events[Math.floor(i * step)])
  }
  return out
}

function payloadFilename(item) {
  if (!item) return ''
  if (typeof item.filename === 'string') return item.filename
  if (Array.isArray(item.file) && item.file[0]?.filename) return item.file[0].filename
  if (Array.isArray(item.file?.filename)) return item.file.filename[0]
  if (typeof item.file?.filename === 'string') return item.file.filename
  return ''
}

function filterPolygonEvents(events, gate, xField, yField) {
  if (!gate?.vertices?.length) return events
  const xs = gate.vertices.map((p) => p.x)
  const ys = gate.vertices.map((p) => p.y)
  const minX = Math.min(...xs)
  const maxX = Math.max(...xs)
  const minY = Math.min(...ys)
  const maxY = Math.max(...ys)
  return events.filter((e) => {
    const x = e[xField]
    const y = e[yField]
    return x >= minX && x <= maxX && y >= minY && y <= maxY && pointInPolygon({ x, y }, gate.vertices)
  })
}

function summarizeGate(events, gate, xField, yField, mode) {
  if (!gate?.vertices?.length) return { count: 0, pct: 0 }
  let count = 0
  if (gate.mode === 'separator') {
    const thresh = gate.vertices[0].x
    count = events.filter((e) => e[xField] >= thresh).length
  } else if (gate.mode === 'positive_1d' || gate.mode === 'negative_1d') {
    const x1 = gate.vertices[0].x
    const x2 = gate.vertices[1].x
    const lo = Math.min(x1, x2)
    const hi = Math.max(x1, x2)
    count = events.filter((e) => e[xField] >= lo && e[xField] <= hi).length
  } else if (mode === 'histogram') {
    const xs = gate.vertices.map((p) => p.x)
    const lo = Math.min(...xs)
    const hi = Math.max(...xs)
    count = events.filter((e) => e[xField] >= lo && e[xField] <= hi).length
  } else {
    count = filterPolygonEvents(events, gate, xField, yField).length
  }
  return { count, pct: events.length ? (count / events.length) * 100 : 0 }
}

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

function validateFileForConfirm(file, payload, gates) {
  const filename = file?.filename || ''
  const issues = []
  const events = Array.isArray(payload?.events) ? payload.events : []
  if (!filename) return issues
  if (!events.length) {
    issues.push({ filename, message: 'events are still loading' })
    return issues
  }

  const cellGate = resolveGateForFile(gates, 'cell', file)
  const singletGate = resolveGateForFile(gates, 'singlet', file)
  if (!gateIsFinalized(cellGate)) issues.push({ filename, message: 'cell gate unfinished' })
  if (!gateIsFinalized(singletGate)) issues.push({ filename, message: 'singlet gate unfinished' })

  const cells = gateIsFinalized(cellGate)
    ? filterPolygonEvents(events, cellGate, cellGate.xChannel, cellGate.yChannel)
    : []
  const singlets = gateIsFinalized(singletGate)
    ? filterPolygonEvents(cells, singletGate, singletGate.xChannel, singletGate.yChannel)
    : []

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

  const positive = summarizeGate(singlets, positiveGate, 'peak', 'count', 'histogram')
  if (positive.count < MIN_CONFIRM_EVENTS) {
    issues.push({ filename, message: `positive gate has only ${positive.count.toLocaleString()} events` })
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
  }
}

function useApi(path, options) {
  return fetch(`${API_BASE}${path}`, options).then((res) => {
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

function normalizeEventCount(value) {
  const numeric = Number(value)
  if (!Number.isFinite(numeric) || numeric <= 0) return DEFAULT_EVENT_COUNT
  const exact = EVENT_COUNT_STEPS.find((step) => step === numeric)
  if (exact) return exact
  return EVENT_COUNT_STEPS.reduce((best, step) => (Math.abs(step - numeric) < Math.abs(best - numeric) ? step : best), DEFAULT_EVENT_COUNT)
}

function stdDev(arr) {
  const n = arr.length
  if (n <= 1) return 0
  const mean = arr.reduce((a, b) => a + b, 0) / n
  const variance = arr.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / (n - 1)
  return Math.sqrt(variance)
}

function iqr(arr) {
  if (arr.length === 0) return 0
  const sorted = [...arr].sort((a, b) => a - b)
  const q1 = sorted[Math.floor(sorted.length * 0.25)]
  const q3 = sorted[Math.floor(sorted.length * 0.75)]
  return q3 - q1
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

function computeDensityBuckets(points, xField, yField) {
  const n = points.length
  const buckets = Array.from({ length: DENSITY_PALETTE.length }, () => [])
  if (n === 0) return buckets

  const xs = points.map((p) => p[xField])
  const ys = points.map((p) => p[yField])
  const minX = Math.min(...xs)
  const maxX = Math.max(...xs)
  const minY = Math.min(...ys)
  const maxY = Math.max(...ys)

  const rx = maxX - minX || 1
  const ry = maxY - minY || 1

  const numBins = 160
  const gridSide = numBins + 1
  const grid = new Uint32Array(gridSide * gridSide)

  points.forEach((p) => {
    const bx = clamp(Math.floor(((p[xField] - minX) / rx) * numBins), 0, numBins)
    const by = clamp(Math.floor(((p[yField] - minY) / ry) * numBins), 0, numBins)
    grid[by * gridSide + bx] += 1
  })

  // Precompute Gaussian weights for an 11x11 window (radius = 5)
  const sigma = 2.0
  const kernel = []
  const radius = 5
  for (let dx = -radius; dx <= radius; dx++) {
    for (let dy = -radius; dy <= radius; dy++) {
      const distSq = dx * dx + dy * dy
      const weight = Math.exp(-distSq / (2 * sigma * sigma))
      kernel.push({ dx, dy, weight })
    }
  }

  const densities = new Float32Array(n)
  let maxD = 1
  for (let pointIndex = 0; pointIndex < n; pointIndex++) {
    const p = points[pointIndex]
    const bx = clamp(Math.floor(((p[xField] - minX) / rx) * numBins), 0, numBins)
    const by = clamp(Math.floor(((p[yField] - minY) / ry) * numBins), 0, numBins)
    let sum = 0
    for (let i = 0; i < kernel.length; i++) {
      const k = kernel[i]
      const nx = bx + k.dx
      const ny = by + k.dy
      const count = nx < 0 || nx > numBins || ny < 0 || ny > numBins ? 0 : grid[ny * gridSide + nx]
      sum += count * k.weight
    }
    densities[pointIndex] = sum
    if (sum > maxD) maxD = sum
  }

  for (let i = 0; i < n; i++) {
    const bucket = clamp(Math.floor((densities[i] / maxD) * (DENSITY_PALETTE.length - 1)), 0, DENSITY_PALETTE.length - 1)
    buckets[bucket].push(i)
  }
  return buckets
}

function GatePlot({
  title,
  subtitle,
  events,
  xField,
  yField,
  xChannel,
  yChannel,
  labels,
  active,
  gate,
  draft,
  onAddPoint,
  onFinish,
  onContextGate,
  onClickPlot,
  onUpdateVertex,
  onUpdateVertices,
  onDragEnd,
  xDomain: xDomainProp,
  yDomain: yDomainProp,
  statsText,
  warningText,
  drawActive,
  pointSize,
  mode = 'scatter',
  isPositivePlot,
  histogramGateType,
  histogramBins = DEFAULT_HISTOGRAM_BINS,
  histogramTransform = DEFAULT_HISTOGRAM_TRANSFORM,
  onHistogramTransformChange,
  onHistogramBinsChange,
  secondaryGates = [],
  onSelectHistogramGate,
  onUpdateAnyGateVertices,
  availableChannels = [],
  onAxisChange,
  isBead,
  onToggleDraw,
  onClear,
  onToggleHistogramGate,
  negativeGateEnabled = true,
}) {
  const plotRef = useRef(null)
  const menuRef = useRef(null)
  const svgRef = useRef(null)
  const canvasRef = useRef(null)
  const [dragTarget, setDragTarget] = useState(null)
  const [lastMousePos, setLastMousePos] = useState(null)
  const [hoverPt, setHoverPt] = useState(null)
  const [dragPreviewVertices, setDragPreviewVertices] = useState(null)
  const [axisMenu, setAxisMenu] = useState(null)

  const transformFns = useMemo(() => histogramTransformFns(histogramTransform), [histogramTransform])
  const toPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.forward(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }
  const fromPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.inverse(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }

  const xDomain = useMemo(() => {
    const eventValues = events.map((e) => e[xField])
    const rawDomain = xDomainProp || extent(eventValues)
    if (mode !== 'histogram') return rawDomain
    const visibleDomain = histogramDomainIncludingGates(
      rawDomain,
      eventValues,
      [gate, ...secondaryGates],
    )
    const transformed = visibleDomain.map((value) => transformFns.forward(value)).filter(Number.isFinite)
    if (transformed.length < 2 || transformed[0] === transformed[1]) return [0, 1]
    return [Math.min(...transformed), Math.max(...transformed)]
  }, [events, xField, xDomainProp, mode, transformFns, gate, secondaryGates])

  const densityCurve = useMemo(() => {
    if (mode !== 'histogram') return []
    return getDensityCurve(events.map((e) => toPlotX(e[xField])), xDomain, histogramBins)
  }, [events, mode, xField, xDomain, histogramBins, histogramTransform])

  const yDomain = useMemo(() => {
    if (yDomainProp) return yDomainProp
    if (mode === 'histogram') {
      if (densityCurve.length === 0) return [0, 1]
      return [0, Math.max(...densityCurve.map((p) => p.y)) * 1.12]
    }
    return extent(events.map((e) => e[yField]))
  }, [events, mode, yDomainProp, densityCurve])

  const xScale = useMemo(() => makeScale(xDomain, [PAD.left, PLOT_WIDTH - PAD.right]), [xDomain])
  const yScale = useMemo(() => makeScale(yDomain, [PLOT_HEIGHT - PAD.bottom, PAD.top]), [yDomain])
  const xToScreen = (value) => xScale.toScreen(toPlotX(value))
  const renderGate = dragPreviewVertices && gate ? { ...gate, vertices: dragPreviewVertices } : gate
  const displayGate = draft?.length ? { vertices: draft } : renderGate
  const densityBuckets = useMemo(() => {
    if (mode !== 'scatter') return []
    return computeDensityBuckets(events, xField, yField)
  }, [events, mode, xField, yField])

  const xTicks = useMemo(() => getTicks(xDomain, 5), [xDomain])
  const yTicks = useMemo(() => getTicks(yDomain, 5), [yDomain])

  function screenToData(evt) {
    const box = svgRef.current.getBoundingClientRect()
    const scaleX = PLOT_WIDTH / box.width
    const scaleY = PLOT_HEIGHT / box.height
    const sx = clamp((evt.clientX - box.left) * scaleX, PAD.left, PLOT_WIDTH - PAD.right)
    const sy = clamp((evt.clientY - box.top) * scaleY, PAD.top, PLOT_HEIGHT - PAD.bottom)
    return {
      x: fromPlotX(xScale.fromScreen(sx)),
      y: yScale.fromScreen(sy),
      plotX: xScale.fromScreen(sx)
    }
  }

  function openAxisMenu(axis, evt) {
    if (!onAxisChange || mode === 'histogram') return
    evt.stopPropagation()
    setAxisMenu(axisMenu === axis ? null : axis)
  }

  function chooseAxis(axis, channel) {
    onAxisChange?.(axis, channel)
    setAxisMenu(null)
  }

  useEffect(() => {
    if (!axisMenu) return
    const closeAxisMenu = (event) => {
      if (menuRef.current && !menuRef.current.contains(event.target)) {
        setAxisMenu(null)
      }
    }
    const timer = setTimeout(() => {
      window.addEventListener('pointerdown', closeAxisMenu)
    }, 0)
    return () => {
      clearTimeout(timer)
      window.removeEventListener('pointerdown', closeAxisMenu)
    }
  }, [axisMenu])

  // Draw plot points onto Canvas
  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const ctx = canvas.getContext('2d')
    
    // Scale canvas dimensions for razor-sharp Retina display resolution
    const dpr = window.devicePixelRatio || 1
    canvas.width = PLOT_WIDTH * dpr
    canvas.height = PLOT_HEIGHT * dpr
    
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    if (mode === 'histogram') return

    ctx.save()
    ctx.scale(dpr, dpr)

    // Clip rendering box to lock points inside axes borders
    ctx.beginPath()
    ctx.rect(PAD.left, PAD.top, PLOT_WIDTH - PAD.left - PAD.right, PLOT_HEIGHT - PAD.top - PAD.bottom)
    ctx.clip()

    // Light theme plot bg area
    ctx.fillStyle = '#f8f7f3'
    ctx.fillRect(PAD.left, PAD.top, PLOT_WIDTH - PAD.left - PAD.right, PLOT_HEIGHT - PAD.top - PAD.bottom)
    ctx.strokeStyle = '#c7c3ba'
    ctx.strokeRect(PAD.left, PAD.top, PLOT_WIDTH - PAD.left - PAD.right, PLOT_HEIGHT - PAD.top - PAD.bottom)

    if (events.length === 0) {
      ctx.restore()
      return
    }

    const side = Math.max(1, pointSize * 1.6)
    const half = side / 2
    densityBuckets.forEach((bucket, colorIndex) => {
      if (!bucket.length) return
      ctx.fillStyle = DENSITY_PALETTE[colorIndex] || 'rgba(0, 0, 255, 0.6)'
      ctx.beginPath()
      bucket.forEach((eventIndex) => {
        const p = events[eventIndex]
        const cx = xToScreen(p[xField])
        const cy = yScale.toScreen(p[yField])
        ctx.rect(cx - half, cy - half, side, side)
      })
      ctx.fill()
    })

    ctx.restore()
  }, [events, densityBuckets, xScale, yScale, mode, pointSize, xField, yField])

  function handleClick(evt) {
    if (!active) {
      onClickPlot()
      return
    }
    if (!drawActive) return

    const clickData = screenToData(evt)
    if (isPositivePlot && mode === 'histogram') {
      if (histogramGateType === 'separator') {
        onAddPoint(clickData)
        return
      }
      if (histogramGateType === 'positive' || histogramGateType === 'negative') {
        onAddPoint(clickData)
        return
      }
    }

    if (evt.detail > 1) {
      onFinish()
      return
    }
    if (draft?.length >= 3) {
      const firstPt = draft[0]
      const screenFirstX = xToScreen(firstPt.x)
      const screenFirstY = yScale.toScreen(firstPt.y)
      
      const box = svgRef.current.getBoundingClientRect()
      const scaleX = PLOT_WIDTH / box.width
      const scaleY = PLOT_HEIGHT / box.height
      const screenClickX = (evt.clientX - box.left) * scaleX
      const screenClickY = (evt.clientY - box.top) * scaleY
      
      const dx = screenClickX - screenFirstX
      const dy = screenClickY - screenFirstY
      const dist = Math.sqrt(dx * dx + dy * dy)
      
      if (dist < 15) {
        onFinish()
        return
      }
    }
    onAddPoint(clickData)
  }

  function handleMouseDown(target, evt, targetGate = gate, targetKey = null) {
    if (drawActive) return
    evt.stopPropagation()
    evt.preventDefault()
    setDragTarget({ target, gate: targetGate, gateKey: targetKey })
    setLastMousePos(screenToData(evt))
    setDragPreviewVertices(targetGate?.vertices ? targetGate.vertices.map((v) => ({ ...v })) : null)
    if (targetKey && targetGate?.type) onSelectHistogramGate?.(targetGate.type)
  }

  function handleMouseMove(evt) {
    const currentPt = screenToData(evt)
    setHoverPt(currentPt)

    if (dragTarget === null || lastMousePos === null || !dragTarget.gate?.vertices?.length) return
    const dx = currentPt.x - lastMousePos.x
    const dy = currentPt.y - lastMousePos.y
    const dPlotX = currentPt.plotX - (lastMousePos.plotX ?? toPlotX(lastMousePos.x))
    setLastMousePos(currentPt)

    const dragMode = dragTarget.gate.mode
    const dragKind = dragTarget.target
    const is1D = dragMode === 'separator' || dragMode === 'positive_1d' || dragMode === 'negative_1d'
    const sourceVertices = dragPreviewVertices || dragTarget.gate.vertices

    if (dragKind === 'gate' || dragKind === 'separator') {
      const newVertices = sourceVertices.map((v) => {
        if (is1D) {
          const newPlotX = toPlotX(v.x) + dPlotX
          return { x: fromPlotX(newPlotX), y: 0 }
        }
        return {
          x: v.x + dx,
          y: v.y + dy,
        }
      })
      setDragPreviewVertices(newVertices)
    } else if (typeof dragKind === 'number') {
      const newVertices = sourceVertices.map((v) => ({ ...v }))
      let newPt
      if (is1D) {
        const newPlotX = toPlotX(sourceVertices[dragKind].x) + dPlotX
        newPt = { x: fromPlotX(newPlotX), y: 0 }
      } else {
        newPt = {
          x: sourceVertices[dragKind].x + dx,
          y: sourceVertices[dragKind].y + dy,
        }
      }
      newVertices[dragKind] = newPt
      setDragPreviewVertices(newVertices)
    }
  }

  function handleMouseUp() {
    if (dragTarget !== null) {
      if (dragPreviewVertices) {
        if (dragTarget.gateKey) {
          onUpdateAnyGateVertices?.(dragTarget.gateKey, dragPreviewVertices)
        } else {
          onUpdateVertices(dragPreviewVertices)
        }
      }
      setDragTarget(null)
      setLastMousePos(null)
      setDragPreviewVertices(null)
      onDragEnd()
    }
  }

  function handleContext(evt) {
    evt.preventDefault()
    const data = screenToData(evt)
    const activeContextGate = renderGate || gate
    if (!activeContextGate?.vertices?.length) return
    if (activeContextGate.mode === 'separator') {
      const x = xToScreen(activeContextGate.vertices[0].x)
      const box = svgRef.current.getBoundingClientRect()
      const scaleX = PLOT_WIDTH / box.width
      const clickX = (evt.clientX - box.left) * scaleX
      if (Math.abs(clickX - x) < 12) onContextGate(evt.clientX, evt.clientY)
      return
    }
    if ((activeContextGate.mode === 'positive_1d' || activeContextGate.mode === 'negative_1d') && activeContextGate.vertices.length >= 2) {
      const lo = Math.min(activeContextGate.vertices[0].x, activeContextGate.vertices[1].x)
      const hi = Math.max(activeContextGate.vertices[0].x, activeContextGate.vertices[1].x)
      if (data.x >= lo && data.x <= hi) onContextGate(evt.clientX, evt.clientY)
      return
    }
    if (pointInPolygon(data, activeContextGate.vertices)) {
      onContextGate(evt.clientX, evt.clientY)
    }
  }

  const path = displayGate?.vertices?.length
    ? displayGate.vertices.map((p, i) => `${i === 0 ? 'M' : 'L'}${xToScreen(p.x)},${yScale.toScreen(p.y)}`).join(' ') + (displayGate.vertices.length > 2 && !draft?.length ? ' Z' : '')
    : ''

  const densityPath = useMemo(() => {
    if (mode !== 'histogram' || densityCurve.length === 0) return ''
    const pointsStr = densityCurve.map((p) => `${xScale.toScreen(p.x)},${yScale.toScreen(p.y)}`)
    const firstX = xScale.toScreen(densityCurve[0].x)
    const lastX = xScale.toScreen(densityCurve[densityCurve.length - 1].x)
    const baseline = PLOT_HEIGHT - PAD.bottom
    return `M${firstX},${baseline} L${pointsStr.join(' L')} L${lastX},${baseline} Z`
  }, [densityCurve, xScale, yScale, mode])

  const gateClass = isBead && !isPositivePlot
    ? 'gate-path-bead'
    : 'gate-path'

  function renderSecondaryGate(extraGate, index) {
    if (!extraGate?.vertices?.length) return null
    if (extraGate.mode === 'separator') {
      return (
        <line
          key={`secondary-separator-${index}`}
          x1={xToScreen(extraGate.vertices[0].x)}
          y1={PAD.top}
          x2={xToScreen(extraGate.vertices[0].x)}
          y2={PLOT_HEIGHT - PAD.bottom}
          className="separator-line is-secondary"
          clipPath="url(#plot-clip)"
        />
      )
    }
    if ((extraGate.mode === 'positive_1d' || extraGate.mode === 'negative_1d') && extraGate.vertices.length >= 2) {
      const gateKeyForDrag = extraGate._gateKey || histogramGateKey(extraGate.type, extraGate.filename || '')
      return (
        <g key={`secondary-interval-${index}`}>
        <rect
          x={Math.min(xToScreen(extraGate.vertices[0].x), xToScreen(extraGate.vertices[1].x))}
          y={PAD.top}
          width={Math.abs(xToScreen(extraGate.vertices[1].x) - xToScreen(extraGate.vertices[0].x))}
          height={PLOT_HEIGHT - PAD.bottom - PAD.top}
          className={`${extraGate.mode === 'positive_1d' ? 'gate-path-pos' : 'gate-path-neg'} is-secondary`}
          clipPath="url(#plot-clip)"
          style={{ cursor: 'move' }}
          onMouseDown={(evt) => handleMouseDown('gate', evt, extraGate, gateKeyForDrag)}
          onContextMenu={(evt) => {
            evt.preventDefault()
            evt.stopPropagation()
            onSelectHistogramGate?.(extraGate.type)
            onContextGate(evt.clientX, evt.clientY, gateKeyForDrag)
          }}
        />
        {extraGate.vertices.map((p, handleIndex) => (
          <circle
            key={`secondary-handle-${index}-${handleIndex}`}
            cx={xToScreen(p.x)}
            cy={PLOT_HEIGHT / 2}
            r="7"
            className="gate-handle is-secondary"
            clipPath="url(#plot-clip)"
            style={{ cursor: 'ew-resize' }}
            onMouseDown={(evt) => handleMouseDown(handleIndex, evt, extraGate, gateKeyForDrag)}
            onContextMenu={(evt) => {
              evt.preventDefault()
              evt.stopPropagation()
              onSelectHistogramGate?.(extraGate.type)
              onContextGate(evt.clientX, evt.clientY, gateKeyForDrag)
            }}
          />
        ))}
        </g>
      )
    }
    return null
  }

  return (
    <section className={`plot-panel ${active ? 'is-active' : ''}`} onClick={() => { if (!active) onClickPlot() }}>
      <header className="plot-head">
        <div>
          <h2>{title}</h2>
          {subtitle ? <p>{subtitle}</p> : null}
        </div>
        <div style={{ textAlign: 'right' }}>
          <span style={{ fontSize: 12, color: 'var(--muted)', fontWeight: 'bold' }}>
            {events.length.toLocaleString()} shown
          </span>
          {statsText && (
            <div style={{ fontSize: 12, fontWeight: 800, color: 'var(--ink)', marginTop: 4, whiteSpace: 'nowrap' }}>
              {statsText}
            </div>
          )}
          {warningText && (
            <div className="plot-warning">
              {warningText}
            </div>
          )}
        </div>
      </header>

      {/* Toolbar buttons inside plot panel headers */}
      <div className="plot-toolbar">
        {!isPositivePlot ? (
          <>
            <button
              className={active && drawActive ? 'on' : ''}
              onClick={onToggleDraw}
            >
              <Hexagon size={14} /> Gate
            </button>
          </>
        ) : (
          <>
            <button
              className={`btn-negative-gate ${active && drawActive && histogramGateType === 'negative' ? 'on' : ''}`}
              disabled={!negativeGateEnabled}
              title={negativeGateEnabled ? 'Draw negative histogram gate' : 'A mapped external negative supplies the background'}
              onClick={(e) => onToggleHistogramGate('negative', e)}
            >
              Neg
            </button>
            <button
              className={`btn-positive-gate ${active && drawActive && histogramGateType === 'positive' ? 'on' : ''}`}
              onClick={(e) => onToggleHistogramGate('positive', e)}
            >
              Pos
            </button>
            <button
              className="btn-clear-histogram"
              onClick={onClear}
              title="Remove the positive and negative histogram gates for this file"
            >
              <Eraser size={13} /> Clear
            </button>
            <label className="histogram-bin-control">
              <span>Bins</span>
              <input
                type="range"
                min={HISTOGRAM_BIN_MIN}
                max={HISTOGRAM_BIN_MAX}
                step="5"
                value={histogramBins}
                onChange={(e) => onHistogramBinsChange?.(Number(e.target.value))}
                onInput={(e) => onHistogramBinsChange?.(Number(e.currentTarget.value))}
              />
              <strong>{histogramBins}</strong>
            </label>
            <div className="histogram-transform-control">
              <TransformDropdown
                value={histogramTransform}
                onChange={(nextValue) => onHistogramTransformChange?.(nextValue)}
              />
            </div>
          </>
        )}
      </div>

      <div ref={plotRef} style={{ position: 'relative', width: '100%', aspectRatio: `${PLOT_WIDTH} / ${PLOT_HEIGHT}` }}>
        {mode === 'scatter' && (
          <canvas
            ref={canvasRef}
            style={{
              position: 'absolute',
              top: 0,
              left: 0,
              width: '100%',
              height: '100%',
              pointerEvents: 'none',
              zIndex: 1,
            }}
          />
        )}
        <svg
          ref={svgRef}
          width="100%"
          viewBox={`0 0 ${PLOT_WIDTH} ${PLOT_HEIGHT}`}
          className="plot-svg"
          style={{
            position: 'absolute',
            top: 0,
            left: 0,
            width: '100%',
            height: '100%',
            zIndex: 2,
            background: 'transparent',
            cursor: drawActive ? 'crosshair' : 'default',
          }}
          onClick={handleClick}
          onDoubleClick={onFinish}
          onContextMenu={handleContext}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseUp}
        >
          <defs>
            <clipPath id="plot-clip">
              <rect x={PAD.left} y={PAD.top} width={PLOT_WIDTH - PAD.left - PAD.right} height={PLOT_HEIGHT - PAD.top - PAD.bottom} />
            </clipPath>
          </defs>
          
          {mode === 'histogram' && (
            <>
              {/* Keep plot areas white/light themed */}
              <rect x={PAD.left} y={PAD.top} width={PLOT_WIDTH - PAD.left - PAD.right} height={PLOT_HEIGHT - PAD.top - PAD.bottom} fill="#f8f7f3" stroke="#c7c3ba" />
              {densityPath && <path d={densityPath} fill="rgba(38, 63, 115, 0.3)" stroke="#263f73" strokeWidth="1.5" clipPath="url(#plot-clip)" />}
            </>
          )}
          <line x1={PAD.left} y1={PLOT_HEIGHT - PAD.bottom} x2={PLOT_WIDTH - PAD.right} y2={PLOT_HEIGHT - PAD.bottom} className="axis" />
          <line x1={PAD.left} y1={PAD.top} x2={PAD.left} y2={PLOT_HEIGHT - PAD.bottom} className="axis" />
          <text
            x={PLOT_WIDTH / 2}
            y={PLOT_HEIGHT - 10}
            className={`axis-label ${mode !== 'histogram' ? 'axis-label-clickable' : ''}`}
            onClick={(evt) => openAxisMenu('x', evt)}
          >
            {axisLabel(labels, xChannel)}
          </text>
          <text
            x="16"
            y={PLOT_HEIGHT / 2}
            transform={`rotate(-90 16 ${PLOT_HEIGHT / 2})`}
            className={`axis-label ${mode !== 'histogram' ? 'axis-label-clickable' : ''}`}
            onClick={(evt) => openAxisMenu('y', evt)}
          >
            {mode === 'histogram' ? 'Events per bin' : axisLabel(labels, yChannel)}
          </text>
          
          {xTicks.map((t) => {
            const x = xScale.toScreen(t)
            return (
              <g key={`xtick-${t}`}>
                <line x1={x} y1={PLOT_HEIGHT - PAD.bottom} x2={x} y2={PLOT_HEIGHT - PAD.bottom + 5} stroke="#6d756f" strokeWidth="1" />
                <text x={x} y={PLOT_HEIGHT - PAD.bottom + 18} textAnchor="middle" fontSize="10" fill="#6d756f" fontWeight="bold">{formatTickValue(mode === 'histogram' ? fromPlotX(t) : t)}</text>
              </g>
            )
          })}
          {yTicks.map((t) => {
            const y = yScale.toScreen(t)
            return (
              <g key={`ytick-${t}`}>
                <line x1={PAD.left - 5} y1={y} x2={PAD.left} y2={y} stroke="#6d756f" strokeWidth="1" />
                <text x={PAD.left - 8} y={y + 3} textAnchor="end" fontSize="10" fill="#6d756f" fontWeight="bold">{formatTickValue(t)}</text>
              </g>
            )
          })}

          {secondaryGates.map(renderSecondaryGate)}

          {/* Separator Line representation */}
          {displayGate?.mode === 'separator' && (
            <line
              x1={xToScreen(displayGate.vertices[0].x)}
              y1={PAD.top}
              x2={xToScreen(displayGate.vertices[0].x)}
              y2={PLOT_HEIGHT - PAD.bottom}
              className="separator-line"
              clipPath="url(#plot-clip)"
              style={{ cursor: 'ew-resize' }}
              onMouseDown={(evt) => handleMouseDown('gate', evt)}
            />
          )}

          {/* 1D intervals: positive or negative */}
          {(displayGate?.mode === 'positive_1d' || displayGate?.mode === 'negative_1d') && (
            <rect
              x={Math.min(xToScreen(displayGate.vertices[0].x), xToScreen(displayGate.vertices[1].x))}
              y={PAD.top}
              width={Math.abs(xToScreen(displayGate.vertices[1].x) - xToScreen(displayGate.vertices[0].x))}
              height={PLOT_HEIGHT - PAD.bottom - PAD.top}
              className={displayGate.mode === 'positive_1d' ? 'gate-path-pos' : 'gate-path-neg'}
              clipPath="url(#plot-clip)"
              style={{ cursor: 'move' }}
              onMouseDown={(evt) => handleMouseDown('gate', evt)}
            />
          )}

          {/* Render boundary handles for 1D gates */}
          {(displayGate?.mode === 'positive_1d' || displayGate?.mode === 'negative_1d') && displayGate.vertices.map((p, index) => (
            <circle
              key={`handle-${index}`}
              cx={xToScreen(p.x)}
              cy={PLOT_HEIGHT / 2}
              r="8"
              className="gate-handle"
              clipPath="url(#plot-clip)"
              style={{ cursor: 'ew-resize' }}
              onMouseDown={(evt) => handleMouseDown(index, evt)}
            />
          ))}

          {/* 1D Gate drawing preview feedback */}
          {draft?.length === 1 && hoverPt && (mode === 'histogram') && (
            <rect
              x={Math.min(xToScreen(draft[0].x), xToScreen(hoverPt.x))}
              y={PAD.top}
              width={Math.abs(xToScreen(hoverPt.x) - xToScreen(draft[0].x))}
              height={PLOT_HEIGHT - PAD.bottom - PAD.top}
              fill={histogramGateType === 'positive' ? 'rgba(214, 82, 56, 0.1)' : 'rgba(38, 63, 115, 0.1)'}
              stroke={histogramGateType === 'positive' ? '#d65238' : '#263f73'}
              strokeWidth="2"
              strokeDasharray="4 4"
              clipPath="url(#plot-clip)"
            />
          )}
          
          {path && !['separator', 'positive_1d', 'negative_1d'].includes(displayGate?.mode) && (
            <path
              d={path}
              className={`${gateClass} ${draft?.length ? 'is-draft' : ''}`}
              clipPath="url(#plot-clip)"
              style={{ cursor: (!drawActive && gate?.vertices?.length) ? 'pointer' : 'default' }}
              onMouseDown={(evt) => handleMouseDown('gate', evt)}
            />
          )}
          {displayGate?.vertices?.length && !['separator', 'positive_1d', 'negative_1d'].includes(displayGate?.mode) && displayGate.vertices.map((p, index) => (
            <circle
              key={`${p.x}-${p.y}-${index}`}
              cx={xToScreen(p.x)}
              cy={yScale.toScreen(p.y)}
              r="8"
              className="gate-handle"
              clipPath="url(#plot-clip)"
              style={{ cursor: 'move' }}
              onMouseDown={(evt) => handleMouseDown(index, evt)}
            />
          ))}
        </svg>
        {axisMenu && (
          <div ref={menuRef} className={`axis-menu axis-menu-${axisMenu}`} onClick={(evt) => evt.stopPropagation()}>
            {availableChannels.map((channel) => (
              <button
                key={`${axisMenu}-${channel}`}
                className={channel === (axisMenu === 'x' ? xChannel : yChannel) ? 'on' : ''}
                onClick={() => chooseAxis(axisMenu, channel)}
              >
                {channelTitle(channel)}
              </button>
            ))}
          </div>
        )}
      </div>
    </section>
  )
}

function App() {
  const [status, setStatus] = useState('Loading controls')
  const [files, setFiles] = useState([])
  const [metadata, setMetadata] = useState({})
  const [selected, setSelected] = useState('')
  const [payload, setPayload] = useState(null)
  const [payloadCache, setPayloadCache] = useState({})
  const [gates, setGates] = useState({})
  const [activeGate, setActiveGate] = useState('cell')
  const [draft, setDraft] = useState([])
  const [histogramGateType, setHistogramGateType] = useState('positive')
  const [contextMenu, setContextMenu] = useState(null)
  const [configs, setConfigs] = useState([])
  const [showConfirmModal, setShowConfirmModal] = useState(false)
  const [showAutogateConfirmModal, setShowAutogateConfirmModal] = useState(false)
  const [showConfirmIssues, setShowConfirmIssues] = useState(false)
  const [showSettingsModal, setShowSettingsModal] = useState(false)
  const [drawActive, setDrawActive] = useState(false)
  const [pointSize, setPointSize] = useState(1.5)
  const [maxPoints, setMaxPoints] = useState(DEFAULT_EVENT_COUNT)
  const [histogramBins, setHistogramBins] = useState(DEFAULT_HISTOGRAM_BINS)
  const [histogramTransform, setHistogramTransform] = useState(DEFAULT_HISTOGRAM_TRANSFORM)
  const [darkMode, setDarkMode] = useState(false)
  const [guiStateLoaded, setGuiStateLoaded] = useState(false)
  const [preloadComplete, setPreloadComplete] = useState(false)
  const [gatesLoaded, setGatesLoaded] = useState(false)
  const [histogramAutogating, setHistogramAutogating] = useState(false)
  const [axisSettings, setAxisSettings] = useState({ cell: {}, singlet: {} })
  const [spectrum, setSpectrum] = useState(null)
  const [spectrumCache, setSpectrumCache] = useState({})
  const loadInputRef = useRef(null)

  // Sync dark mode class to document body
  useEffect(() => {
    document.body.classList.toggle('dark', darkMode)
  }, [darkMode])

  useEffect(() => {
    const handleKeyDown = (e) => {
      if (e.key === 'Escape') {
        setDraft([])
        setDrawActive(false)
      }
    }
    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [])

  useEffect(() => {
    const closeMenu = () => setContextMenu(null)
    window.addEventListener('click', closeMenu)
    return () => window.removeEventListener('click', closeMenu)
  }, [])

  // Initial load config, files, persisted GUI settings, and backend cache
  useEffect(() => {
    Promise.all([
      useApi('/gate_files'),
      useApi('/gate_configs'),
      useApi('/gate_cache'),
      useApi(`/gui_state?module=${encodeURIComponent(GUI_MODULE)}`),
      useApi('/gate_config')
    ])
      .then(([fileData, configData, cacheData, guiState, gateConfig]) => {
        const clean = fileData.files.filter((f) => f.file_exists)
        setFiles(clean)
        setMetadata(fileData.metadata || {})
        setSelected(clean[0]?.filename || '')
        setConfigs(configData.configs || [])

        const compatibleCsv = reconcileGateCsvRows(gateConfig?.rows || [], clean)
        const csvGates = parseConfigRows(compatibleCsv.rows)
        if (Object.keys(csvGates.gates).length > 0) {
          setGates(pruneUnavailableNegativeGates(csvGates.gates, clean))
        } else if (cacheData?.gates && Object.keys(cacheData.gates).length > 0) {
          setGates(pruneUnavailableNegativeGates(cacheData.gates, clean))
        }
        if (typeof csvGates.pointSize === 'number') setPointSize(csvGates.pointSize)
        if (typeof csvGates.maxPoints === 'number') setMaxPoints(normalizeEventCount(csvGates.maxPoints))
        if (typeof csvGates.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(csvGates.histogramBins))
        if (typeof csvGates.histogramTransform === 'string') setHistogramTransform(normalizeHistogramTransform(csvGates.histogramTransform))
        if (cacheData?.pointSize) {
          setPointSize(Number(cacheData.pointSize))
        }
        if (cacheData?.maxPoints && cacheData?.eventCountVersion === EVENT_COUNT_VERSION) {
          setMaxPoints(normalizeEventCount(cacheData.maxPoints))
        }
        if (typeof cacheData?.histogramBins === 'number') {
          setHistogramBins(normalizeHistogramBins(cacheData.histogramBins))
        }
        if (typeof cacheData?.histogramTransform === 'string') {
          setHistogramTransform(normalizeHistogramTransform(cacheData.histogramTransform))
        }
        const persisted = guiState?.config || {}
        if (typeof persisted.pointSize === 'number') setPointSize(persisted.pointSize)
        if (typeof persisted.maxPoints === 'number' && persisted.eventCountVersion === EVENT_COUNT_VERSION) setMaxPoints(normalizeEventCount(persisted.maxPoints))
        if (typeof persisted.darkMode === 'boolean') setDarkMode(persisted.darkMode)
        if (typeof persisted.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(persisted.histogramBins))
        if (typeof persisted.histogramTransform === 'string') setHistogramTransform(normalizeHistogramTransform(persisted.histogramTransform))
        if (persisted.axisSettings && typeof persisted.axisSettings === 'object' && persisted.axisSettingsVersion === AXIS_SETTINGS_VERSION) {
          setAxisSettings(persisted.axisSettings)
        }
        setGatesLoaded(true)
        setGuiStateLoaded(true)
        setStatus('Ready')
      })
      .catch((err) => setStatus(`Could not load controls: ${err.message}`))
  }, [])

  useEffect(() => {
    if (!guiStateLoaded) return
    const timer = setTimeout(() => {
      useApi('/gui_state', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          module: GUI_MODULE,
          config_json: { pointSize, maxPoints: normalizeEventCount(maxPoints), eventCountVersion: EVENT_COUNT_VERSION, histogramBins, histogramTransform, darkMode, axisSettings, axisSettingsVersion: AXIS_SETTINGS_VERSION }
        })
      }).catch(() => {})
    }, 350)
    return () => clearTimeout(timer)
  }, [pointSize, maxPoints, histogramBins, histogramTransform, darkMode, axisSettings, guiStateLoaded])

  // Synchronize frontend gates and settings to backend in-memory cache
  useEffect(() => {
    if (!gatesLoaded) return
    const timer = setTimeout(() => {
      useApi('/gate_cache', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ gates, pointSize, maxPoints: normalizeEventCount(maxPoints), histogramBins, histogramTransform, eventCountVersion: EVENT_COUNT_VERSION })
      }).catch(() => {})
    }, 400)
    return () => clearTimeout(timer)
  }, [gates, pointSize, maxPoints, histogramBins, histogramTransform, gatesLoaded])

  useEffect(() => {
    if (!gatesLoaded || !files.length) return
    const timer = setTimeout(() => {
      useApi('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows: buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform }) }),
      }).catch(() => {})
    }, 700)
    return () => clearTimeout(timer)
  }, [gates, files, pointSize, maxPoints, histogramBins, histogramTransform, gatesLoaded])

  useEffect(() => {
    if (!selected || !gatesLoaded) return
    const selectedFile = files.find((file) => file.filename === selected) || {}
    const selectedControlType = fileControlType(selectedFile)
    const relevantGates = {
      cell: gates[`cell:${selected}`] || gates[`cell:${selectedControlType}`] || null,
      singlet: gates[`singlet:${selected}`] || gates[`singlet:${selectedControlType}`] || null,
      positive: gates[`positive:${selected}`] || null,
    }
    const cacheKey = `${selected}:${JSON.stringify(relevantGates)}`
    if (spectrumCache[cacheKey]) {
      setSpectrum(spectrumCache[cacheKey])
      return
    }
    setSpectrum(null)
    const timer = setTimeout(() => {
      useApi(`/gate_spectrum?filename=${encodeURIComponent(selected)}`)
        .then((data) => {
          const image = data?.spectrum || null
          setSpectrum(image)
          if (image) {
            setSpectrumCache((prev) => ({ ...prev, [cacheKey]: image }))
          }
        })
        .catch(() => setSpectrum(null))
    }, 800)
    return () => clearTimeout(timer)
  }, [selected, gates, gatesLoaded, files, spectrumCache])

  useEffect(() => {
    if (!files.length) return
    setStatus('Preloading controls')
    setPreloadComplete(false)
    useApi(`/gate_preload?max_points=${PRELOAD_POINTS}`)
      .then((data) => {
        const next = {}
        ;(data.payloads || []).forEach((item) => {
          const filename = payloadFilename(item)
          if (filename && !item.error) next[filename] = item
        })
        setPayloadCache(next)
        const active = selected || files[0]?.filename
        if (active && next[active]) setPayload(next[active])
        setDraft([])
        setDrawActive(false)
        setPreloadComplete(true)
        setStatus('Ready')
      })
      .catch((err) => {
        setPreloadComplete(true)
        setStatus(`Could not preload events: ${err.message}`)
      })
  }, [files])

  useEffect(() => {
    if (!selected) return
    if (payloadCache[selected]) {
      setPayload(payloadCache[selected])
      setDraft([])
      setDrawActive(false)
      setStatus('Ready')
      return
    }
    if (!preloadComplete) return
    setStatus(`Loading ${selected}`)
    useApi(`/gate_events?filename=${encodeURIComponent(selected)}&max_points=${PRELOAD_POINTS}`)
      .then((data) => {
        setPayloadCache((prev) => ({ ...prev, [selected]: data }))
        setPayload(data)
        setDraft([])
        setDrawActive(false)
        setStatus('Ready')
      })
      .catch((err) => setStatus(`Could not load events: ${err.message}`))
  }, [selected, payloadCache, preloadComplete])

  const currentFile = payload?.file?.[0] || files.find((f) => f.filename === selected) || {}
  const events = useMemo(() => displayEventsForLimit(payload?.events || [], maxPoints), [payload, maxPoints])
  const labels = payload?.labels || metadata.labels || {}
  const channels = payload?.channels || {}
  const domains = payload?.domains || {}
  const availableScatterChannels = useMemo(() => {
    const fromPayloads = new Set()
    Object.values(payloadCache).forEach((item) => {
      ;(item?.channels?.scatter || []).forEach((channel) => fromPayloads.add(channel))
    })
    ;(channels.scatter || []).forEach((channel) => fromPayloads.add(channel))
    return Array.from(fromPayloads)
  }, [payloadCache, channels])
  const globalScatterDomains = useMemo(() => {
    const out = {}
    Object.values(payloadCache).forEach((item) => {
      const scatterDomains = item?.domains?.scatter || {}
      Object.entries(scatterDomains).forEach(([channel, domain]) => {
        if (!Array.isArray(domain) || domain.length < 2) return
        if (!out[channel]) {
          out[channel] = [Number(domain[0]), Number(domain[1])]
        } else {
          out[channel][0] = Math.min(out[channel][0], Number(domain[0]))
          out[channel][1] = Math.max(out[channel][1], Number(domain[1]))
        }
      })
    })
    return out
  }, [payloadCache])
  const domainForChannel = (channel, fallback) => {
    if (channel && globalScatterDomains[channel]) return globalScatterDomains[channel]
    if (channel && domains.scatter?.[channel]) return domains.scatter[channel]
    return fallback || [0, 1]
  }

  const controlType = fileControlType(currentFile)
  const isBead = controlType === 'beads'
  const activeKey = gateKey(activeGate, selected, controlType)

  const fileCellGate = gates[`cell:${selected}`]
  const fileSingletGate = gates[`singlet:${selected}`]
  const cellGate = fileCellGate?.mode === 'blocked' ? null : (fileCellGate || gates[`cell:${controlType}`])
  const singletGate = fileSingletGate?.mode === 'blocked' ? null : (fileSingletGate || gates[`singlet:${controlType}`])
  const positiveGate = gates[`positive:${selected}`]
  const usesNegativeHistogramGate = fileUsesNegativeHistogramGate(currentFile)
  const negativeGate = usesNegativeHistogramGate ? gates[`negative:${selected}`] : null
  const activeHistogramGate = histogramGateType === 'negative' ? negativeGate : positiveGate
  const secondaryHistogramGates = [
    histogramGateType === 'negative'
      ? (positiveGate ? { ...positiveGate, _gateKey: histogramGateKey('positive', selected) } : null)
      : (negativeGate ? { ...negativeGate, _gateKey: histogramGateKey('negative', selected) } : null)
  ].filter(Boolean)
  useEffect(() => {
    if (!usesNegativeHistogramGate && histogramGateType === 'negative') {
      setHistogramGateType('positive')
      setDraft([])
      setDrawActive(false)
    }
  }, [usesNegativeHistogramGate, histogramGateType])
  const selectedIndex = files.findIndex((file) => file.filename === selected)
  const canGoPrevious = selectedIndex > 0
  const canGoNext = selectedIndex >= 0 && selectedIndex < files.length - 1
  const cellAxes = {
    x: axisSettings.cell?.x || cellGate?.xChannel || channels.fsc_a,
    y: axisSettings.cell?.y || cellGate?.yChannel || channels.ssc_a,
  }
  const singletAxes = {
    x: axisSettings.singlet?.x || singletGate?.xChannel || channels.fsc_h,
    y: axisSettings.singlet?.y || singletGate?.yChannel || channels.fsc_a,
  }

  function updateAxis(plot, axis, channel) {
    setAxisSettings((prev) => ({
      ...prev,
      [plot]: {
        ...(prev[plot] || {}),
        [axis]: channel,
      },
    }))
    setDraft([])
    setDrawActive(false)
  }

  function selectRelativeFile(offset) {
    if (selectedIndex < 0) return
    const next = files[selectedIndex + offset]
    if (next) setSelected(next.filename)
  }

  const cellsFilteredEvents = useMemo(() => {
    if (!cellGate?.vertices?.length) return events
    return filterPolygonEvents(events, cellGate, cellGate.xChannel || cellAxes.x, cellGate.yChannel || cellAxes.y)
  }, [events, cellGate])

  const singletsFilteredEvents = useMemo(() => {
    if (!singletGate?.vertices?.length) return cellsFilteredEvents
    return filterPolygonEvents(cellsFilteredEvents, singletGate, singletGate.xChannel || singletAxes.x, singletGate.yChannel || singletAxes.y)
  }, [cellsFilteredEvents, singletGate])

  const cellSummary = summarizeGate(events, cellGate, cellGate?.xChannel || cellAxes.x, cellGate?.yChannel || cellAxes.y, 'scatter')
  const singletSummary = summarizeGate(cellsFilteredEvents, singletGate, singletGate?.xChannel || singletAxes.x, singletGate?.yChannel || singletAxes.y, 'scatter')
  const positiveSummary = summarizeGate(singletsFilteredEvents, positiveGate, 'peak', 'count', 'histogram')
  const negativeSummary = summarizeGate(singletsFilteredEvents, negativeGate, 'peak', 'count', 'histogram')
  const confirmIssues = useMemo(() => (
    files.flatMap((file) => validateFileForConfirm(file, payloadCache[file.filename], gates))
  ), [files, payloadCache, gates])
  const histogramAutogateTargets = useMemo(() => files.filter(fileUsesHistogramGates), [files])
  const histogramAutogateNegativeTargets = useMemo(() => (
    histogramAutogateTargets.filter(fileUsesNegativeHistogramGate)
  ), [histogramAutogateTargets])
  const histogramAutogateMissing = useMemo(() => (
    histogramAutogateTargets.filter((file) => (
      !gateIsFinalized(resolveGateForFile(gates, 'cell', file)) ||
      !gateIsFinalized(resolveGateForFile(gates, 'singlet', file))
    ))
  ), [histogramAutogateTargets, gates])
  const canAutogateHistograms = gatesLoaded && histogramAutogateTargets.length > 0 && histogramAutogateMissing.length === 0 && !histogramAutogating
  const canConfirm = gatesLoaded && files.length > 0 && confirmIssues.length === 0
  const confirmIssueLines = useMemo(() => formatConfirmIssues(confirmIssues), [confirmIssues])
  const initialLoading = !status.startsWith('Could not') && (!gatesLoaded || (files.length > 0 && (!preloadComplete || !payload)))
  const singletWarningText = !fileUsesHistogramGates(currentFile) && gateIsFinalized(singletGate) && singletsFilteredEvents.length < MIN_CONFIRM_EVENTS
    ? `Only ${singletsFilteredEvents.length.toLocaleString()} events in singlets`
    : ''
  const positiveWarningText = fileUsesHistogramGates(currentFile) && gateIsFinalized(positiveGate) && positiveSummary.count < MIN_CONFIRM_EVENTS
    ? `Only ${positiveSummary.count.toLocaleString()} events in Pos gate`
    : ''

  function handleDragEnd() {
    setStatus('Ready')
  }

  function finishDraft() {
    if (draft.length < 3) return
    const scope = activeGate === 'positive' ? 'file' : controlType
    const filename = activeGate === 'positive' ? selected : ''
    const yChannel = activeGate === 'cell' ? cellAxes.y : activeGate === 'singlet' ? singletAxes.y : ''
    const xChannel = activeGate === 'cell' ? cellAxes.x : activeGate === 'singlet' ? singletAxes.x : channels.peak
    setGates((prev) => ({
      ...prev,
      [activeKey]: {
        type: activeGate,
        scope,
        filename,
        xChannel,
        yChannel,
        mode: activeGate === 'positive' ? 'histogram' : 'scatter',
        vertices: draft,
      },
    }))
    setDraft([])
    setDrawActive(false)
  }

  function updateGateVertex(key, index, newPt) {
    setGates((prev) => {
      if (!prev[key]) return prev
      const newVertices = [...prev[key].vertices]
      newVertices[index] = newPt
      return {
        ...prev,
        [key]: {
          ...prev[key],
          vertices: newVertices
        }
      }
    })
  }

  function updateGateVertices(key, newVertices) {
    setGates((prev) => {
      if (!prev[key]) return prev
      return {
        ...prev,
        [key]: {
          ...prev[key],
          vertices: newVertices
        }
      }
    })
  }

  function clearActiveGate() {
    setDraft([])
    setGates((prev) => {
      const next = { ...prev }
      if (activeGate === 'positive') {
        delete next[histogramGateKey(histogramGateType, selected)]
      } else {
        delete next[activeKey]
      }
      return next
    })
  }

  function clearSelectedHistogramGates() {
    setDraft([])
    setDrawActive(false)
    setGates((prev) => clearHistogramGatesForFile(prev, selected))
  }

  function makeFileSpecific() {
    if (!contextMenu?.gateKey) return
    const source = gates[contextMenu.gateKey]
    if (!source || source.type === 'positive') return
    const targetType = source.scope === 'cells' || source.scope === 'beads' ? source.scope : controlType
    setGates((prev) => {
      const next = { ...prev }
      files
        .filter((file) => fileControlType(file) === targetType)
        .forEach((file) => {
          next[gateKey(source.type, file.filename, targetType, true)] = {
            ...source,
            scope: 'file',
            filename: file.filename,
            vertices: source.vertices?.map((vertex) => ({ ...vertex })) || [],
          }
        })
      delete next[gateKey(source.type, '', targetType, false)]
      return next
    })
    setContextMenu(null)
  }

  function useGlobalGate() {
    if (!contextMenu?.gateKey) return
    const source = gates[contextMenu.gateKey]
    if (!source || source.scope === 'cells' || source.scope === 'beads') return
    const key = gateKey(source.type, selected, controlType, false)
    setGates((prev) => {
      const next = { ...prev }
      files
        .filter((file) => fileControlType(file) === controlType)
        .forEach((file) => {
          delete next[gateKey(source.type, file.filename, controlType, true)]
        })
      next[key] = {
        ...source,
        scope: controlType,
        filename: '',
        vertices: source.vertices?.map((vertex) => ({ ...vertex })) || [],
      }
      return next
    })
    setContextMenu(null)
  }

  function clearContextGate() {
    if (!contextMenu?.gateKey) return
    const source = gates[contextMenu.gateKey]
    setDraft([])
    setGates((prev) => {
      const next = { ...prev }
      if (source && (source.type === 'cell' || source.type === 'singlet') && source.scope === 'file') {
        next[contextMenu.gateKey] = {
          ...source,
          scope: 'file',
          filename: selected,
          mode: 'blocked',
          vertices: [],
        }
      } else {
        delete next[contextMenu.gateKey]
      }
      return next
    })
    setContextMenu(null)
  }

  async function saveConfig(closeAfter = false, choosePath = false) {
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform })
    try {
      const saved = await useApi('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows }),
      })
      if (saved.cancelled) {
        setStatus(saved.message || 'Save cancelled')
        return
      }
      if (saved.success === false) {
        setStatus(saved.message || 'Could not save gate CSV')
        return
      }
      setStatus(`Saved ${saved.rows ?? rows.length} rows to ${saved.path}`)
      const cfg = await useApi('/gate_configs')
      setConfigs(cfg.configs || [])
      if (closeAfter) {
        await useApi('/gate_shutdown', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{}' })
        window.close()
      }
    } catch (err) {
      setStatus(`Could not save gate CSV: ${err.message}`)
    }
  }

  async function saveConfigAsCsv() {
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform })
    const csvText = gateRowsToCsv(rows)
    let filename = CONFIG_NAME
    try {
      filename = await saveCsvWithSystemPicker(CONFIG_NAME, csvText)
    } catch (err) {
      if (err?.name === 'AbortError') {
        setStatus('Save cancelled')
      } else {
        setStatus(`Could not open save picker: ${err.message}`)
      }
      return
    }

    setStatus('Saving gate CSV')
    try {
      const saved = await useApi('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows }),
      }).catch(() => null)
      const cfg = await useApi('/gate_configs')
      setConfigs(cfg.configs || [])
      setStatus(`Saved ${saved?.rows ?? rows.length} rows to ${filename}`)
    } catch (err) {
      setStatus(`Saved ${rows.length} rows to ${filename}; backend sync failed: ${err.message}`)
    }
  }

  async function applyLoadedConfigRows(rows, sourceName, headers = REQUIRED_GATE_CSV_COLUMNS) {
    const fileData = await useApi('/gate_files')
    const currentFiles = (fileData.files || []).filter((file) => file.file_exists)
    setFiles(currentFiles)
    setMetadata(fileData.metadata || {})
    setSelected((current) => (
      currentFiles.some((file) => file.filename === current)
        ? current
        : currentFiles[0]?.filename || ''
    ))

    const compatible = reconcileGateCsvRows(rows, currentFiles)
    validateGateCsvRows(compatible.rows, headers)
    const parsed = parseConfigRows(compatible.rows)
    setGates(pruneUnavailableNegativeGates(parsed.gates, currentFiles))
    if (typeof parsed.pointSize === 'number') setPointSize(parsed.pointSize)
    if (typeof parsed.maxPoints === 'number') setMaxPoints(normalizeEventCount(parsed.maxPoints))
    if (typeof parsed.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(parsed.histogramBins))
    if (typeof parsed.histogramTransform === 'string') setHistogramTransform(normalizeHistogramTransform(parsed.histogramTransform))
    setDraft([])
    await useApi('/gate_config', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filename: CONFIG_NAME, rows: compatible.rows }),
    }).catch(() => null)
    const cfg = await useApi('/gate_configs')
    setConfigs(cfg.configs || [])
    setStatus(`Loaded ${compatible.rows.length} rows from ${sourceName}`)
  }

  async function loadConfigFromPicker() {
    const openPicker = window.showOpenFilePicker
    if (typeof openPicker !== 'function') {
      loadInputRef.current?.click()
      return
    }
    try {
      const [handle] = await openPicker({
        multiple: false,
        types: [{
          description: 'CSV files',
          accept: { 'text/csv': ['.csv'] },
        }],
      })
      if (!handle) return
      const file = await handle.getFile()
      const text = await file.text()
      const parsedCsv = parseCsvDocument(text)
      await applyLoadedConfigRows(parsedCsv.rows, file.name, parsedCsv.headers)
    } catch (err) {
      if (err?.name === 'AbortError') {
        setStatus('Load cancelled')
      } else {
        setStatus(`Could not load gate CSV: ${err.message}`)
      }
    }
  }

  async function loadConfigFromFile(file) {
    if (!file) return
    try {
      const text = await file.text()
      const parsedCsv = parseCsvDocument(text)
      await applyLoadedConfigRows(parsedCsv.rows, file.name, parsedCsv.headers)
    } catch (err) {
      setStatus(`Could not load gate CSV: ${err.message}`)
    } finally {
      if (loadInputRef.current) loadInputRef.current.value = ''
    }
  }

  async function autogateHistograms() {
    if (!canAutogateHistograms) return
    setHistogramAutogating(true)
    setStatus('Auto-generating histogram gates')
    try {
      const result = await useApi('/gate_histogram_autogate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ gates }),
      })
      if (result.success === false) throw new Error(result.error || 'Histogram autogating failed')
      setGates((previous) => pruneUnavailableNegativeGates({ ...previous, ...(result.gates || {}) }, files))
      setSpectrumCache({})
      const generated = Number(result.gates_generated ?? Object.keys(result.gates || {}).length)
      const preserved = Number(result.gates_preserved ?? 0)
      setStatus(
        `Auto-generated ${generated} missing histogram gate${generated === 1 ? '' : 's'} ` +
        `for ${result.files_processed} controls; preserved ${preserved} existing gate${preserved === 1 ? '' : 's'}.`,
      )
    } catch (err) {
      setStatus(`Could not auto-generate histogram gates: ${err.message}`)
    } finally {
      setHistogramAutogating(false)
    }
  }

  return (
    <main className={`app-shell ${initialLoading || histogramAutogating ? 'is-initial-loading' : ''}`}>
      <aside className="sidebar">
        <div className="brand">
          <ScatterChart size={26} />
          <div>
            <h1>Control File Gating</h1>
          </div>
        </div>
        <div className="file-list">
          {files.map((file) => (
            <button
              key={file.filename}
              className={`file-row ${selected === file.filename ? 'selected' : ''} ${fileUsesHistogramGates(file) ? '' : 'is-af'}`}
              onClick={() => setSelected(file.filename)}
            >
              <span>{file.fluorophore}</span>
              <strong>{file.marker}</strong>
              <small>{file.channel}{fileUsesHistogramGates(file) ? '' : ' · AF'}</small>
            </button>
          ))}
        </div>
      </aside>

      <section className="workspace">
        <header className="gating-topbar">
          <div className="control-heading">
            <h1>
              {currentFile.fluorophore} <span>{currentFile.marker}</span>
              {payload?.total_events && (
                <span style={{ fontSize: 15, fontWeight: 500, color: 'var(--muted)', marginLeft: 10 }}>
                  ({payload.total_events.toLocaleString()} events · {controlType})
                </span>
              )}
            </h1>
            <div className="control-navigation" aria-label="Control navigation">
              <button title="Previous file" aria-label="Previous file" className="icon-button" disabled={!canGoPrevious} onClick={() => selectRelativeFile(-1)}>
                <ChevronLeft size={18} />
              </button>
              <button title="Next file" aria-label="Next file" className="icon-button" disabled={!canGoNext} onClick={() => selectRelativeFile(1)}>
                <ChevronRight size={18} />
              </button>
            </div>
          </div>
          <div className="actions">
            <button
              className="auto-histogram-button"
              title={canAutogateHistograms
                ? 'Auto-generate missing required histogram gates'
                : histogramAutogateMissing.length > 0
                  ? 'Complete FSC/SSC and singlet gates for all non-AF controls first'
                  : 'No non-AF single-color controls are available'}
              disabled={!canAutogateHistograms}
              onClick={() => { if (canAutogateHistograms) setShowAutogateConfirmModal(true) }}
            >
              <WandSparkles size={17} /> Auto-gate histograms
            </button>
            <button title="Settings" onClick={() => setShowSettingsModal(true)}>
              <Settings size={17} />
            </button>
            <input
              ref={loadInputRef}
              type="file"
              accept=".csv,text/csv"
              className="hidden-file-input"
              onChange={(event) => loadConfigFromFile(event.target.files?.[0] || null)}
            />
            <button title="Save gates as CSV" onClick={() => saveConfigAsCsv()}><Save size={17} /> Save</button>
            <button title="Load gates from CSV" onClick={() => loadConfigFromPicker()}><FolderOpen size={17} /> Load</button>
            <span
              className="confirm-wrapper"
              onMouseEnter={() => !canConfirm && setShowConfirmIssues(true)}
              onMouseLeave={() => setShowConfirmIssues(false)}
              onFocus={() => !canConfirm && setShowConfirmIssues(true)}
              onBlur={() => setShowConfirmIssues(false)}
            >
              <button
                title={canConfirm ? 'Export and close' : 'Finish required gates before confirming'}
                className="confirm"
                disabled={!canConfirm}
                onClick={() => { if (canConfirm) setShowConfirmModal(true) }}
              >
                <CheckCircle2 size={18} /> Confirm
              </button>
              {!canConfirm && showConfirmIssues && (
                <div className="confirm-tooltip" role="status">
                  <strong>Cannot confirm yet</strong>
                  {confirmIssueLines.map((line) => (
                    <span key={line}>{line}</span>
                  ))}
                </div>
              )}
            </span>
          </div>
        </header>

        {(status.startsWith('Could not') || (gatesLoaded && files.length === 0)) && (
          <div
            className={`gating-status-banner ${status.startsWith('Could not') ? 'is-error' : ''}`}
            role={status.startsWith('Could not') ? 'alert' : 'status'}
          >
            <strong>
              {status.startsWith('Could not')
                ? 'Controls could not be loaded'
                : 'No control files found'}
            </strong>
            <span>
              {status.startsWith('Could not')
                ? status
                : 'Add FCS control files to the configured SCC folder, then reopen or refresh the gating GUI.'}
            </span>
          </div>
        )}

        <div className="plot-grid">
          <GatePlot
            title={isBead ? 'Beads' : 'Cells'}
            subtitle={`${channelTitle(cellAxes.x)} / ${channelTitle(cellAxes.y)} population gate`}
            events={events}
            xField={cellAxes.x}
            yField={cellAxes.y}
            xChannel={cellAxes.x}
            yChannel={cellAxes.y}
            labels={labels}
            active={activeGate === 'cell'}
            gate={cellGate}
            draft={activeGate === 'cell' ? draft : []}
            onAddPoint={(p) => {
              if (!p) {
                setDraft([])
                return
              }
              if (draft.length === 0) {
                setGates((prev) => {
                  const next = { ...prev }
                  delete next[gateKey('cell', selected, controlType)]
                  delete next[gateKey('cell', selected, controlType, true)]
                  return next
                })
              }
              setDraft((d) => [...d, p])
            }}
            onFinish={finishDraft}
            onClickPlot={() => { if (activeGate !== 'cell') { setActiveGate('cell'); setDraft([]); setDrawActive(false) } }}
            onUpdateVertex={(index, pt) => {
              if (activeGate === 'cell' && draft.length) {
                setDraft((d) => {
                  const next = [...d]
                  next[index] = pt
                  return next
                })
              } else {
                updateGateVertex(gates[`cell:${selected}`] ? `cell:${selected}` : `cell:${controlType}`, index, pt)
              }
            }}
            onUpdateVertices={(newVertices) => {
              updateGateVertices(gates[`cell:${selected}`] ? `cell:${selected}` : `cell:${controlType}`, newVertices)
            }}
            onDragEnd={handleDragEnd}
            xDomain={domainForChannel(cellAxes.x, domains.fsc_a)}
            yDomain={domainForChannel(cellAxes.y, domains.ssc_a)}
            statsText={cellGate ? `${cellSummary.count.toLocaleString()} (${cellSummary.pct.toFixed(1)}%)` : null}
            drawActive={drawActive}
            pointSize={pointSize}
            onContextGate={(x, y) => setContextMenu({ x, y, gateKey: gates[`cell:${selected}`] ? `cell:${selected}` : `cell:${controlType}` })}
            availableChannels={availableScatterChannels}
            onAxisChange={(axis, channel) => updateAxis('cell', axis, channel)}
            isBead={isBead}
            onToggleDraw={(e) => {
              e.stopPropagation()
              if (activeGate === 'cell' && drawActive) {
                setDrawActive(false)
              } else {
                setActiveGate('cell')
                setDrawActive(true)
                setDraft([])
              }
            }}
            onClear={(e) => {
              e.stopPropagation()
              setActiveGate('cell')
              clearActiveGate()
            }}
          />
          <GatePlot
            title="Singlets"
            subtitle={`${channelTitle(singletAxes.x)} / ${channelTitle(singletAxes.y)} doublet filter`}
            events={cellsFilteredEvents}
            xField={singletAxes.x}
            yField={singletAxes.y}
            xChannel={singletAxes.x}
            yChannel={singletAxes.y}
            labels={labels}
            active={activeGate === 'singlet'}
            gate={singletGate}
            draft={activeGate === 'singlet' ? draft : []}
            onAddPoint={(p) => {
              if (!p) {
                setDraft([])
                return
              }
              if (draft.length === 0) {
                setGates((prev) => {
                  const next = { ...prev }
                  delete next[gateKey('singlet', selected, controlType)]
                  delete next[gateKey('singlet', selected, controlType, true)]
                  return next
                })
              }
              setDraft((d) => [...d, p])
            }}
            onFinish={finishDraft}
            onClickPlot={() => { if (activeGate !== 'singlet') { setActiveGate('singlet'); setDraft([]); setDrawActive(false) } }}
            onUpdateVertex={(index, pt) => {
              if (activeGate === 'singlet' && draft.length) {
                setDraft((d) => {
                  const next = [...d]
                  next[index] = pt
                  return next
                })
              } else {
                updateGateVertex(gates[`singlet:${selected}`] ? `singlet:${selected}` : `singlet:${controlType}`, index, pt)
              }
            }}
            onUpdateVertices={(newVertices) => {
              updateGateVertices(gates[`singlet:${selected}`] ? `singlet:${selected}` : `singlet:${controlType}`, newVertices)
            }}
            onDragEnd={handleDragEnd}
            xDomain={domainForChannel(singletAxes.x, domains.fsc_h)}
            yDomain={domainForChannel(singletAxes.y, domains.fsc_a)}
            statsText={singletGate ? `${singletSummary.count.toLocaleString()} (${singletSummary.pct.toFixed(1)}%)` : null}
            warningText={singletWarningText}
            drawActive={drawActive}
            pointSize={pointSize}
            onContextGate={(x, y) => setContextMenu({ x, y, gateKey: gates[`singlet:${selected}`] ? `singlet:${selected}` : `singlet:${controlType}` })}
            availableChannels={availableScatterChannels}
            onAxisChange={(axis, channel) => updateAxis('singlet', axis, channel)}
            isBead={isBead}
            onToggleDraw={(e) => {
              e.stopPropagation()
              if (activeGate === 'singlet' && drawActive) {
                setDrawActive(false)
              } else {
                setActiveGate('singlet')
                setDrawActive(true)
                setDraft([])
              }
            }}
            onClear={(e) => {
              e.stopPropagation()
              setActiveGate('singlet')
              clearActiveGate()
            }}
          />
          {fileUsesHistogramGates(currentFile) ? (
            <GatePlot
              title="Histogram"
              subtitle=""
              events={singletsFilteredEvents}
              xField="peak"
              yField="count"
              xChannel={channels.peak}
              yChannel=""
              labels={labels}
              active={activeGate === 'positive'}
              gate={activeHistogramGate}
              draft={activeGate === 'positive' ? draft : []}
              mode="histogram"
              onAddPoint={(p) => {
                if (!p) {
                  setDraft([])
                  return
                }
                if (histogramGateType === 'positive' || histogramGateType === 'negative') {
                  if (draft.length === 0) {
                    setDraft([p])
                  } else {
                    const minX = Math.min(draft[0].x, p.x)
                    const maxX = Math.max(draft[0].x, p.x)
                    const key = histogramGateKey(histogramGateType, selected)
                    setGates((prev) => ({
                      ...prev,
                      [key]: {
                        type: histogramGateType === 'positive' ? 'positive' : 'negative',
                        scope: 'file',
                        filename: selected,
                        xChannel: channels.peak,
                        yChannel: '',
                        mode: histogramGateType === 'positive' ? 'positive_1d' : 'negative_1d',
                        vertices: [{ x: minX, y: 0 }, { x: maxX, y: 0 }]
                      }
                    }))
                    setDraft([])
                    setDrawActive(false)
                  }
                  return
                }
              }}
              onFinish={finishDraft}
              onClickPlot={() => { if (activeGate !== 'positive') { setActiveGate('positive'); setDraft([]); setDrawActive(false) } }}
              onUpdateVertex={(index, pt) => {
                updateGateVertex(histogramGateKey(histogramGateType, selected), index, pt)
              }}
              onUpdateVertices={(newVertices) => {
                updateGateVertices(histogramGateKey(histogramGateType, selected), newVertices)
              }}
              onDragEnd={handleDragEnd}
              xDomain={domains.peak}
              statsText={[
                positiveGate ? `Pos ${positiveSummary.count.toLocaleString()} (${positiveSummary.pct.toFixed(1)}%)` : null,
                negativeGate ? `Neg ${negativeSummary.count.toLocaleString()} (${negativeSummary.pct.toFixed(1)}%)` : null,
              ].filter(Boolean).join(' · ') || null}
              warningText={positiveWarningText}
              drawActive={drawActive}
              pointSize={pointSize}
              onContextGate={(x, y, key) => setContextMenu({ x, y, gateKey: key || histogramGateKey(histogramGateType, selected) })}
              isPositivePlot={true}
              histogramGateType={histogramGateType}
              negativeGateEnabled={usesNegativeHistogramGate}
              histogramBins={histogramBins}
              histogramTransform={histogramTransform}
              onHistogramTransformChange={(value) => setHistogramTransform(normalizeHistogramTransform(value))}
              onHistogramBinsChange={(value) => setHistogramBins(normalizeHistogramBins(value))}
              secondaryGates={secondaryHistogramGates}
              onSelectHistogramGate={(type) => {
                setActiveGate('positive')
                setHistogramGateType(type)
                setDrawActive(false)
                setDraft([])
              }}
              onUpdateAnyGateVertices={(key, newVertices) => updateGateVertices(key, newVertices)}
              onToggleHistogramGate={(type, e) => {
                e.stopPropagation()
                if (type === 'negative' && !usesNegativeHistogramGate) return
                setActiveGate('positive')
                setHistogramGateType(type)
                setDraft([])
                if (activeGate === 'positive' && histogramGateType === type && drawActive) {
                  setDrawActive(false)
                  return
                }
                setDrawActive(true)
              }}
              onClear={(e) => {
                e.stopPropagation()
                setActiveGate('positive')
                clearSelectedHistogramGates()
              }}
            />
          ) : null}
        </div>

        <div className="spectrum-container" style={{ marginTop: 28, display: 'flex', justifyContent: 'center' }}>
          <div className="plot-panel spectrum-panel" style={{ minHeight: 'auto', padding: '16px 20px', width: '100%' }}>
              <h2 style={{ fontSize: 18, margin: '0 0 12px 0', color: 'var(--ink)' }}>Spectrum</h2>
              {spectrum ? (
                <img
                  src={spectrum}
                  alt="Detector Spectrum"
                  style={{ width: '100%', height: 'auto', display: 'block', borderRadius: 6, border: '1px solid var(--line)' }}
                />
              ) : (
                <div className="spectrum-placeholder" aria-label="Spectrum loading placeholder" />
              )}
          </div>
        </div>
      </section>

      {contextMenu && gates[contextMenu.gateKey] && (
        <div className="context-menu" style={{ left: contextMenu.x, top: contextMenu.y }}>
          {gates[contextMenu.gateKey].type !== 'positive' && gates[contextMenu.gateKey].type !== 'negative' && (
            gates[contextMenu.gateKey].scope === 'cells' || gates[contextMenu.gateKey].scope === 'beads' ? (
              <button onClick={makeFileSpecific}><FileDown size={15} /> Make file specific</button>
            ) : (
              <button onClick={useGlobalGate}><Upload size={15} /> Make gate global</button>
            )
          )}
          <button onClick={clearContextGate}><Eraser size={15} /> Clear gate</button>
        </div>
      )}

      {showConfirmModal && (
        <div className="confirm-modal-overlay">
          <div className="confirm-modal">
            <h2>Confirm Gating Configuration?</h2>
            <p>Are you sure you want to save all gate configurations and exit the Gating GUI?</p>
            <div className="confirm-modal-actions">
              <button className="cancel-btn" onClick={() => setShowConfirmModal(false)}>Cancel</button>
              <button className="confirm-btn" disabled={!canConfirm} onClick={() => { if (canConfirm) { setShowConfirmModal(false); saveConfig(true) } }}>Confirm & Exit</button>
            </div>
          </div>
        </div>
      )}

      {showAutogateConfirmModal && (
        <div className="confirm-modal-overlay">
          <div className="confirm-modal">
            <h2>Auto-generate Histogram Gates?</h2>
            <p>
              Generate missing positive gates for {histogramAutogateTargets.length} controls.
              {' '}{histogramAutogateNegativeTargets.length} also require an internal negative gate.
              Existing gates are preserved.
            </p>
            <div className="confirm-modal-actions">
              <button className="cancel-btn" onClick={() => setShowAutogateConfirmModal(false)}>Cancel</button>
              <button
                className="confirm-btn"
                disabled={!canAutogateHistograms}
                onClick={() => {
                  if (!canAutogateHistograms) return
                  setShowAutogateConfirmModal(false)
                  autogateHistograms()
                }}
              >
                Generate Gates
              </button>
            </div>
          </div>
        </div>
      )}

      {initialLoading && (
        <div className="initial-loading-overlay" role="status" aria-live="polite" aria-label="Loading files">
          <div className="initial-loading-card">
            <div className="initial-loading-spinner" aria-hidden="true" />
            <strong>Loading files...</strong>
          </div>
        </div>
      )}

      {histogramAutogating && (
        <div className="initial-loading-overlay" role="status" aria-live="polite" aria-label="Auto-generating histogram gates">
          <div className="initial-loading-card">
            <div className="initial-loading-spinner" aria-hidden="true" />
            <strong>Generating histogram gates...</strong>
            <span>
              Required gates · {histogramAutogateTargets.length} controls
              {' · '}{histogramAutogateNegativeTargets.length} internal negatives
            </span>
          </div>
        </div>
      )}

      {/* Settings Modal */}
      {showSettingsModal && (
        <div className="confirm-modal-overlay" onClick={() => setShowSettingsModal(false)}>
          <div className="confirm-modal" onClick={(e) => e.stopPropagation()}>
            <h2>Gating Settings</h2>
            <div style={{ display: 'grid', gap: 20, margin: '24px 0', textAlign: 'left' }}>
              <div style={{ display: 'grid', gap: 8 }}>
                <label style={{ fontWeight: 'bold', fontSize: 14 }}>Theme:</label>
                <button type="button" onClick={() => setDarkMode(!darkMode)} style={{ justifyContent: 'flex-start' }}>
                  {darkMode ? <Sun size={17} /> : <Moon size={17} />}
                  {darkMode ? 'Light mode' : 'Dark mode'}
                </button>
              </div>
              <div style={{ display: 'grid', gap: 6 }}>
                <label style={{ fontWeight: 'bold', fontSize: 14 }}>Point Size ({pointSize}px):</label>
                <input
                  type="range"
                  min="0.5"
                  max="4.0"
                  step="0.25"
                  value={pointSize}
                  onChange={(e) => setPointSize(Number(e.target.value))}
                  style={{ width: '100%', cursor: 'pointer' }}
                />
              </div>
              <div style={{ display: 'grid', gap: 6 }}>
                <label style={{ fontWeight: 'bold', fontSize: 14 }}>
                  Show Downsampled Events ({eventStepLabel(maxPoints, payload?.total_events)}):
                </label>
                <input
                  type="range"
                  min="0"
                  max={EVENT_COUNT_STEPS.length - 1}
                  step="1"
                  value={eventStepIndex(maxPoints)}
                  onChange={(e) => setMaxPoints(EVENT_COUNT_STEPS[Number(e.target.value)])}
                  style={{ width: '100%', cursor: 'pointer' }}
                />
                <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: 11, color: 'var(--muted)', fontWeight: 800 }}>
                  <span>1K</span>
                  <span>3K</span>
                  <span>10K</span>
                  <span>50K</span>
                  <span>100K</span>
                </div>
              </div>
              <div style={{ display: 'grid', gap: 6 }}>
                <label style={{ fontWeight: 'bold', fontSize: 14 }}>Histogram ({histogramBins} bins):</label>
                <input
                  type="range"
                  min={HISTOGRAM_BIN_MIN}
                  max={HISTOGRAM_BIN_MAX}
                  step="5"
                  value={histogramBins}
                  onChange={(e) => setHistogramBins(normalizeHistogramBins(e.target.value))}
                  onInput={(e) => setHistogramBins(normalizeHistogramBins(e.currentTarget.value))}
                  style={{ width: '100%', cursor: 'pointer' }}
                />
              </div>
              <div style={{ display: 'grid', gap: 6 }}>
                <TransformDropdown
                  value={histogramTransform}
                  onChange={(nextValue) => setHistogramTransform(normalizeHistogramTransform(nextValue))}
                />
              </div>
            </div>
            <div className="confirm-modal-actions">
              <button className="confirm-btn" onClick={() => setShowSettingsModal(false)}>Close</button>
            </div>
          </div>
        </div>
      )}
    </main>
  )
}

export default App
