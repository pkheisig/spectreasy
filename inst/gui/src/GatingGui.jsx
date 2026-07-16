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
  Info,
  PanelLeftClose,
  PanelLeftOpen,
} from 'lucide-react'
import './GatingGui.css'
import { reconcileGateCsvRows } from './gatingCsvCompatibility.js'
import {
  fileUsesHistogramGates,
  fileUsesNegativeHistogramGate,
  pruneUnavailableNegativeGates,
} from './gatingEligibility.js'
import { clearHistogramGatesForFile, histogramDomainIncludingGates } from './histogramGates.js'
import {
  buildViewSettingRows,
  normalizePlotView,
  parseViewSettings,
} from './gatingViewSettings.js'
import SpectrumCanvas from './SpectrumCanvas.jsx'
import { decodeSpectrumData } from './spectrumData.js'
import { resolveApiBase, resolveApiToken } from './apiBase'
import { withGatingApiToken } from './gatingApi.js'

const API_BASE = resolveApiBase()
const API_TOKEN = resolveApiToken()
const CONFIG_NAME = 'ssc_gate_config.csv'
const GUI_MODULE = 'control_gating'

const unboxGuiState = (value) => {
  if (Array.isArray(value)) {
    if (value.length === 1) return unboxGuiState(value[0])
    return value.map(unboxGuiState)
  }
  if (value && typeof value === 'object') {
    return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, unboxGuiState(item)]))
  }
  return value
}
const THEME_STORAGE_KEY = 'spectreasy-theme'
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
const MIN_CONFIRM_EVENTS = 200
const REQUIRED_GATE_CSV_COLUMNS = ['gate_type', 'scope', 'filename', 'x_channel', 'y_channel', 'plot_mode', 'vertex_index', 'x', 'y']
const PLOT_WIDTH = 520
const PLOT_HEIGHT = 420
const PAD = { left: 54, right: 18, top: 18, bottom: 46 }

function HistogramSparkleIcon() {
  return (
    <svg className="histogram-sparkle-icon" viewBox="0 0 34 26" aria-hidden="true">
      <path className="histogram-sparkle-curve" d="M3 21 C7 21 7 17 10 17 C13 17 13 20 16 20 C20 20 20 6 24 6 C28 6 27 21 31 21" />
      <path className="histogram-sparkle-star star-one" d="M8 3 L9 5.4 L11.5 6.3 L9 7.2 L8 9.7 L7 7.2 L4.5 6.3 L7 5.4 Z" />
      <path className="histogram-sparkle-star star-two" d="M29 1 L29.7 2.8 L31.5 3.5 L29.7 4.2 L29 6 L28.3 4.2 L26.5 3.5 L28.3 2.8 Z" />
      <path className="histogram-sparkle-star star-three" d="M17 8 L17.6 9.4 L19 10 L17.6 10.6 L17 12 L16.4 10.6 L15 10 L16.4 9.4 Z" />
    </svg>
  )
}

function axisLabel(labels, channel) {
  const cleanChannel = String(channel || '').replace(/\s*\(\s*\)\s*$/, '').trim()
  const desc = String(labels?.[channel] || '').replace(/^\s*\(\s*\)\s*$/, '').trim()
  return desc && desc !== cleanChannel ? `${cleanChannel} (${desc})` : cleanChannel
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

function pointInPolygonValues(x, y, polygon) {
  if (!polygon || polygon.length < 3) return false
  let inside = false
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const xi = polygon[i].x
    const yi = polygon[i].y
    const xj = polygon[j].x
    const yj = polygon[j].y
    const intersect = yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi
    if (intersect) inside = !inside
  }
  return inside
}

function pointInPolygon(point, polygon) {
  return pointInPolygonValues(point.x, point.y, polygon)
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

function decodeCompactPayload(value) {
  const payload = unboxGuiState(value)
  const compact = payload?.events_compact
  if (compact?.format !== 'float32-column-major' || typeof compact.data !== 'string') return payload
  const binary = window.atob(compact.data.replace(/\s/g, ''))
  const bytes = new Uint8Array(binary.length)
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i)
  return {
    ...payload,
    events_compact: {
      format: compact.format,
      fields: Array.isArray(compact.fields) ? compact.fields : [compact.fields],
      rows: Number(compact.rows) || 0,
      values: new Float32Array(bytes.buffer),
    },
  }
}

function materializePayloadEvents(payload, maxPoints) {
  if (Array.isArray(payload?.events)) return displayEventsForLimit(payload.events, maxPoints)
  const compact = payload?.events_compact
  const rows = Number(compact?.rows) || 0
  const fields = compact?.fields || []
  const values = compact?.values
  if (!rows || !fields.length || !(values instanceof Float32Array)) return []
  const limit = Math.min(rows, normalizeEventCount(maxPoints))
  const step = rows / Math.max(limit, 1)
  const events = new Array(limit)
  for (let outputIndex = 0; outputIndex < limit; outputIndex++) {
    const rowIndex = Math.floor(outputIndex * step)
    const event = {}
    for (let fieldIndex = 0; fieldIndex < fields.length; fieldIndex++) {
      event[fields[fieldIndex]] = values[fieldIndex * rows + rowIndex]
    }
    events[outputIndex] = event
  }
  return events
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

function useApi(path, options = {}) {
  return fetch(`${API_BASE}${path}`, withGatingApiToken(options, API_TOKEN)).then((res) => {
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
const DARK_DENSITY_PALETTE = Array.from({ length: DENSITY_PALETTE.length }, (_, index) => {
  const value = index / Math.max(DENSITY_PALETTE.length - 1, 1)
  const hue = 218 - value * 198
  const lightness = 72 - value * 10
  return `hsla(${hue}, 100%, ${lightness}%, 0.94)`
})

function computeDensityBuckets(points, xField, yField, xDomain, yDomain) {
  const n = points.length
  const buckets = Array.from({ length: DENSITY_PALETTE.length }, () => [])
  if (n === 0) return buckets

  const xs = points.map((p) => p[xField])
  const ys = points.map((p) => p[yField])
  const [minX, maxX] = xDomain || [Math.min(...xs), Math.max(...xs)]
  const [minY, maxY] = yDomain || [Math.min(...ys), Math.max(...ys)]

  const rx = maxX - minX || 1
  const ry = maxY - minY || 1

  // Build a smooth density field, then sample it continuously at every event.
  // Assigning one density to every event in a hard grid cell makes the colour
  // field look tiled when a long scatter tail compresses the main population
  // into a relatively small part of the plot.
  const numBins = 256
  const gridSide = numBins + 1
  const grid = new Float32Array(gridSide * gridSide)

  points.forEach((p) => {
    if (p[xField] < minX || p[xField] > maxX || p[yField] < minY || p[yField] > maxY) return
    const bx = clamp(Math.floor(((p[xField] - minX) / rx) * numBins), 0, numBins)
    const by = clamp(Math.floor(((p[yField] - minY) / ry) * numBins), 0, numBins)
    grid[by * gridSide + bx] += 1
  })

  // A separable Gaussian blur gives every grid node a local density estimate.
  // Bilinear sampling below removes the remaining cell boundaries.
  const sigma = 2.0
  const radius = Math.ceil(sigma * 3)
  const kernel = new Float32Array(radius * 2 + 1)
  let kernelSum = 0
  for (let offset = -radius; offset <= radius; offset++) {
    const weight = Math.exp(-(offset * offset) / (2 * sigma * sigma))
    kernel[offset + radius] = weight
    kernelSum += weight
  }
  for (let i = 0; i < kernel.length; i++) kernel[i] /= kernelSum

  const horizontal = new Float32Array(grid.length)
  const smoothed = new Float32Array(grid.length)
  for (let y = 0; y < gridSide; y++) {
    const row = y * gridSide
    for (let x = 0; x < gridSide; x++) {
      let sum = 0
      for (let offset = -radius; offset <= radius; offset++) {
        const nx = clamp(x + offset, 0, numBins)
        sum += grid[row + nx] * kernel[offset + radius]
      }
      horizontal[row + x] = sum
    }
  }
  for (let y = 0; y < gridSide; y++) {
    for (let x = 0; x < gridSide; x++) {
      let sum = 0
      for (let offset = -radius; offset <= radius; offset++) {
        const ny = clamp(y + offset, 0, numBins)
        sum += horizontal[ny * gridSide + x] * kernel[offset + radius]
      }
      smoothed[y * gridSide + x] = sum
    }
  }

  const densities = new Float32Array(n)
  let maxD = 1
  for (let pointIndex = 0; pointIndex < n; pointIndex++) {
    const p = points[pointIndex]
    if (p[xField] < minX || p[xField] > maxX || p[yField] < minY || p[yField] > maxY) continue
    const gx = clamp(((p[xField] - minX) / rx) * numBins, 0, numBins)
    const gy = clamp(((p[yField] - minY) / ry) * numBins, 0, numBins)
    const x0 = Math.floor(gx)
    const y0 = Math.floor(gy)
    const x1 = Math.min(x0 + 1, numBins)
    const y1 = Math.min(y0 + 1, numBins)
    const tx = gx - x0
    const ty = gy - y0
    const top = smoothed[y0 * gridSide + x0] * (1 - tx) + smoothed[y0 * gridSide + x1] * tx
    const bottom = smoothed[y1 * gridSide + x0] * (1 - tx) + smoothed[y1 * gridSide + x1] * tx
    const density = top * (1 - ty) + bottom * ty
    densities[pointIndex] = density
    if (density > maxD) maxD = density
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
  darkMode = false,
  viewDomain = null,
  onViewDomainChange,
}) {
  const plotRef = useRef(null)
  const menuRef = useRef(null)
  const svgRef = useRef(null)
  const canvasRef = useRef(null)
  const [dragTarget, setDragTarget] = useState(null)
  const [lastMousePos, setLastMousePos] = useState(null)
  const [hoverPt, setHoverPt] = useState(null)
  const [dragPreviewVertices, setDragPreviewVertices] = useState(null)
  const [panTarget, setPanTarget] = useState(null)
  const [axisMenu, setAxisMenu] = useState(null)
  const [zoomXDomain, setZoomXDomain] = useState(null)
  const [zoomYDomain, setZoomYDomain] = useState(null)
  const [plotHovered, setPlotHovered] = useState(false)
  const panMovedRef = useRef(false)
  const pendingViewRef = useRef(null)

  const transformFns = useMemo(() => histogramTransformFns(histogramTransform), [histogramTransform])
  const toPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.forward(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }
  const fromPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.inverse(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }

  const baseXDomain = useMemo(() => {
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
  const xDomain = zoomXDomain || baseXDomain

  const densityCurve = useMemo(() => {
    if (mode !== 'histogram') return []
    return getDensityCurve(events.map((e) => toPlotX(e[xField])), xDomain, histogramBins)
  }, [events, mode, xField, xDomain, histogramBins, histogramTransform])

  const baseYDomain = useMemo(() => {
    if (yDomainProp) return yDomainProp
    if (mode === 'histogram') {
      if (densityCurve.length === 0) return [0, 1]
      return [0, Math.max(...densityCurve.map((p) => p.y)) * 1.12]
    }
    return extent(events.map((e) => e[yField]))
  }, [events, mode, yDomainProp, densityCurve])
  const yDomain = zoomYDomain || baseYDomain

  const xScale = useMemo(() => makeScale(xDomain, [PAD.left, PLOT_WIDTH - PAD.right]), [xDomain])
  const yScale = useMemo(() => makeScale(yDomain, [PLOT_HEIGHT - PAD.bottom, PAD.top]), [yDomain])
  const xToScreen = (value) => xScale.toScreen(toPlotX(value))
  const renderGate = dragPreviewVertices && gate ? { ...gate, vertices: dragPreviewVertices } : gate
  const displayGate = draft?.length ? { vertices: draft } : renderGate
  const densityBuckets = useMemo(() => {
    if (mode !== 'scatter') return []
    const densityXDomain = extent(events.map((event) => event[xField]))
    const densityYDomain = extent(events.map((event) => event[yField]))
    return computeDensityBuckets(events, xField, yField, densityXDomain, densityYDomain)
  }, [events, mode, xField, yField])

  const xTicks = useMemo(() => getTicks(xDomain, 5), [xDomain])
  const yTicks = useMemo(() => getTicks(yDomain, 5), [yDomain])

  useEffect(() => {
    const saved = normalizePlotView(viewDomain, mode !== 'histogram')
    const savedX = mode === 'histogram' && saved?.x
      ? saved.x.map((value) => transformFns.forward(value))
      : saved?.x
    setZoomXDomain(savedX || null)
    setZoomYDomain(mode === 'histogram' ? null : (saved?.y || null))
    pendingViewRef.current = savedX ? { x: savedX, y: mode === 'histogram' ? null : (saved?.y || null) } : null
  }, [
    xField,
    yField,
    mode,
    xDomainProp?.[0],
    xDomainProp?.[1],
    yDomainProp?.[0],
    yDomainProp?.[1],
    histogramTransform,
    viewDomain?.x?.[0],
    viewDomain?.x?.[1],
    viewDomain?.y?.[0],
    viewDomain?.y?.[1],
  ])

  function emitViewChange(view) {
    if (!view) {
      onViewDomainChange?.(null)
      return
    }
    if (mode === 'histogram') {
      onViewDomainChange?.({ x: view.x.map((value) => fromPlotX(value)), y: null })
      return
    }
    onViewDomainChange?.(view)
  }

  function zoomDomainAround(domain, baseDomain, factor, anchor) {
    const [baseMin, baseMax] = baseDomain
    const [currentMin, currentMax] = domain
    const baseSpan = baseMax - baseMin
    const currentSpan = currentMax - currentMin
    if (!Number.isFinite(baseSpan) || baseSpan <= 0 || !Number.isFinite(currentSpan) || currentSpan <= 0) return domain
    const targetSpan = clamp(currentSpan * factor, baseSpan * 0.04, baseSpan)
    const ratio = clamp((anchor - currentMin) / currentSpan, 0, 1)
    let nextMin = anchor - targetSpan * ratio
    let nextMax = nextMin + targetSpan
    if (nextMin < baseMin) {
      nextMax += baseMin - nextMin
      nextMin = baseMin
    }
    if (nextMax > baseMax) {
      nextMin -= nextMax - baseMax
      nextMax = baseMax
    }
    return [Math.max(baseMin, nextMin), Math.min(baseMax, nextMax)]
  }

  function applyPlotZoom(factor, anchorX, anchorY) {
    if (mode === 'histogram') return
    const nextX = zoomDomainAround(xDomain, baseXDomain, factor, anchorX)
    const nextY = zoomDomainAround(yDomain, baseYDomain, factor, anchorY)
    setZoomXDomain(nextX)
    setZoomYDomain(nextY)
    pendingViewRef.current = { x: nextX, y: nextY }
    emitViewChange(pendingViewRef.current)
  }

  function handleWheel(evt) {
    if (mode === 'histogram') return
    evt.preventDefault()
    evt.stopPropagation()
    const point = screenToData(evt)
    const factor = evt.deltaY < 0 ? 0.97 : 1.03
    applyPlotZoom(factor, point.plotX, point.y)
  }

  useEffect(() => {
    const plotNode = svgRef.current
    if (!plotNode) return undefined
    const wheelListener = (event) => handleWheel(event)
    plotNode.addEventListener('wheel', wheelListener, { passive: false })
    return () => plotNode.removeEventListener('wheel', wheelListener)
  }, [xDomain, yDomain, baseXDomain, baseYDomain, mode])

  useEffect(() => {
    if (!plotHovered) return undefined
    const handleZoomKey = (event) => {
      if (!(event.ctrlKey || event.metaKey)) return
      const zoomIn = event.key === '+' || event.key === '='
      const zoomOut = event.key === '-' || event.key === '_'
      const reset = event.key === '0'
      if (!zoomIn && !zoomOut && !reset) return
      if (mode === 'histogram' && !reset) return
      event.preventDefault()
      if (reset) {
        setZoomXDomain(null)
        setZoomYDomain(null)
        pendingViewRef.current = null
        emitViewChange(null)
        return
      }
      applyPlotZoom(
        zoomIn ? 0.95 : 1.05,
        (xDomain[0] + xDomain[1]) / 2,
        (yDomain[0] + yDomain[1]) / 2,
      )
    }
    window.addEventListener('keydown', handleZoomKey)
    return () => window.removeEventListener('keydown', handleZoomKey)
  }, [plotHovered, xDomain, yDomain, baseXDomain, baseYDomain, mode])

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

    ctx.fillStyle = darkMode ? '#0b1110' : '#f8f7f3'
    ctx.fillRect(PAD.left, PAD.top, PLOT_WIDTH - PAD.left - PAD.right, PLOT_HEIGHT - PAD.top - PAD.bottom)
    ctx.strokeStyle = darkMode ? '#52615b' : '#c7c3ba'
    ctx.strokeRect(PAD.left, PAD.top, PLOT_WIDTH - PAD.left - PAD.right, PLOT_HEIGHT - PAD.top - PAD.bottom)

    if (events.length === 0) {
      ctx.restore()
      return
    }

    const side = Math.max(1, pointSize * 1.6)
    const half = side / 2
    densityBuckets.forEach((bucket, colorIndex) => {
      if (!bucket.length) return
      const palette = darkMode ? DARK_DENSITY_PALETTE : DENSITY_PALETTE
      ctx.fillStyle = palette[colorIndex] || (darkMode ? 'rgba(90, 190, 255, 0.9)' : 'rgba(0, 0, 255, 0.6)')
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
  }, [events, densityBuckets, xScale, yScale, mode, pointSize, xField, yField, darkMode])

  function handleClick(evt) {
    if (panMovedRef.current) {
      panMovedRef.current = false
      return
    }
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

  function shiftDomain(domain, amount, limitDomain) {
    const span = domain[1] - domain[0]
    let nextMin = domain[0] + amount
    let nextMax = domain[1] + amount
    if (nextMin < limitDomain[0]) {
      nextMax = limitDomain[0] + span
      nextMin = limitDomain[0]
    }
    if (nextMax > limitDomain[1]) {
      nextMin = limitDomain[1] - span
      nextMax = limitDomain[1]
    }
    return [nextMin, nextMax]
  }

  function handlePlotMouseDown(evt) {
    if (evt.button !== 0 || drawActive || !active) return
    const scatterCanPan = mode === 'scatter' && (zoomXDomain || zoomYDomain)
    if (!scatterCanPan && mode !== 'histogram') return
    const box = svgRef.current.getBoundingClientRect()
    const sx = (evt.clientX - box.left) * (PLOT_WIDTH / box.width)
    const sy = (evt.clientY - box.top) * (PLOT_HEIGHT / box.height)
    if (sx < PAD.left || sx > PLOT_WIDTH - PAD.right || sy < PAD.top || sy > PLOT_HEIGHT - PAD.bottom) return
    evt.preventDefault()
    panMovedRef.current = false
    setPanTarget({
      clientX: evt.clientX,
      clientY: evt.clientY,
      xDomain: [...xDomain],
      yDomain: [...yDomain],
    })
  }

  function handleMouseMove(evt) {
    const currentPt = screenToData(evt)
    setHoverPt(currentPt)

    if (panTarget) {
      const box = svgRef.current.getBoundingClientRect()
      const plotWidth = box.width * (PLOT_WIDTH - PAD.left - PAD.right) / PLOT_WIDTH
      const plotHeight = box.height * (PLOT_HEIGHT - PAD.top - PAD.bottom) / PLOT_HEIGHT
      const dxPixels = evt.clientX - panTarget.clientX
      const dyPixels = evt.clientY - panTarget.clientY
      if (Math.abs(dxPixels) > 2 || Math.abs(dyPixels) > 2) panMovedRef.current = true

      const startXSpan = panTarget.xDomain[1] - panTarget.xDomain[0]
      const xShift = -(dxPixels / plotWidth) * startXSpan
      if (mode === 'histogram') {
        const baseSpan = baseXDomain[1] - baseXDomain[0]
        const extendedLimits = [baseXDomain[0] - baseSpan * 4, baseXDomain[1] + baseSpan * 4]
        const nextX = shiftDomain(panTarget.xDomain, xShift, extendedLimits)
        setZoomXDomain(nextX)
        pendingViewRef.current = { x: nextX, y: null }
      } else {
        const startYSpan = panTarget.yDomain[1] - panTarget.yDomain[0]
        const yShift = (dyPixels / plotHeight) * startYSpan
        const nextX = shiftDomain(panTarget.xDomain, xShift, baseXDomain)
        const nextY = shiftDomain(panTarget.yDomain, yShift, baseYDomain)
        setZoomXDomain(nextX)
        setZoomYDomain(nextY)
        pendingViewRef.current = { x: nextX, y: nextY }
      }
      return
    }

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
    if (panTarget) {
      setPanTarget(null)
      if (panMovedRef.current && pendingViewRef.current) {
        emitViewChange(pendingViewRef.current)
      }
    }
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
            r="6"
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
          data-x-domain={xDomain.join(',')}
          data-y-domain={yDomain.join(',')}
          style={{
            position: 'absolute',
            top: 0,
            left: 0,
            width: '100%',
            height: '100%',
            zIndex: 2,
            background: 'transparent',
            cursor: drawActive
              ? 'crosshair'
              : panTarget
                ? 'grabbing'
                : (mode === 'histogram' || zoomXDomain || zoomYDomain)
                  ? 'grab'
                  : 'default',
          }}
          onMouseDown={handlePlotMouseDown}
          onClick={handleClick}
          onDoubleClick={onFinish}
          onContextMenu={handleContext}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseUp}
          onMouseEnter={() => setPlotHovered(true)}
          onMouseOut={(event) => {
            if (!event.currentTarget.contains(event.relatedTarget)) setPlotHovered(false)
          }}
        >
          <defs>
            <clipPath id="plot-clip">
              <rect x={PAD.left} y={PAD.top} width={PLOT_WIDTH - PAD.left - PAD.right} height={PLOT_HEIGHT - PAD.top - PAD.bottom} />
            </clipPath>
          </defs>
          
          {mode === 'histogram' && (
            <>
              <rect className="plot-bg" x={PAD.left} y={PAD.top} width={PLOT_WIDTH - PAD.left - PAD.right} height={PLOT_HEIGHT - PAD.top - PAD.bottom} />
              {densityPath && <path className="histogram-density" d={densityPath} strokeWidth="1.35" clipPath="url(#plot-clip)" />}
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
                <line className="axis-tick-line" x1={x} y1={PLOT_HEIGHT - PAD.bottom} x2={x} y2={PLOT_HEIGHT - PAD.bottom + 5} />
                <text className="axis-tick-label" x={x} y={PLOT_HEIGHT - PAD.bottom + 18} textAnchor="middle">{formatTickValue(mode === 'histogram' ? fromPlotX(t) : t)}</text>
              </g>
            )
          })}
          {yTicks.map((t) => {
            const y = yScale.toScreen(t)
            return (
              <g key={`ytick-${t}`}>
                <line className="axis-tick-line" x1={PAD.left - 5} y1={y} x2={PAD.left} y2={y} />
                <text className="axis-tick-label" x={PAD.left - 8} y={y + 3} textAnchor="end">{formatTickValue(t)}</text>
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
              r="6"
              className={`gate-handle ${isBead && !isPositivePlot ? 'gate-handle-bead' : ''}`}
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
              strokeWidth="1.6"
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
              r="6"
              className={`gate-handle ${isBead && !isPositivePlot ? 'gate-handle-bead' : ''}`}
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

function App({ embedded = false, cockpitTheme = null, onRequestExit = null }) {
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
  const [darkMode, setDarkMode] = useState(() => {
    if (embedded && (cockpitTheme === 'dark' || cockpitTheme === 'light')) return cockpitTheme === 'dark'
    const stored = window.localStorage.getItem(THEME_STORAGE_KEY)
    if (stored === 'dark' || stored === 'light') return stored === 'dark'
    return window.matchMedia?.('(prefers-color-scheme: dark)').matches || false
  })
  const [guiStateLoaded, setGuiStateLoaded] = useState(false)
  const [sidebarWidth, setSidebarWidth] = useState(192)
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false)
  const [preloadComplete, setPreloadComplete] = useState(false)
  const [gatesLoaded, setGatesLoaded] = useState(false)
  const [histogramAutogating, setHistogramAutogating] = useState(false)
  const [axisSettings, setAxisSettings] = useState({ cell: {}, singlet: {} })
  const [viewSettings, setViewSettings] = useState({ cell: {}, singlet: {}, histogram: {} })
  const [spectrum, setSpectrum] = useState(null)
  const [spectrumCache, setSpectrumCache] = useState({})
  const [spectraPrecomputing, setSpectraPrecomputing] = useState(false)
  const loadInputRef = useRef(null)
  const payloadCacheRef = useRef({})
  const spectrumBatchRef = useRef('')

  useEffect(() => {
    if (embedded) return
    window.localStorage.setItem(THEME_STORAGE_KEY, darkMode ? 'dark' : 'light')
  }, [darkMode, embedded])

  useEffect(() => {
    if (embedded && (cockpitTheme === 'dark' || cockpitTheme === 'light')) setDarkMode(cockpitTheme === 'dark')
  }, [cockpitTheme, embedded])

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
        if (hasViewSettings(csvGates.viewSettings)) {
          setViewSettings(csvGates.viewSettings)
        } else if (hasViewSettings(cacheData?.viewSettings)) {
          setViewSettings(cacheData.viewSettings)
        }
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
        const persisted = unboxGuiState(guiState?.config || {})
        if (typeof persisted.pointSize === 'number') setPointSize(persisted.pointSize)
        if (typeof persisted.maxPoints === 'number' && persisted.eventCountVersion === EVENT_COUNT_VERSION) setMaxPoints(normalizeEventCount(persisted.maxPoints))
        if (!embedded && typeof persisted.darkMode === 'boolean') setDarkMode(persisted.darkMode)
        if (typeof persisted.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(persisted.histogramBins))
        if (typeof persisted.histogramTransform === 'string') setHistogramTransform(normalizeHistogramTransform(persisted.histogramTransform))
        if (typeof persisted.sidebarWidth === 'number' && Number.isFinite(persisted.sidebarWidth)) {
          setSidebarWidth(Math.min(380, Math.max(160, persisted.sidebarWidth)))
        }
        if (typeof persisted.sidebarCollapsed === 'boolean') setSidebarCollapsed(persisted.sidebarCollapsed)
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
          config_json: { pointSize, maxPoints: normalizeEventCount(maxPoints), eventCountVersion: EVENT_COUNT_VERSION, histogramBins, histogramTransform, ...(!embedded ? { darkMode } : {}), axisSettings, axisSettingsVersion: AXIS_SETTINGS_VERSION, sidebarWidth, sidebarCollapsed }
        })
      }).catch(() => {})
    }, 350)
    return () => clearTimeout(timer)
  }, [pointSize, maxPoints, histogramBins, histogramTransform, darkMode, axisSettings, sidebarWidth, sidebarCollapsed, guiStateLoaded, embedded])

  // Synchronize frontend gates and settings to backend in-memory cache
  useEffect(() => {
    if (!gatesLoaded) return
    const timer = setTimeout(() => {
      useApi('/gate_cache', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ gates, pointSize, maxPoints: normalizeEventCount(maxPoints), histogramBins, histogramTransform, viewSettings, eventCountVersion: EVENT_COUNT_VERSION })
      }).catch(() => {})
    }, 400)
    return () => clearTimeout(timer)
  }, [gates, pointSize, maxPoints, histogramBins, histogramTransform, viewSettings, gatesLoaded])

  useEffect(() => {
    if (!gatesLoaded || !files.length) return
    const timer = setTimeout(() => {
      useApi('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows: buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform, viewSettings }) }),
      }).catch(() => {})
    }, 700)
    return () => clearTimeout(timer)
  }, [gates, files, pointSize, maxPoints, histogramBins, histogramTransform, viewSettings, gatesLoaded])

  const selectedSpectrumFile = files.find((file) => file.filename === selected) || {}
  const selectedSpectrumState = spectrumGateState(gates, selectedSpectrumFile)
  const spectrumEligible = Boolean(selected) && selectedSpectrumState.eligible
  const spectrumUsesHistogramGate = selectedSpectrumState.usesHistogram
  const spectrumCacheKey = spectrumCacheKeyForFile(gates, selectedSpectrumFile)
  const cachedSpectrum = spectrumCacheKey ? spectrumCache[spectrumCacheKey] : null
  const allSpectraEligible = gatesLoaded && preloadComplete && files.length > 0 &&
    files.every((file) => spectrumGateState(gates, file).eligible)
  const missingSpectrumFiles = useMemo(() => (
    allSpectraEligible
      ? files.filter((file) => !spectrumCache[spectrumCacheKeyForFile(gates, file)])
      : []
  ), [allSpectraEligible, files, gates, spectrumCache])
  const spectrumBatchKey = allSpectraEligible
    ? files.map((file) => spectrumCacheKeyForFile(gates, file)).join('|')
    : ''

  useEffect(() => {
    if (!gatesLoaded || !preloadComplete || !spectrumEligible) {
      setSpectrum(null)
      return undefined
    }
    if (cachedSpectrum) {
      setSpectrum(cachedSpectrum)
      return
    }
    setSpectrum(null)
    if (allSpectraEligible) return undefined
    const controller = new AbortController()
    let disposed = false
    const timer = setTimeout(() => {
      useApi('/gate_spectra', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filenames: [selected], gates }),
        signal: controller.signal,
      })
        .then((data) => {
          if (disposed) return
          const spectra = unboxGuiState(data?.spectra || {})
          const spectrumValue = spectra?.[selected]
          const spectrumData = decodeSpectrumData(spectrumValue)
          setSpectrum(spectrumData)
          if (spectrumData) {
            setSpectrumCache((prev) => ({ ...prev, [spectrumCacheKey]: spectrumData }))
          }
        })
        .catch((error) => {
          if (!disposed && error?.name !== 'AbortError') setSpectrum(null)
        })
    }, 800)
    return () => {
      disposed = true
      clearTimeout(timer)
      controller.abort()
    }
  }, [selected, gates, gatesLoaded, preloadComplete, spectrumEligible, spectrumCacheKey, cachedSpectrum, allSpectraEligible])

  useEffect(() => {
    if (!allSpectraEligible || !missingSpectrumFiles.length || !spectrumBatchKey) return undefined
    if (spectrumBatchRef.current === spectrumBatchKey) return undefined
    spectrumBatchRef.current = spectrumBatchKey
    const controller = new AbortController()
    let disposed = false
    setSpectraPrecomputing(true)
    setStatus('Precomputing spectra')
    useApi('/gate_spectra', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        filenames: missingSpectrumFiles.map((file) => file.filename),
        gates,
      }),
      signal: controller.signal,
    })
      .then((data) => {
        if (disposed) return
        if (data?.success === false) throw new Error(data.error || 'Spectrum precomputation failed')
        const spectra = unboxGuiState(data?.spectra || {})
        setSpectrumCache((previous) => {
          const next = { ...previous }
          missingSpectrumFiles.forEach((file) => {
            const spectrumData = decodeSpectrumData(spectra?.[file.filename])
            const key = spectrumCacheKeyForFile(gates, file)
            if (key && spectrumData) next[key] = spectrumData
          })
          return next
        })
        setStatus('Ready')
      })
      .catch((error) => {
        if (disposed || error?.name === 'AbortError') return
        spectrumBatchRef.current = ''
        setStatus(`Could not precompute spectra: ${error.message}`)
      })
      .finally(() => {
        if (!disposed) setSpectraPrecomputing(false)
      })
    return () => {
      disposed = true
      controller.abort()
    }
  }, [allSpectraEligible, missingSpectrumFiles, spectrumBatchKey, gates])

  useEffect(() => {
    if (!files.length) return undefined
    const requestedPoints = normalizeEventCount(maxPoints)
    const controller = new AbortController()
    let disposed = false
    setPayload(null)
    setPayloadCache({})
    payloadCacheRef.current = {}
    setPreloadComplete(false)
    setStatus('Preloading controls')

    useApi(`/gate_preload_compact?max_points=${requestedPoints}`, {
      signal: controller.signal,
    })
      .then((data) => {
        if (disposed) return
        const next = {}
        ;(data?.payloads || []).forEach((item) => {
          const decoded = decodeCompactPayload(item)
          const filename = decoded?.file?.filename || decoded?.filename
          if (filename && !decoded?.error) next[filename] = decoded
        })
        payloadCacheRef.current = next
        setPayloadCache(next)
        setPreloadComplete(true)
        setStatus('Ready')
      })
      .catch((error) => {
        if (disposed || error?.name === 'AbortError') return
        setPreloadComplete(true)
        setStatus(`Could not preload controls: ${error.message}`)
      })

    return () => {
      disposed = true
      controller.abort()
    }
  }, [files, maxPoints])

  useEffect(() => {
    if (!preloadComplete || !selected) return
    const nextPayload = payloadCacheRef.current[selected] || null
    setPayload(nextPayload)
    setDraft([])
    setDrawActive(false)
    setStatus(nextPayload ? 'Ready' : `Could not find preloaded events for ${selected}`)
  }, [selected, preloadComplete])

  const currentFile = payload?.file || files.find((f) => f.filename === selected) || {}
  const events = useMemo(() => materializePayloadEvents(payload, maxPoints), [payload, maxPoints])
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
  const updateViewSetting = (plot, target, view) => {
    setViewSettings((previous) => {
      const next = {
        ...previous,
        [plot]: { ...(previous[plot] || {}) },
      }
      const normalized = normalizePlotView(view, plot !== 'histogram')
      if (normalized) next[plot][target] = normalized
      else delete next[plot][target]
      return next
    })
  }
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
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform, viewSettings })
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
        if (embedded) {
          onRequestExit?.()
          return
        }
        await useApi('/gate_shutdown', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{}' })
        window.close()
      }
    } catch (err) {
      setStatus(`Could not save gate CSV: ${err.message}`)
    }
  }

  async function saveConfigAsCsv() {
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform, viewSettings })
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
    setViewSettings(parsed.viewSettings)
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
      const generatedGates = unboxGuiState(result.gates || {})
      const nextGates = pruneUnavailableNegativeGates({ ...gates, ...generatedGates }, files)
      const returnedSpectra = unboxGuiState(result.spectra || {})
      setGates(nextGates)
      setSpectrumCache((previous) => {
        const next = { ...previous }
        files.forEach((file) => {
          const spectrumData = decodeSpectrumData(returnedSpectra?.[file.filename])
          const key = spectrumCacheKeyForFile(nextGates, file)
          if (key && spectrumData) next[key] = spectrumData
        })
        return next
      })
      const generated = Number(result.gates_generated ?? Object.keys(generatedGates).length)
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

  const beginSidebarResize = (event) => {
    if (sidebarCollapsed) return
    event.preventDefault()
    const startX = event.clientX
    const startWidth = sidebarWidth
    const previousCursor = document.body.style.cursor
    const previousUserSelect = document.body.style.userSelect
    document.body.style.cursor = 'col-resize'
    document.body.style.userSelect = 'none'

    const move = (moveEvent) => {
      setSidebarWidth(Math.min(380, Math.max(160, startWidth + moveEvent.clientX - startX)))
    }
    const finish = () => {
      window.removeEventListener('pointermove', move)
      window.removeEventListener('pointerup', finish)
      window.removeEventListener('pointercancel', finish)
      document.body.style.cursor = previousCursor
      document.body.style.userSelect = previousUserSelect
    }
    window.addEventListener('pointermove', move)
    window.addEventListener('pointerup', finish)
    window.addEventListener('pointercancel', finish)
  }

  return (
    <main
      className={`app-shell ${darkMode ? 'dark' : 'light'} ${initialLoading || histogramAutogating ? 'is-initial-loading' : ''}`}
    >
      <aside
        className={`sidebar gating-sidebar ${sidebarCollapsed ? 'is-collapsed' : ''}`}
        style={{ '--gating-sidebar-width': `${sidebarWidth}px` }}
      >
        <button
          type="button"
          className="gating-sidebar-toggle"
          title={sidebarCollapsed ? 'Show control sidebar' : 'Hide control sidebar'}
          aria-label={sidebarCollapsed ? 'Show control sidebar' : 'Hide control sidebar'}
          onClick={() => setSidebarCollapsed((collapsed) => !collapsed)}
        >
          {sidebarCollapsed ? <PanelLeftOpen size={14} /> : <PanelLeftClose size={14} />}
        </button>
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
              <div className="file-row-top">
                <span>{file.fluorophore}</span>
                <small>{file.channel}{fileUsesHistogramGates(file) ? '' : ' · AF'}</small>
              </div>
              <strong>{file.marker}</strong>
            </button>
          ))}
        </div>
      </aside>

      {!sidebarCollapsed && (
        <div
          className="gating-sidebar-resizer"
          role="separator"
          aria-label="Resize control sidebar"
          aria-orientation="vertical"
          onPointerDown={beginSidebarResize}
        />
      )}

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
              <HistogramSparkleIcon />
            </button>
            <div className="gating-help">
              <button type="button" className="icon-button gating-help-button" aria-label="Gating help" aria-describedby="gating-help-popover">
                <Info size={17} />
              </button>
              <div id="gating-help-popover" className="gating-help-popover" role="tooltip">
                <strong>Quick guide</strong>
                <dl>
                  <div><dt>Navigate</dt><dd>Use the arrows or select a control in the sidebar.</dd></div>
                  <div><dt>Scatter gates</dt><dd>Press Gate, click polygon points, then click the first point or double-click to close. Drag the border to move the gate or its nodes to reshape it.</dd></div>
                  <div><dt>Gate scope</dt><dd>Cell and bead gates are global for their control type. Right-click a gate to make it file-specific, make it global, or clear it. A gate remains tied to the axis channels on which it was created. After changing a plot axis, create a new gate; moving the old gate does not convert it to the new channels.</dd></div>
                  <div><dt>Histogram gates</dt><dd>Use Pos or Neg and click twice to define an interval. Neg is disabled when a matched external AF background exists. Clear removes histogram gates only for this SCC file.</dd></div>
                  <div><dt>Auto-gate</dt><dd>The starred histogram button creates only missing required histogram gates and preserves existing gates.</dd></div>
                  <div><dt>Plot view</dt><dd>Use the mouse wheel or Ctrl +/− to zoom scatter plots, then drag to pan. Histograms pan horizontally without zoom. Ctrl+0 resets the current view. Views are saved in the gate CSV.</dd></div>
                  <div><dt>Axes and files</dt><dd>Click an underlined scatter-axis label to change its channel. Save and Load write or read the gate CSV; Confirm validates required gates and {embedded ? 'returns to the cockpit' : 'exits'}.</dd></div>
                </dl>
              </div>
            </div>
            <div className="gating-settings">
              <button
                className="icon-button"
                title="Settings"
                aria-label="Settings"
                aria-expanded={showSettingsModal}
                onClick={() => setShowSettingsModal((open) => !open)}
              >
                <Settings size={17} />
              </button>
              {showSettingsModal && (
                <div className="gating-settings-popover" role="dialog" aria-label="Gating settings">
                  <div className="gating-settings-heading">
                    <strong>Settings</strong>
                    <span>Changes are saved automatically</span>
                  </div>
                  {!embedded && <div className="gating-settings-appearance">
                    <span>Appearance</span>
                    <button type="button" onClick={() => setDarkMode(!darkMode)}>
                      {darkMode ? <Moon size={14} /> : <Sun size={14} />}
                      {darkMode ? 'Dark' : 'Light'}
                    </button>
                  </div>}
                  <label className="gating-settings-row">
                    <span>Point size</span>
                    <input
                      type="range"
                      min="0.5"
                      max="4.0"
                      step="0.25"
                      value={pointSize}
                      onChange={(e) => setPointSize(Number(e.target.value))}
                    />
                    <strong>{pointSize.toFixed(2)}</strong>
                  </label>
                  <label className="gating-settings-row">
                    <span>Events</span>
                    <input
                      type="range"
                      min="0"
                      max={EVENT_COUNT_STEPS.length - 1}
                      step="1"
                      value={eventStepIndex(maxPoints)}
                      onChange={(e) => setMaxPoints(EVENT_COUNT_STEPS[Number(e.target.value)])}
                    />
                    <strong>{eventStepLabel(maxPoints, payload?.total_events)}</strong>
                  </label>
                  <label className="gating-settings-row">
                    <span>Histogram bins</span>
                    <input
                      type="range"
                      min={HISTOGRAM_BIN_MIN}
                      max={HISTOGRAM_BIN_MAX}
                      step="5"
                      value={histogramBins}
                      onChange={(e) => setHistogramBins(normalizeHistogramBins(e.target.value))}
                      onInput={(e) => setHistogramBins(normalizeHistogramBins(e.currentTarget.value))}
                    />
                    <strong>{histogramBins}</strong>
                  </label>
                  <div className="gating-settings-transform">
                    <TransformDropdown
                      value={histogramTransform}
                      onChange={(nextValue) => setHistogramTransform(normalizeHistogramTransform(nextValue))}
                    />
                  </div>
                </div>
              )}
            </div>
            <input
              ref={loadInputRef}
              type="file"
              accept=".csv"
              className="hidden-file-input"
              onChange={(event) => loadConfigFromFile(event.target.files?.[0] || null)}
            />
            <button className="icon-button" aria-label="Save gates as CSV" title="Save gates as CSV" onClick={() => saveConfigAsCsv()}><Save size={17} /></button>
            <button className="icon-button" aria-label="Load gates from CSV" title="Load gates from CSV" onClick={() => loadConfigFromPicker()}><FolderOpen size={17} /></button>
            <span
              className="confirm-wrapper"
              onMouseEnter={() => !canConfirm && setShowConfirmIssues(true)}
              onMouseLeave={() => setShowConfirmIssues(false)}
              onFocus={() => !canConfirm && setShowConfirmIssues(true)}
              onBlur={() => setShowConfirmIssues(false)}
            >
              <button
                title={canConfirm ? (embedded ? 'Save and return to cockpit' : 'Export and close') : 'Finish required gates before confirming'}
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
            viewDomain={viewSettings.cell?.global || null}
            onViewDomainChange={(view) => updateViewSetting('cell', 'global', view)}
            statsText={cellGate ? `${cellSummary.count.toLocaleString()} (${cellSummary.pct.toFixed(1)}%)` : null}
            drawActive={drawActive}
            pointSize={pointSize}
            darkMode={darkMode}
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
            viewDomain={viewSettings.singlet?.global || null}
            onViewDomainChange={(view) => updateViewSetting('singlet', 'global', view)}
            statsText={singletGate ? `${singletSummary.count.toLocaleString()} (${singletSummary.pct.toFixed(1)}%)` : null}
            warningText={singletWarningText}
            drawActive={drawActive}
            pointSize={pointSize}
            darkMode={darkMode}
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
              viewDomain={viewSettings.histogram?.[selected] || null}
              onViewDomainChange={(view) => updateViewSetting('histogram', selected, view)}
              statsText={[
                positiveGate ? `Pos ${positiveSummary.count.toLocaleString()} (${positiveSummary.pct.toFixed(1)}%)` : null,
                negativeGate ? `Neg ${negativeSummary.count.toLocaleString()} (${negativeSummary.pct.toFixed(1)}%)` : null,
              ].filter(Boolean).join(' · ') || null}
              warningText={positiveWarningText}
              drawActive={drawActive}
              pointSize={pointSize}
              darkMode={darkMode}
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
                <SpectrumCanvas spectrum={spectrum} darkMode={darkMode} />
              ) : !spectrumEligible ? (
                <div className="spectrum-gate-required">
                  {spectrumUsesHistogramGate
                    ? 'Complete the cell, singlet, and positive histogram gates to create the spectrum.'
                    : 'Complete the cell and singlet gates to create the AF spectrum.'}
                </div>
              ) : (
                <div
                  className="spectrum-placeholder"
                  aria-label={spectraPrecomputing ? 'Precomputing all spectra' : 'Spectrum loading placeholder'}
                />
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
            <p>Are you sure you want to save all gate configurations and {embedded ? 'return to the cockpit' : 'exit the Gating GUI'}?</p>
            <div className="confirm-modal-actions">
              <button className="cancel-btn" onClick={() => setShowConfirmModal(false)}>Cancel</button>
              <button className="confirm-btn" disabled={!canConfirm} onClick={() => { if (canConfirm) { setShowConfirmModal(false); saveConfig(true) } }}>{embedded ? 'Confirm & Return' : 'Confirm & Exit'}</button>
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

    </main>
  )
}

export default App
