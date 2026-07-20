import { useEffect, useRef, useState } from 'react'
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

function normalizeEventCount(value) {
  const numeric = Number(value)
  if (!Number.isFinite(numeric) || numeric <= 0) return DEFAULT_EVENT_COUNT
  const exact = EVENT_COUNT_STEPS.find((step) => step === numeric)
  if (exact) return exact
  return EVENT_COUNT_STEPS.reduce(
    (best, step) => Math.abs(step - numeric) < Math.abs(best - numeric) ? step : best,
    DEFAULT_EVENT_COUNT,
  )
}

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

export {
  AXIS_SETTINGS_VERSION,
  CONFIG_NAME,
  DEFAULT_EVENT_COUNT,
  DEFAULT_HISTOGRAM_BINS,
  DEFAULT_HISTOGRAM_TRANSFORM,
  EVENT_COUNT_STEPS,
  EVENT_COUNT_VERSION,
  GUI_MODULE,
  HISTOGRAM_BIN_MAX,
  HISTOGRAM_BIN_MIN,
  HistogramSparkleIcon,
  MIN_CONFIRM_EVENTS,
  PAD,
  PLOT_HEIGHT,
  PLOT_WIDTH,
  REQUIRED_GATE_CSV_COLUMNS,
  THEME_STORAGE_KEY,
  TransformDropdown,
  axisLabel,
  channelTitle,
  clamp,
  decodeCompactPayload,
  extent,
  fileControlType,
  filterPolygonEvents,
  formatTickValue,
  gateKey,
  getTicks,
  histogramGateKey,
  histogramTransformFns,
  makeScale,
  materializePayloadEvents,
  normalizeHistogramBins,
  normalizeHistogramTransform,
  normalizeEventCount,
  pointInPolygon,
  pointInPolygonValues,
  summarizeGate,
  unboxGuiState,
}
