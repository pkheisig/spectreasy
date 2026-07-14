const SPECTRUM_PALETTE = ['#0000ff', '#00ffff', '#00ff00', '#ffff00', '#ff0000']

function asArray(value) {
  if (value === undefined || value === null) return []
  return Array.isArray(value) ? value : [value]
}

export function decodeSpectrumData(value) {
  if (!value || value.format !== 'spectrum-histogram-v1') return null
  const compact = value.counts
  if (!compact || compact.format !== 'uint32-column-major') return null
  if (compact.values instanceof Uint32Array) return value
  if (typeof compact.data !== 'string') return null
  const binary = globalThis.atob(compact.data.replace(/\s/g, ''))
  const bytes = new Uint8Array(binary.length)
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i)
  const rows = Number(compact.rows) || 0
  const columns = Number(compact.columns) || 0
  if (!rows || !columns || bytes.byteLength !== rows * columns * 4) return null
  const view = new DataView(bytes.buffer)
  const values = new Uint32Array(rows * columns)
  for (let i = 0; i < values.length; i++) values[i] = view.getUint32(i * 4, true)
  return {
    ...value,
    channels: asArray(value.channels).map(String),
    labels: asArray(value.labels).map(String),
    bin_mid: asArray(value.bin_mid).map(Number),
    fill_limits: asArray(value.fill_limits).map(Number),
    counts: { ...compact, rows, columns, values },
  }
}

function interpolateChannel(start, end, amount) {
  return Math.round(start + (end - start) * amount)
}

export function spectrumColor(amount) {
  const scaled = Math.max(0, Math.min(1, amount)) * (SPECTRUM_PALETTE.length - 1)
  const index = Math.min(SPECTRUM_PALETTE.length - 2, Math.floor(scaled))
  const fraction = scaled - index
  const from = SPECTRUM_PALETTE[index]
  const to = SPECTRUM_PALETTE[index + 1]
  const channel = (hex, offset) => Number.parseInt(hex.slice(offset, offset + 2), 16)
  const red = interpolateChannel(channel(from, 1), channel(to, 1), fraction)
  const green = interpolateChannel(channel(from, 3), channel(to, 3), fraction)
  const blue = interpolateChannel(channel(from, 5), channel(to, 5), fraction)
  return `rgb(${red}, ${green}, ${blue})`
}
