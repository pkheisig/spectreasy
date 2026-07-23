import { DENSITY_PALETTE } from '../gating/GatingPalette.js'
import type { AnalysisRunEvent, AnalysisRunResult } from './types'

export type ContinuousPaletteId =
  | 'control-density'
  | 'viridis'
  | 'sunset'
  | 'plasma'
  | 'inferno'
  | 'magma'
  | 'cividis'
  | 'turbo'
  | 'ice-fire'
  | 'spectral'

export type ResultColorMode = 'density' | 'cluster' | 'pseudotime' | 'identity' | `marker:${string}`

export type ContinuousPalette = {
  id: ContinuousPaletteId
  name: string
  colors: string[]
}
export const CONTINUOUS_PALETTES: ContinuousPalette[] = [
  { id: 'control-density', name: 'Control density', colors: DENSITY_PALETTE },
  { id: 'viridis', name: 'Viridis', colors: ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'] },
  { id: 'sunset', name: 'Sunset', colors: ['#30123b', '#6d2f84', '#b73779', '#eb5a5a', '#f89f5b', '#f9e78a'] },
  { id: 'plasma', name: 'Plasma', colors: ['#0d0887', '#7e03a8', '#cc4778', '#f89540', '#f0f921'] },
  { id: 'inferno', name: 'Inferno', colors: ['#000004', '#420a68', '#932667', '#dd513a', '#fca50a', '#fcffa4'] },
  { id: 'magma', name: 'Magma', colors: ['#000004', '#3b0f70', '#8c2981', '#de4968', '#fe9f6d', '#fcfdbf'] },
  { id: 'cividis', name: 'Cividis', colors: ['#00204c', '#2e4a7d', '#666970', '#a48b58', '#e1c63b', '#ffea46'] },
  { id: 'turbo', name: 'Turbo', colors: ['#30123b', '#4669e8', '#1bcfd4', '#75fe5c', '#f9e721', '#f77b16', '#7a0403'] },
  { id: 'ice-fire', name: 'Ice–fire', colors: ['#00204c', '#2676b8', '#9bd7e4', '#f5f2ea', '#f4a582', '#c83b3b', '#5d0018'] },
  { id: 'spectral', name: 'Spectral', colors: ['#5e4fa2', '#3288bd', '#66c2a5', '#e6f598', '#fee08b', '#f46d43', '#9e0142'] },
]

export const CATEGORICAL_COLORS = [
  '#0c7c86', '#d39a1e', '#cf5148', '#569a43', '#8b55b1', '#ce7131',
  '#416fae', '#bd4f83', '#238b63', '#8c6d31', '#6b6ecf', '#e6550d',
]

function parseColor(color: string): [number, number, number, number] {
  const hex = /^#([0-9a-f]{6})$/i.exec(color)
  if (hex) {
    const value = Number.parseInt(hex[1], 16)
    return [(value >> 16) & 255, (value >> 8) & 255, value & 255, 1]
  }
  const rgba = /^rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)(?:\s*,\s*([\d.]+))?\s*\)$/i.exec(color)
  if (rgba) return [Number(rgba[1]), Number(rgba[2]), Number(rgba[3]), rgba[4] == null ? 1 : Number(rgba[4])]
  return [36, 127, 158, 1]
}

export function paletteColors(id: ContinuousPaletteId, steps = 64): string[] {
  const palette = CONTINUOUS_PALETTES.find((candidate) => candidate.id === id) ?? CONTINUOUS_PALETTES[0]
  if (palette.colors.length === steps) return palette.colors
  if (steps <= 1) return [palette.colors[0]]
  return Array.from({ length: steps }, (_, index) => {
    const position = (index / (steps - 1)) * (palette.colors.length - 1)
    const left = Math.floor(position)
    const right = Math.min(palette.colors.length - 1, left + 1)
    const mix = position - left
    const a = parseColor(palette.colors[left])
    const b = parseColor(palette.colors[right])
    const values = a.map((value, channel) => value + (b[channel] - value) * mix)
    return `rgba(${Math.round(values[0])}, ${Math.round(values[1])}, ${Math.round(values[2])}, ${values[3].toFixed(3)})`
  })
}

export function paletteGradient(id: ContinuousPaletteId): string {
  const palette = CONTINUOUS_PALETTES.find((candidate) => candidate.id === id) ?? CONTINUOUS_PALETTES[0]
  return `linear-gradient(90deg, ${palette.colors.join(', ')})`
}

export function resultColorModes(result: AnalysisRunResult): Array<{ id: ResultColorMode; label: string; continuous: boolean }> {
  const modes: Array<{ id: ResultColorMode; label: string; continuous: boolean }> = [
    { id: 'density', label: 'Density', continuous: true },
  ]
  if (result.events.some((event) => Number.isFinite(event.cluster_id))) {
    modes.push({ id: 'cluster', label: 'Cluster', continuous: false })
  }
  if (result.events.some((event) => Number.isFinite(event.pseudotime))) {
    modes.push({ id: 'pseudotime', label: 'Pseudotime', continuous: true })
  }
  if (result.events.some((event) => Boolean(event.predicted_identity))) {
    modes.push({ id: 'identity', label: 'Predicted identity', continuous: false })
  }
  for (const marker of result.metadata.marker_columns ?? []) {
    modes.push({ id: `marker:${marker.column}`, label: marker.marker, continuous: true })
  }
  return modes
}

export function markerValue(event: AnalysisRunEvent, mode: ResultColorMode): number {
  if (mode === 'pseudotime') return Number.isFinite(event.pseudotime) ? Number(event.pseudotime) : 0
  if (mode.startsWith('marker:')) {
    const key = mode.slice('marker:'.length) as `marker_${number}`
    const value = event[key]
    return Number.isFinite(value) ? Number(value) : 0
  }
  return 0
}

export function identityColor(result: AnalysisRunResult, label: string | undefined): string {
  if (!label || label === 'Unassigned') return '#a8afb1'
  const signatures = result.metadata.identity_annotation?.signatures ?? []
  return signatures.find((signature) => signature.name === label)?.color ?? '#607d8b'
}
