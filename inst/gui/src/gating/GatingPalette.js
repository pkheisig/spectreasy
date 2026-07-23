export const DARK_PLOT_BACKGROUND = '#2c3734'

function densityPalette(size = 64) {
  const jetColors = (v) => {
    const clamp01 = (value) => Math.min(1, Math.max(0, value))
    const r = clamp01(Math.min(4 * v - 1.5, -4 * v + 4.5))
    const g = clamp01(Math.min(4 * v - 0.5, -4 * v + 3.5))
    // Lift only the sparse-event end of the jet scale. Above the first blue
    // segment the original cyan -> green -> yellow -> red progression is unchanged.
    const b = clamp01(Math.min(4 * v + 0.85, -4 * v + 2.5))
    return `rgba(${Math.floor(r * 255)}, ${Math.floor(g * 255)}, ${Math.floor(b * 255)}, 0.7)`
  }
  return Array.from({ length: size }, (_, i) => jetColors(i / Math.max(size - 1, 1)))
}

export const DENSITY_PALETTE = densityPalette()
