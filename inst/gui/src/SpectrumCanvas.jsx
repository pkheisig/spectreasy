import React, { memo, useEffect, useRef } from 'react'
import { spectrumColor } from './spectrumData.js'

function drawSpectrum(canvas, spectrum, darkMode) {
  const width = canvas.clientWidth
  if (!width) return
  const height = Math.max(280, width * (320 / 1150))
  const ratio = Math.min(window.devicePixelRatio || 1, 2)
  canvas.width = Math.round(width * ratio)
  canvas.height = Math.round(height * ratio)
  const context = canvas.getContext('2d')
  if (!context) return
  context.setTransform(ratio, 0, 0, ratio, 0, 0)

  const foreground = darkMode ? '#d9e1de' : '#333333'
  const muted = darkMode ? '#85918d' : '#6c7773'
  const background = darkMode ? '#0b1110' : '#ffffff'
  const grid = darkMode ? 'rgba(217, 225, 222, 0.12)' : 'rgba(51, 51, 51, 0.12)'
  const margins = { left: 64, right: 12, top: 10, bottom: Math.max(92, height * 0.27) }
  const plotWidth = Math.max(1, width - margins.left - margins.right)
  const plotHeight = Math.max(1, height - margins.top - margins.bottom)
  const rows = spectrum.counts.rows
  const columns = spectrum.counts.columns
  const values = spectrum.counts.values
  const binMid = spectrum.bin_mid
  const labels = spectrum.labels
  const yPower = Number(spectrum.y_power) || 1.5
  const maxY = Math.max(1, Number(spectrum.max_y) || 6)
  const transformedMax = Math.pow(maxY + 0.5, yPower)
  const minCount = Math.max(1, Number(spectrum.min_bin_count) || 1)
  const fillLow = Number(spectrum.fill_limits[0]) || 0
  const fillHigh = Math.max(fillLow + Number.EPSILON, Number(spectrum.fill_limits[1]) || 1)
  const binHeight = Number(spectrum.bin_height) || 0.04
  const tileHeight = Math.max(1, (binHeight * 3 / transformedMax) * plotHeight)
  const columnWidth = plotWidth / Math.max(columns, 1)

  context.fillStyle = background
  context.fillRect(0, 0, width, height)
  context.strokeStyle = grid
  context.lineWidth = 1
  for (let tick = 0; tick <= Math.ceil(maxY); tick++) {
    const y = margins.top + plotHeight - (Math.pow(tick, yPower) / transformedMax) * plotHeight
    context.beginPath()
    context.moveTo(margins.left, y)
    context.lineTo(margins.left + plotWidth, y)
    context.stroke()
  }

  for (let column = 0; column < columns; column++) {
    const x = margins.left + column * columnWidth + columnWidth * 0.15
    const tileWidth = Math.max(1, columnWidth * 0.7)
    const offset = column * rows
    for (let row = 0; row < rows; row++) {
      const count = values[offset + row]
      const intensity = binMid[row]
      if (count < minCount || !Number.isFinite(intensity) || intensity < 0) continue
      const transformed = Math.pow(intensity, yPower)
      const y = margins.top + plotHeight - (transformed / transformedMax) * plotHeight
      const fill = Math.log10(count + 1)
      context.fillStyle = spectrumColor((fill - fillLow) / (fillHigh - fillLow))
      context.fillRect(x, y - tileHeight / 2, tileWidth, tileHeight)
    }
  }

  context.strokeStyle = muted
  context.beginPath()
  context.moveTo(margins.left, margins.top)
  context.lineTo(margins.left, margins.top + plotHeight)
  context.lineTo(margins.left + plotWidth, margins.top + plotHeight)
  context.stroke()

  context.fillStyle = foreground
  context.font = '10px system-ui, sans-serif'
  context.textAlign = 'right'
  context.textBaseline = 'middle'
  for (let tick = 0; tick <= Math.ceil(maxY); tick++) {
    const y = margins.top + plotHeight - (Math.pow(tick, yPower) / transformedMax) * plotHeight
    context.fillText(`10^${tick}`, margins.left - 7, y)
  }
  context.save()
  context.translate(16, margins.top + plotHeight / 2)
  context.rotate(-Math.PI / 2)
  context.textAlign = 'center'
  context.font = '12px system-ui, sans-serif'
  context.fillText('Intensity', 0, 0)
  context.restore()

  context.fillStyle = foreground
  context.font = `${Math.max(5, Math.min(7, columnWidth * 0.72))}px system-ui, sans-serif`
  context.textAlign = 'right'
  context.textBaseline = 'middle'
  for (let column = 0; column < columns; column++) {
    const x = margins.left + (column + 0.5) * columnWidth
    const y = margins.top + plotHeight + 5
    context.save()
    context.translate(x, y)
    context.rotate(-Math.PI / 2)
    context.fillText(labels[column] || spectrum.channels[column] || '', 0, 0)
    context.restore()
  }
}

const SpectrumCanvas = memo(function SpectrumCanvas({ spectrum, darkMode }) {
  const canvasRef = useRef(null)

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas || !spectrum) return undefined
    const redraw = () => drawSpectrum(canvas, spectrum, darkMode)
    redraw()
    const observer = new ResizeObserver(redraw)
    observer.observe(canvas)
    return () => observer.disconnect()
  }, [spectrum, darkMode])

  return (
    <canvas
      ref={canvasRef}
      className="spectrum-canvas"
      role="img"
      aria-label={`Detector spectrum from ${Number(spectrum.event_count || 0).toLocaleString()} gated events`}
    >
      Detector spectrum
    </canvas>
  )
})

export default SpectrumCanvas
