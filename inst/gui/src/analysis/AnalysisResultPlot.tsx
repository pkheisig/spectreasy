import { useEffect, useRef } from 'react'
import type { AnalysisRunResult } from './types'
import { linearScale, PLOT_MARGIN, robustExtent } from './geometry'

const clusterPalette = ['#54b9c5', '#f0c45d', '#e66d62', '#83c36a', '#bd7bd2', '#e89a56', '#6a91d4', '#df73a3']

export function AnalysisResultPlot({ result }: { result: AnalysisRunResult }) {
  const canvasRef = useRef<HTMLCanvasElement>(null)

  useEffect(() => {
    const canvas = canvasRef.current
    if (!canvas) return
    const width = 720
    const height = 440
    const ratio = window.devicePixelRatio || 1
    canvas.width = width * ratio
    canvas.height = height * ratio
    canvas.style.width = '100%'
    canvas.style.aspectRatio = `${width}/${height}`
    const context = canvas.getContext('2d')
    if (!context) return
    context.setTransform(ratio, 0, 0, ratio, 0, 0)
    context.fillStyle = '#07131d'
    context.fillRect(0, 0, width, height)
    const xDomain = robustExtent(result.events.map((event) => event.dimension_1))
    const yDomain = robustExtent(result.events.map((event) => event.dimension_2))
    const x = linearScale(xDomain, [PLOT_MARGIN.left, width - PLOT_MARGIN.right])
    const y = linearScale(yDomain, [height - PLOT_MARGIN.bottom, PLOT_MARGIN.top])
    const pseudotimes = result.events.map((event) => event.pseudotime).filter((value): value is number => Number.isFinite(value))
    const low = Math.min(...pseudotimes, 0)
    const high = Math.max(...pseudotimes, 1)
    for (const event of result.events) {
      if (Number.isFinite(event.pseudotime)) {
        const normalized = ((event.pseudotime ?? low) - low) / (high - low || 1)
        context.fillStyle = `hsl(${210 - normalized * 185} 72% ${52 + normalized * 10}%)`
      } else if (Number.isFinite(event.cluster_id)) {
        context.fillStyle = clusterPalette[((event.cluster_id ?? 1) - 1) % clusterPalette.length]
      } else {
        context.fillStyle = '#60b8c8'
      }
      context.globalAlpha = 0.8
      context.beginPath()
      context.arc(x.toScreen(event.dimension_1), y.toScreen(event.dimension_2), 1.7, 0, Math.PI * 2)
      context.fill()
    }
    context.globalAlpha = 1
    context.strokeStyle = '#496475'
    context.strokeRect(PLOT_MARGIN.left, PLOT_MARGIN.top, width - PLOT_MARGIN.left - PLOT_MARGIN.right, height - PLOT_MARGIN.top - PLOT_MARGIN.bottom)
    context.fillStyle = '#a9bdc8'
    context.font = '11px "Avenir Next", sans-serif'
    context.fillText('Dimension 1', width / 2 - 28, height - 15)
    context.save()
    context.translate(16, height / 2 + 25)
    context.rotate(-Math.PI / 2)
    context.fillText('Dimension 2', 0, 0)
    context.restore()
  }, [result])

  return (
    <div className="analysis-result-plot">
      <canvas ref={canvasRef} aria-label={`${result.metadata.method.name} result plot`} />
    </div>
  )
}
