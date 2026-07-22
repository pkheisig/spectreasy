import { useEffect, useMemo, useRef, useState } from 'react'
import {
  dragHistogramGateInDisplaySpace,
  histogramDomainIncludingGates,
  padHistogramDomain,
  panHistogramDomain,
} from '../histogramGates.js'
import { normalizePlotView } from '../gatingViewSettings.js'
import GatePlotView from './GatePlotView.jsx'
import {
  DARK_DENSITY_PALETTE,
  DEFAULT_HISTOGRAM_BINS,
  DEFAULT_HISTOGRAM_TRANSFORM,
  DENSITY_PALETTE,
  HISTOGRAM_BIN_MAX,
  HISTOGRAM_BIN_MIN,
  PAD,
  PLOT_HEIGHT,
  PLOT_WIDTH,
  TransformDropdown,
  axisLabel,
  channelTitle,
  clamp,
  computeDensityBuckets,
  extent,
  formatTickValue,
  gateKey,
  getDensityCurve,
  getTicks,
  histogramGateKey,
  histogramTransformFns,
  makeScale,
  pointInPolygon,
} from './GatingCore.jsx'

export default function GatePlot({
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

  const histogramValues = useMemo(
    () => mode === 'histogram' ? events.map((event) => event[xField]) : [],
    [events, mode, xField],
  )
  const transformFns = useMemo(
    () => histogramTransformFns(histogramTransform, histogramValues),
    [histogramTransform, histogramValues],
  )
  const toPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.forward(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }
  const fromPlotX = (value) => {
    const out = mode === 'histogram' ? transformFns.inverse(value) : Number(value)
    return Number.isFinite(out) ? out : 0
  }

  const baseXDomain = useMemo(() => {
    const eventValues = histogramValues
    const rawDomain = xDomainProp || extent(eventValues)
    if (mode !== 'histogram') return rawDomain
    const visibleDomain = histogramDomainIncludingGates(
      rawDomain,
      eventValues,
      [gate, ...secondaryGates],
    )
    const transformed = visibleDomain.map((value) => transformFns.forward(value)).filter(Number.isFinite)
    if (transformed.length < 2 || transformed[0] === transformed[1]) return [0, 1]
    return padHistogramDomain([Math.min(...transformed), Math.max(...transformed)])
  }, [histogramValues, xDomainProp, mode, transformFns, gate, secondaryGates])
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
      rawXDomain: mode === 'histogram' ? xDomain.map((value) => fromPlotX(value)) : null,
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

      if (mode === 'histogram') {
        const baseRawDomain = baseXDomain.map((value) => fromPlotX(value))
        const baseRawSpan = baseRawDomain[1] - baseRawDomain[0]
        const extendedRawLimits = [
          baseRawDomain[0] - baseRawSpan * 4,
          baseRawDomain[1] + baseRawSpan * 4,
        ]
        const nextRawX = panHistogramDomain(
          panTarget.rawXDomain,
          dxPixels,
          plotWidth,
          extendedRawLimits,
        )
        const nextX = nextRawX.map((value) => toPlotX(value))
        setZoomXDomain(nextX)
        pendingViewRef.current = { x: nextX, y: null }
      } else {
        const startXSpan = panTarget.xDomain[1] - panTarget.xDomain[0]
        const xShift = -(dxPixels / plotWidth) * startXSpan
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
    const displayDx = currentPt.plotX - (lastMousePos.plotX ?? toPlotX(lastMousePos.x))
    setLastMousePos(currentPt)

    const dragMode = dragTarget.gate.mode
    const dragKind = dragTarget.target
    const is1D = dragMode === 'separator' || dragMode === 'positive_1d' || dragMode === 'negative_1d'
    const sourceVertices = dragPreviewVertices || dragTarget.gate.vertices

    if (dragKind === 'gate' || dragKind === 'separator') {
      if (is1D) {
        setDragPreviewVertices(
          dragHistogramGateInDisplaySpace(sourceVertices, displayDx, transformFns),
        )
        return
      }
      const newVertices = sourceVertices.map((v) => {
        return {
          x: v.x + dx,
          y: v.y + dy,
        }
      })
      setDragPreviewVertices(newVertices)
    } else if (typeof dragKind === 'number') {
      if (is1D) {
        setDragPreviewVertices(
          dragHistogramGateInDisplaySpace(sourceVertices, displayDx, transformFns, dragKind),
        )
        return
      }
      const newVertices = sourceVertices.map((v) => ({ ...v }))
      const newPt = {
        x: sourceVertices[dragKind].x + dx,
        y: sourceVertices[dragKind].y + dy,
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

  return <GatePlotView view={{
    active,
    onClickPlot,
    title,
    subtitle,
    events,
    statsText,
    warningText,
    isPositivePlot,
    drawActive,
    onToggleDraw,
    histogramGateType,
    negativeGateEnabled,
    onToggleHistogramGate,
    onClear,
    plotRef,
    mode,
    canvasRef,
    svgRef,
    xDomain,
    yDomain,
    panTarget,
    zoomXDomain,
    zoomYDomain,
    handlePlotMouseDown,
    handleClick,
    onFinish,
    handleContext,
    handleMouseMove,
    handleMouseUp,
    setPlotHovered,
    densityPath,
    labels,
    xChannel,
    yChannel,
    openAxisMenu,
    xTicks,
    xScale,
    fromPlotX,
    yTicks,
    yScale,
    secondaryGates,
    renderSecondaryGate,
    displayGate,
    xToScreen,
    handleMouseDown,
    isBead,
    hoverPt,
    path,
    gateClass,
    draft,
    gate,
    axisMenu,
    menuRef,
    availableChannels,
    chooseAxis
  }} />
}
