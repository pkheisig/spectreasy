import { Eraser, Hexagon } from 'lucide-react'
import {
  HISTOGRAM_BIN_MAX,
  HISTOGRAM_BIN_MIN,
  PAD,
  PLOT_HEIGHT,
  PLOT_WIDTH,
  TransformDropdown,
  axisLabel,
  channelTitle,
  formatTickValue,
} from './GatingCore.jsx'

export default function GatePlotView({ view }) {
  const {
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
    histogramBins,
    onHistogramBinsChange,
    histogramTransform,
    onHistogramTransformChange,
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
  } = view

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
              fill={histogramGateType === 'positive' ? 'rgba(214, 82, 56, 0.1)' : 'var(--negative-gate-preview)'}
              stroke={histogramGateType === 'positive' ? '#d65238' : 'var(--negative-gate)'}
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
