import { useEffect, useRef, useState } from 'react'
import {
  CheckCircle2,
  ChevronLeft,
  ChevronRight,
  Eraser,
  FileDown,
  FolderOpen,
  Info,
  Moon,
  PanelLeftClose,
  PanelLeftOpen,
  Save,
  ScatterChart,
  Settings,
  Sun,
  Upload,
  X,
} from 'lucide-react'
import { fileUsesHistogramGates } from '../gatingEligibility.js'
import SpectrumCanvas from '../SpectrumCanvas.jsx'
import GatePlot from './GatePlot.jsx'
import {
  EVENT_COUNT_STEPS,
  HISTOGRAM_BIN_MAX,
  HISTOGRAM_BIN_MIN,
  HistogramSparkleIcon,
  TransformDropdown,
  channelTitle,
  eventStepIndex,
  eventStepLabel,
  gateKey,
  histogramGateKey,
  normalizeHistogramBins,
  normalizeHistogramTransform,
} from './GatingCore.jsx'

function CommitRange({ value, onCommit, formatValue, ...inputProps }) {
  const [draftValue, setDraftValue] = useState(Number(value))
  const interactingRef = useRef(false)

  useEffect(() => {
    if (!interactingRef.current) setDraftValue(Number(value))
  }, [value])

  const commit = (nextValue) => {
    const numericValue = Number(nextValue)
    interactingRef.current = false
    setDraftValue(numericValue)
    if (numericValue !== Number(value)) onCommit(numericValue)
  }

  return <>
    <input
      {...inputProps}
      type="range"
      value={draftValue}
      onPointerDown={(event) => {
        interactingRef.current = true
        event.currentTarget.setPointerCapture?.(event.pointerId)
      }}
      onChange={(event) => setDraftValue(Number(event.currentTarget.value))}
      onPointerUp={(event) => {
        commit(event.currentTarget.value)
        if (event.currentTarget.hasPointerCapture?.(event.pointerId)) {
          event.currentTarget.releasePointerCapture(event.pointerId)
        }
      }}
      onPointerCancel={(event) => commit(event.currentTarget.value)}
      onKeyUp={(event) => commit(event.currentTarget.value)}
      onBlur={(event) => commit(event.currentTarget.value)}
    />
    <strong>{formatValue(draftValue)}</strong>
  </>
}

export default function GatingWorkspace({ view }) {
  const {
    darkMode,
    initialLoading,
    histogramAutogating,
    sidebarCollapsed,
    sidebarWidth,
    setSidebarCollapsed,
    files,
    selected,
    setSelected,
    beginSidebarResize,
    currentFile,
    payload,
    controlType,
    canGoPrevious,
    canGoNext,
    selectRelativeFile,
    canAutogateHistograms,
    histogramAutogateMissing,
    setShowAutogateConfirmModal,
    embedded,
    onRequestClose,
    showSettingsModal,
    setShowSettingsModal,
    setDarkMode,
    pointSize,
    setPointSize,
    maxPoints,
    setMaxPoints,
    histogramBins,
    setHistogramBins,
    histogramTransform,
    setHistogramTransform,
    loadInputRef,
    loadConfigFromFile,
    saveConfigAsCsv,
    loadConfigFromPicker,
    canConfirm,
    setShowConfirmIssues,
    showConfirmIssues,
    confirmIssueLines,
    setShowConfirmModal,
    status,
    gatesLoaded,
    isBead,
    cellAxes,
    singletAxes,
    events,
    labels,
    activeGate,
    cellGate,
    draft,
    setDraft,
    setGates,
    finishDraft,
    setActiveGate,
    setDrawActive,
    updateGateVertex,
    updateGateVertices,
    handleDragEnd,
    domainForChannel,
    domains,
    viewSettings,
    updateViewSetting,
    cellSummary,
    drawActive,
    setContextMenu,
    availableScatterChannels,
    updateAxis,
    clearActiveGate,
    cellsFilteredEvents,
    singletGate,
    singletSummary,
    singletWarningText,
    singletsFilteredEvents,
    activeHistogramGate,
    channels,
    histogramGateType,
    setHistogramGateType,
    negativeGate,
    positiveGate,
    positiveSummary,
    negativeSummary,
    positiveWarningText,
    usesNegativeHistogramGate,
    secondaryHistogramGates,
    clearSelectedHistogramGates,
    spectrum,
    spectrumEligible,
    spectrumUsesHistogramGate,
    spectraPrecomputing,
    contextMenu,
    gates,
    makeFileSpecific,
    useGlobalGate,
    clearContextGate,
    showConfirmModal,
    saveConfig,
    showAutogateConfirmModal,
    histogramAutogateTargets,
    histogramAutogateNegativeTargets,
    autogateHistograms
  } = view

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
                    <CommitRange
                      min="0.5"
                      max="4.0"
                      step="0.25"
                      value={pointSize}
                      onCommit={setPointSize}
                      formatValue={(value) => value.toFixed(2)}
                    />
                  </label>
                  <label className="gating-settings-row">
                    <span>Events</span>
                    <CommitRange
                      min="0"
                      max={EVENT_COUNT_STEPS.length - 1}
                      step="1"
                      value={eventStepIndex(maxPoints)}
                      onCommit={(index) => setMaxPoints(EVENT_COUNT_STEPS[index])}
                      formatValue={(index) => eventStepLabel(EVENT_COUNT_STEPS[index], payload?.total_events)}
                    />
                  </label>
                  <label className="gating-settings-row">
                    <span>Histogram bins</span>
                    <CommitRange
                      min={HISTOGRAM_BIN_MIN}
                      max={HISTOGRAM_BIN_MAX}
                      step="5"
                      value={histogramBins}
                      onCommit={(value) => setHistogramBins(normalizeHistogramBins(value))}
                      formatValue={normalizeHistogramBins}
                    />
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
            {embedded && onRequestClose && (
              <button
                type="button"
                className="icon-button applet-close-button"
                onClick={onRequestClose}
                aria-label="Close control gating and return to cockpit"
                autoFocus
              >
                <X size={16} /> Close
              </button>
            )}
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
