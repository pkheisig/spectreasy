import { useEffect, useMemo, useRef, useState } from 'react'
import './GatingGui.css'
import { reconcileGateCsvRows } from './gatingCsvCompatibility.js'
import {
  fileUsesHistogramGates,
  fileUsesNegativeHistogramGate,
  pruneUnavailableNegativeGates,
} from './gatingEligibility.js'
import { clearHistogramGatesForFile } from './histogramGates.js'
import { normalizePlotView } from './gatingViewSettings.js'
import { decodeSpectrumData } from './spectrumData.js'
import GatingWorkspace from './gating/GatingWorkspace.jsx'
import { useGatingConfigIO } from './gating/useGatingConfigIO.js'
import {
  GATING_SIDEBAR_MIN_WIDTH,
  autogateHistogramsAction,
  beginGatingSidebarResize,
} from './gating/GatingActions.js'
import { appletCacheKey, loadCachedAppletData } from './appletDataCache.ts'
import {
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
  HISTOGRAM_TRANSFORM_VERSION,
  HistogramSparkleIcon,
  MIN_CONFIRM_EVENTS,
  REQUIRED_GATE_CSV_COLUMNS,
  THEME_STORAGE_KEY,
  TransformDropdown,
  buildRows,
  channelTitle,
  decodeCompactPayload,
  eventStepIndex,
  eventStepLabel,
  fileControlType,
  filterPolygonEvents,
  formatConfirmIssues,
  gateIsFinalized,
  gateKey,
  gateRowsToCsv,
  hasViewSettings,
  histogramGateKey,
  materializePayloadEvents,
  normalizeEventCount,
  normalizeHistogramBins,
  normalizeHistogramTransform,
  parseConfigRows,
  parseCsvDocument,
  resolveGateForFile,
  spectrumCacheKeyForFile,
  spectrumGateState,
  summarizeGate,
  unboxGuiState,
  gatingApiRequest,
  validateFileForConfirm,
  validateGateCsvRows,
} from './gating/GatingCore.jsx'

const DEFAULT_SIDEBAR_WIDTH = 252
const SIDEBAR_WIDTH_VERSION = 2

function App({ embedded = false, cockpitTheme = null, projectPath = '', projectRevision = 'empty', initialFiles = null, initialMetadata = {}, onRequestExit = null, onRequestClose = null }) {
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
  const [sidebarWidth, setSidebarWidth] = useState(DEFAULT_SIDEBAR_WIDTH)
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false)
  const [preloadComplete, setPreloadComplete] = useState(false)
  const [gatesLoaded, setGatesLoaded] = useState(false)
  const [histogramAutogating, setHistogramAutogating] = useState(false)
  const [histogramAutogateNotice, setHistogramAutogateNotice] = useState(null)
  const [axisSettings, setAxisSettings] = useState({ cell: {}, singlet: {} })
  const [viewSettings, setViewSettings] = useState({ cell: {}, singlet: {}, histogram: {} })
  const [spectrum, setSpectrum] = useState(null)
  const [spectrumCache, setSpectrumCache] = useState({})
  const [spectraPrecomputing, setSpectraPrecomputing] = useState(false)
  const loadInputRef = useRef(null)
  const payloadCacheRef = useRef({})
  const spectrumBatchRef = useRef('')
  const bootPromiseRef = useRef(null)
  const activeProjectPath = projectPath || window.sessionStorage.getItem('spectreasy-project-path') || ''
  const fileInventoryKey = files.map((file) => file.filename).join('\u0000')

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
    if (bootPromiseRef.current) return
    const fileInventory = Array.isArray(initialFiles)
      ? Promise.resolve({ files: initialFiles, metadata: initialMetadata })
      : gatingApiRequest('/gate_files')
    bootPromiseRef.current = Promise.all([
      fileInventory,
      gatingApiRequest('/gate_configs'),
      gatingApiRequest('/gate_cache'),
      gatingApiRequest(`/gui_state?module=${encodeURIComponent(GUI_MODULE)}`),
      gatingApiRequest('/gate_config')
    ])
    bootPromiseRef.current
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
        if (typeof csvGates.histogramTransform === 'string' && csvGates.histogramTransformVersion === HISTOGRAM_TRANSFORM_VERSION) {
          setHistogramTransform(normalizeHistogramTransform(csvGates.histogramTransform))
        }
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
        if (typeof cacheData?.histogramTransform === 'string' && cacheData?.histogramTransformVersion === HISTOGRAM_TRANSFORM_VERSION) {
          setHistogramTransform(normalizeHistogramTransform(cacheData.histogramTransform))
        }
        const persisted = unboxGuiState(guiState?.config || {})
        if (typeof persisted.pointSize === 'number') setPointSize(persisted.pointSize)
        if (typeof persisted.maxPoints === 'number' && persisted.eventCountVersion === EVENT_COUNT_VERSION) setMaxPoints(normalizeEventCount(persisted.maxPoints))
        if (!embedded && typeof persisted.darkMode === 'boolean') setDarkMode(persisted.darkMode)
        if (typeof persisted.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(persisted.histogramBins))
        if (typeof persisted.histogramTransform === 'string' && persisted.histogramTransformVersion === HISTOGRAM_TRANSFORM_VERSION) {
          setHistogramTransform(normalizeHistogramTransform(persisted.histogramTransform))
        }
        if (
          persisted.sidebarWidthVersion === SIDEBAR_WIDTH_VERSION &&
          typeof persisted.sidebarWidth === 'number' &&
          Number.isFinite(persisted.sidebarWidth)
        ) {
          setSidebarWidth(Math.min(380, Math.max(GATING_SIDEBAR_MIN_WIDTH, persisted.sidebarWidth)))
        } else {
          setSidebarWidth(DEFAULT_SIDEBAR_WIDTH)
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
      gatingApiRequest('/gui_state', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          module: GUI_MODULE,
          config_json: { pointSize, maxPoints: normalizeEventCount(maxPoints), eventCountVersion: EVENT_COUNT_VERSION, histogramBins, histogramTransform, histogramTransformVersion: HISTOGRAM_TRANSFORM_VERSION, ...(!embedded ? { darkMode } : {}), axisSettings, axisSettingsVersion: AXIS_SETTINGS_VERSION, sidebarWidth, sidebarWidthVersion: SIDEBAR_WIDTH_VERSION, sidebarCollapsed }
        })
      }).catch(() => {})
    }, 350)
    return () => clearTimeout(timer)
  }, [pointSize, maxPoints, histogramBins, histogramTransform, darkMode, axisSettings, sidebarWidth, sidebarCollapsed, guiStateLoaded, embedded])

  // Synchronize frontend gates and settings to backend in-memory cache
  useEffect(() => {
    if (!gatesLoaded) return
    const timer = setTimeout(() => {
      gatingApiRequest('/gate_cache', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ gates, pointSize, maxPoints: normalizeEventCount(maxPoints), histogramBins, histogramTransform, histogramTransformVersion: HISTOGRAM_TRANSFORM_VERSION, viewSettings, eventCountVersion: EVENT_COUNT_VERSION })
      }).catch(() => {})
    }, 400)
    return () => clearTimeout(timer)
  }, [gates, pointSize, maxPoints, histogramBins, histogramTransform, viewSettings, gatesLoaded])

  useEffect(() => {
    if (!gatesLoaded || !files.length) return
    const timer = setTimeout(() => {
      gatingApiRequest('/gate_config', {
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
      gatingApiRequest('/gate_spectra', {
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
    gatingApiRequest('/gate_spectra', {
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
    // Wait for persisted gate settings before choosing the cache key or loading
    // events. Otherwise the default event count starts one preload and the saved
    // event count immediately starts a second, much more expensive preload.
    if (!guiStateLoaded || !files.length) return undefined
    const requestedPoints = normalizeEventCount(maxPoints)
    const preloadKey = appletCacheKey('gating-events', activeProjectPath, projectRevision, requestedPoints)
    const preloadGroup = appletCacheKey('gating-events')
    const controller = new AbortController()
    let disposed = false
    setPayload(null)
    setPayloadCache({})
    payloadCacheRef.current = {}
    setPreloadComplete(false)
    setStatus('Preloading controls')

    loadCachedAppletData(preloadKey, async () => {
      const data = await gatingApiRequest(`/gate_preload_compact?max_points=${requestedPoints}`, {
        signal: controller.signal,
      })
      const next = {}
      ;(data?.payloads || []).forEach((item) => {
        const decoded = decodeCompactPayload(item)
        const filename = decoded?.file?.filename || decoded?.filename
        if (filename && !decoded?.error) next[filename] = decoded
      })
      return next
    }, preloadGroup)
      .then((next) => {
        if (disposed) return
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
  }, [fileInventoryKey, maxPoints, activeProjectPath, projectRevision, guiStateLoaded])

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

  const {
    saveConfig,
    saveConfigAsCsv,
    loadConfigFromPicker,
    loadConfigFromFile,
  } = useGatingConfigIO({
    gates,
    files,
    pointSize,
    maxPoints,
    histogramBins,
    histogramTransform,
    viewSettings,
    embedded,
    onRequestExit,
    loadInputRef,
    setStatus,
    setConfigs,
    setFiles,
    setMetadata,
    setSelected,
    setGates,
    setPointSize,
    setMaxPoints,
    setHistogramBins,
    setHistogramTransform,
    setViewSettings,
    setDraft,
  })

  const autogateHistograms = () => autogateHistogramsAction({
    enabled: canAutogateHistograms,
    gates,
    files,
    setHistogramAutogating,
    setHistogramAutogateNotice,
    setStatus,
    setGates,
    setSpectrumCache,
  })

  const beginSidebarResize = (event) => beginGatingSidebarResize(
    event, sidebarCollapsed, sidebarWidth, setSidebarWidth,
  )

  return <GatingWorkspace view={{
    darkMode,
    initialLoading,
    histogramAutogating,
    histogramAutogateNotice,
    setHistogramAutogateNotice,
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
  }} />
}

export default App
