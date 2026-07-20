import { reconcileGateCsvRows } from '../gatingCsvCompatibility.js'
import { pruneUnavailableNegativeGates } from '../gatingEligibility.js'
import {
  CONFIG_NAME,
  REQUIRED_GATE_CSV_COLUMNS,
  buildRows,
  gateRowsToCsv,
  normalizeEventCount,
  normalizeHistogramBins,
  normalizeHistogramTransform,
  parseConfigRows,
  parseCsvDocument,
  saveCsvWithSystemPicker,
  gatingApiRequest,
  validateGateCsvRows,
} from './GatingCore.jsx'

export function useGatingConfigIO(options) {
  const {
    gates, files, pointSize, maxPoints, histogramBins, histogramTransform, viewSettings,
    embedded, onRequestExit, loadInputRef,
    setStatus, setConfigs, setFiles, setMetadata, setSelected, setGates,
    setPointSize, setMaxPoints, setHistogramBins, setHistogramTransform,
    setViewSettings, setDraft,
  } = options

  async function saveConfig(closeAfter = false) {
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform, viewSettings })
    try {
      const saved = await gatingApiRequest('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows }),
      })
      if (saved.cancelled || saved.success === false) {
        setStatus(saved.message || (saved.cancelled ? 'Save cancelled' : 'Could not save gate CSV'))
        return
      }
      setStatus(`Saved ${saved.rows ?? rows.length} rows to ${saved.path}`)
      const cfg = await gatingApiRequest('/gate_configs')
      setConfigs(cfg.configs || [])
      if (!closeAfter) return
      if (embedded) {
        onRequestExit?.()
        return
      }
      await gatingApiRequest('/gate_shutdown', {
        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: '{}',
      })
      window.close()
    } catch (error) {
      setStatus(`Could not save gate CSV: ${error.message}`)
    }
  }

  async function saveConfigAsCsv() {
    const rows = buildRows(gates, files, { pointSize, maxPoints, histogramBins, histogramTransform, viewSettings })
    const csvText = gateRowsToCsv(rows)
    let filename = CONFIG_NAME
    try {
      filename = await saveCsvWithSystemPicker(CONFIG_NAME, csvText)
    } catch (error) {
      setStatus(error?.name === 'AbortError' ? 'Save cancelled' : `Could not open save picker: ${error.message}`)
      return
    }
    setStatus('Saving gate CSV')
    try {
      const saved = await gatingApiRequest('/gate_config', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ filename: CONFIG_NAME, rows }),
      }).catch(() => null)
      const cfg = await gatingApiRequest('/gate_configs')
      setConfigs(cfg.configs || [])
      setStatus(`Saved ${saved?.rows ?? rows.length} rows to ${filename}`)
    } catch (error) {
      setStatus(`Saved ${rows.length} rows to ${filename}; backend sync failed: ${error.message}`)
    }
  }

  async function applyLoadedConfigRows(rows, sourceName, headers = REQUIRED_GATE_CSV_COLUMNS) {
    const fileData = await gatingApiRequest('/gate_files')
    const currentFiles = (fileData.files || []).filter((file) => file.file_exists)
    setFiles(currentFiles)
    setMetadata(fileData.metadata || {})
    setSelected((current) => currentFiles.some((file) => file.filename === current) ? current : currentFiles[0]?.filename || '')
    const compatible = reconcileGateCsvRows(rows, currentFiles)
    validateGateCsvRows(compatible.rows, headers)
    const parsed = parseConfigRows(compatible.rows)
    setGates(pruneUnavailableNegativeGates(parsed.gates, currentFiles))
    if (typeof parsed.pointSize === 'number') setPointSize(parsed.pointSize)
    if (typeof parsed.maxPoints === 'number') setMaxPoints(normalizeEventCount(parsed.maxPoints))
    if (typeof parsed.histogramBins === 'number') setHistogramBins(normalizeHistogramBins(parsed.histogramBins))
    if (typeof parsed.histogramTransform === 'string') setHistogramTransform(normalizeHistogramTransform(parsed.histogramTransform))
    setViewSettings(parsed.viewSettings)
    setDraft([])
    await gatingApiRequest('/gate_config', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ filename: CONFIG_NAME, rows: compatible.rows }),
    }).catch(() => null)
    const cfg = await gatingApiRequest('/gate_configs')
    setConfigs(cfg.configs || [])
    setStatus(`Loaded ${compatible.rows.length} rows from ${sourceName}`)
  }

  async function loadConfigFromPicker() {
    if (typeof window.showOpenFilePicker !== 'function') {
      loadInputRef.current?.click()
      return
    }
    try {
      const [handle] = await window.showOpenFilePicker({ multiple: false })
      if (!handle) return
      const file = await handle.getFile()
      const parsed = parseCsvDocument(await file.text())
      await applyLoadedConfigRows(parsed.rows, file.name, parsed.headers)
    } catch (error) {
      setStatus(error?.name === 'AbortError' ? 'Load cancelled' : `Could not load gate CSV: ${error.message}`)
    }
  }

  async function loadConfigFromFile(file) {
    if (!file) return
    try {
      const parsed = parseCsvDocument(await file.text())
      await applyLoadedConfigRows(parsed.rows, file.name, parsed.headers)
    } catch (error) {
      setStatus(`Could not load gate CSV: ${error.message}`)
    } finally {
      if (loadInputRef.current) loadInputRef.current.value = ''
    }
  }

  return { saveConfig, saveConfigAsCsv, loadConfigFromPicker, loadConfigFromFile }
}
