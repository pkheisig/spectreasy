import { pruneUnavailableNegativeGates } from '../gatingEligibility.js'
import { decodeSpectrumData } from '../spectrumData.js'
import {
  spectrumCacheKeyForFile,
  unboxGuiState,
  gatingApiRequest,
} from './GatingCore.jsx'

export const GATING_SIDEBAR_MIN_WIDTH = 244

export async function autogateHistogramsAction(options) {
  const {
    enabled, gates, files, setHistogramAutogating, setHistogramAutogateNotice,
    setStatus, setGates, setSpectrumCache,
  } = options
  if (!enabled) return
  setHistogramAutogating(true)
  setHistogramAutogateNotice?.(null)
  setStatus('Auto-generating histogram gates')
  try {
    const result = await gatingApiRequest('/gate_histogram_autogate', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ gates }),
    })
    if (result.success === false) throw new Error(result.error || 'Histogram autogating failed')
    const generatedGates = unboxGuiState(result.gates || {})
    const nextGates = pruneUnavailableNegativeGates({ ...gates, ...generatedGates }, files)
    const returnedSpectra = unboxGuiState(result.spectra || {})
    setGates(nextGates)
    setSpectrumCache((previous) => {
      const next = { ...previous }
      files.forEach((file) => {
        const spectrum = decodeSpectrumData(returnedSpectra?.[file.filename])
        const key = spectrumCacheKeyForFile(nextGates, file)
        if (key && spectrum) next[key] = spectrum
      })
      return next
    })
    const generated = Number(result.gates_generated ?? Object.keys(generatedGates).length)
    const preserved = Number(result.gates_preserved ?? 0)
    const decodedFailures = unboxGuiState(result.failures || [])
    const failures = Array.isArray(decodedFailures)
      ? decodedFailures
      : decodedFailures && typeof decodedFailures === 'object'
        ? [decodedFailures]
        : []
    if (failures.length > 0) {
      setHistogramAutogateNotice?.({
        generated,
        processed: Number(result.files_processed ?? files.length),
        failures: failures.map((failure) => ({
          filename: String(failure?.filename || 'Unknown control'),
          error: String(failure?.error || 'Histogram gates could not be generated.'),
        })),
      })
    }
    setStatus(
      `Auto-generated ${generated} missing histogram gate${generated === 1 ? '' : 's'} ` +
      `for ${Number(result.files_succeeded ?? result.files_processed)} controls; ` +
      `preserved ${preserved} existing gate${preserved === 1 ? '' : 's'}.`,
    )
  } catch (error) {
    setHistogramAutogateNotice?.(null)
    setStatus(`Could not auto-generate histogram gates: ${error.message}`)
  } finally {
    setHistogramAutogating(false)
  }
}

export function beginGatingSidebarResize(event, collapsed, width, setWidth) {
  if (collapsed) return
  event.preventDefault()
  const startX = event.clientX
  const previousCursor = document.body.style.cursor
  const previousUserSelect = document.body.style.userSelect
  document.body.style.cursor = 'col-resize'
  document.body.style.userSelect = 'none'
  const move = (moveEvent) => setWidth(Math.min(380, Math.max(GATING_SIDEBAR_MIN_WIDTH, width + moveEvent.clientX - startX)))
  const finish = () => {
    window.removeEventListener('pointermove', move)
    window.removeEventListener('pointerup', finish)
    window.removeEventListener('pointercancel', finish)
    document.body.style.cursor = previousCursor
    document.body.style.userSelect = previousUserSelect
  }
  window.addEventListener('pointermove', move)
  window.addEventListener('pointerup', finish)
  window.addEventListener('pointercancel', finish)
}
