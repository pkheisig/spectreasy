import { useEffect, useMemo, useRef, useState } from 'react'
import { Check, ChevronDown, Clipboard, Download, ExternalLink, FileDown, LoaderCircle, Sparkles, X } from 'lucide-react'
import {
  downloadProjectReport,
  downloadProjectReportPrompt,
  exportProjectReportPdf,
  loadProjectReportObjectUrl,
  loadProjectReportPrompt,
  loadProjectReports,
} from '../api'
import type { Report } from '../types'

type QcReportAppletProps = {
  kind: 'control' | 'sample'
  theme: 'light' | 'dark'
  projectPath: string
  outputRoot: string
  initialReportPath: string
  onRequestExit: () => void
}

function reportTime(report: Report): number {
  if (typeof report.createdEpoch === 'number' && Number.isFinite(report.createdEpoch)) return report.createdEpoch * 1000
  const parsed = Date.parse(report.created)
  return Number.isFinite(parsed) ? parsed : 0
}

function newestReport(reports: Report[], kind: 'control' | 'sample') {
  const type = kind === 'control' ? 'Control QC' : 'Sample QC'
  return reports
    .filter((report) => report.type === type && report.path)
    .reduce<Report | null>((newest, report) => !newest || reportTime(report) > reportTime(newest) ? report : newest, null)
}

export default function QcReportApplet({ kind, theme, projectPath, outputRoot, initialReportPath, onRequestExit }: QcReportAppletProps) {
  const requestKey = `${projectPath}\n${outputRoot}\n${kind}`
  const [reports, setReports] = useState<Report[]>([])
  const [loading, setLoading] = useState(true)
  const [loadFailed, setLoadFailed] = useState(false)
  const [loadedRequestKey, setLoadedRequestKey] = useState('')
  const [exporting, setExporting] = useState(false)
  const [promptMenuOpen, setPromptMenuOpen] = useState(false)
  const [promptBusy, setPromptBusy] = useState(false)
  const [promptCopied, setPromptCopied] = useState(false)
  const [message, setMessage] = useState('')
  const promptMenuRef = useRef<HTMLDivElement>(null)
  const [selectedReportPath, setSelectedReportPath] = useState(initialReportPath)
  const [reportObject, setReportObject] = useState({ key: '', url: '' })
  const reportsForKind = useMemo(() => {
    const type = kind === 'control' ? 'Control QC' : 'Sample QC'
    return reports
      .filter((candidate): candidate is Report & { path: string } => candidate.type === type && Boolean(candidate.path))
      .sort((left, right) => reportTime(right) - reportTime(left))
  }, [reports, kind])
  const report = useMemo(
    () => reportsForKind.find((candidate) => candidate.path === selectedReportPath) ?? reportsForKind[0] ?? null,
    [reportsForKind, selectedReportPath],
  )
  const reportObjectKey = report?.path ? `${projectPath}\n${report.path}` : ''
  const reportObjectUrl = reportObject.key === reportObjectKey ? reportObject.url : ''
  const promptExportAvailable = Boolean(report?.promptPath)

  useEffect(() => {
    let active = true
    void loadProjectReports(projectPath, outputRoot, kind).then((nextReports) => {
      if (!active) return
      setReports(nextReports)
      setSelectedReportPath((current) => {
        if (initialReportPath && nextReports.some((candidate) => candidate.path === initialReportPath)) return initialReportPath
        if (current && nextReports.some((candidate) => candidate.path === current)) return current
        return newestReport(nextReports, kind)?.path ?? ''
      })
      setLoadFailed(false)
      setLoading(false)
      setLoadedRequestKey(requestKey)
    }).catch(() => {
      if (!active) return
      setLoadFailed(true)
      setLoading(false)
      setLoadedRequestKey(requestKey)
    })
    return () => {
      active = false
    }
  }, [initialReportPath, kind, outputRoot, projectPath, requestKey])

  useEffect(() => {
    let active = true
    let objectUrl = ''
    if (!report?.path) return () => { active = false }
    void loadProjectReportObjectUrl(report.path, projectPath).then((nextUrl) => {
      objectUrl = nextUrl
      if (active) setReportObject({ key: `${projectPath}\n${report.path}`, url: nextUrl })
      else URL.revokeObjectURL(nextUrl)
    }).catch(() => {
      if (active) setMessage(`Could not load the ${report.format} report.`)
    })
    return () => {
      active = false
      if (objectUrl) URL.revokeObjectURL(objectUrl)
    }
  }, [projectPath, report?.format, report?.path])

  useEffect(() => {
    if (!promptMenuOpen) return
    const close = (event: MouseEvent) => {
      if (!promptMenuRef.current?.contains(event.target as Node)) setPromptMenuOpen(false)
    }
    const closeOnEscape = (event: KeyboardEvent) => {
      if (event.key === 'Escape') setPromptMenuOpen(false)
    }
    document.addEventListener('mousedown', close)
    document.addEventListener('keydown', closeOnEscape)
    return () => {
      document.removeEventListener('mousedown', close)
      document.removeEventListener('keydown', closeOnEscape)
    }
  }, [promptMenuOpen])

  const reportLabel = (path: string) => path.replace(/\\/g, '/').split('/').slice(-2).join('/')

  const download = async () => {
    if (!report?.path) return
    setMessage('')
    if (!await downloadProjectReport(report.path, projectPath)) setMessage(`Could not download the ${report.format} report.`)
  }

  const exportPdf = async () => {
    if (!report?.path) return
    setMessage('')
    setExporting(true)
    const success = await exportProjectReportPdf(report.path, projectPath)
    setExporting(false)
    if (!success) setMessage('Could not export the existing HTML report as PDF.')
  }

  const openInNewTab = async () => {
    if (!report?.path) return
    const popup = window.open('about:blank', '_blank')
    if (!popup) {
      setMessage('The browser blocked the new report tab.')
      return
    }
    popup.opener = null
    try {
      const objectUrl = await loadProjectReportObjectUrl(report.path, projectPath)
      popup.location.replace(objectUrl)
      window.setTimeout(() => URL.revokeObjectURL(objectUrl), 60000)
    } catch {
      popup.close()
      setMessage(`Could not open the ${report.format} report.`)
    }
  }

  const copyPrompt = async () => {
    if (!report?.promptPath) return
    setMessage('')
    setPromptBusy(true)
    try {
      const prompt = await loadProjectReportPrompt(report.promptPath, projectPath)
      try {
        await navigator.clipboard.writeText(prompt)
      } catch {
        const textarea = document.createElement('textarea')
        textarea.value = prompt
        textarea.style.position = 'fixed'
        textarea.style.opacity = '0'
        document.body.appendChild(textarea)
        textarea.select()
        document.execCommand('copy')
        textarea.remove()
      }
      setPromptCopied(true)
      window.setTimeout(() => setPromptCopied(false), 1600)
      setPromptMenuOpen(false)
    } catch {
      setMessage('Could not copy the saved AI-ready QC prompt.')
    } finally {
      setPromptBusy(false)
    }
  }

  const downloadPrompt = async () => {
    if (!report?.promptPath) return
    setMessage('')
    setPromptBusy(true)
    const success = await downloadProjectReportPrompt(report.promptPath, projectPath)
    setPromptBusy(false)
    setPromptMenuOpen(false)
    if (!success) setMessage('Could not download the saved AI-ready QC prompt.')
  }

  return (
    <section className={`qc-report-applet theme-${theme}`} aria-label={`${kind} QC report viewer`}>
      <header className="qc-report-toolbar">
        <strong>{report ? reportLabel(report.path) : `${kind === 'control' ? 'Controls' : 'Samples'} QC reports`}</strong>
        <div className="qc-report-actions">
          {report && <>
          <div className="qc-prompt-action" ref={promptMenuRef}>
            <button
              type="button"
              className="qc-prompt-trigger"
              aria-haspopup="menu"
              aria-expanded={promptMenuOpen}
              disabled={!promptExportAvailable || promptBusy}
              title={promptExportAvailable ? 'Copy or download the precomputed AI prompt and its numeric QC data' : 'This report predates automatic AI prompt export. Regenerate the QC report to enable this action.'}
              onClick={() => setPromptMenuOpen((open) => !open)}
            >
              {promptCopied ? <Check size={15} /> : <Sparkles size={15} />}
              {promptCopied ? 'Prompt copied' : 'Export AI prompt'}
              <ChevronDown size={14} />
            </button>
            {promptMenuOpen && promptExportAvailable && <div className="qc-prompt-menu" role="menu" aria-label="AI prompt export actions">
              <button type="button" role="menuitem" onClick={() => void copyPrompt()}>
                <Clipboard size={15} /><span><strong>Copy</strong><small>Prompt + all numeric QC data</small></span>
              </button>
              <button type="button" role="menuitem" onClick={() => void downloadPrompt()}>
                <FileDown size={15} /><span><strong>Download prompt</strong><small>{report.promptBytes ? `${Math.max(1, Math.round(report.promptBytes / 1024)).toLocaleString()} KB text file` : 'Complete text file'}</small></span>
              </button>
            </div>}
          </div>
          <button type="button" onClick={() => void openInNewTab()}>
            <ExternalLink size={15} /> Open in new tab
          </button>
          <button type="button" onClick={() => void download()}>
            <Download size={15} /> Download {report.format}
          </button>
          {report.format === 'HTML' && <button type="button" onClick={() => void exportPdf()} disabled={exporting}>
            {exporting ? <LoaderCircle className="is-spinning" size={15} /> : <FileDown size={15} />}
            {exporting ? 'Exporting…' : 'Export PDF'}
          </button>}
          </>}
          <button
            type="button"
            className="qc-report-close"
            onClick={onRequestExit}
            aria-label={`Close ${kind} QC report and return to cockpit`}
            autoFocus
          >
            <X size={15} /> Close
          </button>
        </div>
      </header>
      {message && <div className="qc-report-message" role="status">{message}</div>}
      {loading || loadedRequestKey !== requestKey ? (
        <div className="qc-report-empty"><LoaderCircle className="is-spinning" size={22} /><span>Loading report…</span></div>
      ) : loadFailed ? (
        <div className="qc-report-empty">
          <strong>Could not load the QC report</strong>
          <span>The local R backend did not return the report list.</span>
        </div>
      ) : reportsForKind.length && report?.path && reportObjectUrl ? (
        <iframe
          className="qc-report-frame"
          src={reportObjectUrl}
          title={`${kind} QC report`}
          {...(report.format === 'HTML' ? { sandbox: 'allow-scripts' } : {})}
        />
      ) : reportsForKind.length && report?.path ? (
        <div className="qc-report-empty"><LoaderCircle className="is-spinning" size={22} /><span>Loading report…</span></div>
      ) : (
        <div className="qc-report-empty">
          <strong>No QC report yet</strong>
          <span>Run {kind} unmixing to create it.</span>
        </div>
      )}
    </section>
  )
}
