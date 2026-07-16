import { useEffect, useMemo, useState } from 'react'
import { Download, ExternalLink, FileDown, LoaderCircle } from 'lucide-react'
import {
  downloadProjectReport,
  exportProjectReportPdf,
  loadProjectReports,
  projectFileUrl,
} from '../api'
import type { Report } from '../types'

type QcReportAppletProps = {
  kind: 'control' | 'sample'
  theme: 'light' | 'dark'
  projectPath: string
  outputRoot: string
  initialReportPath: string
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

export default function QcReportApplet({ kind, theme, projectPath, outputRoot, initialReportPath }: QcReportAppletProps) {
  const requestKey = `${projectPath}\n${outputRoot}\n${kind}`
  const [reports, setReports] = useState<Report[]>([])
  const [loading, setLoading] = useState(true)
  const [loadFailed, setLoadFailed] = useState(false)
  const [loadedRequestKey, setLoadedRequestKey] = useState('')
  const [exporting, setExporting] = useState(false)
  const [message, setMessage] = useState('')
  const [selectedReportPath, setSelectedReportPath] = useState(initialReportPath)
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

  const openInNewTab = () => {
    if (!report?.path) return
    window.open(projectFileUrl(report.path, projectPath), '_blank', 'noopener,noreferrer')
  }

  return (
    <section className={`qc-report-applet theme-${theme}`} aria-label={`${kind} QC report viewer`}>
      <header className="qc-report-toolbar">
        <strong>{report ? reportLabel(report.path) : `${kind === 'control' ? 'Controls' : 'Samples'} QC reports`}</strong>
        {report && <div className="qc-report-actions">
          <button type="button" onClick={openInNewTab}>
            <ExternalLink size={15} /> Open in new tab
          </button>
          <button type="button" onClick={() => void download()}>
            <Download size={15} /> Download {report.format}
          </button>
          {report.format === 'HTML' && <button type="button" onClick={() => void exportPdf()} disabled={exporting}>
            {exporting ? <LoaderCircle className="is-spinning" size={15} /> : <FileDown size={15} />}
            {exporting ? 'Exporting…' : 'Export PDF'}
          </button>}
        </div>}
      </header>
      {message && <div className="qc-report-message" role="status">{message}</div>}
      {loading || loadedRequestKey !== requestKey ? (
        <div className="qc-report-empty"><LoaderCircle className="is-spinning" size={22} /><span>Loading report…</span></div>
      ) : loadFailed ? (
        <div className="qc-report-empty">
          <strong>Could not load the QC report</strong>
          <span>The local R backend did not return the report list.</span>
        </div>
      ) : reportsForKind.length && report?.path ? (
        <iframe
          className="qc-report-frame"
          src={projectFileUrl(report.path, projectPath)}
          title={`${kind} QC report`}
          {...(report.format === 'HTML' ? { sandbox: 'allow-scripts' } : {})}
        />
      ) : (
        <div className="qc-report-empty">
          <strong>No QC report yet</strong>
          <span>Run {kind} unmixing to create it.</span>
        </div>
      )}
    </section>
  )
}
