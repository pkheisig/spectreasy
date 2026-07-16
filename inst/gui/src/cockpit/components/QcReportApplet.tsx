import { useEffect, useMemo, useState } from 'react'
import { Download, FileDown, LoaderCircle } from 'lucide-react'
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

export default function QcReportApplet({ kind, theme, projectPath }: QcReportAppletProps) {
  const [reports, setReports] = useState<Report[]>([])
  const [loading, setLoading] = useState(true)
  const [exporting, setExporting] = useState(false)
  const [message, setMessage] = useState('')
  const report = useMemo(() => newestReport(reports, kind), [reports, kind])

  useEffect(() => {
    let active = true
    void loadProjectReports(projectPath).then((nextReports) => {
      if (!active) return
      setReports(nextReports)
      setLoading(false)
    })
    return () => {
      active = false
    }
  }, [kind, projectPath])

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

  return (
    <section className={`qc-report-applet theme-${theme}`} aria-label={`${kind} QC report viewer`}>
      <header className="qc-report-toolbar">
        <strong>{kind === 'control' ? 'Controls' : 'Samples'} QC report</strong>
        {report && <div className="qc-report-actions">
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
      {loading ? (
        <div className="qc-report-empty"><LoaderCircle className="is-spinning" size={22} /><span>Loading report…</span></div>
      ) : report?.path ? (
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
