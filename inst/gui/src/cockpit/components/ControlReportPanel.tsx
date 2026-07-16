import { ArrowRight, FileText } from 'lucide-react'
import type { ProjectState } from '../types'

type Props = {
  project: ProjectState
  onView: (path: string) => void
}

export function ControlReportPanel({ project, onView }: Props) {
  const reports = project.artifacts.filter((artifact) =>
    artifact.group === 'Reports' && /control|qc/i.test(artifact.type),
  )
  reports.sort((left, right) => (right.updatedEpoch ?? 0) - (left.updatedEpoch ?? 0))

  const reportLabel = (path: string) => path.replace(/\\/g, '/').split('/').slice(-2).join('/')

  if (!reports.length) {
    return (
      <section className="surface-card report-history-empty">
        <h2>No control QC report yet</h2>
        <p>Run control unmixing with report generation enabled.</p>
      </section>
    )
  }

  return (
    <div className="report-history-grid" aria-label="Control QC reports">
      {reports.map((report) => {
        const relativePath = report.relativePath ?? report.path
        return (
          <article className="surface-card report-history-card" key={report.id}>
            <span className="report-history-icon"><FileText size={18} /></span>
            <div className="report-history-copy">
              <h2>{reportLabel(relativePath)}</h2>
              <p>{report.updated}</p>
            </div>
            <button className="button button-ghost" onClick={() => onView(relativePath)}>
              View report <ArrowRight size={14} />
            </button>
          </article>
        )
      })}
    </div>
  )
}
