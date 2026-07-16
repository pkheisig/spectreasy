import { ArrowRight } from 'lucide-react'
import type { ProjectState } from '../types'

type Props = {
  project: ProjectState
  onView: () => void
}

export function ControlReportPanel({ project, onView }: Props) {
  const reports = project.artifacts.filter((artifact) =>
    artifact.group === 'Reports' && /control|qc/i.test(artifact.type),
  )
  const latest = reports.reduce<(typeof reports)[number] | undefined>((newest, report) =>
    !newest || (report.updatedEpoch ?? 0) > (newest.updatedEpoch ?? 0) ? report : newest,
  undefined)

  return (
    <section className="surface-card report-snapshot">
      <div className="card-toolbar">
        <div>
          <h2>{latest ? 'Control QC report available' : 'No control QC report yet'}</h2>
          <p>{latest ? latest.name : 'Run control unmixing with report generation enabled.'}</p>
        </div>
        <div className="toolbar-actions">
          {latest && (
            <button className="button button-ghost" onClick={onView}>
              <ArrowRight size={14} /> View report
            </button>
          )}
        </div>
      </div>
    </section>
  )
}
