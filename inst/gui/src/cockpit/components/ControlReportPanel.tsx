import { ArrowRight, RefreshCcw } from 'lucide-react'
import type { ProjectState } from '../types'

type ControlRun = (
  action: 'control' | 'sample' | 'control-report' | 'sample-report' | 'af',
  label: string,
) => Promise<boolean>

type Props = {
  project: ProjectState
  onRun: ControlRun
  onView: () => void
}

export function ControlReportPanel({ project, onRun, onView }: Props) {
  const reports = project.artifacts.filter((artifact) =>
    artifact.group === 'Reports' && /control|qc/i.test(artifact.type),
  )
  const latest = reports.at(-1)

  return (
    <section className="surface-card report-snapshot">
      <div className="card-toolbar">
        <div>
          <h2>{latest ? 'Control QC report available' : 'No control QC report yet'}</h2>
          <p>{latest ? latest.name : 'Generate a report after the control workflow completes.'}</p>
        </div>
        <div className="toolbar-actions">
          {latest && (
            <button className="button button-ghost" onClick={onView}>
              <ArrowRight size={14} /> View report
            </button>
          )}
          <button
            className="button button-primary"
            onClick={() => onRun('control-report', 'Render control QC report')}
          >
            <RefreshCcw size={14} /> {latest ? 'Regenerate' : 'Generate report'}
          </button>
        </div>
      </div>
    </section>
  )
}
