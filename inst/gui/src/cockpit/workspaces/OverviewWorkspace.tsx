import { WorkflowStageCard } from '../components/WorkflowStageCard'
import type { Job, ProjectState, StepState, WorkflowSettings } from '../types'

type OverviewWorkspaceProps = {
  project: ProjectState
  job: Job
  settings: WorkflowSettings
  onControls: (tab: 'mapping' | 'gating' | 'build' | 'qc') => void
  onSamples: () => void
}

function formatReportDate(report: ProjectState['latestControlReport']): string {
  if (!report) return 'Not available'
  const date = new Date(report.updatedEpoch * 1000)
  if (!Number.isFinite(date.getTime())) return report.updated || 'Not available'
  return new Intl.DateTimeFormat(undefined, {
    day: '2-digit', month: 'short', year: 'numeric', hour: '2-digit', minute: '2-digit',
  }).format(date)
}

export function OverviewWorkspace({ project, job, settings, onControls, onSamples }: OverviewWorkspaceProps) {
  const mappingState: StepState = project.mappingDirty ? 'warning' : project.mappingExists ? 'complete' : 'blocked'
  const gatesState: StepState = project.gatesDirty ? 'warning' : project.scan.gates > 0 ? 'complete' : 'blocked'
  const matrixState: StepState = project.matrixState === 'current' ? 'complete' : project.matrixState === 'stale' ? 'stale' : 'blocked'
  const controlRunning = job.state === 'running' && /control/i.test(job.label)
  const sampleRunning = job.state === 'running' && /sample/i.test(job.label)

  let controlAction = 'Open controls'
  let controlTab: 'mapping' | 'gating' | 'build' | 'qc' = 'mapping'
  if (!project.mappingExists) controlAction = 'Set up controls'
  else if (project.mappingDirty) controlAction = 'Review mapping'
  else if (project.scan.gates === 0 || project.gatesDirty) {
    controlAction = 'Review gates'
    controlTab = 'gating'
  } else if (project.matrixState === 'missing') {
    controlAction = 'Build reference'
    controlTab = 'build'
  } else if (project.matrixState === 'stale') {
    controlAction = 'Rebuild reference'
    controlTab = 'build'
  }

  const sampleDisabled = project.matrixState === 'missing' || sampleRunning
  const sampleAction = project.sampleOutputState === 'stale'
    ? 'Re-run samples'
    : project.unmixedSampleCount > 0
      ? 'Open samples'
      : 'Unmix samples'
  const sampleStatus: StepState = sampleRunning
    ? 'ready'
    : project.matrixState === 'missing'
      ? 'blocked'
      : project.sampleOutputState === 'stale'
        ? 'stale'
        : project.unmixedSampleCount > 0
          ? 'complete'
          : 'ready'

  return (
    <div className="overview-stage-grid" aria-label="Workflow overview">
      <WorkflowStageCard
        title="Controls"
        status={controlRunning ? 'ready' : matrixState}
        rows={[
          { label: 'Control files', value: String(project.scan.controls) },
          { label: 'Mapping', value: project.mappingDirty ? 'Unsaved' : project.mappingExists ? 'Confirmed' : 'Missing', state: mappingState },
          { label: 'Gates', value: project.gatesDirty ? 'Unsaved' : project.scan.gates > 0 ? 'Confirmed' : 'Missing', state: gatesState },
          { label: 'Reference matrix', value: project.matrixState === 'current' ? 'Current' : project.matrixState === 'stale' ? 'Stale' : 'Missing', state: matrixState },
          { label: 'Method', value: settings.control.method },
          { label: 'Latest report', value: formatReportDate(project.latestControlReport) },
        ]}
        actionLabel={controlRunning ? 'Controls running' : controlAction}
        actionDisabled={controlRunning}
        onAction={() => onControls(controlTab)}
      />
      <WorkflowStageCard
        title="Samples"
        status={sampleStatus}
        rows={[
          { label: 'Sample files', value: String(project.scan.samples) },
          { label: 'Reference matrix', value: project.matrixState === 'current' ? 'Current' : project.matrixState === 'stale' ? 'Stale' : 'Missing', state: matrixState },
          { label: 'Unmixed outputs', value: String(project.unmixedSampleCount), state: project.sampleOutputState === 'stale' ? 'stale' : undefined },
          { label: 'Method', value: settings.sample.method },
          { label: 'Latest report', value: formatReportDate(project.latestSampleReport) },
        ]}
        actionLabel={sampleRunning ? 'Samples running' : sampleAction}
        actionDisabled={sampleDisabled}
        actionTooltip={project.matrixState === 'missing' ? 'Build the reference matrix before unmixing samples.' : undefined}
        onAction={onSamples}
      />
    </div>
  )
}
