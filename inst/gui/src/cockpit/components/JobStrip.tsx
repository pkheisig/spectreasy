import { AlertCircle, Check, RefreshCcw } from 'lucide-react'
import type { Job } from '../types'
import { StatusPill } from './StatusPill'

export function JobStrip({ job }: { job: Job }) {
  if (job.state === 'idle') return null
  const state = job.state === 'complete' ? 'complete' : job.state === 'failed' ? 'warning' : 'ready'

  return (
    <div className={`job-strip job-${job.state}`} role="status" aria-live="polite">
      <div className="job-icon">
        {job.state === 'running' ? (
          <RefreshCcw size={16} className="spin" />
        ) : job.state === 'complete' ? (
          <Check size={16} />
        ) : (
          <AlertCircle size={16} />
        )}
      </div>
      <div className="job-copy">
        <strong>{job.label}</strong>
        {job.subtask && <span>{job.subtask}</span>}
      </div>
      <div className="job-progress">
        <div className={`progress-track ${job.state === 'running' ? 'is-indeterminate' : ''}`}>
          <span style={job.state === 'running' ? undefined : { width: `${job.progress}%` }} />
        </div>
        {job.state !== 'running' && <span>{job.progress}%</span>}
      </div>
      <StatusPill
        state={state}
        label={job.state === 'running' ? 'Running' : job.state === 'complete' ? 'Complete' : 'Failed'}
        compact
      />
    </div>
  )
}
