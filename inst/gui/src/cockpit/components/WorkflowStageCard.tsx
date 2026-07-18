import { ArrowRight } from 'lucide-react'
import { StatusPill } from './StatusPill'
import type { StepState } from '../types'

export type WorkflowStageCardProps = {
  title: string
  status: StepState
  rows: Array<{ label: string; value: string; state?: StepState }>
  actionLabel: string
  actionDisabled?: boolean
  actionTooltip?: string
  onAction: () => void
}

export function WorkflowStageCard({
  title,
  status,
  rows,
  actionLabel,
  actionDisabled = false,
  actionTooltip,
  onAction,
}: WorkflowStageCardProps) {
  return (
    <article className="surface-card workflow-stage-card">
      <header>
        <h2>{title}</h2>
        <StatusPill state={status} compact />
      </header>
      <dl>
        {rows.map((row) => (
          <div key={row.label}>
            <dt>{row.label}</dt>
            <dd className={row.state ? `stage-value-${row.state}` : undefined}>{row.value}</dd>
          </div>
        ))}
      </dl>
      <button
        className="button button-primary"
        type="button"
        disabled={actionDisabled}
        title={actionTooltip}
        onClick={onAction}
      >
        {actionLabel}
        <ArrowRight size={15} />
      </button>
    </article>
  )
}
