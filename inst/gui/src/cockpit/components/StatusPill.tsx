import { AlertCircle, Check, Clock3, FlaskConical, Info, LoaderCircle, MinusCircle, XCircle } from 'lucide-react'
import type { StepState, ArtifactStatus } from '../types'

type StatusPillProps = {
  state: StepState | ArtifactStatus | 'connected' | 'offline'
  label?: string
  compact?: boolean
}

const labels: Record<string, string> = {
  complete: 'Complete',
  current: 'Current',
  ready: 'Ready',
  warning: 'Warning',
  stale: 'Stale',
  blocked: 'Blocked',
  idle: 'Not started',
  user: 'User supplied',
  synthetic: 'Synthetic',
  missing: 'Missing',
  unknown: 'Unknown',
  connected: 'Connected',
  offline: 'Preview mode',
}

const icons = {
  complete: Check,
  current: Check,
  connected: Check,
  ready: Clock3,
  warning: AlertCircle,
  stale: Clock3,
  blocked: XCircle,
  idle: MinusCircle,
  user: Info,
  synthetic: FlaskConical,
  missing: XCircle,
  unknown: Info,
  offline: LoaderCircle,
}

export function StatusPill({ state, label, compact = false }: StatusPillProps) {
  const Icon = icons[state]
  return (
    <span className={`status-pill status-${state} ${compact ? 'is-compact' : ''}`}>
      <Icon size={compact ? 11 : 13} strokeWidth={2.3} />
      <span>{label ?? labels[state]}</span>
    </span>
  )
}

