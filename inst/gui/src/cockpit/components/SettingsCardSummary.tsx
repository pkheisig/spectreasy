import { RotateCcw } from 'lucide-react'
import type { ReactNode, MouseEvent } from 'react'

type ResetProps = {
  label: string
  onReset: () => void
}

export function ResetSettingsButton({ label, onReset }: ResetProps) {
  function reset(event: MouseEvent<HTMLButtonElement>) {
    event.preventDefault()
    event.stopPropagation()
    onReset()
  }

  return (
    <button className="settings-reset-button" type="button" aria-label={`Reset ${label} to default`} onClick={reset}>
      <RotateCcw size={13} />
      <span>Reset to default</span>
    </button>
  )
}

type SummaryProps = Omit<ResetProps, 'label'> & {
  icon?: ReactNode
  title: string
  detail?: ReactNode
}

export function SettingsCardSummary({ icon, title, detail, onReset }: SummaryProps) {
  return (
    <>
      <span className="settings-summary-label">{icon}{title}</span>
      {detail && <span className="settings-summary-detail">{detail}</span>}
      <ResetSettingsButton label={title} onReset={onReset} />
    </>
  )
}
