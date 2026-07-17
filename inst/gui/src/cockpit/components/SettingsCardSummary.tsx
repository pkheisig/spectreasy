import { RotateCcw } from 'lucide-react'
import { createPortal } from 'react-dom'
import { useState } from 'react'
import type { ReactNode, MouseEvent } from 'react'

type ResetProps = {
  label: string
  onReset: () => void
}

export function ResetSettingsButton({ label, onReset }: ResetProps) {
  const [confirming, setConfirming] = useState(false)

  function reset(event: MouseEvent<HTMLButtonElement>) {
    event.preventDefault()
    event.stopPropagation()
    setConfirming(true)
  }

  return <>
    <button className="settings-reset-button" type="button" aria-label={`Reset ${label} to default`} onClick={reset}>
      <RotateCcw size={13} />
      <span>Reset to default</span>
    </button>
    {confirming && createPortal(
      <div className="cockpit-confirm-overlay" role="presentation" onMouseDown={() => setConfirming(false)}>
        <div className="cockpit-confirm" role="dialog" aria-modal="true" aria-labelledby="reset-settings-title" onMouseDown={(event) => event.stopPropagation()}>
          <h2 id="reset-settings-title">Reset {label}?</h2>
          <p>Your current {label} will be replaced with the default values.</p>
          <div>
            <button className="button button-ghost" type="button" onClick={() => setConfirming(false)}>Cancel</button>
            <button className="button button-danger" type="button" onClick={() => { setConfirming(false); onReset() }}>Reset to default</button>
          </div>
        </div>
      </div>,
      document.body,
    )}
  </>
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
