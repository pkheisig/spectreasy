import { CopyPlus, X } from 'lucide-react'
import { useEffect, useRef } from 'react'
import { createPortal } from 'react-dom'

type Props = {
  workflow: 'control' | 'sample'
  onConfirm: () => void
  onCancel: () => void
}

export function VersionRunDialog({ workflow, onConfirm, onCancel }: Props) {
  const confirmRef = useRef<HTMLButtonElement>(null)
  useEffect(() => {
    confirmRef.current?.focus()
    const close = (event: KeyboardEvent) => {
      if (event.key === 'Escape') onCancel()
    }
    document.addEventListener('keydown', close)
    return () => document.removeEventListener('keydown', close)
  }, [onCancel])
  const name = workflow === 'control' ? 'control' : 'sample'
  return createPortal(
    <div className="cockpit-confirm-overlay" role="presentation">
      <div className="cockpit-confirm" role="dialog" aria-modal="true" aria-labelledby="version-run-title" aria-describedby="version-run-description">
        <CopyPlus size={22} aria-hidden="true" />
        <button className="project-initialize-close" type="button" onClick={onCancel} aria-label="Cancel rerun"><X size={15} /></button>
        <h2 id="version-run-title">Existing {name} results found</h2>
        <p id="version-run-description">
          The current results will remain intact. Continue to write this rerun to a new versioned output directory?
        </p>
        <div>
          <button className="button button-ghost" type="button" onClick={onCancel}>Cancel</button>
          <button ref={confirmRef} className="button button-primary" type="button" onClick={onConfirm}>Create new version</button>
        </div>
      </div>
    </div>,
    document.body,
  )
}
