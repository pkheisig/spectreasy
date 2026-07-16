import { FolderPlus, X } from 'lucide-react'
import { useEffect, useRef } from 'react'
import { createPortal } from 'react-dom'

type Props = {
  projectName: string
  missingFolders: string[]
  busy: boolean
  message: string
  onConfirm: () => void
  onCancel: () => void
}

export function ProjectInitializationDialog({ projectName, missingFolders, busy, message, onConfirm, onCancel }: Props) {
  const dialogRef = useRef<HTMLDivElement>(null)
  const confirmRef = useRef<HTMLButtonElement>(null)

  useEffect(() => {
    confirmRef.current?.focus()
    const handleKeyDown = (event: KeyboardEvent) => {
      if (event.key === 'Escape' && !busy) {
        event.preventDefault()
        onCancel()
        return
      }
      if (event.key !== 'Tab') return
      const focusable = Array.from(dialogRef.current?.querySelectorAll<HTMLElement>('button:not([disabled])') ?? [])
      if (!focusable.length) return
      const first = focusable[0]
      const last = focusable[focusable.length - 1]
      if (event.shiftKey && document.activeElement === first) {
        event.preventDefault()
        last.focus()
      } else if (!event.shiftKey && document.activeElement === last) {
        event.preventDefault()
        first.focus()
      }
    }
    document.addEventListener('keydown', handleKeyDown)
    return () => document.removeEventListener('keydown', handleKeyDown)
  }, [busy, onCancel])

  const folderCopy = missingFolders.map((folder) => `${folder}/`).join(' and ')
  return createPortal(
    <div className="cockpit-confirm-overlay" role="presentation">
      <div ref={dialogRef} className="cockpit-confirm project-initialize-confirm" role="dialog" aria-modal="true" aria-labelledby="project-initialize-title" aria-describedby="project-initialize-description">
        <FolderPlus size={22} aria-hidden="true" />
        <button className="project-initialize-close" type="button" onClick={onCancel} disabled={busy} aria-label="Not now">
          <X size={15} />
        </button>
        <h2 id="project-initialize-title">Prepare {projectName}</h2>
        <p id="project-initialize-description">
          This folder does not yet contain {folderCopy}. Create the missing folders to use it as a Spectreasy project?
        </p>
        {message ? <p className="project-initialize-message" role="status" aria-live="polite">{message}</p> : null}
        <div>
          <button className="button button-ghost" type="button" onClick={onCancel} disabled={busy}>Not now</button>
          <button ref={confirmRef} className="button button-primary" type="button" onClick={onConfirm} disabled={busy}>
            {busy ? 'Creating folders…' : 'Create folders'}
          </button>
        </div>
      </div>
    </div>,
    document.body,
  )
}
