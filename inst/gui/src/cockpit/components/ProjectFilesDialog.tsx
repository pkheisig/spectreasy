import { useCallback, useEffect, useRef, useState } from 'react'
import { File, FolderOpen, RefreshCcw, Trash2, Upload, X } from 'lucide-react'
import { createPortal } from 'react-dom'
import { deleteProjectFile, listProjectFiles, uploadProjectFile } from '../api'
import type { DragEvent } from 'react'
import type { ProjectFileEntry, ProjectFileKind } from '../api'

type Props = {
  projectName: string
  onClose: () => void
  onChanged: () => void | Promise<void>
}

function formatSize(bytes: number) {
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / 1024 ** 2).toFixed(bytes < 10 * 1024 ** 2 ? 1 : 0)} MB`
}

function formatModified(value: string) {
  const parsed = new Date(value)
  return Number.isNaN(parsed.getTime())
    ? ''
    : new Intl.DateTimeFormat(undefined, { dateStyle: 'medium', timeStyle: 'short' }).format(parsed)
}

export function ProjectFilesDialog({ projectName, onClose, onChanged }: Props) {
  const [kind, setKind] = useState<ProjectFileKind>('controls')
  const [files, setFiles] = useState<ProjectFileEntry[]>([])
  const [loading, setLoading] = useState(true)
  const [busy, setBusy] = useState(false)
  const [dragActive, setDragActive] = useState(false)
  const [message, setMessage] = useState('')
  const [deleteCandidate, setDeleteCandidate] = useState('')
  const inputRef = useRef<HTMLInputElement>(null)

  const refresh = useCallback(async (selectedKind = kind) => {
    setLoading(true)
    const result = await listProjectFiles(selectedKind)
    setFiles(result.files)
    setMessage(result.success ? '' : (result.message ?? 'Project files could not be loaded.'))
    setLoading(false)
  }, [kind])

  useEffect(() => {
    let cancelled = false
    void listProjectFiles(kind).then((result) => {
      if (cancelled) return
      setFiles(result.files)
      setMessage(result.success ? '' : (result.message ?? 'Project files could not be loaded.'))
      setLoading(false)
    })
    return () => { cancelled = true }
  }, [kind])

  async function addFiles(selected: FileList | File[]) {
    const incoming = Array.from(selected)
    if (!incoming.length || busy) return
    setBusy(true)
    setMessage('')
    const failures: string[] = []
    let added = 0
    for (const file of incoming) {
      const result = await uploadProjectFile(kind, file)
      if (result.success) added += 1
      else failures.push(result.message)
    }
    await refresh(kind)
    if (added) await onChanged()
    setMessage(failures.length
      ? `${added ? `${added} file${added === 1 ? '' : 's'} added. ` : ''}${failures.join(' ')}`
      : `${added} file${added === 1 ? '' : 's'} added.`)
    setBusy(false)
    if (inputRef.current) inputRef.current.value = ''
  }

  async function removeFile(filename: string) {
    if (busy) return
    setBusy(true)
    const result = await deleteProjectFile(kind, filename)
    setDeleteCandidate('')
    if (result.success) {
      await refresh(kind)
      await onChanged()
    }
    setMessage(result.message)
    setBusy(false)
  }

  function receiveDrop(event: DragEvent<HTMLDivElement>) {
    event.preventDefault()
    setDragActive(false)
    void addFiles(event.dataTransfer.files)
  }

  return createPortal(
    <div className="project-files-overlay" role="presentation" onMouseDown={onClose}>
      <section className="project-files-dialog" role="dialog" aria-modal="true" aria-labelledby="project-files-title" onMouseDown={(event) => event.stopPropagation()}>
        <header className="project-files-header">
          <div className="project-files-heading">
            <FolderOpen size={19} />
            <div>
              <h2 id="project-files-title">Files</h2>
              <span>{projectName}</span>
            </div>
          </div>
          <button className="project-files-icon" type="button" onClick={onClose} aria-label="Close files"><X size={18} /></button>
        </header>

        <div className="project-files-tabs" role="tablist" aria-label="Project file folders">
          <button type="button" role="tab" aria-selected={kind === 'controls'} className={kind === 'controls' ? 'is-active' : ''} onClick={() => { setLoading(true); setKind('controls') }}>Controls</button>
          <button type="button" role="tab" aria-selected={kind === 'samples'} className={kind === 'samples' ? 'is-active' : ''} onClick={() => { setLoading(true); setKind('samples') }}>Samples</button>
        </div>

        <div
          className={`project-files-dropzone ${dragActive ? 'is-dragging' : ''}`}
          onDragEnter={(event) => { event.preventDefault(); setDragActive(true) }}
          onDragOver={(event) => event.preventDefault()}
          onDragLeave={(event) => { if (event.currentTarget === event.target) setDragActive(false) }}
          onDrop={receiveDrop}
        >
          <div className="project-files-toolbar">
            <span>{kind === 'controls' ? 'scc' : 'samples'}/</span>
            <button className="button button-ghost" type="button" disabled={busy} onClick={() => void refresh(kind)}><RefreshCcw size={14} /> Refresh</button>
            <button className="button" type="button" disabled={busy} onClick={() => inputRef.current?.click()}><Upload size={14} /> Add FCS files</button>
            <input ref={inputRef} type="file" accept=".fcs" multiple hidden onChange={(event) => void addFiles(event.target.files ?? [])} />
          </div>

          <div className="project-files-list" aria-busy={loading || busy}>
            {loading ? (
              <div className="project-files-empty">Loading files…</div>
            ) : files.length === 0 ? (
              <div className="project-files-empty"><Upload size={22} /><strong>Drop FCS files here</strong><span>or use Add FCS files</span></div>
            ) : files.map((entry) => (
              <div className="project-file-row" key={entry.name}>
                <File size={17} />
                <div className="project-file-copy">
                  <strong title={entry.name}>{entry.name}</strong>
                  <span>{formatSize(entry.size)}{entry.modified ? ` · ${formatModified(entry.modified)}` : ''}</span>
                </div>
                {deleteCandidate === entry.name ? (
                  <div className="project-file-delete-confirm">
                    <span>Delete?</span>
                    <button type="button" disabled={busy} onClick={() => void removeFile(entry.name)}>Yes</button>
                    <button type="button" disabled={busy} onClick={() => setDeleteCandidate('')}>No</button>
                  </div>
                ) : (
                  <button className="project-files-icon project-file-delete" type="button" disabled={busy} onClick={() => setDeleteCandidate(entry.name)} aria-label={`Delete ${entry.name}`}><Trash2 size={15} /></button>
                )}
              </div>
            ))}
          </div>
        </div>
        <footer className="project-files-footer">
          <span className={message && /could not|not an|already exists|did not/i.test(message) ? 'is-error' : ''}>{message || `${files.length} FCS file${files.length === 1 ? '' : 's'}`}</span>
          <span>Drag and drop adds files to the selected folder.</span>
        </footer>
      </section>
    </div>,
    document.body,
  )
}
