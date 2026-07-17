import { useCallback, useEffect, useState } from 'react'
import { File, Trash2 } from 'lucide-react'
import { deleteProjectFile, listProjectFiles } from '../api'
import type { ProjectFileEntry, ProjectFileKind } from '../api'

type RowsProps = {
  files: ProjectFileEntry[]
  loading: boolean
  busy?: boolean
  emptyLabel?: string
  onDelete: (filename: string) => void | Promise<void>
}

function formatProjectFileSize(bytes: number) {
  if (bytes < 1024) return `${bytes} B`
  if (bytes < 1024 ** 2) return `${(bytes / 1024).toFixed(1)} KB`
  return `${(bytes / 1024 ** 2).toFixed(bytes < 10 * 1024 ** 2 ? 1 : 0)} MB`
}

function formatProjectFileModified(value: string) {
  const parsed = new Date(value)
  return Number.isNaN(parsed.getTime())
    ? ''
    : new Intl.DateTimeFormat(undefined, { dateStyle: 'medium', timeStyle: 'short' }).format(parsed)
}

export function ProjectFileRows({ files, loading, busy = false, emptyLabel = 'No FCS files', onDelete }: RowsProps) {
  const [deleteCandidate, setDeleteCandidate] = useState('')
  return (
    <div className="project-files-list" aria-busy={loading || busy}>
      {loading ? (
        <div className="project-files-empty">Loading files…</div>
      ) : files.length === 0 ? (
        <div className="project-files-empty"><strong>{emptyLabel}</strong></div>
      ) : files.map((entry) => (
        <div className="project-file-row" key={entry.name}>
          <File size={17} />
          <div className="project-file-copy">
            <strong title={entry.name}>{entry.name}</strong>
            <span>{formatProjectFileSize(entry.size)}{entry.modified ? ` · ${formatProjectFileModified(entry.modified)}` : ''}</span>
          </div>
          {deleteCandidate === entry.name ? (
            <div className="project-file-delete-confirm">
              <span>Delete?</span>
              <button type="button" disabled={busy} onClick={() => void onDelete(entry.name)}>Yes</button>
              <button type="button" disabled={busy} onClick={() => setDeleteCandidate('')}>No</button>
            </div>
          ) : (
            <button className="project-files-icon project-file-delete" type="button" disabled={busy} onClick={() => setDeleteCandidate(entry.name)} aria-label={`Delete ${entry.name}`}>
              <Trash2 size={15} />
            </button>
          )}
        </div>
      ))}
    </div>
  )
}

type InlineProps = {
  kind: ProjectFileKind
  projectPath: string
  directory: string
  refreshKey: string
  onChanged: () => void | Promise<void>
}

export function InlineProjectFiles({ kind, projectPath, directory, refreshKey, onChanged }: InlineProps) {
  const [files, setFiles] = useState<ProjectFileEntry[]>([])
  const [loading, setLoading] = useState(true)
  const [busy, setBusy] = useState(false)
  const [message, setMessage] = useState('')

  const refresh = useCallback(async () => {
    setLoading(true)
    const result = await listProjectFiles(kind, projectPath)
    setFiles(result.files)
    setMessage(result.success ? '' : (result.message ?? 'Files could not be loaded.'))
    setLoading(false)
  }, [kind, projectPath])

  useEffect(() => {
    let cancelled = false
    void listProjectFiles(kind, projectPath).then((result) => {
      if (cancelled) return
      setFiles(result.files)
      setMessage(result.success ? '' : (result.message ?? 'Files could not be loaded.'))
      setLoading(false)
    })
    return () => { cancelled = true }
  }, [kind, projectPath, refreshKey])

  async function remove(filename: string) {
    setBusy(true)
    const result = await deleteProjectFile(kind, filename, projectPath)
    if (result.success) {
      await refresh()
      await onChanged()
    }
    setMessage(result.success ? '' : result.message)
    setBusy(false)
  }

  return (
    <section className="surface-card project-inline-files" aria-label={`${kind === 'controls' ? 'Control' : 'Sample'} files`}>
      <div className="project-inline-files-path">{directory}/</div>
      <ProjectFileRows files={files} loading={loading} busy={busy} onDelete={remove} />
      {message && <div className="project-inline-files-message">{message}</div>}
    </section>
  )
}
