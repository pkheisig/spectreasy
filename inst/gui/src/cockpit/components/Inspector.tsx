import { Clipboard, Download, ExternalLink, FileWarning, MoreHorizontal, Pencil, ShieldCheck, Trash2 } from 'lucide-react'
import type { Artifact } from '../types'
import { StatusPill } from './StatusPill'

type InspectorProps = { artifact: Artifact | null; onDownload: (artifact: Artifact) => void }

export function Inspector({ artifact, onDownload }: InspectorProps) {
  if (!artifact) {
    return <aside className="inspector inspector-empty"><div className="empty-inspector-mark"><MoreHorizontal size={18} /></div><span className="eyebrow">Inspector</span><h3>Select an artifact</h3><p>Choose a file from the project index to see metadata, lineage, and safe actions.</p></aside>
  }
  return (
    <aside className="inspector">
      <div className="inspector-topline"><span className="eyebrow">Inspector</span><button className="icon-button" aria-label="More artifact actions"><MoreHorizontal size={17} /></button></div>
      <div className="inspector-file-icon"><FileWarning size={21} /></div>
      <h3>{artifact.name}</h3>
      <p className="inspector-type">{artifact.type} · {artifact.group}</p>
      <StatusPill state={artifact.status} />
      <div className="inspector-section"><span className="eyebrow">Lineage</span><dl><div><dt>Run</dt><dd>{artifact.run ?? 'User supplied'}</dd></div><div><dt>Updated</dt><dd>{artifact.updated}</dd></div><div><dt>Location</dt><dd className="path-value">{artifact.path}</dd></div></dl></div>
      {artifact.status === 'stale' && <div className="warning-box"><FileWarning size={15} /><span>This artifact is stale because an upstream input changed.</span></div>}
      <div className="inspector-actions"><button className="button button-primary" onClick={() => onDownload(artifact)}><Download size={15} /> Download</button><button className="button button-ghost"><ExternalLink size={15} /> Open</button></div>
      <div className="inspector-actions secondary"><button className="text-action"><Clipboard size={14} /> Copy path</button><button className="text-action"><Pencil size={14} /> Rename</button><button className="text-action danger"><Trash2 size={14} /> Delete</button></div>
      <div className="local-note"><ShieldCheck size={15} /><span>Protected local artifact. Deleting generated files requires confirmation.</span></div>
    </aside>
  )
}

