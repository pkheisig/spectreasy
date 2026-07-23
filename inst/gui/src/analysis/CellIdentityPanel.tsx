import { useMemo, useState, type CSSProperties } from 'react'
import { CircleHelp, ListFilter, Plus, Settings, Sparkles, Trash2 } from 'lucide-react'
import { analysisRequest } from './api'
import { commonImmuneIdentityTemplate } from './identityTemplates'
import type {
  AnalysisIdentityAnnotation,
  AnalysisIdentitySignature,
  AnalysisRunEvent,
  AnalysisRunResult,
} from './types'

const IDENTITY_COLORS = ['#197783', '#d06d32', '#7d58a6', '#4f8f45', '#bc4d65', '#416fae', '#9c7427', '#607d8b']

function newIdentity(index: number): AnalysisIdentitySignature {
  return {
    id: `identity-${crypto.randomUUID()}`,
    name: `Identity ${index}`,
    color: IDENTITY_COLORS[(index - 1) % IDENTITY_COLORS.length],
    positive_markers: [],
    negative_markers: [],
  }
}

export function CellIdentityPanel({
  projectPath,
  result,
  onAnnotated,
}: {
  projectPath: string
  result: AnalysisRunResult
  onAnnotated: (annotation: AnalysisIdentityAnnotation, events: Array<Pick<AnalysisRunEvent, 'event_id' | 'predicted_identity' | 'identity_score' | 'identity_margin'>>) => void
}) {
  const markers = result.metadata.marker_columns?.map((item) => item.marker) ?? result.metadata.markers
  const [signatures, setSignatures] = useState<AnalysisIdentitySignature[]>(() => [newIdentity(1), newIdentity(2)])
  const [minScore, setMinScore] = useState(0.55)
  const [minMargin, setMinMargin] = useState(0.08)
  const [evidenceSensitivity, setEvidenceSensitivity] = useState(1)
  const [running, setRunning] = useState(false)
  const [message, setMessage] = useState('')
  const [guideOpen, setGuideOpen] = useState(false)
  const [settingsOpen, setSettingsOpen] = useState(false)
  const [topMarkerCount, setTopMarkerCount] = useState(8)
  const [minimumMarkerAuc, setMinimumMarkerAuc] = useState(0.55)
  const [markerRows, setMarkerRows] = useState<Array<{
    cluster_id: number
    marker: string
    auc: number
    robust_effect: number
    p_adjusted: number
  }>>([])
  const hasClusters = result.events.some((event) => Number.isFinite(event.cluster_id))
  const immuneTemplate = useMemo(() => commonImmuneIdentityTemplate(markers), [markers])

  const issue = useMemo(() => {
    if (signatures.length < 2) return 'Define at least two identities.'
    const names = signatures.map((signature) => signature.name.trim().toLowerCase())
    if (names.some((name) => !name)) return 'Name every identity.'
    if (new Set(names).size !== names.length) return 'Identity names must be unique.'
    const empty = signatures.find((signature) => !signature.positive_markers.length && !signature.negative_markers.length)
    if (empty) return `Choose markers for ${empty.name}.`
    return ''
  }, [signatures])

  function updateSignature(id: string, patch: Partial<AnalysisIdentitySignature>) {
    setSignatures((current) => current.map((signature) => signature.id === id ? { ...signature, ...patch } : signature))
  }

  function toggleMarker(id: string, marker: string, role: 'positive' | 'negative') {
    setSignatures((current) => current.map((signature) => {
      if (signature.id !== id) return signature
      const target = role === 'positive' ? signature.positive_markers : signature.negative_markers
      const selected = target.includes(marker)
      return {
        ...signature,
        positive_markers: role === 'positive'
          ? (selected ? target.filter((value) => value !== marker) : [...target, marker])
          : signature.positive_markers.filter((value) => value !== marker),
        negative_markers: role === 'negative'
          ? (selected ? target.filter((value) => value !== marker) : [...target, marker])
          : signature.negative_markers.filter((value) => value !== marker),
      }
    }))
  }

  function applyImmuneTemplate() {
    if (immuneTemplate.length < 2) return
    setSignatures(immuneTemplate)
    setMessage(`${immuneTemplate.length} editable immune identity patterns loaded`)
  }

  async function annotate() {
    if (issue || running) return
    setRunning(true)
    setMessage('')
    try {
      const response = await analysisRequest<{
        success: boolean
        result: {
          metadata: AnalysisIdentityAnnotation
          events: Array<Pick<AnalysisRunEvent, 'event_id' | 'predicted_identity' | 'identity_score' | 'identity_margin'>>
        }
      }>('/analysis/annotate', projectPath, {
        method: 'POST',
        body: JSON.stringify({
          analysisId: result.metadata.analysis_id,
          signatures: signatures.map(({ name, color, positive_markers, negative_markers }) => ({
            name: name.trim(), color, positive_markers, negative_markers,
          })),
          minScore,
          minMargin,
          evidenceSensitivity,
        }),
      })
      onAnnotated(response.result.metadata, response.result.events)
      setMessage(`${response.result.metadata.assigned_count.toLocaleString()} assigned · ${response.result.metadata.unassigned_count.toLocaleString()} unassigned`)
    } catch (reason) {
      setMessage(reason instanceof Error ? reason.message : String(reason))
    } finally {
      setRunning(false)
    }
  }

  async function discoverMarkers() {
    if (!hasClusters || running) return
    setRunning(true)
    setMessage('')
    try {
      const response = await analysisRequest<{
        success: boolean
        result: { markers: typeof markerRows }
      }>('/analysis/markers', projectPath, {
        method: 'POST',
        body: JSON.stringify({
          analysisId: result.metadata.analysis_id,
          topN: topMarkerCount,
          minimumAuc: minimumMarkerAuc,
        }),
      })
      setMarkerRows(response.result.markers)
      setMessage(`${response.result.markers.length.toLocaleString()} distinguishing cluster markers found`)
    } catch (reason) {
      setMessage(reason instanceof Error ? reason.message : String(reason))
    } finally {
      setRunning(false)
    }
  }

  return (
    <aside className="analysis-identity-panel" aria-label="Cell identity annotation">
      <header>
        <div><strong>Cell identities</strong><small>Marker-pattern matching</small></div>
        <div className="analysis-identity-header-actions">
          <button type="button" className="analysis-icon-button" aria-label="Cell identity guide" onClick={() => setGuideOpen((open) => !open)}><CircleHelp size={15} /></button>
          <button type="button" className="analysis-icon-button" aria-label="Cell identity settings" aria-expanded={settingsOpen} onClick={() => setSettingsOpen((open) => !open)}><Settings size={15} /></button>
        </div>
        {settingsOpen ? (
          <div className="analysis-identity-settings" role="dialog" aria-label="Cell identity settings menu">
            <strong>How strict should assignment be?</strong>
            <p>A cell must pass both confidence checks. Otherwise it stays Unassigned.</p>
            <label>
              <span>Minimum identity match<small>How well the cell must fit its best marker pattern. Raise this to accept only clearer matches.</small></span>
              <input aria-label="Minimum identity match" type="number" min="0" max="1" step="0.01" value={minScore} onChange={(event) => setMinScore(Math.max(0, Math.min(1, Number(event.target.value))))} />
            </label>
            <label>
              <span>Minimum lead over second choice<small>How far the best identity must beat the next best. Raise this to leave more ambiguous cells Unassigned.</small></span>
              <input aria-label="Minimum lead over second choice" type="number" min="0" max="1" step="0.01" value={minMargin} onChange={(event) => setMinMargin(Math.max(0, Math.min(1, Number(event.target.value))))} />
            </label>
            <label>
              <span>Marker sensitivity<small>How strongly expression above or below the population’s robust center counts. 1 is balanced; higher values emphasize smaller differences.</small></span>
              <input aria-label="Marker sensitivity" type="number" min="0.25" max="4" step="0.05" value={evidenceSensitivity} onChange={(event) => setEvidenceSensitivity(Math.max(0.25, Math.min(4, Number(event.target.value))))} />
            </label>
            <strong>Cluster marker discovery</strong>
            <label>
              <span>Markers shown per cluster<small>Maximum number of the strongest positive cluster-versus-rest markers to show.</small></span>
              <input aria-label="Markers shown per cluster" type="number" min="1" max="100" step="1" value={topMarkerCount} onChange={(event) => setTopMarkerCount(Math.max(1, Math.min(100, Number(event.target.value))))} />
            </label>
            <label>
              <span>Minimum separation score<small>Rank AUC: 0.5 means no separation and 1 means perfect separation. Raise this to keep only clearer markers.</small></span>
              <input aria-label="Minimum marker separation score" type="number" min="0.5" max="1" step="0.01" value={minimumMarkerAuc} onChange={(event) => setMinimumMarkerAuc(Math.max(0.5, Math.min(1, Number(event.target.value))))} />
            </label>
          </div>
        ) : null}
      </header>

      {guideOpen ? (
        <div className="analysis-identity-guide">
          <strong>How it works</strong>
          <p>For each identity, mark proteins expected to be high as <strong>POS</strong> and exclusion proteins expected to be low as <strong>NEG</strong>. Spectreasy compares every cell with those patterns after robustly centering each marker in this population.</p>
          <p>The best-matching identity is assigned only when it is both strong enough and clearly ahead of the second choice. Coordinates never determine the label. Verify assignments by coloring the result by identity and by its defining markers.</p>
          <p><strong>Discover cluster markers</strong> compares each cluster with every other cell using rank AUC and a robust effect size. It explains clusters; it does not automatically claim a biological identity.</p>
          <p><strong>Load immune template</strong> is a panel-aware starting point for CD4 T, CD8 T, B, NK, and monocyte patterns. Only identities supported by the markers in this result are offered. Review and edit every pattern for the tissue and panel before annotation.</p>
        </div>
      ) : null}

      <div className="analysis-identity-list">
        {signatures.map((signature) => (
          <section key={signature.id} className="analysis-identity-card" style={{ '--identity-color': signature.color } as CSSProperties}>
            <div className="analysis-identity-name">
              <input type="color" aria-label={`${signature.name} color`} value={signature.color} onChange={(event) => updateSignature(signature.id, { color: event.target.value })} />
              <input aria-label="Identity name" value={signature.name} onChange={(event) => updateSignature(signature.id, { name: event.target.value })} />
              <button type="button" aria-label={`Remove ${signature.name}`} disabled={signatures.length <= 2} onClick={() => setSignatures((current) => current.filter((item) => item.id !== signature.id))}><Trash2 size={13} /></button>
            </div>
            <div className="analysis-identity-marker-row">
              <span>POS</span>
              <div>{markers.map((marker) => <button type="button" key={marker} className={signature.positive_markers.includes(marker) ? 'is-positive' : ''} onClick={() => toggleMarker(signature.id, marker, 'positive')}>{marker}</button>)}</div>
            </div>
            <div className="analysis-identity-marker-row">
              <span>NEG</span>
              <div>{markers.map((marker) => <button type="button" key={marker} className={signature.negative_markers.includes(marker) ? 'is-negative' : ''} onClick={() => toggleMarker(signature.id, marker, 'negative')}>{marker}</button>)}</div>
            </div>
          </section>
        ))}
      </div>

      <button
        type="button"
        className="analysis-identity-discover"
        disabled={immuneTemplate.length < 2 || running}
        title={immuneTemplate.length >= 2 ? 'Load editable identities supported by this marker panel' : 'At least two supported immune identities are required'}
        onClick={applyImmuneTemplate}
      ><Sparkles size={13} /> Load immune template ({immuneTemplate.length})</button>
      <button type="button" className="analysis-identity-discover" disabled={!hasClusters || running} title={hasClusters ? 'Rank markers for each cluster' : 'Run clustering first'} onClick={() => void discoverMarkers()}><ListFilter size={13} /> Discover cluster markers</button>
      {markerRows.length ? (
        <div className="analysis-cluster-marker-results">
          <table>
            <thead><tr><th>Cluster</th><th>Marker</th><th>AUC</th><th>Effect</th></tr></thead>
            <tbody>{markerRows.map((row) => <tr key={`${row.cluster_id}-${row.marker}`}><td>{row.cluster_id}</td><td>{row.marker}</td><td>{Number(row.auc).toFixed(2)}</td><td>{Number(row.robust_effect).toFixed(2)}</td></tr>)}</tbody>
          </table>
        </div>
      ) : null}

      <button type="button" className="analysis-identity-add" onClick={() => setSignatures((current) => [...current, newIdentity(current.length + 1)])}><Plus size={13} /> Add identity</button>

      <div className={`analysis-identity-status ${issue ? 'is-blocked' : ''}`} role="status">{issue || message || 'Ready to annotate'}</div>
      <button type="button" className="analysis-primary" disabled={Boolean(issue) || running} onClick={() => void annotate()}><Sparkles size={14} /> {running ? 'Annotating…' : 'Annotate cells'}</button>
    </aside>
  )
}
