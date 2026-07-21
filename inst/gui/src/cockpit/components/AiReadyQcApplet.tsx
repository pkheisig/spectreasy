import { useCallback, useEffect, useMemo, useState } from 'react'
import { AlertTriangle, Clipboard, Download, ExternalLink, RefreshCw, ShieldCheck, Sparkles } from 'lucide-react'
import { downloadTextFile, generateAiQc, loadAiQc, loadProjectReportObjectUrl } from '../api'
import type { AiQcDetail, AiQcGrade, AiQcPrivacy, AiQcResponse, AiQcScope } from '../types'

const emptyGrades: Record<AiQcGrade, number> = { good: 0, review: 0, poor: 0, not_graded: 0 }

export default function AiReadyQcApplet({ theme, projectPath, outputRoot }: { theme: 'light' | 'dark'; projectPath: string; outputRoot: string }) {
  const [tab, setTab] = useState<'overview' | 'findings' | 'settings' | 'prompt'>('overview')
  const [data, setData] = useState<AiQcResponse | null>(null)
  const [loading, setLoading] = useState(true)
  const [scope, setScope] = useState<AiQcScope>('combined')
  const [privacy, setPrivacy] = useState<AiQcPrivacy>('standard')
  const [detail, setDetail] = useState<AiQcDetail>('standard')
  const [reference, setReference] = useState('auto')
  const [context, setContext] = useState('')
  const [gradeFilter, setGradeFilter] = useState<AiQcGrade | 'all'>('all')
  const [copied, setCopied] = useState(false)

  const refresh = useCallback(async () => {
    setLoading(true)
    const next = await loadAiQc(projectPath, outputRoot, detail)
    setData(next)
    setScope(next.defaultScope)
    setPrivacy(next.privacyMode)
    setLoading(false)
  }, [detail, outputRoot, projectPath])

  useEffect(() => {
    let active = true
    void loadAiQc(projectPath, outputRoot, detail).then((next) => {
      if (!active) return
      setData(next)
      setScope(next.defaultScope)
      setPrivacy(next.privacyMode)
      setLoading(false)
    })
    return () => { active = false }
  }, [detail, outputRoot, projectPath])

  async function generate() {
    setLoading(true)
    const next = await generateAiQc({ projectPath, outputRoot, scope, privacy, detail, reference, context })
    setData(next)
    setLoading(false)
    if (!next.error) setTab('overview')
  }

  async function copyPrompt() {
    const prompt = data?.prompt ?? ''
    if (!prompt) return
    try {
      await navigator.clipboard.writeText(prompt)
    } catch {
      const area = document.createElement('textarea'); area.value = prompt; area.style.position = 'fixed'; area.style.opacity = '0'; document.body.appendChild(area); area.select(); document.execCommand('copy'); area.remove()
    }
    setCopied(true); window.setTimeout(() => setCopied(false), 1400)
  }

  async function openArtifact(path: string) {
    const url = await loadProjectReportObjectUrl(path, projectPath)
    window.open(url, '_blank', 'noopener,noreferrer')
    window.setTimeout(() => URL.revokeObjectURL(url), 60000)
  }

  const findings = useMemo(() => (data?.findings ?? []).filter((finding) => gradeFilter === 'all' || finding.grade === gradeFilter), [data?.findings, gradeFilter])
  const counts = data?.gradeCounts ?? emptyGrades
  const largePrompt = (data?.estimatedTokens ?? 0) > 12000

  return <main className={`ai-qc-applet theme-${theme}`}>
    <header className="ai-qc-header"><div><span className="eyebrow">Local scientific review</span><h1>AI-ready QC</h1><p>Deterministic evidence, explicit grading provenance, and a paste-ready prompt. Nothing is sent by Spectreasy.</p></div><div className="ai-qc-header-actions"><span className={`ai-qc-status is-${data?.status ?? 'not_generated'}`}>{loading ? 'generating' : data?.status?.replace('_', ' ') ?? 'not generated'}</span><button className="button button-ghost" onClick={() => void refresh()} disabled={loading}><RefreshCw size={14} /> Refresh</button><button className="button button-primary" onClick={() => void generate()} disabled={loading}><Sparkles size={14} /> Generate</button></div></header>
    <div className="ai-qc-local-notice"><ShieldCheck size={16} /><strong>{privacy}</strong> privacy active · local files only · no event-level data · no provider traffic</div>
    {data?.stale && <div className="ai-qc-warning"><AlertTriangle size={16} /> Source artifacts changed after this export. Refresh before interpretation.</div>}
    {data?.error && <div className="ai-qc-warning"><AlertTriangle size={16} /> {data.error}</div>}
    <nav className="ai-qc-tabs" aria-label="AI-ready QC views">{(['overview', 'findings', 'settings', 'prompt'] as const).map((item) => <button key={item} className={tab === item ? 'is-active' : ''} onClick={() => setTab(item)}>{item === 'settings' ? 'Data settings' : item}</button>)}</nav>
    <section className="ai-qc-body">
      {tab === 'overview' && <><div className="ai-qc-grade-grid">{(Object.keys(counts) as AiQcGrade[]).map((grade) => <button key={grade} className={`ai-qc-grade-card is-${grade}`} onClick={() => { setGradeFilter(grade); setTab('findings') }}><span>{grade.replace('_', ' ')}</span><strong>{counts[grade]}</strong></button>)}</div><div className="ai-qc-overview-grid"><article><h2>Readiness</h2><dl><div><dt>Scope</dt><dd>{data?.scope ?? scope}</dd></div><div><dt>Status</dt><dd>{data?.status ?? 'not generated'}</dd></div><div><dt>Schema</dt><dd>{data?.schemaVersion ? `${data.schemaName} ${data.schemaVersion}` : 'not generated'}</dd></div><div><dt>Profile</dt><dd>{data?.profileName ? `${data.profileName} ${data.profileVersion} (n=${data.referenceN})` : 'none compatible'}</dd></div><div><dt>Prompt</dt><dd>{(data?.promptCharacters ?? 0).toLocaleString()} chars · ~{(data?.estimatedTokens ?? 0).toLocaleString()} tokens</dd></div></dl>{data?.missingSections.length ? <p className="ai-qc-muted">Partial evidence: {data.missingSections.join(', ')}</p> : <p className="ai-qc-muted">All collected sections are available.</p>}</article><article><h2>Artifacts</h2>{data?.artifactPaths.length ? <ul className="ai-qc-artifacts">{data.artifactPaths.map((path) => <li key={path}><button onClick={() => void openArtifact(path)}>{path}<ExternalLink size={13} /></button></li>)}</ul> : <p className="ai-qc-muted">Generate an export to create JSON, text, Markdown, prompt, and manifest files.</p>}</article></div></>}
      {tab === 'findings' && <><div className="ai-qc-filter-row"><label>Grade<select value={gradeFilter} onChange={(event) => setGradeFilter(event.target.value as AiQcGrade | 'all')}><option value="all">All grades</option>{(Object.keys(emptyGrades) as AiQcGrade[]).map((grade) => <option key={grade} value={grade}>{grade.replace('_', ' ')}</option>)}</select></label></div><div className="ai-qc-findings">{findings.length ? findings.map((finding, index) => <details key={`${finding.metricId}-${finding.entity}-${index}`}><summary><span className={`ai-qc-dot is-${finding.grade}`} /> <strong>{finding.metricId}</strong><span>{finding.entity}</span><em>{finding.grade.replace('_', ' ')}</em></summary><p>{finding.explanation || 'No additional explanation was recorded.'}</p><small>{finding.stage}</small></details>) : <p className="ai-qc-muted">No findings match this filter.</p>}</div></>}
      {tab === 'settings' && <div className="ai-qc-settings"><label>Scope<select value={scope} onChange={(event) => setScope(event.target.value as AiQcScope)}>{(['control', 'sample', 'combined'] as AiQcScope[]).map((item) => <option key={item} value={item} disabled={Boolean(data?.availableScopes.length) && !data?.availableScopes.includes(item)}>{item}</option>)}</select></label><label>Privacy<select value={privacy} onChange={(event) => setPrivacy(event.target.value as AiQcPrivacy)}><option value="standard">Standard</option><option value="strict">Strict</option><option value="none">None — explicit identifiers</option></select></label><label>Detail<select value={detail} onChange={(event) => setDetail(event.target.value as AiQcDetail)}><option value="compact">Compact</option><option value="standard">Standard</option><option value="full">Full</option></select></label><label>Reference profile<input value={reference} onChange={(event) => setReference(event.target.value)} placeholder="auto" /></label><label className="ai-qc-context">Optional analyst context<textarea value={context} onChange={(event) => setContext(event.target.value)} rows={5} placeholder="Context is included in the prompt, not used to alter measurements or grades." /></label><p className="ai-qc-muted">Standard aliases samples and removes absolute paths. Strict also aliases controls and removes project/operator/date/free-text metadata. None must be selected explicitly.</p></div>}
      {tab === 'prompt' && <div className="ai-qc-prompt-view">{largePrompt && <div className="ai-qc-warning"><AlertTriangle size={16} /> Large prompt: approximately {(data?.estimatedTokens ?? 0).toLocaleString()} tokens. Use compact detail if needed.</div>}<div className="ai-qc-prompt-toolbar"><span>{(data?.promptCharacters ?? 0).toLocaleString()} characters · ~{(data?.estimatedTokens ?? 0).toLocaleString()} tokens</span><button className="button button-ghost" onClick={() => void copyPrompt()} disabled={!data?.prompt}><Clipboard size={14} /> {copied ? 'Copied' : 'Copy'}</button><button className="button button-ghost" onClick={() => downloadTextFile('spectreasy_ai_qc_prompt.txt', data?.prompt ?? '')} disabled={!data?.prompt}><Download size={14} /> Save</button></div><pre>{data?.prompt || 'Generate AI-ready QC to preview the prompt.'}</pre></div>}
    </section>
  </main>
}
