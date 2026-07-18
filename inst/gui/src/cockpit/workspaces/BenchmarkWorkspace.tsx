import { useEffect, useMemo, useState } from 'react'
import { Play, Plus } from 'lucide-react'
import { BenchmarkEntryCard } from '../components/BenchmarkEntryCard'
import { createBenchmark, listBenchmarks, runBenchmarkStage, saveBenchmark } from '../api'
import type { BenchmarkEntry, BenchmarkProject } from '../types'

type Props = { projectPath: string; onOpenReport: (path: string, stage: 'controls' | 'samples') => void }

function newEntry(index: number): BenchmarkEntry {
  return {
    entryId: `entry-${crypto.randomUUID()}`,
    label: `Entry ${index}`,
    enabled: true,
    method: 'Spectreasy',
    controlSettings: {}, sampleSettings: {},
    controlStatus: { state: 'not_run', message: '' }, sampleStatus: { state: 'not_run', message: '' },
    controlOutputs: { matrixFile: '', detectorNoiseFile: '', variantLibraryFile: '', reportFile: '', outputDir: '' },
    sampleOutputs: { reportFile: '', outputDir: '' },
  }
}

export function BenchmarkWorkspace({ projectPath, onOpenReport }: Props) {
  const [benchmark, setBenchmark] = useState<BenchmarkProject | null>(null)
  const [stage, setStage] = useState<'controls' | 'samples'>('controls')
  const [busy, setBusy] = useState(false)
  const [message, setMessage] = useState('')
  const [confirmation, setConfirmation] = useState<'run' | 'delete' | null>(null)
  const [deleteId, setDeleteId] = useState('')

  useEffect(() => {
    let cancelled = false
    void listBenchmarks(projectPath).then((items) => { if (!cancelled) setBenchmark(items[0] ?? null) })
    return () => { cancelled = true }
  }, [projectPath])

  const enabledCount = benchmark?.entries.filter((entry) => entry.enabled).length ?? 0
  const replaceLabels = useMemo(() => benchmark?.entries.filter((entry) => entry.enabled && (stage === 'controls' ? entry.controlOutputs.outputDir : entry.sampleOutputs.outputDir)).map((entry) => entry.label) ?? [], [benchmark, stage])

  async function persist(next: BenchmarkProject) {
    setBenchmark(next)
    const saved = await saveBenchmark(next, projectPath)
    if (saved) setBenchmark(saved)
    else setMessage('Benchmark configuration could not be saved.')
  }

  async function create() {
    setBusy(true)
    const created = await createBenchmark('Method comparison', projectPath)
    setBenchmark(created)
    setMessage(created ? '' : 'Benchmark configuration could not be created.')
    setBusy(false)
  }

  function updateEntry(entryId: string, patch: Partial<BenchmarkEntry>) {
    if (!benchmark) return
    const current = benchmark.entries.find((entry) => entry.entryId === entryId)
    if (!current) return
    const methodChanged = patch.method && patch.method !== current.method
    const controlSettingsChanged = patch.controlSettings !== undefined
    const sampleSettingsChanged = patch.sampleSettings !== undefined
    const next = {
      ...benchmark,
      entries: benchmark.entries.map((entry) => entry.entryId === entryId ? {
        ...entry, ...patch,
        controlStatus: (methodChanged || controlSettingsChanged) && entry.controlStatus.state === 'complete'
          ? { state: 'stale' as const, message: 'Configuration changed.' }
          : entry.controlStatus,
        sampleStatus: (methodChanged || controlSettingsChanged || sampleSettingsChanged) && entry.sampleStatus.state === 'complete'
          ? { state: 'stale' as const, message: 'Configuration changed.' }
          : entry.sampleStatus,
      } : entry),
    }
    void persist(next)
  }

  function requestRun() {
    if (!benchmark || enabledCount < 2) return
    if (replaceLabels.length) setConfirmation('run')
    else void run()
  }

  async function run() {
    if (!benchmark) return
    setConfirmation(null)
    setBusy(true)
    setBenchmark({ ...benchmark, entries: benchmark.entries.map((entry) => entry.enabled ? { ...entry, [stage === 'controls' ? 'controlStatus' : 'sampleStatus']: { state: 'queued', message: '' } } : entry) })
    const result = await runBenchmarkStage(stage, benchmark.benchmarkId, projectPath)
    if (result.benchmark) setBenchmark(result.benchmark)
    else {
      const items = await listBenchmarks(projectPath)
      setBenchmark(items.find((item) => item.benchmarkId === benchmark.benchmarkId) ?? benchmark)
    }
    setMessage(result.message)
    setBusy(false)
  }

  if (!benchmark) return (
    <section className="surface-card benchmark-empty">
      <button className="button button-primary" type="button" disabled={busy || !projectPath} onClick={() => void create()}><Plus size={15} />Create benchmark</button>
    </section>
  )

  return (
    <div className="benchmark-workspace">
      <header className="benchmark-header">
        <input aria-label="Benchmark name" value={benchmark.name} onChange={(event) => setBenchmark({ ...benchmark, name: event.target.value })} onBlur={() => void persist(benchmark)} />
        <button className="button" type="button" onClick={() => void persist({ ...benchmark, entries: [...benchmark.entries, newEntry(benchmark.entries.length + 1)] })}><Plus size={15} />Add entry</button>
        <button className="button button-primary" type="button" disabled={busy || enabledCount < 2} title={enabledCount < 2 ? 'Enable at least two entries.' : undefined} onClick={requestRun}><Play size={14} />Run {stage}</button>
      </header>
      <div className="subnav benchmark-tabs" role="tablist">
        <button className={stage === 'controls' ? 'is-active' : ''} role="tab" onClick={() => setStage('controls')}>01 Controls</button>
        <button className={stage === 'samples' ? 'is-active' : ''} role="tab" onClick={() => setStage('samples')}>02 Samples</button>
      </div>
      {message && <p className="benchmark-message">{message}</p>}
      <div className="benchmark-entry-grid">
        {benchmark.entries.map((entry) => (
          <BenchmarkEntryCard
            key={entry.entryId}
            entry={entry}
            stage={stage}
            onChange={(patch) => updateEntry(entry.entryId, patch)}
            onDelete={() => { setDeleteId(entry.entryId); setConfirmation('delete') }}
            onReport={(path) => onOpenReport(path, stage)}
          />
        ))}
      </div>
      {confirmation && (
        <div className="dialog-backdrop" role="presentation">
          <div className="dialog-card benchmark-confirmation" role="dialog" aria-modal="true" aria-labelledby="benchmark-confirmation-title">
            <h2 id="benchmark-confirmation-title">{confirmation === 'run' ? `Run benchmark ${stage}?` : 'Delete benchmark entry?'}</h2>
            {confirmation === 'run' ? <><p>Existing outputs will be replaced for:</p><ul>{replaceLabels.map((label) => <li key={label}>{label}</li>)}</ul></> : <p>The linked control and sample entry will be removed.</p>}
            <div className="dialog-actions">
              <button className="button" type="button" onClick={() => setConfirmation(null)}>Cancel</button>
              <button className="button button-primary" type="button" onClick={() => {
                if (confirmation === 'run') void run()
                else {
                  void persist({ ...benchmark, entries: benchmark.entries.filter((entry) => entry.entryId !== deleteId) })
                  setConfirmation(null)
                }
              }}>{confirmation === 'run' ? 'Run benchmark' : 'Delete entry'}</button>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}
