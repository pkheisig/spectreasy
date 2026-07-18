import { FileText, Trash2 } from 'lucide-react'
import { GuiSelect } from './GuiSelect'
import { StatusPill } from './StatusPill'
import type { BenchmarkEntry, BenchmarkSettingValue, StepState } from '../types'

type Definition = { key: string; label: string; type?: 'number' | 'boolean'; step?: number; defaultValue: BenchmarkSettingValue }

const commonControl: Definition[] = [
  { key: 'af_n_bands', label: 'AF band count', type: 'number', defaultValue: 100 },
  { key: 'scc_background_method', label: 'SCC background method', defaultValue: 'scatter_knn' },
  { key: 'scc_background_k', label: 'Background neighbours', type: 'number', defaultValue: 2 },
  { key: 'spectral_variant_som_nodes', label: 'Variant SOM nodes', type: 'number', defaultValue: 16 },
  { key: 'spectral_variant_top_k', label: 'Variant top-K', type: 'number', defaultValue: 3 },
  { key: 'spectral_variant_cosine_threshold', label: 'Variant cosine threshold', type: 'number', step: 0.01, defaultValue: 0.98 },
  { key: 'spectral_variant_max_variants', label: 'Maximum variants', type: 'number', defaultValue: 8 },
  { key: 'spectral_variant_min_events', label: 'Minimum variant events', type: 'number', defaultValue: 50 },
]
const spectreasyControl: Definition[] = [{ key: 'spectreasy_weight_quantile', label: 'Spectreasy weight quantile', type: 'number', step: 0.01, defaultValue: 0.65 }, ...commonControl]
const autospectralControl: Definition[] = [...commonControl,
  { key: 'autospectral_n_candidates', label: 'Candidate events', type: 'number', defaultValue: 1000 },
  { key: 'autospectral_n_spectral', label: 'Spectrum events', type: 'number', defaultValue: 200 },
  { key: 'autospectral_min_events', label: 'Minimum selector events', type: 'number', defaultValue: 10 },
  { key: 'autospectral_refine', label: 'Refine AF bank', type: 'boolean', defaultValue: false },
]
const sampleFields: Definition[] = [
  { key: 'spectreasy_weight_quantile', label: 'Spectreasy weight quantile', type: 'number', step: 0.01, defaultValue: 0.65 },
  { key: 'spectral_variant_top_k', label: 'Variant top-K', type: 'number', defaultValue: 3 },
  { key: 'spectral_variant_min_abundance', label: 'Variant minimum abundance', type: 'number', defaultValue: 1 },
  { key: 'spectral_variant_positive_fraction', label: 'Variant positive fraction', type: 'number', step: 0.01, defaultValue: 0.02 },
  { key: 'spectral_variant_min_improvement', label: 'Variant minimum improvement', type: 'number', step: 0.01, defaultValue: 0.01 },
  { key: 'estimate_af', label: 'Estimate AF', type: 'boolean', defaultValue: false },
  { key: 'write_fcs', label: 'Write FCS', type: 'boolean', defaultValue: true },
]

function statusState(state: BenchmarkEntry['controlStatus']['state']): StepState {
  if (state === 'complete') return 'complete'
  if (state === 'failed') return 'blocked'
  if (state === 'stale') return 'stale'
  if (state === 'running' || state === 'queued') return 'ready'
  return 'idle'
}

type Props = {
  entry: BenchmarkEntry
  stage: 'controls' | 'samples'
  onChange: (patch: Partial<BenchmarkEntry>) => void
  onDelete: () => void
  onReport: (path: string) => void
}

export function BenchmarkEntryCard({ entry, stage, onChange, onDelete, onReport }: Props) {
  const settingsKey = stage === 'controls' ? 'controlSettings' : 'sampleSettings'
  const settings = entry[settingsKey]
  const definitions = stage === 'controls'
    ? entry.method === 'Spectreasy' ? spectreasyControl : entry.method === 'AutoSpectral' ? autospectralControl : []
    : sampleFields.filter((field) => entry.method === 'Spectreasy' || field.key !== 'spectreasy_weight_quantile')
  const advanced = stage === 'controls'
    ? [{ key: 'n_threads', label: 'Threads', type: 'number' as const, defaultValue: 1 }, { key: 'rwls_max_iter', label: 'RWLS iterations', type: 'number' as const, defaultValue: 1 }, { key: 'seed', label: 'Seed', type: 'number' as const, defaultValue: 1 }]
    : [{ key: 'n_threads', label: 'Threads', type: 'number' as const, defaultValue: 1 }, { key: 'chunk_size', label: 'Chunk size', type: 'number' as const, defaultValue: 50000 }, { key: 'seed', label: 'Seed', type: 'number' as const, defaultValue: 1 }]
  const status = stage === 'controls' ? entry.controlStatus : entry.sampleStatus
  const outputs = stage === 'controls' ? entry.controlOutputs : entry.sampleOutputs
  const reportFile = outputs.reportFile

  const updateSetting = (key: string, value: BenchmarkSettingValue) => onChange({ [settingsKey]: { ...settings, [key]: value } })
  const renderField = (field: Definition) => field.type === 'boolean' ? (
    <label className="toggle-label benchmark-toggle" key={field.key}>
      <input type="checkbox" checked={Boolean(settings[field.key] ?? field.defaultValue)} onChange={(event) => updateSetting(field.key, event.target.checked)} />
      <span className="toggle-ui" />
      <span>{field.label}</span>
    </label>
  ) : (
    <label key={field.key}>
      {field.label}
      <input
        type={field.type === 'number' ? 'number' : 'text'}
        step={field.step}
        value={String(settings[field.key] ?? field.defaultValue)}
        onChange={(event) => updateSetting(field.key, field.type === 'number' ? Number(event.target.value) : event.target.value)}
      />
    </label>
  )

  return (
    <article className={`surface-card benchmark-entry-card ${entry.enabled ? '' : 'is-disabled'}`}>
      <header>
        <input className="benchmark-entry-name" aria-label="Entry name" value={entry.label} onChange={(event) => onChange({ label: event.target.value })} />
        <StatusPill state={statusState(status.state)} label={status.state.replaceAll('_', ' ')} compact />
      </header>
      <div className="benchmark-entry-toolbar">
        <GuiSelect value={entry.method} onChange={(event) => onChange({ method: event.target.value })} aria-label="Method">
          <option>Spectreasy</option><option>AutoSpectral</option><option>OLS</option><option>WLS</option><option>RWLS</option><option>NNLS</option>
        </GuiSelect>
        <label className="toggle-label"><input type="checkbox" checked={entry.enabled} onChange={(event) => onChange({ enabled: event.target.checked })} /><span className="toggle-ui" /><span>Enabled</span></label>
        <button className="icon-button" type="button" aria-label={`Delete ${entry.label}`} title="Delete entry" onClick={onDelete}><Trash2 size={14} /></button>
      </div>
      {stage === 'samples' && (
        <dl className="benchmark-inherited">
          <div><dt>Control matrix</dt><dd>{entry.controlOutputs.matrixFile || 'Blocked until controls complete'}</dd></div>
          <div><dt>Variant library</dt><dd>{entry.controlOutputs.variantLibraryFile || 'Not available'}</dd></div>
        </dl>
      )}
      <div className="benchmark-fields">{definitions.map(renderField)}</div>
      <details className="benchmark-advanced"><summary>Advanced</summary><div className="benchmark-fields">{advanced.map(renderField)}</div></details>
      <footer>
        <span>{status.message || status.state.replaceAll('_', ' ')}</span>
        {reportFile && <button className="text-action" type="button" onClick={() => onReport(reportFile)}><FileText size={14} />Report</button>}
        {outputs.outputDir && <code title={outputs.outputDir}>{outputs.outputDir}</code>}
      </footer>
    </article>
  )
}
