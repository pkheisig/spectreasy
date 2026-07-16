import { lazy, Suspense, useCallback, useEffect, useMemo, useRef, useState } from 'react'
import {
  Check,
  Clipboard,
  FlaskConical,
  MonitorCog,
  RefreshCw,
  RotateCcw,
  TerminalSquare,
} from 'lucide-react'
import { resolveApiBase } from '../apiBase'
import { installPackageCommand, launchCommand } from './installCommands'
import './setup.css'
import './sandbox.css'

const CockpitApp = lazy(() => import('../cockpit/CockpitApp'))
const sandboxApi = resolveApiBase()
const installCommand = installPackageCommand

type RState = 'missing' | 'outdated' | 'broken' | 'ready'
type PackageState = 'missing' | 'unknown' | 'broken' | 'ready'
type PreviewKey = 'missing-r' | 'outdated-r' | 'r-only' | 'broken-package' | 'ready'
type Environment = {
  r: RState
  rVersion: string | null
  spectreasy: PackageState
  spectreasyVersion: string | null
  detail: string
}
type SandboxJob = {
  busy: boolean
  action: string
  phase: string
  error: string
  log: string[]
}
type SandboxSnapshot = {
  environment: Environment
  job: SandboxJob
  paths: { sandboxRoot: string; environmentPrefix: string; repositoryRoot: string }
  minimumRVersion: string
  supported: boolean
}

const previews: Record<PreviewKey, { label: string; environment: Environment }> = {
  'missing-r': {
    label: 'No R',
    environment: { r: 'missing', rVersion: null, spectreasy: 'missing', spectreasyVersion: null, detail: 'No R installation was detected.' },
  },
  'outdated-r': {
    label: 'Old R',
    environment: { r: 'outdated', rVersion: '4.4.3', spectreasy: 'missing', spectreasyVersion: null, detail: 'Spectreasy requires R 4.5.0 or newer.' },
  },
  'r-only': {
    label: 'R only',
    environment: { r: 'ready', rVersion: '4.5.3', spectreasy: 'missing', spectreasyVersion: null, detail: 'R is compatible; Spectreasy is not installed.' },
  },
  'broken-package': {
    label: 'Dependency issue',
    environment: { r: 'ready', rVersion: '4.5.3', spectreasy: 'broken', spectreasyVersion: null, detail: 'Spectreasy is present, but one or more dependencies cannot be loaded.' },
  },
  ready: {
    label: 'Ready',
    environment: { r: 'ready', rVersion: '4.5.3', spectreasy: 'ready', spectreasyVersion: '1.0.0', detail: 'R and the Spectreasy namespace are ready.' },
  },
}

function CopyButton({ value, label }: { value: string; label: string }) {
  const [copied, setCopied] = useState(false)
  const timer = useRef<number | null>(null)

  useEffect(() => () => {
    if (timer.current !== null) window.clearTimeout(timer.current)
  }, [])

  const copy = async () => {
    await navigator.clipboard.writeText(value)
    setCopied(true)
    if (timer.current !== null) window.clearTimeout(timer.current)
    timer.current = window.setTimeout(() => setCopied(false), 1600)
  }

  return (
    <button className="setup-copy-button" type="button" onClick={() => void copy()}>
      {copied ? <Check size={14} /> : <Clipboard size={14} />}
      {copied ? 'Copied' : label}
    </button>
  )
}

function stageFor(environment: Environment): PreviewKey {
  if (environment.r === 'missing' || environment.r === 'broken') return 'missing-r'
  if (environment.r === 'outdated') return 'outdated-r'
  if (environment.spectreasy === 'broken') return 'broken-package'
  if (environment.spectreasy !== 'ready') return 'r-only'
  return 'ready'
}

function SandboxControls({
  source,
  preview,
  snapshot,
  holdAtReady,
  onSource,
  onPreview,
  onAction,
  onHoldAtReady,
}: {
  source: 'real' | 'preview'
  preview: PreviewKey
  snapshot: SandboxSnapshot | null
  holdAtReady: boolean
  onSource: (source: 'real' | 'preview') => void
  onPreview: (preview: PreviewKey) => void
  onAction: (action: string) => void
  onHoldAtReady: (value: boolean) => void
}) {
  const job = snapshot?.job
  return (
    <aside className="sandbox-controls" aria-label="Sandbox controls">
      <div className="sandbox-controls-heading">
        <FlaskConical size={17} />
        <div><strong>Onboarding sandbox</strong><span>Host R is never inspected</span></div>
      </div>

      <div className="sandbox-segmented" role="group" aria-label="Sandbox source">
        <button className={source === 'real' ? 'is-active' : ''} type="button" onClick={() => onSource('real')}>Real sandbox</button>
        <button className={source === 'preview' ? 'is-active' : ''} type="button" onClick={() => onSource('preview')}>Preview states</button>
      </div>

      {source === 'real' ? (
        <div className="sandbox-control-section">
          <span className="sandbox-control-label">Isolated environment</span>
          <p>Install a real R runtime and package into the isolated cache. A first full install can take several minutes and use about 2 GB.</p>
          <div className="sandbox-presets">
            <button type="button" disabled={job?.busy} onClick={() => onAction('reset')}>Empty</button>
            <button type="button" disabled={job?.busy} onClick={() => onAction('install-old-r')}>Old R</button>
            <button type="button" disabled={job?.busy} onClick={() => onAction('prepare-r-only')}>R only</button>
            <button type="button" disabled={job?.busy} onClick={() => onAction('break-spectreasy')}>Dependency issue</button>
            <button type="button" disabled={job?.busy} onClick={() => onAction('prepare-ready')}>Ready</button>
          </div>
          <div className="sandbox-path"><span>Prefix</span><code>{snapshot?.paths.environmentPrefix ?? 'Connecting…'}</code></div>
          <button className="sandbox-purge-button" type="button" disabled={job?.busy} onClick={() => onAction('purge')}>
            <RotateCcw size={12} /> Delete all sandbox files
          </button>
        </div>
      ) : (
        <div className="sandbox-control-section">
          <span className="sandbox-control-label">Instant UI state</span>
          <p>No downloads or filesystem changes.</p>
          <div className="sandbox-preview-list">
            {(Object.entries(previews) as Array<[PreviewKey, (typeof previews)[PreviewKey]]>).map(([key, item]) => (
              <button className={preview === key ? 'is-active' : ''} type="button" key={key} onClick={() => onPreview(key)}>
                <span>{item.label}</span><small>{key === 'ready' ? 'R + package' : item.environment.detail}</small>
              </button>
            ))}
          </div>
        </div>
      )}

      <label className="sandbox-ready-toggle">
        <input type="checkbox" checked={holdAtReady} onChange={(event) => onHoldAtReady(event.target.checked)} />
        <span><strong>Hold at ready screen</strong><small>Otherwise readiness opens the cockpit automatically.</small></span>
      </label>

      {job && (job.busy || job.log.length > 0 || job.error) ? (
        <div className={`sandbox-job ${job.error ? 'has-error' : ''}`}>
          <div><strong>{job.busy ? 'Working' : job.error ? 'Action failed' : 'Last action complete'}</strong><span>{job.phase}</span></div>
          <pre aria-label="Sandbox installation log">{job.log.slice(-18).join('\n')}</pre>
        </div>
      ) : null}
    </aside>
  )
}

function OnboardingStage({
  environment,
  job,
  source,
  onAction,
  onPreview,
  onOpenCockpit,
}: {
  environment: Environment
  job: SandboxJob | undefined
  source: 'real' | 'preview'
  onAction: (action: string) => void
  onPreview: (preview: PreviewKey) => void
  onOpenCockpit: () => void
}) {
  const stage = stageFor(environment)
  const busy = Boolean(job?.busy)
  const action = (realAction: string, previewStage: PreviewKey) => source === 'real' ? onAction(realAction) : onPreview(previewStage)
  const copyLabel = stage === 'broken-package' ? 'Copy repair code' : 'Copy installation code'

  const content = {
    'missing-r': {
      title: 'R is not available',
      description: 'Install a compatible R runtime before Spectreasy can be checked or started.',
      action: 'Install isolated R 4.5.3',
      run: () => action('install-r', 'r-only'),
    },
    'outdated-r': {
      title: `R ${environment.rVersion ?? ''} is too old`,
      description: 'This R installation cannot satisfy Spectreasy. Choose or install R 4.5.0 or newer.',
      action: 'Replace with isolated R 4.5.3',
      run: () => action('install-r', 'r-only'),
    },
    'r-only': {
      title: 'Install Spectreasy',
      description: `R ${environment.rVersion} is compatible. Install Spectreasy and its dependencies into this same R library.`,
      action: 'Run installation in sandbox',
      run: () => action('install-spectreasy', 'ready'),
    },
    'broken-package': {
      title: 'Repair the Spectreasy installation',
      description: 'The package was found, but its namespace or a dependency could not be loaded.',
      action: 'Repair sandbox package',
      run: () => action('install-spectreasy', 'ready'),
    },
    ready: {
      title: 'The cockpit can start',
      description: `R ${environment.rVersion} and Spectreasy ${environment.spectreasyVersion} are available in the selected runtime.`,
      action: 'Open cockpit',
      run: onOpenCockpit,
    },
  }[stage]

  return (
    <div className="setup-shell sandbox-setup-shell">
      <header className="setup-header">
        <div className="setup-brand">
          <span className="setup-brand-mark" aria-hidden="true"><i /><i /><i /><i /></span>
          <span><strong>spectreasy</strong></span>
        </div>
        <span className="setup-header-note">Sandboxed onboarding inspection</span>
      </header>

      <main className="setup-main sandbox-main">
        <section className="setup-intro sandbox-intro">
          <div>
            <h1>{stage === 'ready' ? 'Runtime requirements satisfied' : 'Prepare the local runtime'}</h1>
          </div>
          <div className="setup-spectrum-rule" aria-hidden="true">{Array.from({ length: 22 }, (_, index) => <i key={index} />)}</div>
        </section>

        <div className={`sandbox-diagnostic state-${stage}`} role="status" aria-live="polite">
          <span className="setup-status-light" />
          <div><strong>{content.title}</strong><span>{content.description}</span></div>
        </div>

        <div className="sandbox-onboarding-grid">
          <section className="setup-panel setup-panel-primary sandbox-action-panel">
            {stage === 'r-only' || stage === 'broken-package' ? (
              <div className="setup-code setup-code-block sandbox-install-code">
                <pre><code>{installCommand}</code></pre>
                <CopyButton value={installCommand} label={copyLabel} />
              </div>
            ) : stage === 'ready' ? (
              <div className="setup-code setup-code-single"><code>{launchCommand}</code><CopyButton value={launchCommand} label="Copy command" /></div>
            ) : (
              <div className="sandbox-r-source"><strong>Install R 4.5.0 or newer</strong><a href="https://cran.r-project.org/" target="_blank" rel="noreferrer">Open CRAN</a></div>
            )}

            <button className="sandbox-primary-action" type="button" disabled={busy} onClick={content.run}>
              {busy ? <RefreshCw className="is-spinning" size={15} /> : stage === 'ready' ? <TerminalSquare size={15} /> : <MonitorCog size={15} />}
              {busy ? `Running ${job?.phase || 'sandbox action'}…` : content.action}
            </button>

            <dl className="sandbox-runtime-summary">
              <div><dt>R</dt><dd>{environment.rVersion ? environment.rVersion : 'Not installed'}</dd></div>
              <div><dt>Spectreasy</dt><dd>{environment.spectreasyVersion ?? environment.spectreasy}</dd></div>
              {environment.spectreasy === 'broken' || environment.spectreasy === 'ready' ? (
                <div><dt>Dependencies</dt><dd>{environment.spectreasy === 'broken' ? 'Load failure' : 'Ready'}</dd></div>
              ) : null}
            </dl>
          </section>
        </div>
      </main>
    </div>
  )
}

export default function SandboxExperience() {
  const [snapshot, setSnapshot] = useState<SandboxSnapshot | null>(null)
  const [error, setError] = useState('')
  const [source, setSource] = useState<'real' | 'preview'>('real')
  const [preview, setPreview] = useState<PreviewKey>('missing-r')
  const [holdAtReady, setHoldAtReady] = useState(false)

  const refresh = useCallback(async () => {
    try {
      const response = await fetch(`${sandboxApi}/sandbox/status`, { cache: 'no-store' })
      if (!response.ok) throw new Error(`Sandbox controller returned HTTP ${response.status}.`)
      setSnapshot(await response.json() as SandboxSnapshot)
      setError('')
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught))
    }
  }, [])

  useEffect(() => {
    void refresh()
    const interval = window.setInterval(() => void refresh(), 1200)
    return () => window.clearInterval(interval)
  }, [refresh])

  const runAction = useCallback(async (action: string) => {
    try {
      const response = await fetch(`${sandboxApi}/sandbox/action`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ action }),
      })
      const result = await response.json() as { error?: string }
      if (!response.ok) throw new Error(result.error || `Sandbox controller returned HTTP ${response.status}.`)
      await refresh()
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : String(caught))
    }
  }, [refresh])

  const environment = useMemo(
    () => source === 'preview' ? previews[preview].environment : snapshot?.environment ?? previews['missing-r'].environment,
    [preview, snapshot, source],
  )
  const ready = stageFor(environment) === 'ready'
  const showCockpit = ready && !holdAtReady

  const controls = (
    <SandboxControls
      source={source}
      preview={preview}
      snapshot={snapshot}
      holdAtReady={holdAtReady}
      onSource={setSource}
      onPreview={setPreview}
      onAction={(action) => void runAction(action)}
      onHoldAtReady={setHoldAtReady}
    />
  )

  if (error && !snapshot) {
    return (
      <div className="sandbox-controller-error">
        <FlaskConical size={26} />
        <h1>Sandbox controller is offline</h1>
        <p>{error}</p>
        <code>npm run sandbox</code>
        <button type="button" onClick={() => void refresh()}><RefreshCw size={14} /> Retry</button>
      </div>
    )
  }

  if (showCockpit) {
    return (
      <div className="sandbox-cockpit-view">
        <Suspense fallback={<div className="app-loading">Loading sandbox cockpit…</div>}><CockpitApp /></Suspense>
        <div className="sandbox-cockpit-banner"><FlaskConical size={14} /><span>Sandbox cockpit · empty isolated project</span><button type="button" onClick={() => setHoldAtReady(true)}>Inspect ready screen</button></div>
        {controls}
      </div>
    )
  }

  return (
    <div className="sandbox-layout">
      {controls}
      <OnboardingStage
        environment={environment}
        job={source === 'real' ? snapshot?.job : undefined}
        source={source}
        onAction={(action) => void runAction(action)}
        onPreview={setPreview}
        onOpenCockpit={() => setHoldAtReady(false)}
      />
    </div>
  )
}
