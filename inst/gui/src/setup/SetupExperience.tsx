import { useCallback, useEffect, useRef, useState } from 'react'
import {
  Check,
  Clipboard,
  Download,
  ExternalLink,
  MonitorCog,
  RefreshCw,
  TerminalSquare,
} from 'lucide-react'
import { resolveApiBase } from '../apiBase'
import { installAndLaunchCommand, launchCommand } from './installCommands'
import './setup.css'

const installCommand = installAndLaunchCommand

type BackendState = 'checking' | 'offline' | 'detected'
const apiBase = resolveApiBase()

function CopyButton({ value, label }: { value: string; label: string }) {
  const [copied, setCopied] = useState(false)
  const resetTimer = useRef<number | null>(null)

  useEffect(() => () => {
    if (resetTimer.current !== null) window.clearTimeout(resetTimer.current)
  }, [])

  const copy = async () => {
    try {
      await navigator.clipboard.writeText(value)
    } catch {
      const textarea = document.createElement('textarea')
      textarea.value = value
      textarea.style.position = 'fixed'
      textarea.style.opacity = '0'
      document.body.appendChild(textarea)
      textarea.select()
      document.execCommand('copy')
      textarea.remove()
    }
    setCopied(true)
    if (resetTimer.current !== null) window.clearTimeout(resetTimer.current)
    resetTimer.current = window.setTimeout(() => setCopied(false), 1800)
  }

  return (
    <button type="button" className="setup-copy-button" onClick={() => void copy()}>
      {copied ? <Check size={15} /> : <Clipboard size={15} />}
      {copied ? 'Copied' : label}
    </button>
  )
}

function ConnectionState({ state, onCheck }: { state: BackendState; onCheck: () => void }) {
  const detected = state === 'detected'
  return (
    <div className={`setup-connection-state state-${state}`} role="status" aria-live="polite">
      <span className="setup-status-light" />
      <div>
        <strong>{detected ? 'Local backend detected' : state === 'checking' ? 'Checking localhost' : 'No local session detected'}</strong>
        <span>{detected ? 'Use the authenticated cockpit tab opened by R.' : 'Your data remains on this computer.'}</span>
      </div>
      <button type="button" onClick={onCheck} disabled={state === 'checking'}>
        <RefreshCw className={state === 'checking' ? 'is-spinning' : ''} size={14} />
        Check again
      </button>
    </div>
  )
}

export function SetupExperience({ forceOffline = false }: { forceOffline?: boolean }) {
  const [backendState, setBackendState] = useState<BackendState>(forceOffline ? 'offline' : 'checking')
  const launcherHref = `${import.meta.env.BASE_URL}spectreasy-launcher.R`

  const checkBackend = useCallback(async () => {
    if (forceOffline) {
      setBackendState('offline')
      return
    }
    setBackendState('checking')
    const controller = new AbortController()
    const timeout = window.setTimeout(() => controller.abort(), 1400)
    try {
      const response = await fetch(`${apiBase}/status`, { cache: 'no-store', signal: controller.signal })
      setBackendState(response.ok ? 'detected' : 'offline')
    } catch {
      setBackendState('offline')
    } finally {
      window.clearTimeout(timeout)
    }
  }, [forceOffline])

  useEffect(() => {
    if (forceOffline) return
    void checkBackend()
    const interval = window.setInterval(() => void checkBackend(), 4000)
    return () => window.clearInterval(interval)
  }, [checkBackend, forceOffline])

  return (
    <div className="setup-shell">
      <header className="setup-header">
        <div className="setup-brand">
          <span className="setup-brand-mark" aria-hidden="true"><i /><i /><i /><i /></span>
          <span><strong>spectreasy</strong><small>Local spectral analysis</small></span>
        </div>
        <span className="setup-header-note">Hosted interface · local computation</span>
      </header>

      <main className="setup-main">
        <section className="setup-intro">
          <div>
            <h1>Connect a local R session</h1>
            <p>
              The cockpit is ready, but scientific processing runs through R on this computer.
              Start Spectreasy from the R installation you intend to use.
            </p>
          </div>
          <div className="setup-spectrum-rule" aria-hidden="true">
            {Array.from({ length: 22 }, (_, index) => <i key={index} />)}
          </div>
        </section>

        <ConnectionState state={backendState} onCheck={() => void checkBackend()} />

        <div className="setup-grid">
          <section className="setup-panel setup-panel-primary">
            <div className="setup-panel-heading">
              <TerminalSquare size={19} />
              <div><h2>Already have R and Spectreasy?</h2><p>Run this in the R version you want the cockpit to use.</p></div>
            </div>
            <div className="setup-code setup-code-single">
              <code>{launchCommand}</code>
              <CopyButton value={launchCommand} label="Copy command" />
            </div>
            <p className="setup-panel-note">
              R will start the authenticated localhost backend and open a connected cockpit in a new tab.
            </p>
          </section>

          <aside className="setup-panel setup-requirements">
            <ol>
              <li>
                <span>1</span>
                <div><strong>Install R 4.5.0 or newer</strong><p>Use the official installer for your operating system.</p></div>
                <a href="https://cran.r-project.org/" target="_blank" rel="noreferrer">Open CRAN <ExternalLink size={13} /></a>
              </li>
              <li>
                <span>2</span>
                <div><strong>Install Spectreasy</strong><p>Copy the code below into that same R installation.</p></div>
              </li>
              <li>
                <span>3</span>
                <div><strong>Start the cockpit</strong><p>The installation code finishes by launching it.</p></div>
              </li>
            </ol>
          </aside>

          <section className="setup-panel setup-install-panel">
            <div className="setup-panel-heading">
              <MonitorCog size={19} />
              <div><h2>Install or repair Spectreasy</h2><p>No command is run by this website.</p></div>
            </div>
            <div className="setup-code setup-code-block">
              <pre><code>{installCommand}</code></pre>
              <CopyButton value={installCommand} label="Copy installation code" />
            </div>
          </section>

          <aside className="setup-panel setup-helper-panel">
            <Download size={22} />
            <h2>Portable setup helper</h2>
            <p>
              A readable R script that reports the active R version and library, checks Spectreasy and its dependencies,
              and launches only after confirmation. It installs nothing automatically.
            </p>
            <a className="setup-download-button" href={launcherHref} download="spectreasy-launcher.R">
              <Download size={15} /> Download .R helper
            </a>
            <small>One script · macOS, Windows and Linux · no administrator access</small>
          </aside>
        </div>
      </main>

      <footer className="setup-footer">
        <span>Spectreasy 1.0</span>
        <span>No project files or cytometry data are uploaded by this interface.</span>
        <a href="https://github.com/pkheisig/spectreasy" target="_blank" rel="noreferrer">Source code <ExternalLink size={12} /></a>
      </footer>
    </div>
  )
}
