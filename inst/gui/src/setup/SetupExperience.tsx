import { useEffect, useRef, useState } from 'react'
import {
  Check,
  Clipboard,
  ExternalLink,
  PackageOpen,
  TerminalSquare,
} from 'lucide-react'
import { installAnalysisCommand, installPackageCommand, launchCommand } from './installCommands'
import './setup.css'

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

function StepNumber({ children }: { children: string }) {
  return <span className="setup-step-number" aria-hidden="true">{children}</span>
}

export function SetupExperience() {
  return (
    <div className="setup-shell">
      <header className="setup-header">
        <div className="setup-brand">
          <div className="setup-brand-mark" aria-hidden="true">
            <span />
            <span />
            <span />
          </div>
          <strong>spectreasy</strong>
        </div>
      </header>

      <main className="setup-main">
        <section className="setup-intro">
          <h1>Setup</h1>
        </section>

        <div className="setup-grid">
          <section className="setup-panel setup-panel-primary">
            <div className="setup-card-head">
              <StepNumber>1</StepNumber>
              <div>
                <span className="setup-card-state">R + Spectreasy installed?</span>
                <h2>Start the app</h2>
              </div>
              <TerminalSquare size={20} aria-hidden="true" />
            </div>
            <div className="setup-code setup-code-single">
              <code>{launchCommand}</code>
              <CopyButton value={launchCommand} label="Copy command" />
            </div>
          </section>

          <section className="setup-panel setup-install-panel">
            <div className="setup-card-head">
              <StepNumber>2</StepNumber>
              <div>
                <span className="setup-card-state">R installed?</span>
                <h2>Install Spectreasy</h2>
              </div>
              <PackageOpen size={20} aria-hidden="true" />
            </div>
            <div className="setup-code setup-code-block">
              <pre><code>{installPackageCommand}</code></pre>
              <CopyButton value={installPackageCommand} label="Copy installation code" />
            </div>
          </section>

          <section className="setup-panel setup-analysis-panel">
            <div className="setup-card-head">
              <StepNumber>3</StepNumber>
              <div>
                <span className="setup-card-state">Using population analysis?</span>
                <h2>Install analysis methods</h2>
              </div>
              <PackageOpen size={20} aria-hidden="true" />
            </div>
            <p>Installs the optional R methods and a pinned private Python environment. It never changes system Python; a Python 3.11 or newer <code>python3</code> command is required for the Python-backed methods.</p>
            <div className="setup-code setup-code-single">
              <code>{installAnalysisCommand}</code>
              <CopyButton value={installAnalysisCommand} label="Copy command" />
            </div>
          </section>

          <section className="setup-panel setup-r-panel">
            <div className="setup-card-head">
              <StepNumber>4</StepNumber>
              <div>
                <span className="setup-card-state">No R or R below 4.5?</span>
                <h2>Install R 4.5.0 or newer</h2>
              </div>
              <a className="setup-cran-link" href="https://cran.r-project.org/" target="_blank" rel="noreferrer">
                Open CRAN <ExternalLink size={15} />
              </a>
            </div>
            <div className="setup-r-next">
              <span>Then</span>
              <strong>Install Spectreasy</strong>
              <span>using step 2.</span>
            </div>
          </section>
        </div>
      </main>

      <footer className="setup-footer">
        <span>Spectreasy 1.0</span>
        <a href="https://github.com/pkheisig/spectreasy" target="_blank" rel="noreferrer">
          Source code <ExternalLink size={12} />
        </a>
      </footer>
    </div>
  )
}
