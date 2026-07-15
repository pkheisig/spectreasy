import { useEffect, useRef, useState } from 'react'
import type { FormEvent, PointerEvent as ReactPointerEvent } from 'react'
import { Check, Copy, Minus, Power, TerminalSquare, X } from 'lucide-react'
import { createPortal } from 'react-dom'
import { runTerminalCommand } from '../api'

type Props = {
  connected: boolean
  projectPath: string
  widthPct: number
  heightPct: number
  onClose: () => void
  onRefresh: () => void | Promise<void>
  onTerminate: () => Promise<boolean>
  onSizeChange: (size: { terminalWidthPct: number; terminalHeightPct: number }) => void
}

type Entry = { command?: string; output: string; ok?: boolean }

const launchCommand = `R -q -e 'spectreasy::spectreasy_gui()'`

function clampWidth(width: number) {
  const maximum = Math.max(320, window.innerWidth - 40)
  return Math.max(Math.min(420, maximum), Math.min(maximum, width))
}

function clampHeight(height: number) {
  const maximum = Math.max(150, window.innerHeight - 90)
  return Math.max(Math.min(180, maximum), Math.min(maximum, height))
}

export function TerminalPanel({ connected, projectPath, widthPct, heightPct, onClose, onRefresh, onTerminate, onSizeChange }: Props) {
  const [width, setWidth] = useState(() => clampWidth(Math.round(window.innerWidth * widthPct / 100)))
  const [height, setHeight] = useState(() => clampHeight(Math.round(window.innerHeight * heightPct / 100)))
  const [minimized, setMinimized] = useState(false)
  const [command, setCommand] = useState('')
  const [cwd, setCwd] = useState(projectPath)
  const [busy, setBusy] = useState(false)
  const [copied, setCopied] = useState(false)
  const [confirmClose, setConfirmClose] = useState(false)
  const [entries, setEntries] = useState<Entry[]>(() => connected
    ? [{ output: 'Spectreasy R console connected.' }]
    : [{ output: 'The GitHub Pages cockpit cannot launch local processes. Run the command below in Terminal, then refresh the cockpit.' }])
  const inputRef = useRef<HTMLInputElement>(null)
  const outputRef = useRef<HTMLDivElement>(null)
  const panelRef = useRef<HTMLElement>(null)
  const resizeCleanupRef = useRef<(() => void) | null>(null)

  useEffect(() => {
    if (!minimized) inputRef.current?.focus()
  }, [minimized])

  useEffect(() => {
    outputRef.current?.scrollTo({ top: outputRef.current.scrollHeight })
  }, [entries])

  useEffect(() => () => resizeCleanupRef.current?.(), [])

  function beginResize(event: ReactPointerEvent<HTMLDivElement>, axis: 'width-left' | 'width-right' | 'height') {
    if (minimized) return
    event.preventDefault()
    const panel = panelRef.current
    if (!panel) return
    resizeCleanupRef.current?.()
    const startX = event.clientX
    const startY = event.clientY
    const startWidth = width
    const startHeight = height
    let nextWidth = startWidth
    let nextHeight = startHeight
    let frame = 0
    panel.style.transition = 'none'

    const render = () => {
      panel.style.width = `${nextWidth}px`
      panel.style.height = `${nextHeight}px`
      frame = 0
    }
    const move = (moveEvent: PointerEvent) => {
      if (axis === 'width-right') nextWidth = clampWidth(startWidth + (moveEvent.clientX - startX) * 2)
      else if (axis === 'width-left') nextWidth = clampWidth(startWidth - (moveEvent.clientX - startX) * 2)
      else nextHeight = clampHeight(startHeight - moveEvent.clientY + startY)
      if (!frame) frame = window.requestAnimationFrame(render)
    }
    const cleanup = () => {
      if (frame) window.cancelAnimationFrame(frame)
      panel.style.transition = ''
      document.removeEventListener('pointermove', move)
      document.removeEventListener('pointerup', stop)
      resizeCleanupRef.current = null
    }
    const stop = () => {
      cleanup()
      render()
      setWidth(Math.round(nextWidth))
      setHeight(Math.round(nextHeight))
      onSizeChange({
        terminalWidthPct: Math.round((nextWidth / window.innerWidth) * 1000) / 10,
        terminalHeightPct: Math.round((nextHeight / window.innerHeight) * 1000) / 10,
      })
    }
    resizeCleanupRef.current = cleanup
    document.addEventListener('pointermove', move)
    document.addEventListener('pointerup', stop)
  }

  async function submit(event: FormEvent) {
    event.preventDefault()
    const next = command.trim()
    if (!next || busy) return
    setCommand('')
    setBusy(true)
    setEntries((current) => [...current, { command: next, output: '' }])
    const result = await runTerminalCommand(next, cwd)
    setCwd(result.cwd || cwd)
    setEntries((current) => [...current, { output: result.output || (result.success ? 'Done.' : 'Command failed.'), ok: result.success }])
    setBusy(false)
    if (result.refresh) await onRefresh()
    if (result.shutdownRequested) await terminateSession()
  }

  async function copyLaunchCommand() {
    await navigator.clipboard.writeText(launchCommand)
    setCopied(true)
    window.setTimeout(() => setCopied(false), 1400)
  }

  async function terminateSession() {
    setConfirmClose(false)
    setBusy(true)
    const terminated = await onTerminate()
    if (terminated) {
      onClose()
      return
    }
    setBusy(false)
    setEntries((current) => [...current, { output: 'The R session could not be terminated.', ok: false }])
  }

  function requestClose() {
    if (connected) setConfirmClose(true)
    else onClose()
  }

  return createPortal(
    <section
      ref={panelRef}
      className={`cockpit-terminal ${minimized ? 'is-minimized' : ''}`}
      style={{ width, height: minimized ? 42 : height }}
      aria-label="Terminal"
    >
      <div className="terminal-resize-height" onPointerDown={(event) => beginResize(event, 'height')} />
      <div className="terminal-resize-width terminal-resize-left" onPointerDown={(event) => beginResize(event, 'width-left')} />
      <div className="terminal-resize-width terminal-resize-right" onPointerDown={(event) => beginResize(event, 'width-right')} />
      <header className="terminal-header">
        <TerminalSquare size={15} />
        <strong>Terminal</strong>
        <span className={`terminal-connection ${connected ? 'is-connected' : ''}`}>{connected ? 'R console' : 'offline'}</span>
        <span className="terminal-cwd" title={cwd}>{cwd || '~'}</span>
        <button type="button" onClick={() => setMinimized((value) => !value)} aria-label={minimized ? 'Restore terminal' : 'Minimize terminal'}><Minus size={15} /></button>
        <button type="button" onClick={requestClose} aria-label="Close terminal"><X size={15} /></button>
      </header>
      {!minimized && <>
        <div className="terminal-output" ref={outputRef} role="log" aria-live="polite">
          {entries.map((entry, index) => <div key={index} className={entry.ok === false ? 'is-error' : ''}>
            {entry.command && <div className="terminal-command"><span>&gt;</span> {entry.command}</div>}
            {entry.output && <pre>{entry.output}</pre>}
          </div>)}
          {!connected && <button type="button" className="terminal-launch-command" onClick={() => void copyLaunchCommand()}>
            <code>{launchCommand}</code>{copied ? <Check size={14} /> : <Copy size={14} />}
          </button>}
        </div>
        <form className="terminal-input" onSubmit={(event) => void submit(event)}>
          <span>&gt;</span>
          <input ref={inputRef} value={command} onChange={(event) => setCommand(event.target.value)} disabled={busy || !connected} aria-label="R console command" autoComplete="off" spellCheck={false} placeholder={connected ? 'Enter R code' : 'Start the local R backend to enable commands'} />
        </form>
      </>}
      {confirmClose && createPortal(
        <div className="cockpit-confirm-overlay" role="presentation" onMouseDown={() => setConfirmClose(false)}>
          <div className="cockpit-confirm terminal-close-confirm" role="dialog" aria-modal="true" aria-labelledby="terminal-close-title" onMouseDown={(event) => event.stopPropagation()}>
            <Power size={21} />
            <h2 id="terminal-close-title">Close the terminal?</h2>
            <p>You can keep the R session running in the background, or terminate it completely. Terminating stops the cockpit backend.</p>
            <div>
              <button className="button button-ghost" type="button" onClick={() => setConfirmClose(false)}>Cancel</button>
              <button className="button" type="button" onClick={onClose}>Keep R running</button>
              <button className="button button-danger" type="button" onClick={() => void terminateSession()}>Terminate R session</button>
            </div>
          </div>
        </div>,
        document.body,
      )}
    </section>,
    document.body,
  )
}
