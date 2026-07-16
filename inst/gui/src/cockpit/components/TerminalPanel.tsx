import { Activity, CircleAlert, CircleCheck, Clock3, ScrollText, WifiOff } from 'lucide-react'
import { useEffect, useRef } from 'react'
import type { ExecutionLogEntry } from '../types'

type Props = {
  connected: boolean
  entries: ExecutionLogEntry[]
}

function EntryIcon({ kind }: { kind: ExecutionLogEntry['kind'] }) {
  if (kind === 'success') return <CircleCheck size={14} />
  if (kind === 'warning' || kind === 'error') return <CircleAlert size={14} />
  if (kind === 'command') return <Activity size={14} />
  return <Clock3 size={14} />
}

export function TerminalPanel({ connected, entries }: Props) {
  const outputRef = useRef<HTMLDivElement>(null)

  useEffect(() => {
    outputRef.current?.scrollTo({ top: outputRef.current.scrollHeight })
  }, [entries])

  return (
    <section className="cockpit-logs-popout" aria-label="Execution logs">
      <header className="logs-header">
        <div>
          <ScrollText size={16} />
          <strong>Execution logs</strong>
        </div>
        <span className={`logs-connection ${connected ? 'is-connected' : ''}`}>
          {connected ? <Activity size={13} /> : <WifiOff size={13} />}
          {connected ? 'R connected' : 'R offline'}
        </span>
      </header>
      <div className="logs-output" ref={outputRef} role="log" aria-live="polite">
        {entries.length === 0 ? (
          <div className="logs-empty">
            <ScrollText size={22} />
            <strong>No workflow activity yet</strong>
            <span>Cockpit operations and R output will appear here.</span>
          </div>
        ) : entries.map((entry) => (
          <article className={`log-entry log-${entry.kind}`} key={entry.id}>
            <span className="log-entry-icon"><EntryIcon kind={entry.kind} /></span>
            <div>
              <time>{entry.time}</time>
              {entry.kind === 'command' ? <pre>{entry.text}</pre> : <p>{entry.text}</p>}
            </div>
          </article>
        ))}
      </div>
    </section>
  )
}
