import { lazy, Suspense, useEffect, useState } from 'react'
import { probeLocalBackend, resolveApiToken } from './apiBase'
import { SetupExperience } from './setup/SetupExperience'

const CockpitApp = lazy(() => import('./cockpit/CockpitApp'))
const SandboxExperience = lazy(() => import('./setup/SandboxExperience'))

type CockpitConnection = 'checking' | 'connected' | 'offline'

function ConnectedCockpit() {
  const [connection, setConnection] = useState<CockpitConnection>('checking')

  useEffect(() => {
    let cancelled = false
    const connect = async () => {
      for (let attempt = 0; attempt < 2; attempt += 1) {
        if (await probeLocalBackend()) {
          if (!cancelled) setConnection('connected')
          return
        }
        if (attempt === 0) {
          await new Promise((resolve) => window.setTimeout(resolve, 300))
        }
      }
      if (!cancelled) setConnection('offline')
    }

    const timer = window.setTimeout(() => void connect(), 0)
    return () => {
      cancelled = true
      window.clearTimeout(timer)
    }
  }, [])

  useEffect(() => {
    if (connection !== 'offline') return
    let cancelled = false
    const retry = async () => {
      if (await probeLocalBackend()) {
        if (!cancelled) setConnection('connected')
      }
    }
    const interval = window.setInterval(() => void retry(), 1500)
    return () => {
      cancelled = true
      window.clearInterval(interval)
    }
  }, [connection])

  if (connection === 'checking') {
    return <div className="app-loading">Connecting to the local R session…</div>
  }

  if (connection === 'offline') {
    return <SetupExperience />
  }

  return (
    <Suspense fallback={<div className="app-loading">Loading Spectreasy cockpit…</div>}>
      <CockpitApp />
    </Suspense>
  )
}

export default function App() {
  const params = new URLSearchParams(window.location.search)
  if (params.get('sandbox') === '1') {
    return (
      <Suspense fallback={<div className="app-loading">Loading onboarding sandbox…</div>}>
        <SandboxExperience />
      </Suspense>
    )
  }
  const showSetup = params.get('setup') === '1' || !resolveApiToken()

  if (showSetup) return <SetupExperience />

  return <ConnectedCockpit />
}
