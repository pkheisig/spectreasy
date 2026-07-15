import { lazy, Suspense } from 'react'
import { resolveApiToken } from './apiBase'
import { SetupExperience } from './setup/SetupExperience'

const CockpitApp = lazy(() => import('./cockpit/CockpitApp'))
const SandboxExperience = lazy(() => import('./setup/SandboxExperience'))

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

  if (showSetup) return <SetupExperience forceOffline={params.get('emulate') === 'offline'} />

  return (
    <Suspense fallback={<div className="app-loading">Loading Spectreasy cockpit…</div>}>
      <CockpitApp />
    </Suspense>
  )
}
