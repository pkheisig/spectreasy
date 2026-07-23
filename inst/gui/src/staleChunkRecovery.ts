const STALE_CHUNK_RELOAD_KEY = 'spectreasy-stale-chunk-reload'
export const STALE_CHUNK_RELOAD_COOLDOWN_MS = 15_000

export function isStaleChunkLoadError(error: unknown): boolean {
  const message = error instanceof Error ? error.message : String(error ?? '')
  return /failed to fetch dynamically imported module|importing a module script failed|loading (?:css )?chunk .+ failed/i.test(message)
}

export function shouldReloadStaleChunk(lastAttempt: number, now: number): boolean {
  return !Number.isFinite(lastAttempt) || lastAttempt <= 0 || now - lastAttempt >= STALE_CHUNK_RELOAD_COOLDOWN_MS
}

export function installStaleChunkRecovery() {
  const recover = (event: Event) => {
    const now = Date.now()
    let lastAttempt = 0
    try {
      lastAttempt = Number(window.sessionStorage.getItem(STALE_CHUNK_RELOAD_KEY) ?? 0)
    } catch {
      // A restricted session store must not prevent recovery.
    }
    if (!shouldReloadStaleChunk(lastAttempt, now)) return

    event.preventDefault()
    try {
      window.sessionStorage.setItem(STALE_CHUNK_RELOAD_KEY, String(now))
    } catch {
      // Reloading still repairs the stale page even without the loop guard.
    }
    window.location.reload()
  }

  window.addEventListener('vite:preloadError', recover)
  const resetTimer = window.setTimeout(() => {
    try {
      const lastAttempt = Number(window.sessionStorage.getItem(STALE_CHUNK_RELOAD_KEY) ?? 0)
      if (lastAttempt > 0 && Date.now() - lastAttempt >= STALE_CHUNK_RELOAD_COOLDOWN_MS) {
        window.sessionStorage.removeItem(STALE_CHUNK_RELOAD_KEY)
      }
    } catch {
      // Nothing to reset when session storage is unavailable.
    }
  }, STALE_CHUNK_RELOAD_COOLDOWN_MS)

  return () => {
    window.removeEventListener('vite:preloadError', recover)
    window.clearTimeout(resetTimer)
  }
}
