export function resolveApiBase(): string {
  if (typeof window !== 'undefined') {
    const queryBase = new URLSearchParams(window.location.search).get('api')?.trim()
    if (queryBase) return queryBase.replace(/\/$/, '')
  }
  const envBase = (import.meta.env.VITE_API_BASE as string | undefined)?.trim()
  if (envBase) return envBase.replace(/\/$/, '')
  if (typeof window !== 'undefined' && window.location.hostname.endsWith('github.io')) {
    return 'http://127.0.0.1:8000'
  }
  if (typeof window !== 'undefined' && window.location.port === '5174') return 'http://127.0.0.1:8000'
  if (typeof window !== 'undefined') return window.location.origin.replace(/\/$/, '')
  return 'http://127.0.0.1:8000'
}

export function resolveApiToken(): string {
  if (typeof window === 'undefined') return ''
  return new URLSearchParams(window.location.search).get('token')?.trim() ?? ''
}

export async function probeLocalBackend(timeoutMs = 1800): Promise<boolean> {
  const controller = new AbortController()
  const timeout = globalThis.setTimeout(() => controller.abort(), timeoutMs)
  const token = resolveApiToken()

  try {
    const response = await fetch(`${resolveApiBase()}/status`, {
      cache: 'no-store',
      headers: token ? { 'X-Spectreasy-Token': token } : undefined,
      signal: controller.signal,
    })
    if (!response.ok) return false
    const payload = await response.json() as { status?: unknown }
    const status = Array.isArray(payload.status) ? payload.status[0] : payload.status
    return status === 'ok'
  } catch {
    return false
  } finally {
    globalThis.clearTimeout(timeout)
  }
}
