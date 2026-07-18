function locationToken(): string {
  if (typeof window === 'undefined') return ''
  const queryToken = new URLSearchParams(window.location.search).get('token')?.trim()
  if (queryToken) return queryToken
  return new URLSearchParams(window.location.hash.replace(/^#/, '')).get('token')?.trim() ?? ''
}

const sessionApiToken = locationToken()

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
  return sessionApiToken
}

export function removeApiTokenFromLocation(): void {
  if (typeof window === 'undefined' || !sessionApiToken) return
  const query = new URLSearchParams(window.location.search)
  query.delete('token')
  const fragment = new URLSearchParams(window.location.hash.replace(/^#/, ''))
  fragment.delete('token')
  const nextQuery = query.toString()
  const nextFragment = fragment.toString()
  const cleanUrl = `${window.location.pathname}${nextQuery ? `?${nextQuery}` : ''}${nextFragment ? `#${nextFragment}` : ''}`
  window.history.replaceState(window.history.state, '', cleanUrl)
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
