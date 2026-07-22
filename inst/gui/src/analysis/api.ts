import { resolveApiBase, resolveApiToken } from '../apiBase'

const API_BASE = resolveApiBase()
const API_TOKEN = resolveApiToken()

type RequestOptions = RequestInit & { signal?: AbortSignal }

function scalar(value: unknown): unknown {
  return Array.isArray(value) && value.length === 1 ? scalar(value[0]) : value
}

export async function analysisRequest<T>(
  path: string,
  projectPath: string,
  options: RequestOptions = {},
): Promise<T> {
  const method = String(options.method ?? 'GET').toUpperCase()
  let requestPath = path
  const next: RequestInit = { ...options, cache: 'no-store' }
  const headers = new Headers(options.headers)
  if (API_TOKEN) headers.set('X-Spectreasy-Token', API_TOKEN)
  if (method === 'GET') {
    requestPath += `${requestPath.includes('?') ? '&' : '?'}project_path=${encodeURIComponent(projectPath)}`
  } else {
    headers.set('Content-Type', 'application/json')
    const supplied = typeof options.body === 'string' ? JSON.parse(options.body) as Record<string, unknown> : {}
    next.body = JSON.stringify({ ...supplied, projectPath })
  }
  next.headers = headers
  const response = await fetch(`${API_BASE}${requestPath}`, next)
  if (!response.ok) throw new Error(`${response.status} ${response.statusText}`)
  const payload = await response.json() as T & { success?: boolean; error?: string }
  if (scalar(payload.success) === false) throw new Error(String(scalar(payload.error) || 'Analysis request failed.'))
  return payload
}
