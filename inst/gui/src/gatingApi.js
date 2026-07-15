export function withGatingApiToken(options = {}, token = '') {
  const headers = new Headers(options.headers)
  if (token) headers.set('X-Spectreasy-Token', token)
  return { ...options, headers }
}
