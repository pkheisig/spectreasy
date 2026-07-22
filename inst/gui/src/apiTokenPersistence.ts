export const API_TOKEN_STORAGE_KEY = 'spectreasy-api-token'

export function selectApiToken(locationToken: string, environmentToken: string, storedToken: string): string {
  return [locationToken, environmentToken, storedToken]
    .map((token) => token.trim())
    .find(Boolean) ?? ''
}
