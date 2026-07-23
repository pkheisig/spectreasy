type CacheEntry = {
  promise: Promise<unknown>
  group: string
}

const MAX_CACHE_ENTRIES = 12
const appletDataCache = new Map<string, CacheEntry>()

export function appletCacheKey(...parts: unknown[]): string {
  return JSON.stringify(parts)
}

export function loadCachedAppletData<T>(key: string, loader: () => Promise<T>, group = key): Promise<T> {
  const cached = appletDataCache.get(key)
  if (cached) {
    appletDataCache.delete(key)
    appletDataCache.set(key, cached)
    return (cached.promise as Promise<T>).catch((error) => {
      // A cached request may have belonged to an effect instance that was
      // discarded by React StrictMode. The next consumer supplies its own
      // loader, so retry an aborted shared entry once instead of inheriting
      // the previous consumer's cancellation.
      if (error?.name !== 'AbortError') throw error
      if (appletDataCache.get(key)?.promise === cached.promise) appletDataCache.delete(key)
      return loadCachedAppletData(key, loader, group)
    })
  }

  for (const [cachedKey, entry] of appletDataCache) {
    if (entry.group === group && cachedKey !== key) appletDataCache.delete(cachedKey)
  }

  const promise = loader().catch((error) => {
    appletDataCache.delete(key)
    throw error
  })
  appletDataCache.set(key, { promise, group })
  while (appletDataCache.size > MAX_CACHE_ENTRIES) {
    const oldest = appletDataCache.keys().next().value
    if (typeof oldest !== 'string') break
    appletDataCache.delete(oldest)
  }
  return promise
}

export function setCachedAppletData<T>(key: string, value: T, group = key): void {
  for (const [cachedKey, entry] of appletDataCache) {
    if (entry.group === group && cachedKey !== key) appletDataCache.delete(cachedKey)
  }
  appletDataCache.set(key, { promise: Promise.resolve(value), group })
}

export function clearAppletDataCache(): void {
  appletDataCache.clear()
}
