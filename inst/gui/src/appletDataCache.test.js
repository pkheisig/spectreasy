import assert from 'node:assert/strict'
import test from 'node:test'
import {
  appletCacheKey,
  clearAppletDataCache,
  loadCachedAppletData,
  setCachedAppletData,
} from './appletDataCache.ts'

test('deduplicates applet loads for the same project revision', async () => {
  clearAppletDataCache()
  const key = appletCacheKey('matrix', '/project', 'revision-1', 'matrix.csv')
  let calls = 0
  const loader = async () => {
    calls += 1
    return { rows: 4 }
  }
  const [first, second] = await Promise.all([
    loadCachedAppletData(key, loader),
    loadCachedAppletData(key, loader),
  ])
  assert.equal(calls, 1)
  assert.equal(first, second)
})

test('replaces stale entries within a project data group', async () => {
  clearAppletDataCache()
  const oldKey = appletCacheKey('gating', '/project', 'revision-1', 10000)
  const newKey = appletCacheKey('gating', '/project', 'revision-2', 10000)
  const group = appletCacheKey('gating', '/project')
  await loadCachedAppletData(oldKey, async () => 'old', group)
  await loadCachedAppletData(newKey, async () => 'new', group)
  let reloaded = false
  await loadCachedAppletData(oldKey, async () => {
    reloaded = true
    return 'old-again'
  }, group)
  assert.equal(reloaded, true)
})

test('allows a saved payload to replace a cached response', async () => {
  clearAppletDataCache()
  const key = appletCacheKey('matrix', '/project', 'matrix.csv')
  setCachedAppletData(key, ['saved'])
  assert.deepEqual(await loadCachedAppletData(key, async () => ['stale']), ['saved'])
})
