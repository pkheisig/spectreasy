import assert from 'node:assert/strict'
import test from 'node:test'
import {
  STALE_CHUNK_RELOAD_COOLDOWN_MS,
  isStaleChunkLoadError,
  shouldReloadStaleChunk,
} from './staleChunkRecovery.ts'

test('recognizes stale Vite lazy-module failures', () => {
  assert.equal(isStaleChunkLoadError(new TypeError('Failed to fetch dynamically imported module: http://127.0.0.1/assets/GatingGui-old.js')), true)
  assert.equal(isStaleChunkLoadError(new Error('Importing a module script failed.')), true)
  assert.equal(isStaleChunkLoadError(new Error('The gating API returned 500.')), false)
})

test('allows one stale-chunk reload per cooldown window', () => {
  const now = 50000
  assert.equal(shouldReloadStaleChunk(0, now), true)
  assert.equal(shouldReloadStaleChunk(now - 1000, now), false)
  assert.equal(shouldReloadStaleChunk(now - STALE_CHUNK_RELOAD_COOLDOWN_MS, now), true)
})
