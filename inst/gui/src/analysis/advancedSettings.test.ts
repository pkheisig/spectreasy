import assert from 'node:assert/strict'
import test from 'node:test'
import { advancedSettingsIssues, resolvedAdvancedSettings } from './advancedSettings.ts'
import type { AnalysisMethod } from './types.ts'

const method = {
  id: 'example',
  name: 'Example',
  family: 'reduction',
  parameters: [
    { id: 'neighbors', label: 'Neighbors', type: 'integer', default: 15 },
    { id: 'metric', label: 'Metric', type: 'select', default: 'euclidean' },
    { id: 'normalize', label: 'Normalize', type: 'boolean', default: true },
  ],
} as AnalysisMethod

test('advanced settings preserve defaults and apply only explicit overrides', () => {
  assert.deepEqual(resolvedAdvancedSettings([method], {}), {
    example: { neighbors: 15, metric: 'euclidean', normalize: true },
  })
  assert.deepEqual(resolvedAdvancedSettings([method], {
    example: { neighbors: 22 },
  }), {
    example: { neighbors: 22, metric: 'euclidean', normalize: true },
  })
})

test('cross-parameter settings fail before execution', () => {
  const wishbone = {
    ...method,
    id: 'wishbone',
    name: 'Wishbone',
    family: 'trajectory',
    parameters: [
      { id: 'neighbors', label: 'Retained neighbors', type: 'integer', default: 15 },
      { id: 'candidate_neighbors', label: 'Candidate neighbors', type: 'integer', default: 20 },
    ],
  } as AnalysisMethod
  assert.deepEqual(
    advancedSettingsIssues([wishbone], { wishbone: { candidate_neighbors: 12 } }),
    ['Wishbone: candidate neighbors must exceed retained neighbors.'],
  )
})
