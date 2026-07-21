import assert from 'node:assert/strict'
import test from 'node:test'
import { normalizeAiQcResponse } from './aiQcResponse.ts'

test('AI-QC response normalization is defensive and preserves conservative defaults', () => {
  const normalized = normalizeAiQcResponse({
    status: ['stale'],
    stale: [true],
    default_scope: ['combined'],
    available_scopes: ['control', 'combined', 'invalid'],
    privacy_mode: ['strict'],
    grade_counts: { good: [4], review: [2], poor: [1], not_graded: [3] },
    findings: [{ metric_id: 'QC-X', entity: 'sample_001', stage: 'samples', grade: 'review', explanation: 'Review evidence.' }],
    prompt: ['prompt'], prompt_characters: [6], estimated_tokens: [2],
    schema: { name: ['spectreasy-ai-qc'], version: ['1.0.0'] },
    profile: { name: ['reviewed'], version: ['2'], reference_n: [8] },
    source_hashes: [{ path: 'spectreasy_outputs/matrix.csv', sha256: 'abc' }],
  })
  assert.equal(normalized.status, 'stale')
  assert.equal(normalized.stale, true)
  assert.deepEqual(normalized.availableScopes, ['control', 'combined'])
  assert.equal(normalized.privacyMode, 'strict')
  assert.deepEqual(normalized.gradeCounts, { good: 4, review: 2, poor: 1, not_graded: 3 })
  assert.equal(normalized.findings[0]?.metricId, 'QC-X')
  assert.equal(normalized.schemaVersion, '1.0.0')
  assert.equal(normalized.referenceN, 8)
  assert.equal(normalized.sourceHashes[0]?.sha256, 'abc')
})

test('AI-QC invalid status and grades become failed and not graded', () => {
  const normalized = normalizeAiQcResponse({ status: 'mystery', findings: [{ grade: 'excellent' }] })
  assert.equal(normalized.status, 'failed')
  assert.equal(normalized.findings[0]?.grade, 'not_graded')
  assert.equal(normalized.prompt, '')
})
