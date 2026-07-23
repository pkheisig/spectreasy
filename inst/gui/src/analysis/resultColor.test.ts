import assert from 'node:assert/strict'
import test from 'node:test'
import { identityColor, markerValue, paletteColors, resultColorModes } from './resultColor.ts'
import type { AnalysisRunResult } from './types.ts'

function fixture(): AnalysisRunResult {
  return {
    metadata: {
      analysis_id: 'analysis',
      method: {} as AnalysisRunResult['metadata']['method'],
      source_file: 'sample.fcs',
      population_id: 'root',
      population_path: 'All events',
      markers: ['CD3', 'CD19'],
      marker_columns: [{ marker: 'CD3', column: 'marker_1' }, { marker: 'CD19', column: 'marker_2' }],
      coordinate_count: 2,
      coordinate_labels: ['UMAP 1', 'UMAP 2'],
      seed: 1,
      event_count: 1,
      root_event_id: null,
      runtime_seconds: 0,
      events_file: 'events.csv',
      identity_annotation: {
        identity_id: 'identity',
        analysis_id: 'analysis',
        method: 'robust-signed-marker-score',
        signatures: [{ name: 'T cell', color: '#123456', positive_markers: ['CD3'], negative_markers: ['CD19'] }],
        thresholds: { min_score: 0.55, min_margin: 0.08 },
        counts: { 'T cell': 1, Unassigned: 0 },
        event_count: 1,
        assigned_count: 1,
        unassigned_count: 0,
        model: { path: 'model.rds', sha256: 'a' },
        scores: { path: 'scores.csv', sha256: 'b' },
        metadata: { path: 'metadata.json', sha256: 'c' },
        created_at: 'now',
      },
    },
    events: [{
      event_id: 1,
      dimension_1: 1,
      dimension_2: 2,
      marker_1: 3.5,
      marker_2: -1,
      cluster_id: 2,
      pseudotime: 0.4,
      predicted_identity: 'T cell',
    }],
    metadata_file: 'metadata.json',
  }
}

test('result color modes expose only data that exists', () => {
  assert.deepEqual(
    resultColorModes(fixture()).map((mode) => mode.id),
    ['density', 'cluster', 'pseudotime', 'identity', 'marker:marker_1', 'marker:marker_2'],
  )
})

test('marker values and identity colors are resolved safely', () => {
  const result = fixture()
  assert.equal(markerValue(result.events[0], 'marker:marker_1'), 3.5)
  assert.equal(identityColor(result, 'T cell'), '#123456')
  assert.equal(identityColor(result, 'Unassigned'), '#a8afb1')
})

test('continuous palettes are interpolated to a high-resolution scale', () => {
  const colors = paletteColors('sunset', 64)
  assert.equal(colors.length, 64)
  assert.match(colors[0], /^rgba\(/)
  assert.notEqual(colors[0], colors[63])
})
