import assert from 'node:assert/strict'
import test from 'node:test'
import { resultTableCsv } from './resultExport.ts'
import type { AnalysisRunResult } from './types.ts'

const result = {
  metadata: {
    analysis_id: 'example',
    method: { name: 'UMAP' },
    coordinate_count: 2,
    coordinate_labels: ['UMAP 1', 'UMAP 2'],
    marker_columns: [
      { marker: 'CD3', channel: 'CD3-A', column: 'marker_1' },
      { marker: 'CD19', channel: 'CD19-A', column: 'marker_2' },
    ],
  },
  events: [{
    event_id: 7,
    dimension_1: 1.5,
    dimension_2: -2,
    marker_1: 0.4,
    marker_2: 1.8,
    cluster_id: 3,
    predicted_identity: 'T cell',
  }],
} as unknown as AnalysisRunResult

test('plot table export contains display marker names, coordinates, clusters, and identities', () => {
  const csv = resultTableCsv(result)
  const [header, row] = csv.split('\r\n')
  assert.equal(header, 'event_id,UMAP 1,UMAP 2,CD3,CD19,cluster_id,predicted_identity')
  assert.equal(row, '7,1.5,-2,0.4,1.8,3,T cell')
  assert.equal(header.includes('CD3-A'), false)
})

test('pooled plot table export preserves source file and original FCS row', () => {
  const pooled = structuredClone(result)
  pooled.events[0].source_file = 'samples/a.fcs'
  pooled.events[0].source_event_id = 42
  pooled.events[0].sample_id = 2
  const csv = resultTableCsv(pooled)
  assert.match(csv.split('\r\n')[0], /source_file,source_event_id,sample_id/)
  assert.match(csv, /samples\/a\.fcs,42,2/)
})
