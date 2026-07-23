import assert from 'node:assert/strict'
import test from 'node:test'
import { coordinateLabel, eventValue, resultCoordinates } from './resultPlot.ts'
import type { AnalysisRunResult } from './types.ts'

function result(coordinateCount: number, labels: string[], dimensions: Record<string, number>): AnalysisRunResult {
  return {
    metadata: {
      analysis_id: 'test',
      method: {} as AnalysisRunResult['metadata']['method'],
      source_file: 'sample.fcs',
      population_id: 'root',
      population_path: 'All events',
      markers: ['A', 'B', 'C'],
      coordinate_count: coordinateCount,
      coordinate_labels: labels,
      seed: 1,
      event_count: 1,
      root_event_id: null,
      runtime_seconds: 0,
      events_file: 'events.csv',
    },
    events: [{ event_id: 1, dimension_1: 1, dimension_2: 2, ...dimensions }],
    metadata_file: 'metadata.json',
  }
}

test('3D eligibility follows returned coordinates rather than the method name', () => {
  assert.deepEqual(resultCoordinates(result(2, ['PC 1', 'PC 2'], {})), ['dimension_1', 'dimension_2'])
  assert.deepEqual(
    resultCoordinates(result(3, ['PC 1', 'PC 2', 'PC 3'], { dimension_3: 3 })),
    ['dimension_1', 'dimension_2', 'dimension_3'],
  )
})

test('coordinate labels and non-finite values have safe fallbacks', () => {
  const value = result(3, ['DC 1', 'DC 2', 'DC 3'], { dimension_3: Number.NaN })
  assert.equal(coordinateLabel(value, 'dimension_3'), 'DC 3')
  assert.equal(eventValue(value.events[0], 'dimension_3'), 0)
})
