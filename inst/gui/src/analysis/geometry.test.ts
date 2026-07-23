/// <reference types="node" />
import test from 'node:test'
import assert from 'node:assert/strict'
import {
  densityContourSegments,
  densityBuckets,
  histogramCounts,
  inverseTransformValue,
  linearScale,
  nearestEvent,
  robustExtent,
  transformValue,
} from './geometry.ts'

test('asinh plot transforms are numerically reversible', () => {
  for (const value of [-10000, -150, 0, 150, 10000]) {
    assert.ok(Math.abs(inverseTransformValue(transformValue(value, 'asinh'), 'asinh') - value) < 1e-9)
  }
})

test('biexponential plot transforms match control gating and are reversible', () => {
  for (const value of [-100000, -500, -1, 0, 1, 500, 100000]) {
    const transformed = transformValue(value, 'biexponential')
    assert.ok(Number.isFinite(transformed))
    assert.ok(Math.abs(inverseTransformValue(transformed, 'biexponential') - value) < Math.max(1e-8, Math.abs(value) * 1e-10))
  }
  assert.equal(transformValue(-500, 'biexponential'), -transformValue(500, 'biexponential'))
})

test('histogram and contour geometry preserve compact multimodal structure', () => {
  const values = [0, 0.1, 0.2, 5, 5.1, 5.2]
  const counts = histogramCounts(values, [0, 6], 6)
  assert.equal(counts.reduce((sum, value) => sum + value, 0), values.length)
  assert.ok(counts[0] > 0)
  assert.ok(counts[5] > 0)

  const events = Array.from({ length: 80 }, (_, index) => ({
    x: index < 40 ? 1 + (index % 5) * 0.03 : 4 + (index % 5) * 0.03,
    y: index < 40 ? 1 + (index % 7) * 0.03 : 4 + (index % 7) * 0.03,
  }))
  const contours = densityContourSegments(events, [0, 5], [0, 5], 32)
  assert.ok(contours.length > 20)
  assert.ok(new Set(contours.map((segment) => segment.level)).size >= 3)
})

test('robust extent ignores non-finite values and pads constant data', () => {
  assert.deepEqual(robustExtent([Number.NaN, Number.POSITIVE_INFINITY]), [0, 1])
  assert.deepEqual(robustExtent([5, 5, 5]), [4.5, 5.5])
  const extent = robustExtent([0, 1, 2, 3, 100000])
  assert.ok(extent[0] < 0)
  assert.ok(extent[1] > 3)
})

test('density buckets and root-event hit testing preserve event identity', () => {
  const events = [
    { event_id: 11, x: 1, y: 1 },
    { event_id: 12, x: 1.1, y: 1.1 },
    { event_id: 13, x: 9, y: 9 },
  ]
  const density = densityBuckets(events, [0, 10], [0, 10], 5)
  assert.equal(density[0], density[1])
  assert.ok(density[0] > density[2])

  const scale = linearScale([0, 10], [0, 100])
  const found = nearestEvent(events, { x: 90, y: 90 }, scale, scale, 'linear', 'linear', 5)
  assert.equal(found?.event_id, 13)
})
