/// <reference types="node" />
import test from 'node:test'
import assert from 'node:assert/strict'
import {
  densityBuckets,
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
