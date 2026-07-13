import assert from 'node:assert/strict'
import test from 'node:test'

import {
  fileUsesHistogramGates,
  fileUsesNegativeHistogramGate,
  pruneUnavailableNegativeGates,
} from './gatingEligibility.js'

test('only controls with histogram plots participate in histogram autogating', () => {
  assert.equal(fileUsesHistogramGates({ is_af: true, uses_histogram_gates: false }), false)
  assert.equal(fileUsesHistogramGates({ is_af: false, uses_histogram_gates: true }), true)
  assert.equal(fileUsesHistogramGates({ is_af: true, uses_histogram_gates: 'FALSE' }), false)
  assert.equal(fileUsesHistogramGates({ is_af: false, uses_histogram_gates: 'TRUE' }), true)
})

test('SCCs with external negatives do not expose or retain negative histogram gates', () => {
  const sccWithExternalNegative = {
    filename: 'PE control.fcs',
    uses_histogram_gates: true,
    uses_negative_histogram_gate: false,
  }
  assert.equal(fileUsesNegativeHistogramGate(sccWithExternalNegative), false)
  assert.equal(fileUsesNegativeHistogramGate({ uses_histogram_gates: true }), true)

  const gates = {
    'positive:PE control.fcs': { type: 'positive' },
    'negative:PE control.fcs': { type: 'negative' },
    'negative:FITC cells.fcs': { type: 'negative' },
  }
  assert.deepEqual(pruneUnavailableNegativeGates(gates, [sccWithExternalNegative]), {
    'positive:PE control.fcs': { type: 'positive' },
    'negative:FITC cells.fcs': { type: 'negative' },
  })
})
