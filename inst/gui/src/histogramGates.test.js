import assert from 'node:assert/strict'
import test from 'node:test'

import { clearHistogramGatesForFile, histogramDomainIncludingGates } from './histogramGates.js'

test('histogram domain expands to show gates inside the observed event range', () => {
  const domain = histogramDomainIncludingGates(
    [0, 25000],
    [-442, 12, 25000, 99319],
    [{ vertices: [{ x: 73604 }, { x: 81688 }] }],
  )

  assert.deepEqual(domain, [0, 81688])
})

test('histogram domain limits corrupt gates to the observed event range', () => {
  const domain = histogramDomainIncludingGates(
    [0, 25000],
    [0, 100000],
    [{ vertices: [{ x: 1e20 }, { x: 2e20 }] }],
  )

  assert.deepEqual(domain, [0, 100000])
})

test('clear removes both histogram gates for only the selected file', () => {
  const gates = {
    'positive:a.fcs': { type: 'positive' },
    'negative:a.fcs': { type: 'negative' },
    'positive:b.fcs': { type: 'positive' },
    'negative:b.fcs': { type: 'negative' },
    'cell:cells': { type: 'cell' },
  }

  assert.deepEqual(clearHistogramGatesForFile(gates, 'a.fcs'), {
    'positive:b.fcs': { type: 'positive' },
    'negative:b.fcs': { type: 'negative' },
    'cell:cells': { type: 'cell' },
  })
})
