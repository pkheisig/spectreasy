import assert from 'node:assert/strict'
import test from 'node:test'

import {
  clearHistogramGatesForFile,
  histogramDomainIncludingGates,
  padHistogramDomain,
  panHistogramDomain,
} from './histogramGates.js'

test('histogram domain includes the complete negative and positive event tails', () => {
  const domain = histogramDomainIncludingGates(
    [0, 25000],
    [-442, 12, 25000, 99319],
    [{ vertices: [{ x: 73604 }, { x: 81688 }] }],
  )

  assert.deepEqual(domain, [-442, 99319])
})

test('histogram domain limits corrupt gates to the observed event range', () => {
  const domain = histogramDomainIncludingGates(
    [0, 25000],
    [0, 100000],
    [{ vertices: [{ x: 1e20 }, { x: 2e20 }] }],
  )

  assert.deepEqual(domain, [0, 100000])
})

test('histogram domain adds balanced display-space margin beyond both tails', () => {
  assert.deepEqual(padHistogramDomain([-2, 8], 0.05), [-2.5, 8.5])
})

test('histogram panning preserves the control-specific raw viewport width', () => {
  const initial = [69000, 95000]
  const shifted = panHistogramDomain(initial, 500, 500, [-100000, 200000])

  assert.deepEqual(shifted, [43000, 69000])
  assert.equal(shifted[1] - shifted[0], initial[1] - initial[0])
})

test('histogram panning clamps without compressing the raw viewport', () => {
  const initial = [3000, 29000]
  const shifted = panHistogramDomain(initial, 1000, 500, [-5000, 200000])

  assert.deepEqual(shifted, [-5000, 21000])
  assert.equal(shifted[1] - shifted[0], initial[1] - initial[0])
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
