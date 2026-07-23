import assert from 'node:assert/strict'
import test from 'node:test'

import {
  adaptiveHistogramCofactor,
  adaptiveHistogramTransform,
  clearHistogramGatesForFile,
  dragHistogramGateInDisplaySpace,
  histogramDomainIncludingGates,
  padHistogramDomain,
  panHistogramDomain,
} from './histogramGates.js'

test('adaptive histogram transform calibrates independently for each control range', () => {
  const narrow = [-100, -60, -20, 0, 200, 1000, 5000, 10000]
  const wide = [-1000, -600, -200, 0, 20000, 100000, 500000, 1000000]

  assert.ok(adaptiveHistogramCofactor(wide) > adaptiveHistogramCofactor(narrow))
})

test('adaptive histogram transform compresses wide controls and remains invertible', () => {
  const values = [-1000, -500, 0, 5000, 200000, 1000000]
  const transform = adaptiveHistogramTransform(values)
  const transformed = values.map(transform.forward)

  assert.ok(transformed.at(-1) - transformed[0] < 20)
  values.forEach((value, index) => {
    assert.ok(Math.abs(transform.inverse(transformed[index]) - value) < 1e-6)
  })
})

test('adaptive histogram transform keeps a minority lower population visible', () => {
  const lowerPopulation = Array.from({ length: 15 }, (_, index) => 1000 + index * 140)
  const upperPopulation = Array.from({ length: 85 }, (_, index) => 500000 + index * 6000)
  const transform = adaptiveHistogramTransform([...lowerPopulation, ...upperPopulation])

  assert.ok(transform.forward(1000) > 1)
  assert.ok(transform.forward(1000000) < 10)
})

test('moving a histogram gate preserves its visible width', () => {
  const transform = adaptiveHistogramTransform([-1000, -500, 0, 5000, 200000, 1000000])
  const gate = [{ x: -500, y: 0 }, { x: 5000, y: 0 }]
  const initialDisplay = gate.map((vertex) => transform.forward(vertex.x))
  const moved = dragHistogramGateInDisplaySpace(gate, 1.75, transform)
  const movedDisplay = moved.map((vertex) => transform.forward(vertex.x))

  assert.ok(Math.abs((movedDisplay[1] - movedDisplay[0]) - (initialDisplay[1] - initialDisplay[0])) < 1e-10)
  assert.ok(Math.abs(movedDisplay[0] - initialDisplay[0] - 1.75) < 1e-10)
  assert.ok(Math.abs(movedDisplay[1] - initialDisplay[1] - 1.75) < 1e-10)
})

test('resizing a histogram gate moves only the dragged border in display space', () => {
  const transform = adaptiveHistogramTransform([-1000, -500, 0, 5000, 200000, 1000000])
  const gate = [{ x: -500, y: 0 }, { x: 5000, y: 0 }]
  const initialDisplay = gate.map((vertex) => transform.forward(vertex.x))
  const resized = dragHistogramGateInDisplaySpace(gate, -0.8, transform, 1)
  const resizedDisplay = resized.map((vertex) => transform.forward(vertex.x))

  assert.equal(resized[0].x, gate[0].x)
  assert.ok(Math.abs(resizedDisplay[1] - initialDisplay[1] + 0.8) < 1e-10)
})

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
