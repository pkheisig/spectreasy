import test from 'node:test'
import assert from 'node:assert/strict'
import { buildViewSettingRows, parseViewSettings } from './gatingViewSettings.js'

test('viewport settings round-trip independently for cell and bead controls', () => {
  const views = {
    cell: {
      cells: { x: [10, 200], y: [30, 400] },
      beads: { x: [100, 2000], y: [300, 4000] },
    },
    singlet: {
      cells: { x: [11, 201], y: [31, 401] },
      beads: { x: [101, 2001], y: [301, 4001] },
    },
    histogram: {
      'one.fcs': { x: [-5, 55], y: null },
      'removed.fcs': { x: [0, 10], y: null },
    },
  }
  const rows = buildViewSettingRows(views, [{ filename: 'one.fcs' }])
  const parsed = parseViewSettings(rows)

  assert.deepEqual(parsed.cell, views.cell)
  assert.deepEqual(parsed.singlet, views.singlet)
  assert.deepEqual(parsed.histogram, {})
  assert.equal(rows.some((row) => row.x_channel === 'view_histogram'), false)
  assert.deepEqual(new Set(rows.map((row) => row.scope)), new Set(['cells', 'beads']))
})

test('legacy saved histogram view rows are ignored', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'file', filename: 'one.fcs', x_channel: 'view_histogram', y_channel: 'x_domain', x: -5, y: 55 },
  ])

  assert.deepEqual(parsed.histogram, {})
})

test('legacy global scatter views migrate into both control-type scopes', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'global', x_channel: 'view_cell', y_channel: 'x_domain', x: 10, y: 20 },
    { gate_type: 'setting', scope: 'global', x_channel: 'view_cell', y_channel: 'y_domain', x: 30, y: 40 },
  ])
  assert.deepEqual(parsed.cell, {
    cells: { x: [10, 20], y: [30, 40] },
    beads: { x: [10, 20], y: [30, 40] },
  })
})

test('control-specific scatter views remain separate', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'cells', x_channel: 'view_cell', y_channel: 'x_domain', x: 10, y: 20 },
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'x_domain', x: 100, y: 200 },
  ])

  assert.deepEqual(parsed.cell.cells.x, [10, 20])
  assert.deepEqual(parsed.cell.beads.x, [100, 200])
})

test('invalid or incomplete viewport rows are ignored', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'x_domain', x: 5, y: 5 },
    { gate_type: 'setting', scope: 'global', x_channel: 'point_size', y_channel: '', x: 2, y: '' },
  ])
  assert.deepEqual(parsed, { cell: {}, singlet: {}, histogram: {} })
})
