import test from 'node:test'
import assert from 'node:assert/strict'
import { buildViewSettingRows, parseViewSettings } from './gatingViewSettings.js'

test('viewport settings round-trip with global scatter and file histogram scopes', () => {
  const views = {
    cell: {
      global: { x: [10, 200], y: [30, 400] },
    },
    singlet: {
      global: { x: [11, 201], y: [31, 401] },
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
  assert.deepEqual(parsed.histogram, { 'one.fcs': { x: [-5, 55], y: null } })
  assert.equal(rows.some((row) => row.filename === 'removed.fcs'), false)
})

test('legacy control-type scatter views migrate to one fixed global envelope', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'cells', x_channel: 'view_cell', y_channel: 'x_domain', x: 10, y: 20 },
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'x_domain', x: 100, y: 200 },
    { gate_type: 'setting', scope: 'cells', x_channel: 'view_cell', y_channel: 'y_domain', x: 30, y: 40 },
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'y_domain', x: 300, y: 400 },
  ])
  assert.deepEqual(parsed.cell, { global: { x: [10, 200], y: [30, 400] } })
})

test('invalid or incomplete viewport rows are ignored', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'x_domain', x: 5, y: 5 },
    { gate_type: 'setting', scope: 'global', x_channel: 'point_size', y_channel: '', x: 2, y: '' },
  ])
  assert.deepEqual(parsed, { cell: {}, singlet: {}, histogram: {} })
})
