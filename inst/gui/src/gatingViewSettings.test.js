import test from 'node:test'
import assert from 'node:assert/strict'
import { buildViewSettingRows, parseViewSettings } from './gatingViewSettings.js'

test('viewport settings round-trip with control-type and file scopes', () => {
  const views = {
    cell: {
      cells: { x: [10, 20], y: [30, 40] },
      beads: { x: [100, 200], y: [300, 400] },
    },
    singlet: {
      cells: { x: [11, 21], y: [31, 41] },
      beads: { x: [101, 201], y: [301, 401] },
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

test('invalid or incomplete viewport rows are ignored', () => {
  const parsed = parseViewSettings([
    { gate_type: 'setting', scope: 'beads', x_channel: 'view_cell', y_channel: 'x_domain', x: 5, y: 5 },
    { gate_type: 'setting', scope: 'global', x_channel: 'point_size', y_channel: '', x: 2, y: '' },
  ])
  assert.deepEqual(parsed, { cell: {}, singlet: {}, histogram: {} })
})
