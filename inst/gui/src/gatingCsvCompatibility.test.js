import test from 'node:test'
import assert from 'node:assert/strict'
import { reconcileGateCsvRows } from './gatingCsvCompatibility.js'

test('keeps group gates and current overrides while dropping missing-file rows', () => {
  const rows = [
    { gate_type: 'setting', filename: '', x_channel: 'point_size', x: 2 },
    { gate_type: 'cell', scope: 'cells', filename: '', vertex_index: 1, x: 10, y: 10 },
    { gate_type: 'singlet', scope: 'cells', filename: '', vertex_index: 1, x: 20, y: 20 },
    { gate_type: 'cell', scope: 'beads', filename: '', vertex_index: 1, x: 30, y: 30 },
    { gate_type: 'singlet', scope: 'beads', filename: '', vertex_index: 1, x: 40, y: 40 },
    { gate_type: 'cell', scope: 'file', filename: 'removed.fcs', vertex_index: 1, x: 'invalid', y: 'invalid' },
    { gate_type: 'positive', scope: 'file', filename: 'removed.fcs', vertex_index: 1, x: 100, y: 0 },
    { gate_type: 'positive', scope: 'file', filename: 'replacement.fcs', vertex_index: 1, x: 200, y: 0 },
  ]

  const result = reconcileGateCsvRows(rows, [{ filename: 'replacement.fcs' }])

  assert.deepEqual(result.rows, [rows[0], rows[1], rows[2], rows[3], rows[4], rows[7]])
  assert.deepEqual(result.ignoredFiles, ['removed.fcs'])
  assert.equal(result.ignoredRowCount, 2)
  assert.ok(result.rows.some((row) => row.scope === 'cells' && row.gate_type === 'cell'))
  assert.ok(result.rows.some((row) => row.scope === 'cells' && row.gate_type === 'singlet'))
  assert.ok(result.rows.some((row) => row.scope === 'beads' && row.gate_type === 'cell'))
  assert.ok(result.rows.some((row) => row.scope === 'beads' && row.gate_type === 'singlet'))
})
