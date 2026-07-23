import test from 'node:test'
import assert from 'node:assert/strict'
import { scalarNullableNumber } from './workspaceNormalization.ts'

test('nullable workspace numbers preserve API-wrapped null values', () => {
  assert.equal(scalarNullableNumber(null), null)
  assert.equal(scalarNullableNumber([null]), null)
  assert.equal(scalarNullableNumber([]), null)
  assert.equal(scalarNullableNumber(''), null)
  assert.equal(scalarNullableNumber(['12.5']), 12.5)
  assert.equal(scalarNullableNumber(['not-a-number']), null)
})
