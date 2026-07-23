import assert from 'node:assert/strict'
import test from 'node:test'
import { commonImmuneIdentityTemplate } from './identityTemplates.ts'

test('common immune templates use display names and ignore area suffixes', () => {
  const signatures = commonImmuneIdentityTemplate(['CD3-A', 'CD4-A', 'CD8-A', 'CD19-A', 'CD56-A', 'CD14-A'])
  assert.deepEqual(signatures.map((signature) => signature.name), [
    'CD4 T cell', 'CD8 T cell', 'B cell', 'NK cell', 'Monocyte',
  ])
  assert.deepEqual(signatures[0].positive_markers, ['CD3-A', 'CD4-A'])
})

test('common immune templates omit identities without required positive evidence', () => {
  const signatures = commonImmuneIdentityTemplate(['CD3', 'CD8', 'CD56'])
  assert.deepEqual(signatures.map((signature) => signature.name), ['CD8 T cell', 'NK cell'])
  assert.ok(signatures[0].negative_markers.includes('CD56'))
})
