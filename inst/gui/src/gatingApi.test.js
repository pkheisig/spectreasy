import assert from 'node:assert/strict'
import test from 'node:test'

import { withGatingApiToken } from './gatingApi.js'

test('standalone gating requests include the backend token', () => {
  const options = withGatingApiToken({
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: '{}',
  }, 'secret-token')

  assert.equal(options.method, 'POST')
  assert.equal(options.body, '{}')
  assert.equal(options.headers.get('Content-Type'), 'application/json')
  assert.equal(options.headers.get('X-Spectreasy-Token'), 'secret-token')
})

test('gating requests omit the token header when the backend has no token', () => {
  const options = withGatingApiToken({}, '')

  assert.equal(options.headers.has('X-Spectreasy-Token'), false)
})
