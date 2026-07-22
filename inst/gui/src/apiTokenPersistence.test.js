import assert from 'node:assert/strict'
import test from 'node:test'
import { selectApiToken } from './apiTokenPersistence.ts'

test('location token wins and can seed persistent tab authentication', () => {
  assert.equal(selectApiToken(' url-token ', 'dev-token', 'stored-token'), 'url-token')
})

test('live dev token repairs a clean URL before stale tab storage', () => {
  assert.equal(selectApiToken('', ' current-dev-token ', 'stale-token'), 'current-dev-token')
})

test('stored token keeps clean production URLs authenticated across reloads', () => {
  assert.equal(selectApiToken('', '', ' stored-token '), 'stored-token')
  assert.equal(selectApiToken('', '', ''), '')
})
