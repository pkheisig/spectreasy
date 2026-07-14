import assert from 'node:assert/strict'
import test from 'node:test'

import { decodeSpectrumData, spectrumColor } from './spectrumData.js'

test('decodes compact column-major spectrum counts', () => {
  const values = new Uint32Array([1, 2, 3, 4000])
  const data = Buffer.from(values.buffer).toString('base64')
  const decoded = decodeSpectrumData({
    format: 'spectrum-histogram-v1',
    channels: 'B1-A',
    labels: 'B1-A',
    bin_mid: [0.5, 1.5],
    fill_limits: [0, 4],
    counts: {
      format: 'uint32-column-major',
      rows: 2,
      columns: 2,
      data,
    },
  })

  assert.deepEqual(decoded.channels, ['B1-A'])
  assert.deepEqual(decoded.labels, ['B1-A'])
  assert.deepEqual([...decoded.counts.values], [1, 2, 3, 4000])
})

test('rejects malformed compact spectrum counts', () => {
  assert.equal(decodeSpectrumData({ format: 'other' }), null)
  assert.equal(decodeSpectrumData({
    format: 'spectrum-histogram-v1',
    counts: { format: 'uint32-column-major', rows: 2, columns: 2, data: 'AAAA' },
  }), null)
})

test('spectrum palette clamps at blue and red endpoints', () => {
  assert.equal(spectrumColor(-1), 'rgb(0, 0, 255)')
  assert.equal(spectrumColor(1), 'rgb(255, 0, 0)')
  assert.equal(spectrumColor(2), 'rgb(255, 0, 0)')
})
