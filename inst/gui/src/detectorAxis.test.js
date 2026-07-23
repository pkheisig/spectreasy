import assert from 'node:assert/strict'
import test from 'node:test'
import {
  DETECTOR_LABEL_ROTATION,
  compactDetectorLabel,
  detectorLaserKey,
  detectorLaserMeta,
  detectorLaserSegments,
  detectorSpectralColors,
} from './detectorAxis.ts'

test('groups standard spectral detector names into contiguous laser segments', () => {
  const detectors = ['UV1-A', 'UV2-A', 'V1-A', 'B1-A', 'B2-A', 'YG1-A', 'R1-A']
  assert.deepEqual(
    detectorLaserSegments(detectors).map(({ key, startIndex, endIndex }) => ({ key, startIndex, endIndex })),
    [
      { key: 'UV', startIndex: 0, endIndex: 1 },
      { key: 'V', startIndex: 2, endIndex: 2 },
      { key: 'B', startIndex: 3, endIndex: 4 },
      { key: 'YG', startIndex: 5, endIndex: 5 },
      { key: 'R', startIndex: 6, endIndex: 6 },
    ],
  )
})

test('builds a distinct spectral color sequence inside each laser group', () => {
  const colors = detectorSpectralColors(['UV1-A', 'UV2-A', 'UV3-A', 'B1-A', 'B2-A', 'B3-A'])
  assert.equal(colors.length, 6)
  assert.notEqual(colors[0], colors[1])
  assert.notEqual(colors[1], colors[2])
  assert.notEqual(colors[2], colors[3])
})

test('normalizes detector group and display labels', () => {
  assert.equal(DETECTOR_LABEL_ROTATION, -90)
  assert.equal(detectorLaserKey('YG10-A'), 'YG')
  assert.equal(detectorLaserKey('V4-A'), 'V')
  assert.equal(compactDetectorLabel('UV12-A'), 'UV12')
})

test('uses backend metadata and supported cytometer naming schemes', () => {
  assert.equal(detectorLaserKey({ detector: 'FL16-A', label: '405nm - 420/10-A', laser: 'Violet', emission: 420 }), 'V')
  assert.equal(detectorLaserKey({ detector: '355CH1-A', laser: 'UV', emission: 370 }), 'UV')
  assert.equal(detectorLaserKey({ detector: '637CH17-A', laser: 'Red', emission: 660 }), 'R')
  assert.equal(detectorLaserKey({ detector: 'UV1 (375)-A', laser: 'UV', emission: 375 }), 'UV')
  assert.equal(detectorLaserKey('808CH36-A'), 'IR')
})

test('keeps laser sequence and spectrum colors complete across panel configurations', () => {
  const configurations = [
    {
      name: 'aurora 4L UV',
      entries: ['UV1-A', 'V1-A', 'B1-A', 'R1-A'],
      expected: ['UV', 'V', 'B', 'R'],
    },
    {
      name: 'aurora 4L YG',
      entries: ['V1-A', 'B1-A', 'YG1-A', 'R1-A'],
      expected: ['V', 'B', 'YG', 'R'],
    },
    {
      name: 'id7000 5L',
      entries: ['355CH1-A', '405CH1-A', '488CH4-A', '561CH10-A', '637CH17-A'],
      expected: ['UV', 'V', 'B', 'YG', 'R'],
    },
    {
      name: 'xenith',
      entries: [
        { detector: 'FL16-A', label: '405nm - 420/10-A', laser: 'Violet', emission: 420 },
        { detector: 'FL10-A', label: '488nm - 530/30-A', laser: 'Blue', emission: 530 },
        { detector: 'FL6-A', label: '561nm - 586/15-A', laser: 'YellowGreen', emission: 586 },
        { detector: 'FL2-A', label: '637nm - 670/30-A', laser: 'Red', emission: 670 },
      ],
      expected: ['V', 'B', 'YG', 'R'],
    },
  ]

  configurations.forEach(({ name, entries, expected }) => {
    const segments = detectorLaserSegments(entries)
    assert.deepEqual(segments.map(segment => segment.key), expected, name)
    assert.equal(detectorSpectralColors(entries).length, entries.length, name)
    assert.equal(segments.some(segment => segment.key === 'Other'), false, name)
  })
})

test('uses cytometer-specific excitation labels', () => {
  assert.equal(detectorLaserMeta('UV', 'xenith').wavelength, '349 nm')
  assert.equal(detectorLaserMeta('R', 'discover').wavelength, '637 nm')
  assert.equal(detectorLaserMeta('R', 'id7000').wavelength, '637 nm')
  assert.equal(detectorLaserMeta('R', 'aurora').wavelength, '640 nm')
  assert.equal(detectorLaserMeta('IR', 'xenith').wavelength, '781 nm')
})

test('accepts R scalar wrappers for cytometer names', () => {
  assert.equal(detectorLaserMeta('UV', ['xenith']).wavelength, '349 nm')
  assert.equal(detectorLaserSegments(['R1-A'], ['discover'])[0].wavelength, '637 nm')
})
