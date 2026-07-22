export type DetectorLaserKey = 'DeepUV' | 'UV' | 'V' | 'B' | 'YG' | 'R' | 'IR' | 'Other'

export type DetectorAxisEntry = {
  detector: string
  label?: string
  laser?: string
  emission?: number
}

export const DETECTOR_AXIS_FOOTER_HEIGHT = 132
export const DETECTOR_LABEL_ROTATION = -90
export const detectorAxisChartWidth = (detectorCount: number) => Math.max(1040, detectorCount * 27)

type DetectorAxisInput = string | DetectorAxisEntry

export type DetectorLaserSegment = {
  key: DetectorLaserKey
  label: string
  wavelength: string
  color: string
  textColor: string
  startIndex: number
  endIndex: number
}

const LASER_META: Record<DetectorLaserKey, Omit<DetectorLaserSegment, 'key' | 'startIndex' | 'endIndex'>> = {
  DeepUV: { label: 'DUV', wavelength: '320 nm', color: '#4c1d95', textColor: '#ffffff' },
  UV: { label: 'UV', wavelength: '355 nm', color: '#7c3aed', textColor: '#ffffff' },
  V: { label: 'V', wavelength: '405 nm', color: '#a21caf', textColor: '#ffffff' },
  B: { label: 'B', wavelength: '488 nm', color: '#2563eb', textColor: '#ffffff' },
  YG: { label: 'YG', wavelength: '561 nm', color: '#9acb30', textColor: '#17220b' },
  R: { label: 'R', wavelength: '640 nm', color: '#ef3e36', textColor: '#ffffff' },
  IR: { label: 'IR', wavelength: '808 nm', color: '#7f1d1d', textColor: '#ffffff' },
  Other: { label: 'Other', wavelength: '', color: '#64748b', textColor: '#ffffff' },
}

function laserMetadataKey(laser = ''): DetectorLaserKey | null {
  const key = laser.toLowerCase().replace(/[^a-z0-9]/g, '')
  if (key === 'deepuv' || key === 'duv' || key === '320') return 'DeepUV'
  if (key === 'uv' || key === 'ultraviolet' || key === '349' || key === '355') return 'UV'
  if (key === 'v' || key === 'violet' || key === '405') return 'V'
  if (key === 'b' || key === 'blue' || key === '488') return 'B'
  if (key === 'yg' || key === 'yellowgreen' || key === '561') return 'YG'
  if (key === 'r' || key === 'red' || key === '637' || key === '640') return 'R'
  if (key === 'ir' || key === 'infrared' || key === '781' || key === '808') return 'IR'
  if (key === 'other') return 'Other'
  return null
}

function numericLaserKey(value: string): DetectorLaserKey | null {
  const match = value.match(/^(320|349|355|405|488|561|637|640|781|808)(?:NM|CH|\b)/i)
  if (!match) return null
  return laserMetadataKey(match[1])
}

export function detectorLaserKey(input: DetectorAxisInput): DetectorLaserKey {
  const entry = typeof input === 'string' ? { detector: input } : input
  const metadataKey = laserMetadataKey(entry.laser)
  if (metadataKey) return metadataKey
  const key = entry.detector.toUpperCase().replace(/\s+/g, '').replace(/-A$/i, '')
  const numericKey = numericLaserKey(key)
  if (numericKey) return numericKey
  if (/^(?:D?UV)/.test(key)) return key.startsWith('DUV') ? 'DeepUV' : 'UV'
  if (/^IR/.test(key)) return 'IR'
  if (/^YG/.test(key) || /^[YG]\d/.test(key)) return 'YG'
  if (/^V\d/.test(key)) return 'V'
  if (/^B\d/.test(key)) return 'B'
  if (/^R\d/.test(key)) return 'R'
  const labelKey = numericLaserKey((entry.label || '').trim().toUpperCase())
  if (labelKey) return labelKey
  return 'Other'
}

export function detectorLaserMeta(key: DetectorLaserKey, cytometer: unknown = ''): Omit<DetectorLaserSegment, 'key' | 'startIndex' | 'endIndex'> {
  const metadata = { ...LASER_META[key] }
  const cytometerValue = Array.isArray(cytometer) ? cytometer[0] : cytometer
  const cytometerKey = String(cytometerValue ?? '').toLowerCase().replace(/[^a-z0-9]/g, '')
  if (key === 'UV' && cytometerKey.includes('xenith')) metadata.wavelength = '349 nm'
  if (key === 'R' && /(discover|id7000|xenith)/.test(cytometerKey)) metadata.wavelength = '637 nm'
  if (key === 'IR' && cytometerKey.includes('xenith')) metadata.wavelength = '781 nm'
  return metadata
}

export function detectorLaserSegments(detectors: DetectorAxisInput[], cytometer: unknown = ''): DetectorLaserSegment[] {
  const segments: DetectorLaserSegment[] = []
  detectors.forEach((detector, index) => {
    const key = detectorLaserKey(detector)
    const previous = segments[segments.length - 1]
    if (previous?.key === key) {
      previous.endIndex = index
      return
    }
    segments.push({ key, ...detectorLaserMeta(key, cytometer), startIndex: index, endIndex: index })
  })
  return segments
}

export function wavelengthToColor(wavelength: number): string {
  let r = 0
  let g = 0
  let b = 0
  if (wavelength >= 350 && wavelength < 440) {
    r = -(wavelength - 440) / (440 - 350)
    b = 1
  } else if (wavelength >= 440 && wavelength < 490) {
    g = (wavelength - 440) / (490 - 440)
    b = 1
  } else if (wavelength >= 490 && wavelength < 510) {
    g = 1
    b = -(wavelength - 510) / (510 - 490)
  } else if (wavelength >= 510 && wavelength < 580) {
    r = (wavelength - 510) / (580 - 510)
    g = 1
  } else if (wavelength >= 580 && wavelength < 645) {
    r = 1
    g = -(wavelength - 645) / (645 - 580)
  } else if (wavelength >= 645 && wavelength <= 780) {
    r = 1
  } else if (wavelength > 780) {
    r = 0.5
    b = 0.2
  }

  let factor = 1
  if (wavelength >= 350 && wavelength < 420) {
    factor = 0.3 + 0.7 * (wavelength - 350) / (420 - 350)
  } else if (wavelength > 700 && wavelength <= 780) {
    factor = 0.3 + 0.7 * (780 - wavelength) / (780 - 700)
  } else if (wavelength > 780) {
    factor = 0.3
  }

  return `rgb(${Math.round(r * factor * 255)}, ${Math.round(g * factor * 255)}, ${Math.round(b * factor * 255)})`
}

export function mapDetectorToEmission(detectorName: string): number {
  const clean = detectorName.replace(/-A$/i, '').toUpperCase()
  const parenthetical = clean.match(/\((\d{3})\)/)
  if (parenthetical) return Number(parenthetical[1])
  const filterCenter = clean.match(/-\s*(\d{3})(?=\/)/)
  if (filterCenter) return Number(filterCenter[1])
  const id7000 = clean.match(/^(320|355|405|488|561|637|808)CH(\d+)$/)
  if (id7000) {
    const startChannel: Record<string, number> = { 320: 1, 355: 1, 405: 1, 488: 4, 561: 10, 637: 17, 808: 36 }
    const startEmission: Record<string, number> = { 320: 350, 355: 370, 405: 420, 488: 500, 561: 570, 637: 660, 808: 810 }
    return startEmission[id7000[1]] + (Number(id7000[2]) - startChannel[id7000[1]]) * 15
  }
  const group = detectorLaserKey(clean)
  const index = Number.parseInt(clean.replace(/^[A-Z]+/, ''), 10)
  const wavelengths: Record<DetectorLaserKey, Record<number, number>> = {
    UV: { 1: 370, 2: 395, 3: 420, 4: 440, 5: 450, 6: 480, 7: 480, 8: 500, 9: 520, 10: 550, 11: 570, 12: 580, 13: 600, 14: 660, 15: 750, 16: 800 },
    V: { 1: 420, 2: 440, 3: 450, 4: 480, 5: 480, 6: 500, 7: 550, 8: 570, 9: 580, 10: 600, 11: 660, 12: 680, 13: 690, 14: 700, 15: 730, 16: 780 },
    B: { 1: 500, 2: 520, 3: 550, 4: 550, 5: 570, 6: 580, 7: 600, 8: 600, 9: 660, 10: 680, 11: 690, 12: 700, 13: 750, 14: 780 },
    YG: { 1: 570, 2: 580, 3: 600, 4: 600, 5: 660, 6: 680, 7: 700, 8: 730, 9: 750, 10: 780 },
    R: { 1: 660, 2: 680, 3: 700, 4: 730, 5: 730, 6: 750, 7: 780, 8: 800 },
    DeepUV: {},
    IR: { 1: 810, 2: 825, 3: 840, 4: 855, 5: 870, 6: 885 },
    Other: {},
  }
  return wavelengths[group][index] ?? 0
}

export function detectorSpectralColors(detectors: DetectorAxisInput[]): string[] {
  const segments = detectorLaserSegments(detectors)
  const colors = Array(detectors.length).fill('#64748b') as string[]
  segments.forEach((segment) => {
    const count = segment.endIndex - segment.startIndex + 1
    for (let index = segment.startIndex; index <= segment.endIndex; index += 1) {
      const position = count <= 1 ? 0.5 : (index - segment.startIndex) / (count - 1)
      const entry = detectors[index]
      const mappedWavelength = typeof entry === 'string'
        ? mapDetectorToEmission(entry)
        : (Number.isFinite(Number(entry.emission)) && Number(entry.emission) > 0
            ? Number(entry.emission)
            : mapDetectorToEmission(entry.label || entry.detector))
      const fallbackWavelength = 420 + position * 300
      colors[index] = wavelengthToColor(mappedWavelength || fallbackWavelength)
    }
  })
  return colors
}

export function compactDetectorLabel(label: string): string {
  return label.trim().replace(/-A$/i, '')
}
