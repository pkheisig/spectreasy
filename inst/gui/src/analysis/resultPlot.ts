import type { AnalysisRunEvent, AnalysisRunResult } from './types'

export type CoordinateKey = `dimension_${number}`

export function resultCoordinates(result: AnalysisRunResult): CoordinateKey[] {
  const reported = Number(result.metadata.coordinate_count)
  const count = Number.isFinite(reported) && reported >= 2
    ? Math.floor(reported)
    : Math.max(2, ...result.events.flatMap((event) => Object.keys(event)
      .map((key) => /^dimension_(\d+)$/.exec(key)?.[1])
      .map(Number)
      .filter(Number.isFinite)))
  return Array.from({ length: count }, (_, index) => `dimension_${index + 1}` as CoordinateKey)
}

export function eventValue(event: AnalysisRunEvent, key: CoordinateKey): number {
  const value = event[key]
  return Number.isFinite(value) ? Number(value) : 0
}

export function coordinateLabel(result: AnalysisRunResult, key: CoordinateKey): string {
  const index = Number(key.replace('dimension_', '')) - 1
  return result.metadata.coordinate_labels?.[index] || `Dimension ${index + 1}`
}
