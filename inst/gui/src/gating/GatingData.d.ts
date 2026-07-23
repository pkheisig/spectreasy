export function computeDensityBuckets<T>(
  points: T[],
  xField: keyof T,
  yField: keyof T,
  xDomain?: [number, number],
  yDomain?: [number, number],
): number[][]
