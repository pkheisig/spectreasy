export function scalarNullableNumber(value: unknown): number | null {
  const unboxed = Array.isArray(value) ? value[0] : value
  if (unboxed == null || unboxed === '') return null
  const parsed = Number(unboxed)
  return Number.isFinite(parsed) ? parsed : null
}
