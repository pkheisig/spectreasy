import type { AiQcGrade, AiQcPrivacy, AiQcResponse, AiQcScope, AiQcStatus } from './types.ts'

function scalar(value: unknown, fallback = ''): string {
  if (Array.isArray(value)) return scalar(value[0], fallback)
  if (value == null) return fallback
  return String(value)
}

function strings(value: unknown): string[] {
  if (Array.isArray(value)) return value.map((item) => scalar(item)).filter(Boolean)
  const valueString = scalar(value)
  return valueString ? [valueString] : []
}

function rows(value: unknown): Array<Record<string, unknown>> {
  if (Array.isArray(value)) return value.filter((row): row is Record<string, unknown> => Boolean(row && typeof row === 'object'))
  if (!value || typeof value !== 'object') return []
  const columns = value as Record<string, unknown>
  const length = Math.max(0, ...Object.values(columns).filter(Array.isArray).map((column) => column.length))
  return Array.from({ length }, (_, index) => Object.fromEntries(Object.entries(columns).map(([key, column]) => [key, Array.isArray(column) ? column[index] : column])))
}

const statuses = new Set<AiQcStatus>(['not_generated', 'generating', 'ready', 'stale', 'partial', 'failed'])
const scopes = new Set<AiQcScope>(['control', 'sample', 'combined'])
const grades = new Set<AiQcGrade>(['good', 'review', 'poor', 'not_graded'])

export function normalizeAiQcResponse(value: unknown): AiQcResponse {
  const raw = value && typeof value === 'object' ? value as Record<string, unknown> : {}
  const statusValue = scalar(raw.status, raw.error ? 'failed' : 'not_generated') as AiQcStatus
  const scopeValue = scalar(raw.scope, scalar(raw.default_scope, 'combined')) as AiQcScope
  const defaultScopeValue = scalar(raw.default_scope, scopeValue) as AiQcScope
  const gradeRows = (raw.grade_counts && typeof raw.grade_counts === 'object' ? raw.grade_counts : {}) as Record<string, unknown>
  const schema = (raw.schema && typeof raw.schema === 'object' ? raw.schema : {}) as Record<string, unknown>
  const profile = (raw.profile && typeof raw.profile === 'object' ? raw.profile : {}) as Record<string, unknown>
  const gradeCounts = Object.fromEntries([...grades].map((grade) => [grade, Math.max(0, Number(scalar(gradeRows[grade], '0')) || 0)])) as Record<AiQcGrade, number>
  return {
    status: statuses.has(statusValue) ? statusValue : 'failed',
    stale: scalar(raw.stale, 'false') === 'true',
    scope: scopes.has(scopeValue) ? scopeValue : 'combined',
    defaultScope: scopes.has(defaultScopeValue) ? defaultScopeValue : 'combined',
    availableScopes: strings(raw.available_scopes).filter((item): item is AiQcScope => scopes.has(item as AiQcScope)),
    privacyMode: ['standard', 'strict', 'none'].includes(scalar(raw.privacy_mode)) ? scalar(raw.privacy_mode) as AiQcPrivacy : 'standard',
    gradeCounts,
    findings: rows(raw.findings).map((item) => {
      const grade = scalar(item.grade, 'not_graded') as AiQcGrade
      return { metricId: scalar(item.metric_id ?? item.metricId), entity: scalar(item.entity), stage: scalar(item.stage), grade: grades.has(grade) ? grade : 'not_graded', explanation: scalar(item.explanation) }
    }),
    artifactPaths: strings(raw.artifact_paths), missingSections: strings(raw.missing_sections), warnings: strings(raw.warnings),
    prompt: scalar(raw.prompt), promptCharacters: Math.max(0, Number(scalar(raw.prompt_characters, '0')) || 0), estimatedTokens: Math.max(0, Number(scalar(raw.estimated_tokens, '0')) || 0),
    schemaName: scalar(schema.name), schemaVersion: scalar(schema.version), profileName: scalar(profile.name), profileVersion: scalar(profile.version),
    referenceN: Math.max(0, Number(scalar(profile.reference_n ?? profile.n_clean, '0')) || 0),
    sourceHashes: rows(raw.source_hashes).map((item) => ({ path: scalar(item.path), sha256: scalar(item.sha256) })).filter((item) => item.path && item.sha256),
    error: scalar(raw.error),
  }
}
