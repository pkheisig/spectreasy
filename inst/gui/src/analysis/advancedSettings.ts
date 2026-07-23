import type { AnalysisAdvancedSettings, AnalysisMethod } from './types'

export function resolvedAdvancedSettings(
  methods: AnalysisMethod[],
  values: AnalysisAdvancedSettings,
): AnalysisAdvancedSettings {
  return Object.fromEntries(methods.map((method) => [
    method.id,
    Object.fromEntries(method.parameters.map((parameter) => [
      parameter.id,
      values[method.id]?.[parameter.id] ?? parameter.default,
    ])),
  ]))
}

export function advancedSettingsIssues(
  methods: AnalysisMethod[],
  values: AnalysisAdvancedSettings,
): string[] {
  const resolved = resolvedAdvancedSettings(methods, values)
  const issues: string[] = []
  for (const method of methods) {
    const settings = resolved[method.id] ?? {}
    for (const parameter of method.parameters) {
      const value = settings[parameter.id]
      if (parameter.type === 'number' || parameter.type === 'integer') {
        const numeric = Number(value)
        if (!Number.isFinite(numeric)) issues.push(`${method.name}: ${parameter.label} must be numeric.`)
        else if (parameter.minimum != null && numeric < parameter.minimum) issues.push(`${method.name}: ${parameter.label} is below its minimum.`)
        else if (parameter.maximum != null && numeric > parameter.maximum) issues.push(`${method.name}: ${parameter.label} is above its maximum.`)
      }
    }
    if (method.id === 'flowsom' && Number(settings.alpha_end) > Number(settings.alpha_start)) {
      issues.push('FlowSOM: final learning rate cannot exceed the initial rate.')
    }
    if (
      (method.id === 'wanderlust' || method.id === 'wishbone')
      && Number(settings.candidate_neighbors) <= Number(settings.neighbors)
    ) {
      issues.push(`${method.name}: candidate neighbors must exceed retained neighbors.`)
    }
  }
  return issues
}
