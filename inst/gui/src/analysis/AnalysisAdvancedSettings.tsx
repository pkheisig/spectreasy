import { ChevronDown, RotateCcw } from 'lucide-react'
import type {
  AnalysisAdvancedSettings,
  AnalysisMethod,
  AnalysisParameterValue,
} from './types'

type Props = {
  methods: AnalysisMethod[]
  values: AnalysisAdvancedSettings
  onChange: (methodId: string, parameterId: string, value: AnalysisParameterValue) => void
  onReset: (methodId: string) => void
}

function parameterValue(
  values: AnalysisAdvancedSettings,
  method: AnalysisMethod,
  parameterId: string,
): AnalysisParameterValue {
  const parameter = method.parameters.find((candidate) => candidate.id === parameterId)
  return values[method.id]?.[parameterId] ?? parameter?.default ?? ''
}

export function AnalysisAdvancedSettingsPanel({ methods, values, onChange, onReset }: Props) {
  const configurable = methods.filter((method, index) => (
    method.parameters.length > 0
    && methods.findIndex((candidate) => candidate.id === method.id) === index
  ))
  if (!configurable.length) return null

  return (
    <details className="analysis-advanced-settings">
      <summary>
        <span><ChevronDown size={12} /> Advanced settings</span>
        <small>{configurable.map((method) => method.name).join(' · ')}</small>
      </summary>
      <div className="analysis-advanced-body">
        {configurable.map((method) => (
          <section key={method.id} className="analysis-advanced-method">
            <header>
              <div><strong>{method.name}</strong><small>{method.family}</small></div>
              <button type="button" onClick={() => onReset(method.id)} title={`Reset ${method.name} settings`}>
                <RotateCcw size={11} /> Reset
              </button>
            </header>
            <div className="analysis-advanced-grid">
              {method.parameters.map((parameter) => {
                const value = parameterValue(values, method, parameter.id)
                if (parameter.type === 'boolean') {
                  return (
                    <label key={parameter.id} className="analysis-advanced-toggle" title={parameter.description}>
                      <input
                        type="checkbox"
                        checked={Boolean(value)}
                        onChange={(event) => onChange(method.id, parameter.id, event.target.checked)}
                      />
                      <span>{parameter.label}<small>{parameter.description}</small></span>
                    </label>
                  )
                }
                if (parameter.type === 'select') {
                  return (
                    <label key={parameter.id} title={parameter.description}>
                      <span>{parameter.label}</span>
                      <select
                        value={String(value)}
                        onChange={(event) => onChange(method.id, parameter.id, event.target.value)}
                      >
                        {(parameter.choices ?? []).map((choice) => (
                          <option key={choice.value} value={choice.value}>{choice.label}</option>
                        ))}
                      </select>
                      {parameter.description ? <small>{parameter.description}</small> : null}
                    </label>
                  )
                }
                return (
                  <label key={parameter.id} title={parameter.description}>
                    <span>{parameter.label}</span>
                    <input
                      type="number"
                      value={Number(value)}
                      min={parameter.minimum ?? undefined}
                      max={parameter.maximum ?? undefined}
                      step={parameter.step ?? (parameter.type === 'integer' ? 1 : 'any')}
                      onChange={(event) => onChange(method.id, parameter.id, Number(event.target.value))}
                    />
                    {parameter.description ? <small>{parameter.description}</small> : null}
                  </label>
                )
              })}
            </div>
          </section>
        ))}
      </div>
    </details>
  )
}
