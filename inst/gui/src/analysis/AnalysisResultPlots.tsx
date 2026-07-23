import { Download, Plus } from 'lucide-react'
import { useState } from 'react'
import { AnalysisResultPlot } from './AnalysisResultPlot'
import { downloadText, resultTableCsv } from './resultExport'
import type { AnalysisRunResult } from './types'

export function AnalysisResultPlots({ result }: { result: AnalysisRunResult }) {
  const [plots, setPlots] = useState(() => [crypto.randomUUID()])

  return (
    <div className="analysis-result-plots">
      <div className="analysis-result-plot-actions">
        <button type="button" onClick={() => setPlots((current) => [...current, crypto.randomUUID()])}>
          <Plus size={13} /> Add plot
        </button>
        <button
          type="button"
          onClick={() => downloadText(
            `${result.metadata.analysis_id}-plot-data.csv`,
            resultTableCsv(result),
            'text/csv;charset=utf-8',
          )}
        >
          <Download size={13} /> Export plot data
        </button>
      </div>
      <div className="analysis-result-plot-grid">
        {plots.map((plotId) => (
          <section key={plotId} className="analysis-result-plot-card">
            <AnalysisResultPlot
              result={result}
              canRemove
              onRemove={() => setPlots((current) => current.filter((id) => id !== plotId))}
            />
          </section>
        ))}
        {plots.length === 0 ? (
          <div className="analysis-empty-result-plots">
            <strong>No result plots</strong>
            <span>The analysis object is still available. Add a plot to visualize it again.</span>
          </div>
        ) : null}
      </div>
    </div>
  )
}
