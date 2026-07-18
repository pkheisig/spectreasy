import { useEffect, useState } from 'react'
import { loadAfProfileSimilarity } from '../api'
import { GuiSelect } from './GuiSelect'
import type { AfProfileSimilarityPayload } from '../api'

type Props = { primaryName: string; profileNames: string[] }

export function AfSimilarityHeatmap({ primaryName, profileNames }: Props) {
  const [compareName, setCompareName] = useState('')
  const [payload, setPayload] = useState<AfProfileSimilarityPayload | null>(null)
  const effectiveCompareName = compareName === primaryName ? '' : compareName

  useEffect(() => {
    let cancelled = false
    void loadAfProfileSimilarity(primaryName, effectiveCompareName).then((result) => { if (!cancelled) setPayload(result) })
    return () => { cancelled = true }
  }, [effectiveCompareName, primaryName])

  const labels = payload?.labels ?? []
  const memberships = payload?.profileMembership ?? []
  const showValues = labels.length <= 20
  const cellSize = showValues ? 44 : 26
  const groupBoundary = memberships.findIndex((membership) => membership !== memberships[0])
  return (
    <section className="surface-card af-similarity-card">
      <div className="card-toolbar">
        <div><span className="eyebrow">AF-band similarity</span><h2>{primaryName}</h2></div>
        <label className="af-compare-select">Compare with
          <GuiSelect value={effectiveCompareName} onChange={(event) => setCompareName(event.target.value)}>
            <option value="">None</option>
            {profileNames.filter((name) => name !== primaryName).map((name) => <option key={name} value={name}>{name}</option>)}
          </GuiSelect>
        </label>
      </div>
      {!payload?.success ? <p className="profile-status">{payload?.message || 'Loading similarity…'}</p> : (
        <div className="af-heatmap-scroll">
          <div className="af-heatmap" style={{ gridTemplateColumns: `minmax(140px,auto) repeat(${labels.length},${cellSize}px)` }}>
            <span />
            {labels.map((label) => <span className="af-heatmap-column" key={`column-${label}`} title={label}>{label}</span>)}
            {labels.map((rowLabel, row) => [
              <span className="af-heatmap-row" key={`row-${rowLabel}`} title={rowLabel}>{rowLabel}</span>,
              ...labels.map((columnLabel, column) => {
                const value = payload.similarity[row]?.[column] ?? 0
                return <span
                  className={`af-heatmap-tile${row === groupBoundary ? ' is-group-row' : ''}${column === groupBoundary ? ' is-group-column' : ''}`}
                  key={`${rowLabel}-${columnLabel}`}
                  title={`${rowLabel} × ${columnLabel}: ${value.toFixed(3)}`}
                  style={{ height: cellSize, background: `color-mix(in srgb, var(--paper) ${(1 - value) * 90}%, var(--cockpit-accent))` }}
                >{showValues ? value.toFixed(2) : ''}</span>
              }),
            ])}
          </div>
        </div>
      )}
    </section>
  )
}
