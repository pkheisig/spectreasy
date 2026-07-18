import { createPortal } from 'react-dom'
import type { AfProfileData } from '../api'

type Props = { profile: AfProfileData; onClose: () => void }

function recorded(value: string | number | undefined): string {
  return value == null || value === '' ? 'Not recorded' : String(value)
}

function formattedDate(value: string | undefined): string {
  if (!value) return 'Not recorded'
  const parsed = new Date(value)
  return Number.isFinite(parsed.getTime()) ? new Intl.DateTimeFormat(undefined, { dateStyle: 'medium', timeStyle: value.includes(':') ? 'short' : undefined }).format(parsed) : value
}

export function AfProfileInfoDialog({ profile, onClose }: Props) {
  const rows = [
    ['Profile', profile.name],
    ['Cytometer', profile.metadata.cytometer],
    ['Acquisition date', profile.metadata.acquisitionDate ? formattedDate(`${profile.metadata.acquisitionDate}T00:00:00`) : undefined],
    ['Tissue', profile.metadata.tissue],
    ['Sample type', profile.metadata.sampleType],
    ['Preprocessing', profile.metadata.preprocessing],
    ['Source FCS', profile.metadata.sourceFcs],
    ['Created', formattedDate(profile.metadata.createdAt)],
    ['Bands', profile.metadata.bandCount ?? profile.spectra.length],
    ['Detectors', profile.detectors.length],
    ['Spectreasy', profile.metadata.spectreasyVersion],
  ] as const
  return createPortal(
    <div className="cockpit-confirm-overlay" role="presentation" onMouseDown={onClose}>
      <div className="cockpit-confirm af-profile-info-dialog" role="dialog" aria-modal="true" aria-labelledby="af-profile-info-title" onMouseDown={(event) => event.stopPropagation()}>
        <h2 id="af-profile-info-title">AF profile information</h2>
        <dl>{rows.map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{recorded(value)}</dd></div>)}</dl>
        <div><button className="button button-primary" type="button" onClick={onClose}>Close</button></div>
      </div>
    </div>,
    document.body,
  )
}
