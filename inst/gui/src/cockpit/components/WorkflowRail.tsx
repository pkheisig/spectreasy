import { BarChart3, Beaker, BookOpen, Boxes, CircleDot, FileText, FlaskConical, FolderKanban, Gauge, Settings2, SlidersHorizontal } from 'lucide-react'
import type { LucideIcon } from 'lucide-react'
import type { SectionId } from '../types'

type RailProps = { activeSection: SectionId; onChange: (section: SectionId) => void }

const sections: Array<{ id: SectionId; number: string; title: string; icon: LucideIcon; state: 'complete' | 'warning' | 'ready' | 'idle' | 'blocked'; note: string }> = [
  { id: 'overview', number: '01', title: 'Project setup', icon: FolderKanban, state: 'complete', note: 'Scanned & ready' },
  { id: 'controls', number: '02', title: 'Controls', icon: CircleDot, state: 'complete', note: 'Reference current' },
  { id: 'samples', number: '03', title: 'Samples', icon: Beaker, state: 'warning', note: '2 files changed' },
  { id: 'reports', number: '04', title: 'Reports', icon: FileText, state: 'warning', note: '1 report stale' },
  { id: 'matrix', number: '05', title: 'Matrix review', icon: SlidersHorizontal, state: 'ready', note: '3 matrices found' },
  { id: 'panel', number: '06', title: 'Panel builder', icon: Boxes, state: 'ready', note: 'Available' },
  { id: 'af', number: '07', title: 'AF library', icon: FlaskConical, state: 'ready', note: '1 saved profile' },
  { id: 'comparison', number: '08', title: 'Method comparison', icon: BarChart3, state: 'idle', note: 'Experimental' },
  { id: 'simulator', number: '09', title: 'Synthetic SCC', icon: Gauge, state: 'blocked', note: 'Needs trusted matrix' },
  { id: 'settings', number: '10', title: 'Settings & logs', icon: Settings2, state: 'idle', note: 'Project tools' },
]

export function WorkflowRail({ activeSection, onChange }: RailProps) {
  return (
    <nav className="workflow-rail" aria-label="Workflow sections">
      <div className="rail-intro"><span className="eyebrow">The cockpit</span><span>10 workspaces</span></div>
      <div className="rail-list">
        {sections.map(({ id, number, title, icon: Icon, state, note }) => (
          <button key={id} className={`rail-item ${activeSection === id ? 'is-active' : ''}`} onClick={() => onChange(id)}>
            <span className="rail-number">{number}</span>
            <span className="rail-icon"><Icon size={17} /></span>
            <span className="rail-copy"><strong>{title}</strong><small>{note}</small></span>
            <span className={`rail-state rail-state-${state}`} />
          </button>
        ))}
      </div>
      <div className="rail-help"><BookOpen size={15} /><span>Workflow guide</span><span className="shortcut-key">?</span></div>
    </nav>
  )
}
