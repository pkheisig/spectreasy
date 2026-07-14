import {
  ChevronRight,
  FileOutput,
  FolderCog,
  Wrench,
} from 'lucide-react'
import { useState } from 'react'
import type { PointerEvent as ReactPointerEvent } from 'react'
import type { LucideIcon } from 'lucide-react'
import type { ProjectState, SectionId } from '../types'

type RailProps = {
  activeSection: SectionId
  project: ProjectState
  showCounts: boolean
  width: number
  onWidthChange: (width: number) => void
  onChange: (section: SectionId) => void
}

type NavigationItem = { id: SectionId; title: string; detail: string }
type NavigationGroup = { id: string; title: string; icon: LucideIcon; items: NavigationItem[] }

const groups: NavigationGroup[] = [
  {
    id: 'workflow',
    title: 'Workflow',
    icon: FolderCog,
    items: [
      { id: 'controls', title: 'Controls', detail: 'Map, gate, build reference' },
      { id: 'samples', title: 'Samples', detail: 'Unmix sample files' },
    ],
  },
  {
    id: 'outputs',
    title: 'Results',
    icon: FileOutput,
    items: [
      { id: 'control-reports', title: 'Controls QC report', detail: 'Review control QC output' },
      { id: 'sample-reports', title: 'Samples QC report', detail: 'Review sample QC output' },
    ],
  },
  {
    id: 'tools',
    title: 'Other tools',
    icon: Wrench,
    items: [
      { id: 'panel', title: 'Panel builder', detail: 'Explore fluorophores' },
      { id: 'af', title: 'AF profiles', detail: 'Extract and manage AF' },
      { id: 'matrix', title: 'Adjust matrix', detail: 'Inspect residual crosstalk' },
    ],
  },
]

function groupForSection(section: SectionId) {
  return groups.find((group) => group.items.some((item) => item.id === section))?.id
}

function sectionCount(section: SectionId, project: ProjectState) {
  const counts: Partial<Record<SectionId, number>> = {
    controls: project.scan.controls,
    samples: project.scan.samples,
    'control-reports': project.artifacts.filter((artifact) => artifact.type === 'Control QC report' || artifact.type === 'QC report').length,
    'sample-reports': project.artifacts.filter((artifact) => artifact.type === 'Sample QC report').length,
    matrix: project.scan.matrices,
    af: project.artifacts.filter((artifact) => artifact.group === 'AF Profiles').length,
  }
  return counts[section]
}

export function WorkflowRail({ activeSection, project, showCounts, width, onWidthChange, onChange }: RailProps) {
  const activeGroup = groupForSection(activeSection)
  const [openGroups, setOpenGroups] = useState<Record<string, boolean>>({})

  function beginResize(event: ReactPointerEvent<HTMLDivElement>) {
    event.preventDefault()
    const startX = event.clientX
    const startWidth = width
    const move = (moveEvent: PointerEvent) => onWidthChange(Math.max(170, Math.min(340, startWidth + moveEvent.clientX - startX)))
    const stop = () => {
      document.removeEventListener('pointermove', move)
      document.removeEventListener('pointerup', stop)
    }
    document.addEventListener('pointermove', move)
    document.addEventListener('pointerup', stop)
  }

  return (
    <nav className="workflow-rail" aria-label="Main navigation" style={{ width }}>
      <div className="rail-intro">
        <span className="eyebrow">Project</span>
        <strong title={project.projectName}>{project.projectName}</strong>
      </div>
      <div className="rail-list">
        {groups.map((group) => {
          const containsActive = group.id === activeGroup
          const isOpen = Boolean(openGroups[group.id]) || containsActive
          const GroupIcon = group.icon
          return (
            <div className={`rail-group ${containsActive ? 'contains-active' : ''}`} key={group.id}>
              <button
                className="rail-group-toggle"
                type="button"
                onClick={() => setOpenGroups((current) => ({ ...current, [group.id]: !current[group.id] }))}
                aria-expanded={isOpen}
                aria-controls={`rail-group-${group.id}`}
              >
                <GroupIcon size={16} />
                <span>{group.title}</span>
                <ChevronRight className="rail-chevron" size={15} />
              </button>
              {isOpen && (
                <div className="rail-subsections" id={`rail-group-${group.id}`}>
                  {group.items.map((item) => {
                    const count = sectionCount(item.id, project)
                    return (
                      <button
                        key={item.id}
                        className={`rail-subitem ${activeSection === item.id ? 'is-active' : ''}`}
                        onClick={() => onChange(item.id)}
                        aria-current={activeSection === item.id ? 'page' : undefined}
                      >
                        <span className="rail-subitem-copy">
                          <strong>{item.title}</strong>
                          <small>{item.detail}</small>
                        </span>
                        {showCounts && count !== undefined && <span className="rail-count">{count}</span>}
                      </button>
                    )
                  })}
                </div>
              )}
            </div>
          )
        })}
      </div>
      <div className="workflow-rail-resizer" role="separator" aria-label="Resize navigation" aria-orientation="vertical" onPointerDown={beginResize} />
    </nav>
  )
}
