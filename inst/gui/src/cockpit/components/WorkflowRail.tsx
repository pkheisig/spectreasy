import {
  ChevronRight,
  FolderCog,
  Wrench,
} from 'lucide-react'
import { useEffect, useState } from 'react'
import type { PointerEvent as ReactPointerEvent } from 'react'
import type { LucideIcon } from 'lucide-react'
import type { ProjectState, SectionId } from '../types'

type RailProps = {
  activeSection: SectionId
  project: ProjectState
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
    id: 'tools',
    title: 'Other tools',
    icon: Wrench,
    items: [
      { id: 'af', title: 'AF profiles', detail: 'Extract and manage AF' },
      { id: 'panel', title: 'Panel builder', detail: 'Explore fluorophores' },
      { id: 'matrix', title: 'Adjust matrix', detail: 'Inspect residual crosstalk' },
    ],
  },
]
const NAVIGATION_GROUPS_STORAGE_KEY = 'spectreasy-navigation-groups'

function loadOpenGroups(): Record<string, boolean> {
  try {
    const stored = window.localStorage.getItem(NAVIGATION_GROUPS_STORAGE_KEY)
    if (!stored) return {}
    const parsed = JSON.parse(stored)
    if (!parsed || typeof parsed !== 'object') return {}
    return Object.fromEntries(
      Object.entries(parsed)
        .filter(([key, value]) => (key === 'workflow' || key === 'tools') && typeof value === 'boolean'),
    ) as Record<string, boolean>
  } catch {
    return {}
  }
}

function groupForSection(section: SectionId) {
  return groups.find((group) => group.items.some((item) => item.id === section))?.id
}

export function WorkflowRail({ activeSection, project, width, onWidthChange, onChange }: RailProps) {
  const activeGroup = groupForSection(activeSection)
  const [openGroups, setOpenGroups] = useState<Record<string, boolean>>(loadOpenGroups)

  useEffect(() => {
    window.localStorage.setItem(NAVIGATION_GROUPS_STORAGE_KEY, JSON.stringify(openGroups))
  }, [openGroups])

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
          const isOpen = openGroups[group.id] ?? containsActive
          const GroupIcon = group.icon
          return (
            <div className={`rail-group ${containsActive ? 'contains-active' : ''}`} key={group.id}>
              <button
                className="rail-group-toggle"
                type="button"
                onClick={() => setOpenGroups((current) => ({
                  ...current,
                  [group.id]: !(current[group.id] ?? containsActive),
                }))}
                aria-expanded={isOpen}
                aria-controls={`rail-group-${group.id}`}
              >
                <GroupIcon size={16} />
                <span>{group.title}</span>
                <ChevronRight className="rail-chevron" size={15} />
              </button>
                <div className={`rail-subsections ${isOpen ? 'is-open' : ''}`} id={`rail-group-${group.id}`} aria-hidden={!isOpen}>
                  {group.items.map((item) => {
                    return (
                      <button
                        key={item.id}
                        className={`rail-subitem ${activeSection === item.id ? 'is-active' : ''}`}
                        onClick={() => onChange(item.id)}
                        tabIndex={isOpen ? 0 : -1}
                        aria-current={activeSection === item.id ? 'page' : undefined}
                      >
                        <span className="rail-subitem-copy">
                          <strong>{item.title}</strong>
                          <small>{item.detail}</small>
                        </span>
                      </button>
                    )
                  })}
                </div>
            </div>
          )
        })}
      </div>
      <p className="rail-privacy">No data is uploaded. Everything stays on your computer.</p>
      <div className="workflow-rail-resizer" role="separator" aria-label="Resize navigation" aria-orientation="vertical" onPointerDown={beginResize} />
    </nav>
  )
}
