import {
  ChevronRight,
  CircleGauge,
  FileOutput,
  FolderCog,
  SlidersHorizontal,
  Wrench,
} from 'lucide-react'
import { useState } from 'react'
import type { LucideIcon } from 'lucide-react'
import type { ProjectState, SectionId } from '../types'

type RailProps = {
  activeSection: SectionId
  project: ProjectState
  showCounts: boolean
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
    title: 'Results and matrices',
    icon: FileOutput,
    items: [
      { id: 'reports', title: 'Reports', detail: 'Review QC output' },
      { id: 'matrix', title: 'Adjust matrix', detail: 'Inspect residual crosstalk' },
    ],
  },
  {
    id: 'tools',
    title: 'Analysis tools',
    icon: Wrench,
    items: [
      { id: 'panel', title: 'Panel builder', detail: 'Explore fluorophores' },
      { id: 'af', title: 'AF profiles', detail: 'Extract and manage AF' },
      { id: 'comparison', title: 'Compare methods', detail: 'Run diagnostics' },
      { id: 'simulator', title: 'Synthetic SCC', detail: 'Generate test controls' },
    ],
  },
  {
    id: 'system',
    title: 'System',
    icon: SlidersHorizontal,
    items: [
      { id: 'settings', title: 'Settings and logs', detail: 'Workflow and appearance' },
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
    reports: project.scan.reports,
    matrix: project.scan.matrices,
    af: project.artifacts.filter((artifact) => artifact.group === 'AF Profiles').length,
  }
  return counts[section]
}

export function WorkflowRail({ activeSection, project, showCounts, onChange }: RailProps) {
  const activeGroup = groupForSection(activeSection)
  const [openGroups, setOpenGroups] = useState<Record<string, boolean>>({})

  return (
    <nav className="workflow-rail" aria-label="Main navigation">
      <div className="rail-intro">
        <span className="eyebrow">Project</span>
        <strong title={project.projectName}>{project.projectName}</strong>
      </div>
      <div className="rail-list">
        <button
          className={`rail-overview ${activeSection === 'overview' ? 'is-active' : ''}`}
          onClick={() => onChange('overview')}
          aria-current={activeSection === 'overview' ? 'page' : undefined}
        >
          <CircleGauge size={17} />
          <span>Overview</span>
        </button>
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
    </nav>
  )
}
