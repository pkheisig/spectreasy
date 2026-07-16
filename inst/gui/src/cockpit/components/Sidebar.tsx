import { ChevronDown, ChevronRight, CircleDot, Folder, Gauge, Layers3, PanelTop, Search, Table2, TerminalSquare } from 'lucide-react'
import type { LucideIcon } from 'lucide-react'
import { useMemo, useState } from 'react'
import type { Artifact, ProjectState, SectionId } from '../types'

type SidebarProps = {
  project: ProjectState
  activeSection: SectionId
  onSectionChange: (section: SectionId) => void
  selectedArtifact: Artifact | null
  onSelectArtifact: (artifact: Artifact) => void
}

const groupMeta: Array<{ name: string; icon: LucideIcon; section?: SectionId }> = [
  { name: 'Project', icon: Folder },
  { name: 'Controls', icon: CircleDot, section: 'controls' },
  { name: 'Samples', icon: Layers3, section: 'samples' },
  { name: 'Gates', icon: Gauge, section: 'controls' },
  { name: 'Matrices', icon: Table2, section: 'matrix' },
  { name: 'AF Profiles', icon: CircleDot, section: 'af' },
  { name: 'Panel Builder', icon: PanelTop, section: 'panel' },
  { name: 'Logs', icon: TerminalSquare, section: 'settings' },
]

function groupFiles(artifacts: Artifact[], group: string) {
  return artifacts.filter((artifact) => artifact.group === group)
}

export function Sidebar({ project, activeSection, onSectionChange, selectedArtifact, onSelectArtifact }: SidebarProps) {
  const [query, setQuery] = useState('')
  const [openGroups, setOpenGroups] = useState<Record<string, boolean>>({ Project: true, Controls: true, Matrices: true, Reports: true })
  const filteredArtifacts = useMemo(() => {
    const needle = query.trim().toLowerCase()
    return needle ? project.artifacts.filter((artifact) => `${artifact.name} ${artifact.detail} ${artifact.group}`.toLowerCase().includes(needle)) : project.artifacts
  }, [project.artifacts, query])

  const toggleGroup = (group: string) => setOpenGroups((current) => ({ ...current, [group]: !current[group] }))

  return (
    <aside className="sidebar">
      <div className="sidebar-heading">
        <div>
          <span className="eyebrow">Project files</span>
          <h2>Artifact index</h2>
        </div>
        <span className="sidebar-count">{project.artifacts.length}</span>
      </div>
      <label className="search-field">
        <Search size={15} />
        <input value={query} onChange={(event) => setQuery(event.target.value)} placeholder="Find a file or output" />
        <kbd>⌘ K</kbd>
      </label>
      <div className="sidebar-groups">
        {groupMeta.map(({ name, icon: Icon, section }) => {
          const files = groupFiles(filteredArtifacts, name)
          const isOpen = Boolean(openGroups[name])
          const isActive = section === activeSection && name !== 'Project'
          const count = project.artifacts.filter((artifact) => artifact.group === name).length
          return (
            <div className={`file-group ${isOpen ? 'is-open' : ''} ${isActive ? 'is-active' : ''}`} key={name}>
              <button className="file-group-header" onClick={() => { toggleGroup(name); if (section) onSectionChange(section) }} aria-expanded={isOpen}>
                <span className="file-group-title"><Icon size={15} /><span>{name}</span></span>
                {count > 0 && <span className="file-group-count">{count}</span>}
                {isOpen ? <ChevronDown size={14} /> : <ChevronRight size={14} />}
              </button>
              {isOpen && files.length > 0 && (
                <div className="file-group-items">
                  {files.slice(0, 5).map((artifact) => (
                    <button className={`artifact-row ${selectedArtifact?.id === artifact.id ? 'is-selected' : ''}`} key={artifact.id} onClick={() => onSelectArtifact(artifact)}>
                      <span className={`artifact-dot dot-${artifact.status}`} />
                      <span className="artifact-copy"><span className="artifact-name">{artifact.name}</span><span className="artifact-detail">{artifact.detail}</span></span>
                      {artifact.status === 'stale' && <span className="artifact-mark">!</span>}
                    </button>
                  ))}
                  {files.length > 5 && <button className="show-more">+ {files.length - 5} more artifacts</button>}
                </div>
              )}
              {isOpen && files.length === 0 && <div className="empty-group">No matching artifacts</div>}
            </div>
          )
        })}
      </div>
    </aside>
  )
}
