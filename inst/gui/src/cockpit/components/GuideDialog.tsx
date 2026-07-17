import { useEffect, useRef, useState } from 'react'
import { BookOpen, ChevronRight, X } from 'lucide-react'
import { createPortal } from 'react-dom'

const sections = [
  {
    id: 'start', label: 'Getting started', title: 'Set up a project',
    content: <>
      <p>Click on <strong>Project</strong> to open a project directory that contains your experiment or create a new one.</p>
      <ol><li>Place single-colour controls in the configured <strong>Control FCS folder</strong>.</li><li>Place experimental files in the configured <strong>Sample FCS folder</strong>.</li></ol>
      <p>Use <strong>Files</strong> to inspect, add, or remove FCS files.</p>
      <p className="guide-callout"><strong>Projects are instance-specific</strong>, meaning each instance of the Spectreasy GUI (i.e., browser tab) keeps its own selected project. All instances still use the same local R backend. Always check the project name before starting a workflow.</p>
      <p>For code and documentation, visit the <a href="https://github.com/pkheisig/spectreasy" target="_blank" rel="noopener noreferrer"><strong><u>Github repository</u></strong> </a>.</p>
    </>,
  },
  {
    id: 'controls', label: 'Controls', title: 'Unmix control files',
    content: <>
      <p>The Controls workflow follows four ordered stages: <strong>Mapping, Gating, Unmixing, and Quality Control</strong>.</p>
      <ol><li><strong>Mapping:</strong> create or review <code>fcs_mapping.csv</code>. This table tells Spectrasy how to process each control file in the backend, like assigning fluorophore names and background controls. Confirming saves the table and opens gating.</li><li><strong>Gating:</strong> select the cell or bead population, singlets, positive events, and any negative population (when no background control is present in the mapping table). Histogram gates can be automatically generated after singlet gating and then manually adjusted. Gates are applied across all bead- or cell-based controls but can be made file-specific. <strong>Confirm</strong> saves <code>ssc_gate_config.csv</code> and advances to Unmixing.</li><li><strong>Unmixing:</strong> choose the unmixing method (algorithm) and report options. Advanced settings expose method-specific parameters and output controls.</li><li><strong>Quality Control:</strong> open and inspect SCC unmixing quality control reports. Repeated runs remain available as separate entries.</li></ol>
    </>,
  },
  {
    id: 'samples', label: 'Samples', title: 'Unmix sample files',
    content: <>
      <p>Samples has two stages: <strong>Unmixing</strong> and <strong>Quality Control</strong>.</p>
      <ol><li><strong>Unmixing:</strong> choose the unmixing method (algorithm) and report options. Advanced settings expose method-specific parameters and output controls.It is recommended to use the same method as in the Controls workflow.</li><li><strong>Quality Control:</strong> open and inspect sample unmixing quality control reports. Repeated runs remain available as separate entries.</li></ol>
    </>,
  },
  {
    id: 'library', label: 'AF library', title: 'Manage autofluorescence profiles',
    content: <>
      <p>Extract an AF profile from an unstained FCS file, give it a reusable name, and save it to the local Spectreasy library.</p>
      <ul><li><strong>Extract & save</strong> creates the profile using the selected band count and event limit. The profile is saved to the package root directory.</li><li><strong>Link to dataset</strong> links a library profile to the current project and uses it as the unstained reference.</li><li> Modify profiles further using the <strong>Rename</strong> and <strong>Delete</strong> options.</li></ul>
    </>,
  },
  {
    id: 'matrix', label: 'Adjust matrix', title: 'Adjust reference and unmixing matrices',
    content: <>
      <p><ol><li><strong>Load</strong> a matrix file and apply it to a sample with the selected unmixing method.</li><li><strong>Preview</strong> residual relationships and adjust population shapes by dragging on the events. A new modified copy will be saved automatically.</li></ol></p>

    </>,
  },
  {
    id: 'panel', label: 'Panel builder', title: 'Design a spectral panel',
    content: <>
      <p><strong>Panel Builder</strong> compares fluorophore spectra for the selected cytometer and configuration.</p>
      <ol><li><strong>Design</strong> a spectral panel from scratch before starting the experiment</li><li><strong>Optimize</strong> for spectral overlap and panel complexity.</li><li><strong>Validate</strong> experimental data from the QC reports against theoretical indices.</li></ol>

    </>,
  },
  {
    id: 'settings', label: 'Settings & logs', title: 'Settings and session logs',
    content: <>
      <ul><li><strong>Cytometer</strong> supplies detector metadata and panel configuration.</li><li><strong>Method</strong> sets the default unmixing method; each workflow card can show method-specific advanced settings.</li><li><strong>Logs</strong> shows the generated operation and grouped R output.</li><li><strong>Settings</strong> page controls global default values, app aesthetics, and more.</li></ul>
      <p>No data is uploaded to the internet. All processes run in the local R environment started by <code>spectreasy::spectreasy_gui()</code>.</p>
    </>,
  },
] as const

type Props = { onClose: () => void }

export function GuideDialog({ onClose }: Props) {
  const [activeId, setActiveId] = useState<(typeof sections)[number]['id']>('start')
  const dialogRef = useRef<HTMLElement>(null)
  const active = sections.find((section) => section.id === activeId) ?? sections[0]

  useEffect(() => {
    dialogRef.current?.querySelector<HTMLButtonElement>('[role="tab"]')?.focus()
    const keydown = (event: KeyboardEvent) => {
      if (event.key === 'Escape') { event.preventDefault(); onClose(); return }
      if (event.key !== 'Tab') return
      const focusable = Array.from(dialogRef.current?.querySelectorAll<HTMLElement>('button:not([disabled]), [href], [tabindex]:not([tabindex="-1"])') ?? [])
      if (!focusable.length) return
      const first = focusable[0]; const last = focusable[focusable.length - 1]
      if (event.shiftKey && document.activeElement === first) { event.preventDefault(); last.focus() }
      else if (!event.shiftKey && document.activeElement === last) { event.preventDefault(); first.focus() }
    }
    document.addEventListener('keydown', keydown)
    return () => document.removeEventListener('keydown', keydown)
  }, [onClose])

  return createPortal(
    <div className="guide-overlay" role="presentation" onMouseDown={onClose}>
      <section ref={dialogRef} className="guide-dialog" role="dialog" aria-modal="true" aria-labelledby="guide-title" onMouseDown={(event) => event.stopPropagation()}>
        <header className="guide-header">
          <div><BookOpen size={20} /><div><span>Spectreasy cockpit</span><h2 id="guide-title">Guide</h2></div></div>
          <button type="button" onClick={onClose} aria-label="Close guide"><X size={19} /></button>
        </header>
        <div className="guide-layout">
          <nav className="guide-tabs" role="tablist" aria-label="Guide sections">
            {sections.map((section) => <button key={section.id} type="button" role="tab" aria-selected={section.id === activeId} tabIndex={section.id === activeId ? 0 : -1} className={section.id === activeId ? 'is-active' : ''} onClick={() => setActiveId(section.id)}><span>{section.label}</span><ChevronRight size={15} /></button>)}
          </nav>
          <article className="guide-content" role="tabpanel" tabIndex={0}>
            <span className="guide-section-number">{String(sections.findIndex((section) => section.id === active.id) + 1).padStart(2, '0')}</span>
            <h3>{active.title}</h3>
            {active.content}
          </article>
        </div>
      </section>
    </div>, document.body,
  )
}
