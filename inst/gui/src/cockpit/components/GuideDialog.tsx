import { useEffect, useRef, useState } from 'react'
import { BookOpen, ChevronRight, X } from 'lucide-react'
import { createPortal } from 'react-dom'

const sections = [
  {
    id: 'start', label: 'Getting started', title: 'Set up a project',
    content: <>
      <p><strong>Choose Project → Open project</strong> and select the folder that contains your experiment. Spectreasy uses that folder for every path shown in the cockpit.</p>
      <ol><li>Place single-colour controls in <strong>scc/</strong>.</li><li>Place experimental files in <strong>samples/</strong>.</li><li>Use <strong>Files</strong> to inspect, add, or remove FCS files without leaving the cockpit.</li></ol>
      <p className="guide-callout"><strong>Project selection is tab-specific.</strong> Each tab keeps its own selected project, while both tabs still use the same local R backend. Always check the project name before starting a long workflow.</p>
    </>,
  },
  {
    id: 'controls', label: 'Controls', title: 'Build the reference matrix',
    content: <>
      <p>The Controls workflow follows four ordered stages: <strong>Mapping, Gating, Unmixing, and Quality Control</strong>.</p>
      <ol><li><strong>Mapping:</strong> create or review <code>fcs_mapping.csv</code>. Confirming saves the table and opens gating.</li><li><strong>Gating:</strong> select the cell or bead population, singlets, positive events, and any negative population. <strong>Confirm</strong> saves <code>ssc_gate_config.csv</code> and advances to Unmixing.</li><li><strong>Unmixing:</strong> choose the method and report options. Advanced settings expose method-specific and output controls.</li><li><strong>Quality Control:</strong> open any report created in the selected output folder. Repeated runs remain available as separate numbered report folders.</li></ol>
    </>,
  },
  {
    id: 'gating', label: 'Gating', title: 'Review control gates',
    content: <>
      <p>Choose a control in the left sidebar, then work from scatter gates toward the histogram gate. The highlighted sidebar entry is the file currently shown.</p>
      <ul><li><strong>Neg / Pos</strong> selects the population being edited.</li><li><strong>Magic wand</strong> proposes gates; inspect every proposal before confirming.</li><li><strong>Settings</strong> changes point count, point size, histogram bins, transforms, and view ranges.</li><li><strong>Save / Load</strong> writes or opens a gate CSV. Confirm uses the project gate file for the cockpit workflow.</li></ul>
      <p className="guide-callout">Closing with <strong>Exit</strong> returns without advancing. Closing with <strong>Confirm</strong> saves the reviewed gates and opens Controls step 03.</p>
    </>,
  },
  {
    id: 'samples', label: 'Samples', title: 'Unmix sample files',
    content: <>
      <p>Samples has two stages: <strong>Unmixing</strong> and <strong>Quality Control</strong>. The selected project’s <code>samples/</code> files are listed below the settings card.</p>
      <ul><li>Select the reference matrix and, when applicable, detector-noise and spectral-variant files produced by Controls.</li><li><strong>Generate report</strong> enables report format and related options.</li><li><strong>One report per sample</strong> creates separate PDFs; HTML uses a sample selector in one interactive report.</li><li><strong>Write FCS outputs</strong> writes unmixed files to a numbered output folder so earlier runs are preserved.</li></ul>
    </>,
  },
  {
    id: 'reports', label: 'Reports', title: 'Inspect quality control',
    content: <>
      <p>Quality Control lists reports found only inside the current project’s selected output root. The folder prefix identifies the run, for example <code>qc_controls_2/qc_controls_report.html</code>.</p>
      <ul><li><strong>HTML</strong> reports are interactive and may include plot controls and sample selectors.</li><li><strong>PDF</strong> reports are fixed documents intended for review and sharing.</li><li>Use the card actions to open, download, or export an existing HTML report to PDF.</li><li>A <strong>stale</strong> label means an upstream input is newer than that report.</li></ul>
    </>,
  },
  {
    id: 'library', label: 'AF library', title: 'Manage autofluorescence profiles',
    content: <>
      <p>Extract an AF profile from an unstained FCS file, give it a reusable name, and save it to the local Spectreasy library.</p>
      <ul><li><strong>Extract</strong> creates the profile using the selected band count and event limit.</li><li><strong>Activate</strong> links a library profile to the current project and uses it as the unstained reference.</li><li><strong>Apply to matrix</strong> writes the AF rows into a selected matrix.</li><li>Deleting a profile removes it from the local library, not from reports already created.</li></ul>
    </>,
  },
  {
    id: 'tools', label: 'Matrix & panel', title: 'Use the analysis tools',
    content: <>
      <p><strong>Matrix Adjustment</strong> loads a project matrix and sample, previews residual relationships, and saves an adjusted copy. It never overwrites the selected source matrix.</p>
      <p><strong>Panel Builder</strong> compares fluorophores for the selected cytometer and configuration. Use similarity and spread information to inspect combinations before exporting the overview.</p>
      <p className="guide-callout">Both tools inherit the project selected in this browser tab. Their saved interface settings are stored with that project.</p>
    </>,
  },
  {
    id: 'settings', label: 'Settings & logs', title: 'Understand the cockpit controls',
    content: <>
      <ul><li><strong>Cytometer</strong> supplies detector metadata and panel configuration.</li><li><strong>Method</strong> sets the default unmixing method; each workflow card can show method-specific advanced settings.</li><li><strong>Logs</strong> shows the generated operation and grouped R output. Warnings and errors remain separate for visibility.</li><li><strong>Settings</strong> controls theme, typography, scale, sidebar width, and count badges.</li></ul>
      <p>No data is uploaded. Scientific work runs in the local R process started by <code>spectreasy_gui()</code>.</p>
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
