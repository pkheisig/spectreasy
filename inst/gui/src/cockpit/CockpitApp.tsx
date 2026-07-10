import { useCallback, useEffect, useMemo, useState } from 'react'
import { Bell, ChevronDown, CircleHelp, Command, Menu, PanelLeft, RefreshCcw, Settings2, Sparkles, Wifi, X } from 'lucide-react'
import { attemptWorkflowAction, downloadTextFile, initialBackendStatus, loadPanelPayload, loadProjectSnapshot, persistControlMapping, persistGuiState, setProjectContext } from './api'
import { demoProject } from './mockData'
import { Inspector } from './components/Inspector'
import { Sidebar } from './components/Sidebar'
import { WorkflowRail } from './components/WorkflowRail'
import { WorkflowWorkspace } from './workspaces/WorkflowWorkspace'
import { defaultWorkflowSettings } from './types'
import type { Artifact, BackendStatus, Job, MappingRow, PanelPayload, ProjectState, SectionId, WorkflowSettings } from './types'
import './cockpit.css'

const emptyJob: Job = { label: '', state: 'idle', progress: 0, subtask: '' }

function TopBar({ project, backend, onRefresh, onSettings }: { project: ProjectState; backend: BackendStatus; onRefresh: () => void; onSettings: () => void }) {
  return <header className="topbar"><div className="brand-lockup"><div className="brand-mark"><span /><span /><span /></div><div><strong>spectreasy</strong><small>spectral cockpit</small></div></div><div className="project-switcher"><div className="project-avatar">LP</div><div><span className="eyebrow">Active project</span><strong>{project.projectName}</strong></div><ChevronDown size={14} /></div><div className="topbar-spacer" /><div className="context-chip"><span className="chip-label">Cytometer</span><strong>{project.cytometer}</strong><ChevronDown size={13} /></div><div className="context-chip method-chip"><span className="chip-label">Method</span><strong>{project.method}</strong></div><div className="backend-chip"><span className={`backend-dot ${backend.connected ? 'is-connected' : ''}`} /><div><span className="chip-label">R backend</span><strong>{backend.connected ? 'Connected' : 'Preview mode'}</strong></div><Wifi size={13} /></div><button className="topbar-icon" onClick={onRefresh} aria-label="Refresh project"><RefreshCcw size={17} /></button><button className="topbar-icon has-alert" aria-label="Notifications"><Bell size={17} /><span /></button><button className="topbar-icon" onClick={onSettings} aria-label="Settings"><Settings2 size={17} /></button><div className="user-mark">PK</div></header>
}

function UtilityBar({ project, onToggleSidebar }: { project: ProjectState; onToggleSidebar: () => void }) {
  return <div className="utilitybar"><button className="utility-button" onClick={onToggleSidebar}><PanelLeft size={14} /> Project index</button><span className="utility-divider" /><span className="breadcrumb">{project.projectName} <strong>/</strong> cockpit</span><div className="utility-spacer" /><span className="last-saved"><span className="saved-dot" /> Last saved today at 09:48</span><button className="utility-button"><CircleHelp size={14} /> Help</button><button className="utility-button"><Command size={13} /> Shortcuts</button></div>
}

export default function CockpitApp() {
  const [project, setProject] = useState<ProjectState>(() => structuredClone(demoProject))
  const [backend, setBackend] = useState<BackendStatus>(initialBackendStatus)
  const [activeSection, setActiveSection] = useState<SectionId>('overview')
  const [mappingTab, setMappingTab] = useState<'mapping' | 'gating' | 'build' | 'qc'>('mapping')
  const [selectedArtifact, setSelectedArtifact] = useState<Artifact | null>(null)
  const [job, setJob] = useState<Job>(emptyJob)
  const [panelPayload, setPanelPayload] = useState<PanelPayload | null>(null)
  const [sidebarOpen, setSidebarOpen] = useState(true)
  const [toast, setToast] = useState<string | null>(null)
  const [settings, setSettings] = useState<WorkflowSettings>(() => defaultWorkflowSettings(demoProject.projectPath))

  const refreshProject = useCallback(async (initial = false) => {
    const snapshot = await loadProjectSnapshot()
    setProject(snapshot.project)
    setBackend(snapshot.backend)
    if (initial) {
      const saved = snapshot.savedSettings ?? {}
      setSettings((current) => ({
        ...current,
        ...saved,
        projectPath: snapshot.project.projectPath || current.projectPath,
        control: { ...current.control, ...(saved.control ?? {}), method: snapshot.project.method || current.control.method, cytometer: snapshot.project.cytometer || current.control.cytometer },
        sample: { ...current.sample, ...(saved.sample ?? {}), method: snapshot.project.method || current.sample.method },
        af: { ...current.af, ...(saved.af ?? {}) },
      }))
    }
    if (!initial) setToast(snapshot.backend.connected ? 'Project rescanned from the R backend.' : 'Preview state restored. Start the R backend for live artifacts.')
  }, [])

  useEffect(() => {
    const timer = window.setTimeout(() => void refreshProject(true), 0)
    return () => window.clearTimeout(timer)
  }, [refreshProject])

  useEffect(() => {
    if (job.state !== 'running') return
    const timer = window.setTimeout(() => {
      setJob((current) => {
        if (current.progress >= 100) return { ...current, state: 'complete', finishedAt: 'now', subtask: 'Outputs written to a new run folder', output: '9 artifacts' }
        const nextProgress = Math.min(100, current.progress + 20)
        const subtasks = ['Validating inputs and detector sets', 'Applying gates and building background', 'Unmixing events in R', 'Caching QC metrics and plots', 'Writing manifest and reports']
        return { ...current, progress: nextProgress, subtask: subtasks[Math.min(4, Math.floor(nextProgress / 20))] }
      })
    }, 650)
    return () => window.clearTimeout(timer)
  }, [job])

  useEffect(() => {
    if (!toast) return
    const timer = window.setTimeout(() => setToast(null), 4200)
    return () => window.clearTimeout(timer)
  }, [toast])

  function updateMapping(id: string, patch: Partial<MappingRow>) {
    setProject((current) => ({ ...current, mappingDirty: true, mapping: current.mapping.map((row) => row.id === id ? { ...row, ...patch } : row), artifacts: current.artifacts.map((artifact) => ['reference', 'unmixing', 'control-report'].includes(artifact.id) ? { ...artifact, status: 'stale' } : artifact) }))
  }

  function updateSettings(section: 'projectPath' | 'control' | 'sample' | 'af', patch: Partial<WorkflowSettings['control']> | Partial<WorkflowSettings['sample']> | Partial<WorkflowSettings['af']> | { projectPath: string }) {
    if (section === 'projectPath') {
      setSettings((current) => ({ ...current, projectPath: (patch as { projectPath: string }).projectPath }))
      return
    }
    setSettings((current) => ({ ...current, [section]: { ...current[section], ...patch } }))
  }

  async function runAction(action: 'control' | 'sample' | 'report' | 'af', label: string) {
    setJob({ label, state: 'running', progress: 0, subtask: 'Preparing inputs and checking prerequisites', startedAt: 'now' })
    const control = settings.control
    const sample = settings.sample
    const af = settings.af
    const payload = action === 'control'
      ? {
          projectPath: settings.projectPath,
          scc_dir: control.sccDir,
          control_file: control.controlFile,
          output_dir: control.outputDir,
          method: control.method,
          cytometer: control.cytometer,
          auto_create_mapping: control.autoCreateMapping,
          auto_unknown_fluor_policy: control.autoUnknownFluorPolicy,
          gate_file: control.gateFile,
          af_n_bands: control.afNBands,
          af_max_cells: control.afMaxCells,
          af_min_cluster_events: control.afMinClusterEvents,
          af_min_cluster_proportion: control.afMinClusterProportion,
          default_sample_type: control.defaultSampleType,
          histogram_pct_beads: control.histogramPctBeads,
          histogram_direction_beads: control.histogramDirectionBeads,
          histogram_pct_cells: control.histogramPctCells,
          histogram_direction_cells: control.histogramDirectionCells,
          outlier_percentile: control.outlierPercentile,
          debris_percentile: control.debrisPercentile,
          bead_gate_scale: control.beadGateScale,
          max_clusters: control.maxClusters,
          min_cluster_proportion: control.minClusterProportion,
          gate_contour_beads: control.gateContourBeads,
          gate_contour_cells: control.gateContourCells,
          subsample_n: control.subsampleN,
          unmix_scatter_panel_size_mm: control.unmixScatterPanelSizeMm,
          rwls_max_iter: control.rwlsMaxIter,
          unmix_threads: control.unmixThreads,
          seed: control.seed,
          save_qc_plots: control.saveQcPlots,
          save_report: control.saveReport,
          use_scatter_gating: control.useScatterGating,
          clean_scc_with_unstained: control.cleanSccWithUnstained,
          scc_background_method: control.sccBackgroundMethod,
          scc_background_k: control.sccBackgroundK,
          spectral_variant_som_nodes: control.spectralVariantSomNodes,
          spectral_variant_top_k: control.spectralVariantTopK,
          spectral_variant_cosine_threshold: control.spectralVariantCosineThreshold,
          spectral_variant_max_variants: control.spectralVariantMaxVariants,
          spectral_variant_min_events: control.spectralVariantMinEvents,
          spectreasy_weight_quantile: control.spectreasyWeightQuantile,
          autospectral_n_candidates: control.autospectralNCandidates,
          autospectral_n_spectral: control.autospectralNSpectral,
          autospectral_min_events: control.autospectralMinEvents,
          refine: control.refine,
        }
      : action === 'sample'
        ? {
            projectPath: settings.projectPath,
            sample_dir: sample.sampleDir,
            matrix_file: sample.matrixFile,
            detector_noise_file: sample.detectorNoiseFile,
            output_dir: sample.outputDir,
            method: sample.method,
            rwls_max_iter: sample.rwlsMaxIter,
            n_threads: sample.nThreads,
            spectral_variant_top_k: sample.spectralVariantTopK,
            spectral_variant_min_abundance: sample.spectralVariantMinAbundance,
            spectral_variant_positive_fraction: sample.spectralVariantPositiveFraction,
            spectral_variant_min_improvement: sample.spectralVariantMinImprovement,
            spectral_variant_library_file: sample.spectralVariantLibraryFile,
            spectreasy_weight_quantile: sample.spectreasyWeightQuantile,
            estimate_af: sample.estimateAf,
            write_fcs: sample.writeFcs,
            save_report: sample.saveReport,
            save_qc_plots: sample.saveQcPlots,
            plot_n_events: sample.plotNEvents,
            chunk_size: sample.chunkSize,
            seed: sample.seed,
            return_type: sample.returnType,
          }
        : action === 'af'
          ? {
              projectPath: settings.projectPath,
              fcs_file: af.fcsFile,
              save_name: af.saveName,
              save_overwrite: af.saveOverwrite,
              af_n_bands: af.afNBands,
              af_max_cells: af.afMaxCells,
              af_min_cluster_events: af.afMinClusterEvents,
              af_min_cluster_proportion: af.afMinClusterProportion,
              seed: af.seed,
            }
          : { projectPath: settings.projectPath, report_type: 'control' }
    const response = await attemptWorkflowAction(action, payload)
    setToast(response.message)
    if (response.connected) {
      window.setTimeout(() => void refreshProject(), 500)
    }
  }

  async function loadPanel() {
    const payload = await loadPanelPayload('aurora')
    if (payload) {
      setPanelPayload(payload)
      setToast('Panel library refreshed from the R backend.')
    } else {
      setToast('Showing the bundled panel library preview. The R backend is not connected.')
    }
  }

  async function saveProject() {
    const mapping = await persistControlMapping(project.mapping)
    const saved = await persistGuiState(project, settings)
    if (mapping.success) setProject((current) => ({ ...current, mappingDirty: false }))
    setToast(saved && mapping.success ? 'Project preferences and control mapping saved to the local R config.' : 'Project preferences saved for this browser session. Connect R to persist them.')
  }

  async function openProject(path: string) {
    const response = await setProjectContext(path)
    setToast(response.message)
    if (response.success) {
      setSettings((current) => ({ ...current, projectPath: path }))
      await refreshProject()
    }
  }

  function downloadArtifact(artifact: Artifact) {
    const manifest = `Spectreasy artifact\n\nName: ${artifact.name}\nType: ${artifact.type}\nStatus: ${artifact.status}\nPath: ${artifact.path}\nUpdated: ${artifact.updated}\nRun: ${artifact.run ?? 'user supplied'}\n`
    downloadTextFile(`${artifact.name.replaceAll('/', '_')}.txt`, manifest)
    setToast(`Prepared a local download for ${artifact.name}.`)
  }

  function selectArtifact(artifact: Artifact) {
    setSelectedArtifact(artifact)
    const sectionMap: Record<string, SectionId> = { Controls: 'controls', Gates: 'controls', Samples: 'samples', Matrices: 'matrix', Reports: 'reports', 'QC Metrics': 'reports', 'AF Profiles': 'af', 'Panel Builder': 'panel', Logs: 'settings' }
    const section = sectionMap[artifact.group]
    if (section) setActiveSection(section)
  }

  const activeTitle = useMemo(() => ({ overview: 'Project setup', controls: 'Controls', samples: 'Samples', reports: 'Reports', matrix: 'Matrix review', panel: 'Panel builder', af: 'AF library', comparison: 'Method comparison', simulator: 'Synthetic SCC', settings: 'Settings & logs' })[activeSection], [activeSection])

  return <div className={`cockpit-app ${sidebarOpen ? '' : 'sidebar-collapsed'} ${activeSection === 'controls' && mappingTab === 'gating' ? 'is-gating-editor' : ''}`}><TopBar project={project} backend={backend} onRefresh={() => void refreshProject()} onSettings={() => setActiveSection('settings')} /><UtilityBar project={project} onToggleSidebar={() => setSidebarOpen((open) => !open)} /><div className="app-body">{sidebarOpen && <Sidebar project={project} activeSection={activeSection} onSectionChange={setActiveSection} selectedArtifact={selectedArtifact} onSelectArtifact={selectArtifact} />}<main className="main-area"><WorkflowRail activeSection={activeSection} onChange={setActiveSection} /><div className="main-canvas"><div className="mobile-canvas-title"><Menu size={15} /><span>{activeTitle}</span></div><WorkflowWorkspace project={project} backend={backend} job={job} activeSection={activeSection} mappingTab={mappingTab} setMappingTab={setMappingTab} onUpdateMapping={updateMapping} onRun={runAction} onRefresh={() => void refreshProject()} onDownload={downloadArtifact} onSelectArtifact={selectArtifact} onLoadPanel={() => void loadPanel()} panelPayload={panelPayload} onSave={() => void saveProject()} onSectionChange={setActiveSection} settings={settings} onSettingsChange={updateSettings} onOpenProject={openProject} /></div></main><Inspector artifact={selectedArtifact} onDownload={downloadArtifact} /></div>{toast && <div className="toast"><Sparkles size={15} /><span>{toast}</span><button onClick={() => setToast(null)} aria-label="Dismiss notification"><X size={14} /></button></div>}</div>
}
