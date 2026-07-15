import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  Menu,
  Moon,
  FolderOpen,
  Files as FilesIcon,
  Settings,
  Sparkles,
  Sun,
  TerminalSquare,
  X,
} from "lucide-react";
import {
  attemptWorkflowAction,
  createControlMapping,
  createProjectFolder,
  downloadTextFile,
  initialBackendStatus,
  loadPanelPayload,
  loadProjectSnapshot,
  persistControlMapping,
  persistGuiState,
  selectProjectFolder,
  terminateRSession,
} from "./api";
import { emptyProject } from "./mockData";
import { WorkflowRail } from "./components/WorkflowRail";
import { GuiSelect } from "./components/GuiSelect";
import { CockpitApplet } from "./components/CockpitApplet";
import { TerminalPanel } from "./components/TerminalPanel";
import { ProjectFilesDialog } from "./components/ProjectFilesDialog";
import { WorkflowWorkspace } from "./workspaces/WorkflowWorkspace";
import { defaultWorkflowSettings, normalizeInterfaceScale } from "./types";
import type {
  Artifact,
  BackendStatus,
  CockpitAppletId,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  SectionId,
  WorkflowSettings,
} from "./types";
import "./cockpit.css";

const emptyJob: Job = { label: "", state: "idle", progress: 0, subtask: "" };

const cytometerOptions = [
  ["auto", "Auto"],
  ["aurora", "Cytek Aurora"],
  ["northern_lights", "Cytek Northern Lights"],
  ["id7000", "Sony ID7000"],
  ["discover_s8", "BD FACSDiscover S8"],
  ["discover_a8", "BD FACSDiscover A8"],
  ["a5se", "BD FACSymphony A5 SE"],
  ["opteon", "Agilent NovoCyte Opteon"],
  ["mosaic", "Beckman Coulter CytoFLEX Mosaic"],
  ["xenith", "Thermo Fisher Attune Xenith"],
] as const;

const cytometerLabels = Object.fromEntries(cytometerOptions) as Record<string, string>;

function normalizeCockpitCytometer(value: unknown): string {
  const token = String(value ?? "auto")
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_|_$/g, "");
  const legacyAliases: Record<string, string> = {
    aurora_5l: "aurora",
    aurora_4l: "aurora",
    cytek_aurora: "aurora",
    discover: "auto",
    cytek_aurora_discover: "auto",
    bd_facsdiscover_xenith: "xenith",
    thermo_fisher_attune_xenith: "xenith",
  };
  const normalized = legacyAliases[token] ?? token;
  return normalized in cytometerLabels ? normalized : "auto";
}

type TopBarProps = {
  project: ProjectState;
  cytometer: string;
  method: string;
  darkMode: boolean;
  settingsActive: boolean;
  terminalActive: boolean;
  backendConnected: boolean;
  onCytometerChange: (value: string) => void;
  onMethodChange: (value: string) => void;
  onToggleTheme: () => void;
  onFiles: () => void;
  onCreateProject: () => void;
  onOpenProject: () => void;
  onSettings: () => void;
  onTerminal: () => void;
};

function TopBar({
  project,
  cytometer,
  method,
  darkMode,
  settingsActive,
  terminalActive,
  backendConnected,
  onCytometerChange,
  onMethodChange,
  onToggleTheme,
  onFiles,
  onCreateProject,
  onOpenProject,
  onSettings,
  onTerminal,
}: TopBarProps) {
  const [projectMenuOpen, setProjectMenuOpen] = useState(false);
  const projectMenuRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!projectMenuOpen) return;
    const closeMenu = (event: MouseEvent | KeyboardEvent) => {
      if (event instanceof KeyboardEvent && event.key !== "Escape") return;
      if (event instanceof MouseEvent && projectMenuRef.current?.contains(event.target as Node)) return;
      setProjectMenuOpen(false);
    };
    document.addEventListener("mousedown", closeMenu);
    document.addEventListener("keydown", closeMenu);
    return () => {
      document.removeEventListener("mousedown", closeMenu);
      document.removeEventListener("keydown", closeMenu);
    };
  }, [projectMenuOpen]);

  const chooseProjectAction = (action: () => void) => {
    setProjectMenuOpen(false);
    action();
  };

  return (
    <header className="topbar">
      <div className="brand-lockup">
        <div className="brand-mark">
          <span />
          <span />
          <span />
        </div>
        <div>
          <strong>spectreasy</strong>
        </div>
      </div>
      <div className="project-actions">
        <div className="project-menu-shell" ref={projectMenuRef}>
          <button
            className="project-switcher"
            type="button"
            aria-expanded={projectMenuOpen}
            aria-haspopup="menu"
            onClick={() => setProjectMenuOpen((open) => !open)}
          >
            <FolderOpen size={16} />
            <strong>Project</strong>
          </button>
          <div
            className={`project-action-flyout ${projectMenuOpen ? "is-open" : ""}`}
            role="menu"
            aria-hidden={!projectMenuOpen}
          >
            <button type="button" role="menuitem" tabIndex={projectMenuOpen ? 0 : -1} onClick={() => chooseProjectAction(onCreateProject)}>
              Create new project
            </button>
            <button type="button" role="menuitem" tabIndex={projectMenuOpen ? 0 : -1} onClick={() => chooseProjectAction(onOpenProject)}>
              Open project
            </button>
          </div>
        </div>
        <button className="project-files-button" type="button" onClick={onFiles} disabled={!project.projectPath}>
          <FilesIcon size={15} />
          <span>Files</span>
        </button>
      </div>
      <div className="topbar-spacer" />
      <label className="context-select cytometer-select">
        <span className="chip-label">Cytometer</span>
        <GuiSelect
          aria-label="Cytometer"
          className="cytometer-picker"
          value={cytometer}
          onChange={(event) => onCytometerChange(event.target.value)}
        >
          {cytometerOptions.map(([value, label]) => (
            <option value={value} key={value}>{label}</option>
          ))}
        </GuiSelect>
      </label>
      <label className="context-select method-select">
        <span className="chip-label">Method</span>
        <GuiSelect
          aria-label="Method"
          value={method}
          onChange={(event) => onMethodChange(event.target.value)}
        >
          <option>Spectreasy</option>
          <option>AutoSpectral</option>
          <option>OLS</option>
          <option>WLS</option>
          <option>RWLS</option>
          <option>NNLS</option>
        </GuiSelect>
      </label>
      <button
        className="topbar-icon"
        onClick={onToggleTheme}
        aria-label={darkMode ? "Use light mode" : "Use dark mode"}
        title={darkMode ? "Use light mode" : "Use dark mode"}
      >
        {darkMode ? <Sun size={17} /> : <Moon size={17} />}
      </button>
      <button
        className={`topbar-icon terminal-trigger ${terminalActive ? "is-active" : ""}`}
        onClick={onTerminal}
        aria-label={terminalActive ? "Close terminal" : "Open terminal"}
        aria-pressed={terminalActive}
        title={backendConnected ? "Open terminal" : "Terminal · R backend offline"}
      >
        <TerminalSquare size={18} />
        <span className={`terminal-status-dot ${backendConnected ? "is-connected" : ""}`} />
      </button>
      <button
        className={`topbar-icon ${settingsActive ? "is-active" : ""}`}
        onClick={onSettings}
        aria-label={settingsActive ? "Return to workflow" : "Settings"}
        aria-pressed={settingsActive}
        title={settingsActive ? "Return to workflow" : "Settings"}
      >
        <Settings size={18} />
      </button>
    </header>
  );
}

export default function CockpitApp() {
  const [project, setProject] = useState<ProjectState>(() =>
    structuredClone(emptyProject),
  );
  const [backend, setBackend] = useState<BackendStatus>(initialBackendStatus);
  const [activeSection, setActiveSection] = useState<SectionId>("controls");
  const [activeApplet, setActiveApplet] = useState<CockpitAppletId | null>(null);
  const [sectionBeforeApplet, setSectionBeforeApplet] =
    useState<SectionId>("controls");
  const [mappingTab, setMappingTab] = useState<
    "mapping" | "gating" | "build" | "qc"
  >("mapping");
  const [job, setJob] = useState<Job>(emptyJob);
  const [panelPayload, setPanelPayload] = useState<PanelPayload | null>(null);
  const [toast, setToast] = useState<string | null>(null);
  const [terminalOpen, setTerminalOpen] = useState(false);
  const [filesOpen, setFilesOpen] = useState(false);
  const [settings, setSettings] = useState<WorkflowSettings>(() => {
    const initial = defaultWorkflowSettings('');
    const savedTheme = window.localStorage.getItem("spectreasy-theme");
    if (savedTheme === "light" || savedTheme === "dark") {
      initial.appearance.theme = savedTheme;
    }
    return initial;
  });
  const [settingsReady, setSettingsReady] = useState(false);
  const [sectionBeforeSettings, setSectionBeforeSettings] =
    useState<SectionId>("controls");
  const darkMode = settings.appearance.theme === "dark";

  const refreshProject = useCallback(async (initial = false) => {
    const snapshot = await loadProjectSnapshot();
    setProject(snapshot.project);
    setBackend(snapshot.backend);
    if (!initial && snapshot.project.projectPath) {
      setSettings((current) => ({
        ...current,
        projectPath: snapshot.project.projectPath,
      }));
    }
    if (initial) {
      const saved = snapshot.savedSettings ?? {};
      setSettings((current) => {
        const savedControl: Partial<WorkflowSettings["control"]> =
          saved.control ?? {};
        const savedSample: Partial<WorkflowSettings["sample"]> =
          saved.sample ?? {};
        const explicitCytometer = window.localStorage.getItem(
          "spectreasy-cytometer",
        );
        const afFiles = snapshot.project.mapping
          .filter((row) => row.marker.trim().toLowerCase() === "autofluorescence" || /^af(?:$|_|\b)/i.test(row.fluorophore.trim()))
          .map((row) => `scc/${row.file}`);
        const savedAf: Partial<WorkflowSettings["af"]> = saved.af ?? {};
        const requestedAfFile = savedAf.fcsFile ?? current.af.fcsFile;
        return {
          ...current,
          ...saved,
          projectPath: snapshot.project.projectPath || current.projectPath,
          control: {
            ...current.control,
            ...savedControl,
            sccDir: "scc",
            controlFile: "fcs_mapping.csv",
            outputDir: "spectreasy_outputs/unmix_controls",
            gateFile: "ssc_gate_config.csv",
            method:
              savedControl.method ??
              snapshot.project.method ??
              current.control.method,
            cytometer: normalizeCockpitCytometer(
              explicitCytometer ?? savedControl.cytometer ?? snapshot.project.cytometer,
            ),
          },
          sample: {
            ...current.sample,
            ...savedSample,
            sampleDir: "samples",
            matrixFile: "spectreasy_outputs/unmix_controls/scc_reference_matrix.csv",
            detectorNoiseFile: "spectreasy_outputs/unmix_controls/scc_detector_noise.csv",
            outputDir: "spectreasy_outputs/unmix_samples/unmixed_fcs",
            method:
              savedSample.method ??
              savedControl.method ??
              snapshot.project.method ??
              current.sample.method,
          },
          af: {
            ...current.af,
            ...savedAf,
            fcsFile: afFiles.includes(requestedAfFile) ? requestedAfFile : (afFiles[0] ?? ""),
          },
          appearance: {
            ...current.appearance,
            ...(saved.appearance ?? {}),
            fontFamily:
              String(saved.appearance?.fontFamily ?? "") === "source-sans"
                ? "atkinson"
                : (saved.appearance?.fontFamily ?? current.appearance.fontFamily),
            fontScale: normalizeInterfaceScale(saved.appearance?.fontScale ?? current.appearance.fontScale),
            sidebarWidth:
              saved.appearance?.sidebarWidth === 220
                ? 242
                : (saved.appearance?.sidebarWidth ?? current.appearance.sidebarWidth),
          },
        };
      });
    }
  }, []);

  useEffect(() => {
    let cancelled = false;
    const timer = window.setTimeout(() => {
      void refreshProject(true).finally(() => {
        if (!cancelled) setSettingsReady(true);
      });
    }, 0);
    return () => {
      cancelled = true;
      window.clearTimeout(timer);
    };
  }, [refreshProject]);

  useEffect(() => {
    if (backend.connected) return;
    const timer = window.setInterval(() => {
      void refreshProject(false);
    }, 1500);
    return () => window.clearInterval(timer);
  }, [backend.connected, refreshProject]);

  useEffect(() => {
    if (!settingsReady) return;
    const timer = window.setTimeout(() => void persistGuiState(project, settings), 450);
    return () => window.clearTimeout(timer);
  }, [project, settings, settingsReady]);

  useEffect(() => {
    const exitOverlay = (event: KeyboardEvent) => {
      if (event.key !== "Escape") return;
      if (document.querySelector(".gui-select-menu")) return;
      if (document.querySelector(".project-action-flyout.is-open")) return;
      if (filesOpen) {
        event.preventDefault();
        setFilesOpen(false);
        return;
      }
      if (activeSection === "settings") {
        event.preventDefault();
        setActiveSection(sectionBeforeSettings);
      }
    };
    document.addEventListener("keydown", exitOverlay);
    return () => document.removeEventListener("keydown", exitOverlay);
  }, [activeSection, filesOpen, sectionBeforeSettings]);

  useEffect(() => {
    if (job.state !== "running") return;
    const timer = window.setTimeout(() => {
      setJob((current) => {
        if (current.progress >= 100)
          return {
            ...current,
            state: "complete",
            finishedAt: "now",
            subtask: "Outputs written to a new run folder",
            output: "9 artifacts",
          };
        const nextProgress = Math.min(100, current.progress + 20);
        const subtasks = [
          "Validating inputs and detector sets",
          "Applying gates and building background",
          "Unmixing events in R",
          "Caching QC metrics and plots",
          "Writing manifest and reports",
        ];
        return {
          ...current,
          progress: nextProgress,
          subtask: subtasks[Math.min(4, Math.floor(nextProgress / 20))],
        };
      });
    }, 650);
    return () => window.clearTimeout(timer);
  }, [job]);

  useEffect(() => {
    if (!toast) return;
    const timer = window.setTimeout(() => setToast(null), 4200);
    return () => window.clearTimeout(timer);
  }, [toast]);

  useEffect(() => {
    const root = document.documentElement;
    const appearance = settings.appearance;
    const scale = normalizeInterfaceScale(appearance.fontScale) * 1.4 / 100;
    root.dataset.theme = appearance.theme;
    root.dataset.density = appearance.density;
    root.dataset.shadows = appearance.shadows;
    root.dataset.contrast = appearance.highContrast ? "high" : "normal";
    root.dataset.texture = appearance.backgroundTexture ? "on" : "off";
    root.dataset.motion = appearance.reduceMotion ? "reduced" : "full";
    root.dataset.stickyHeader = appearance.stickyHeader ? "on" : "off";
    root.style.setProperty("--ui-zoom", String(scale));
    root.style.setProperty("--scaled-viewport-height", `${100 / scale}vh`);
    root.style.setProperty("--corner-radius", `${appearance.cornerRadius}px`);
    const fonts = {
      avenir: '"Avenir Next", Avenir, sans-serif',
      futura: 'Futura, "Century Gothic", "Avenir Next", sans-serif',
      atkinson: '"Atkinson Hyperlegible", "Arial Nova", sans-serif',
      charter: 'Charter, "Bitstream Charter", Georgia, serif',
      palatino: 'Palatino, "Palatino Linotype", "Book Antiqua", serif',
      monaco: 'Monaco, Menlo, Consolas, monospace',
      system: '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif',
    };
    root.style.setProperty("--ui-font", fonts[appearance.fontFamily]);
    window.localStorage.setItem("spectreasy-theme", appearance.theme);
  }, [settings.appearance]);

  function updateMapping(id: string, patch: Partial<MappingRow>) {
    setProject((current) => ({
      ...current,
      mappingDirty: true,
      mapping: current.mapping.map((row) =>
        row.id === id ? { ...row, ...patch } : row,
      ),
      artifacts: current.artifacts.map((artifact) =>
        ["reference", "unmixing", "control-report"].includes(artifact.id)
          ? { ...artifact, status: "stale" }
          : artifact,
      ),
    }));
  }

  function updateSettings(
    section: "projectPath" | "control" | "sample" | "af" | "appearance",
    patch:
      | Partial<WorkflowSettings["control"]>
      | Partial<WorkflowSettings["sample"]>
      | Partial<WorkflowSettings["af"]>
      | Partial<WorkflowSettings["appearance"]>
      | { projectPath: string },
  ) {
    if (section === "projectPath") {
      setSettings((current) => ({
        ...current,
        projectPath: (patch as { projectPath: string }).projectPath,
      }));
      return;
    }
    setSettings((current) => ({
      ...current,
      [section]: {
        ...current[section],
        ...patch,
        ...(section === "control" &&
        "method" in patch &&
        patch.method !== "AutoSpectral"
          ? { refine: false }
          : {}),
      },
    }));
  }

  function changeCytometer(cytometer: string) {
    cytometer = normalizeCockpitCytometer(cytometer);
    window.localStorage.setItem("spectreasy-cytometer", cytometer);
    setSettings((current) => ({
      ...current,
      control: { ...current.control, cytometer },
    }));
    setProject((current) => ({
      ...current,
      cytometer: cytometerLabels[cytometer] ?? cytometer,
    }));
  }

  function changeMethod(method: string) {
    setSettings((current) => ({
      ...current,
      control: {
        ...current.control,
        method,
        refine: method === "AutoSpectral" ? current.control.refine : false,
      },
      sample: { ...current.sample, method },
    }));
    setProject((current) => ({ ...current, method }));
  }

  function toggleSettings() {
    if (activeSection === "settings") {
      setActiveSection(sectionBeforeSettings);
      return;
    }
    setSectionBeforeSettings(activeSection);
    setActiveSection("settings");
  }

  function openApplet(applet: CockpitAppletId, section = activeSection) {
    setSectionBeforeApplet(activeSection);
    setActiveSection(section);
    setActiveApplet(applet);
  }

  function navigateToSection(section: SectionId) {
    if (section === "panel") {
      openApplet("panel-builder", section);
      return;
    }
    if (section === "matrix") {
      openApplet("matrix-adjustment", section);
      return;
    }
    if (section === "control-reports") {
      openApplet("control-qc-report", section);
      return;
    }
    if (section === "sample-reports") {
      openApplet("sample-qc-report", section);
      return;
    }
    setActiveSection(section);
  }

  function exitApplet() {
    setActiveApplet(null);
    setActiveSection(sectionBeforeApplet);
    void refreshProject();
  }

  async function runAction(
    action: "control" | "sample" | "control-report" | "sample-report" | "af",
    label: string,
  ) {
    setJob({
      label,
      state: "running",
      progress: 0,
      subtask: "Preparing inputs and checking prerequisites",
      startedAt: "now",
    });
    const control = settings.control;
    const sample = settings.sample;
    const af = settings.af;
    const payload =
      action === "control"
        ? {
            projectPath: settings.projectPath,
            scc_dir: "scc",
            control_file: "fcs_mapping.csv",
            output_dir: "spectreasy_outputs/unmix_controls",
            method: control.method,
            cytometer: control.cytometer,
            auto_create_mapping: control.autoCreateMapping,
            auto_unknown_fluor_policy: control.autoUnknownFluorPolicy,
            gate_file: "ssc_gate_config.csv",
            gating_mode: "reuse",
            af_n_bands: control.afNBands,
            af_max_cells: control.afMaxCells,
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
            n_threads: control.unmixThreads,
            seed: control.seed,
            save_qc_png: control.saveQcPlots,
            save_report: control.saveReport,
            report_format: control.outputFormat,
            scc_background_method: control.sccBackgroundMethod,
            scc_background_k: control.sccBackgroundK,
            spectral_variant_som_nodes: control.spectralVariantSomNodes,
            spectral_variant_top_k: control.spectralVariantTopK,
            spectral_variant_cosine_threshold:
              control.spectralVariantCosineThreshold,
            spectral_variant_max_variants: control.spectralVariantMaxVariants,
            spectral_variant_min_events: control.spectralVariantMinEvents,
            spectreasy_weight_quantile: control.spectreasyWeightQuantile,
            autospectral_n_candidates: control.autospectralNCandidates,
            autospectral_n_spectral: control.autospectralNSpectral,
            autospectral_min_events: control.autospectralMinEvents,
            autospectral_refine: control.refine,
          }
        : action === "sample"
          ? {
              projectPath: settings.projectPath,
              sample_dir: "samples",
              matrix_file: "spectreasy_outputs/unmix_controls/scc_reference_matrix.csv",
              detector_noise_file: "spectreasy_outputs/unmix_controls/scc_detector_noise.csv",
              output_dir: "spectreasy_outputs/unmix_samples/unmixed_fcs",
              method: sample.method,
              rwls_max_iter: sample.rwlsMaxIter,
              n_threads: sample.nThreads,
              spectral_variant_top_k: sample.spectralVariantTopK,
              spectral_variant_min_abundance:
                sample.spectralVariantMinAbundance,
              spectral_variant_positive_fraction:
                sample.spectralVariantPositiveFraction,
              spectral_variant_min_improvement:
                sample.spectralVariantMinImprovement,
              spectral_variant_library_file: sample.spectralVariantLibraryFile,
              spectreasy_weight_quantile: sample.spectreasyWeightQuantile,
              estimate_af: sample.estimateAf,
              write_fcs: sample.writeFcs,
              save_report: sample.saveReport,
              report_format: sample.outputFormat,
              save_qc_plots: sample.saveQcPlots,
              plot_n_events: sample.plotNEvents,
              chunk_size: sample.chunkSize,
              seed: sample.seed,
              return_type: sample.returnType,
            }
          : action === "af"
            ? {
                projectPath: settings.projectPath,
                fcs_file: af.fcsFile,
                save_name: af.saveName,
                save_overwrite: af.saveOverwrite,
                af_n_bands: af.afNBands,
                af_max_cells: af.afMaxCells,
                seed: af.seed,
              }
            : {
                projectPath: settings.projectPath,
                report_type: action === "sample-report" ? "sample" : "control",
                report_format: action === "sample-report" ? sample.outputFormat : control.outputFormat,
                overwrite: window.prompt(
                  "Existing report behavior: version (recommended), overwrite, or cancel",
                  "version",
                ) ?? "cancel",
              };
    if (
      (action === "control-report" || action === "sample-report") &&
      (payload as Record<string, unknown>).overwrite === "cancel"
    ) {
      setJob(emptyJob);
      setToast("Report generation cancelled. No file was changed.");
      return false;
    }
    const response = await attemptWorkflowAction(action, payload);
    setToast(response.message);
    if (response.connected) {
      window.setTimeout(() => void refreshProject(), 500);
    }
    return response.connected;
  }

  async function loadPanel() {
    const payload = await loadPanelPayload("aurora");
    if (payload) {
      setPanelPayload(payload);
      setToast("Panel library refreshed from the R backend.");
    } else {
      setToast(
        "Showing the bundled panel library preview. The R backend is not connected.",
      );
    }
  }

  async function saveMapping() {
    const mapping = await persistControlMapping(project.mapping);
    if (mapping.success)
      setProject((current) => ({ ...current, mappingDirty: false }));
    setToast(mapping.message);
  }

  async function chooseProject(mode: "create" | "open") {
    const response = mode === "create"
      ? await createProjectFolder()
      : await selectProjectFolder();
    if (!response.cancelled && !response.success) setToast(response.message);
    if (response.success) {
      await refreshProject(true);
    }
  }

  async function createMapping() {
    const response = await createControlMapping();
    setToast(response.message);
    if (response.success) await refreshProject();
  }

  async function terminateBackend() {
    const terminated = await terminateRSession();
    if (terminated) {
      setBackend({
        ...initialBackendStatus,
        message: "The local R session was terminated from the cockpit.",
      });
      setToast("R session terminated.");
    } else {
      setToast("The R session could not be terminated.");
    }
    return terminated;
  }

  function downloadArtifact(artifact: Artifact) {
    const manifest = `Spectreasy artifact\n\nName: ${artifact.name}\nType: ${artifact.type}\nStatus: ${artifact.status}\nPath: ${artifact.path}\nUpdated: ${artifact.updated}\nRun: ${artifact.run ?? "user supplied"}\n`;
    downloadTextFile(`${artifact.name.replaceAll("/", "_")}.txt`, manifest);
    setToast(`Prepared a local download for ${artifact.name}.`);
  }

  function selectArtifact(artifact: Artifact) {
    const sectionMap: Record<string, SectionId> = {
      Controls: "controls",
      Gates: "controls",
      Samples: "samples",
      Matrices: "matrix",
      Reports: "control-reports",
      "QC Metrics": "control-reports",
      "AF Profiles": "af",
      "Panel Builder": "panel",
      Logs: "settings",
    };
    const section = sectionMap[artifact.group];
    if (section) navigateToSection(section);
  }

  const activeTitle = useMemo(
    () =>
      ({
        controls: "Controls",
        samples: "Samples",
        "control-reports": "Controls QC report",
        "sample-reports": "Samples QC report",
        matrix: "Matrix review",
        panel: "Panel builder",
        af: "AF library",
        settings: "Settings & logs",
      })[activeSection],
    [activeSection],
  );

  return (
    <div className="cockpit-app">
      <TopBar
        project={project}
        cytometer={settings.control.cytometer}
        method={settings.control.method}
        darkMode={darkMode}
        settingsActive={activeSection === "settings"}
        terminalActive={terminalOpen}
        backendConnected={backend.connected}
        onCytometerChange={changeCytometer}
        onMethodChange={changeMethod}
        onToggleTheme={() =>
          updateSettings("appearance", { theme: darkMode ? "light" : "dark" })
        }
        onFiles={() => setFilesOpen(true)}
        onCreateProject={() => void chooseProject("create")}
        onOpenProject={() => void chooseProject("open")}
        onSettings={toggleSettings}
        onTerminal={() => setTerminalOpen((value) => !value)}
      />
      <div className="app-body">
        <main className="main-area">
          <WorkflowRail
            activeSection={activeSection}
            project={project}
            showCounts={settings.appearance.showSectionCounts}
            width={settings.appearance.sidebarWidth}
            onWidthChange={(sidebarWidth) => updateSettings("appearance", { sidebarWidth })}
            onChange={navigateToSection}
          />
          <div className="main-canvas">
            <div className="mobile-canvas-title">
              <Menu size={15} />
              <span>{activeTitle}</span>
            </div>
            <WorkflowWorkspace
              project={project}
              backend={backend}
              job={job}
              activeSection={activeSection}
              mappingTab={mappingTab}
              setMappingTab={setMappingTab}
              onUpdateMapping={updateMapping}
              onRun={runAction}
              onRefresh={() => void refreshProject()}
              onDownload={downloadArtifact}
              onSelectArtifact={selectArtifact}
              onLoadPanel={() => void loadPanel()}
              panelPayload={panelPayload}
              onSave={() => void saveMapping()}
              onCreateMapping={() => void createMapping()}
              onSectionChange={navigateToSection}
              settings={settings}
              onSettingsChange={updateSettings}
              onOpenApplet={(applet) => openApplet(applet)}
            />
          </div>
        </main>
      </div>
      {toast && (
        <div className="toast">
          <Sparkles size={15} />
          <span>{toast}</span>
          <button onClick={() => setToast(null)} aria-label="Dismiss message">
            <X size={14} />
          </button>
        </div>
      )}
      {activeApplet && (
        <CockpitApplet applet={activeApplet} theme={settings.appearance.theme} onExit={exitApplet} />
      )}
      {filesOpen && (
        <ProjectFilesDialog
          projectName={project.projectName}
          onClose={() => setFilesOpen(false)}
          onChanged={() => refreshProject(false)}
        />
      )}
      {terminalOpen && <TerminalPanel
        connected={backend.connected}
        projectPath={project.projectPath}
        widthPct={settings.appearance.terminalWidthPct}
        heightPct={settings.appearance.terminalHeightPct}
        onClose={() => setTerminalOpen(false)}
        onRefresh={() => refreshProject(false)}
        onTerminate={terminateBackend}
        onSizeChange={(size) => updateSettings("appearance", size)}
      />}
    </div>
  );
}
