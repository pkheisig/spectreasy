import { useCallback, useEffect, useMemo, useState } from "react";
import {
  Menu,
  Moon,
  RefreshCcw,
  Settings2,
  Sparkles,
  Sun,
  X,
} from "lucide-react";
import {
  attemptWorkflowAction,
  downloadTextFile,
  initialBackendStatus,
  loadPanelPayload,
  loadProjectSnapshot,
  persistControlMapping,
  persistGuiState,
  setProjectContext,
} from "./api";
import { demoProject } from "./mockData";
import { WorkflowRail } from "./components/WorkflowRail";
import { WorkflowWorkspace } from "./workspaces/WorkflowWorkspace";
import { defaultWorkflowSettings } from "./types";
import type {
  Artifact,
  BackendStatus,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  SectionId,
  WorkflowSettings,
} from "./types";
import "./cockpit.css";

const emptyJob: Job = { label: "", state: "idle", progress: 0, subtask: "" };

type TopBarProps = {
  project: ProjectState;
  backend: BackendStatus;
  cytometer: string;
  method: string;
  darkMode: boolean;
  settingsActive: boolean;
  onCytometerChange: (value: string) => void;
  onMethodChange: (value: string) => void;
  onToggleTheme: () => void;
  onRefresh: () => void;
  onSettings: () => void;
};

function TopBar({
  project,
  backend,
  cytometer,
  method,
  darkMode,
  settingsActive,
  onCytometerChange,
  onMethodChange,
  onToggleTheme,
  onRefresh,
  onSettings,
}: TopBarProps) {
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
          <small>Spectral analysis</small>
        </div>
      </div>
      <div className="project-switcher">
        <div>
          <span className="eyebrow">Active project</span>
          <strong>{project.projectName}</strong>
        </div>
      </div>
      <div className="topbar-spacer" />
      <label className="context-select">
        <span className="chip-label">Cytometer</span>
        <select
          aria-label="Cytometer"
          value={cytometer}
          onChange={(event) => onCytometerChange(event.target.value)}
        >
          <option value="auto">Auto</option>
          <option value="aurora">Cytek Aurora</option>
          <option value="aurora_5l">Cytek Aurora 5L</option>
          <option value="aurora_4l">Cytek Aurora 4L</option>
          <option value="discover">Cytek Aurora Discover</option>
          <option value="id7000">Sony ID7000</option>
          <option value="xenith">BD FACSDiscover Xenith</option>
        </select>
      </label>
      <label className="context-select method-select">
        <span className="chip-label">Method</span>
        <select
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
        </select>
      </label>
      <div
        className="backend-chip"
        title="R performs the calculations locally. This indicator only reports whether the browser can reach it."
      >
        <span
          className={`backend-dot ${backend.connected ? "is-connected" : ""}`}
        />
        <div>
          <span className="chip-label">Local R</span>
          <strong>{backend.connected ? "Connected" : "Not connected"}</strong>
        </div>
      </div>
      <button
        className="topbar-icon"
        onClick={onRefresh}
        aria-label="Reload project"
        title="Reload project"
      >
        <RefreshCcw size={17} />
      </button>
      <button
        className="topbar-icon"
        onClick={onToggleTheme}
        aria-label={darkMode ? "Use light mode" : "Use dark mode"}
        title={darkMode ? "Use light mode" : "Use dark mode"}
      >
        {darkMode ? <Sun size={17} /> : <Moon size={17} />}
      </button>
      <button
        className={`topbar-icon ${settingsActive ? "is-active" : ""}`}
        onClick={onSettings}
        aria-label={settingsActive ? "Return to workflow" : "Settings"}
        aria-pressed={settingsActive}
        title={settingsActive ? "Return to workflow" : "Settings"}
      >
        <Settings2 size={17} />
      </button>
    </header>
  );
}

export default function CockpitApp() {
  const [project, setProject] = useState<ProjectState>(() =>
    structuredClone(demoProject),
  );
  const [backend, setBackend] = useState<BackendStatus>(initialBackendStatus);
  const [activeSection, setActiveSection] = useState<SectionId>("overview");
  const [mappingTab, setMappingTab] = useState<
    "mapping" | "gating" | "build" | "qc"
  >("mapping");
  const [job, setJob] = useState<Job>(emptyJob);
  const [panelPayload, setPanelPayload] = useState<PanelPayload | null>(null);
  const [toast, setToast] = useState<string | null>(null);
  const [settings, setSettings] = useState<WorkflowSettings>(() => {
    const initial = defaultWorkflowSettings(demoProject.projectPath);
    const savedTheme = window.localStorage.getItem("spectreasy-theme");
    if (savedTheme === "light" || savedTheme === "dark") {
      initial.appearance.theme = savedTheme;
    }
    return initial;
  });
  const [sectionBeforeSettings, setSectionBeforeSettings] =
    useState<SectionId>("overview");
  const darkMode = settings.appearance.theme === "dark";

  const refreshProject = useCallback(async (initial = false) => {
    const snapshot = await loadProjectSnapshot();
    setProject(snapshot.project);
    setBackend(snapshot.backend);
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
        return {
          ...current,
          ...saved,
          projectPath: snapshot.project.projectPath || current.projectPath,
          control: {
            ...current.control,
            ...savedControl,
            method:
              savedControl.method ??
              snapshot.project.method ??
              current.control.method,
            cytometer: explicitCytometer ?? "auto",
          },
          sample: {
            ...current.sample,
            ...savedSample,
            method:
              savedSample.method ??
              savedControl.method ??
              snapshot.project.method ??
              current.sample.method,
          },
          af: { ...current.af, ...(saved.af ?? {}) },
          appearance: { ...current.appearance, ...(saved.appearance ?? {}) },
        };
      });
    }
    if (!initial)
      setToast(
        snapshot.backend.connected
          ? "Project rescanned from the R backend."
          : "Preview state restored. Start the R backend for live artifacts.",
      );
  }, []);

  useEffect(() => {
    const timer = window.setTimeout(() => void refreshProject(true), 0);
    return () => window.clearTimeout(timer);
  }, [refreshProject]);

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
    const scale = appearance.fontScale / 100;
    root.dataset.theme = appearance.theme;
    root.dataset.density = appearance.density;
    root.dataset.sidebar = appearance.sidebarWidth;
    root.dataset.shadows = appearance.shadows;
    root.dataset.contrast = appearance.highContrast ? "high" : "normal";
    root.dataset.texture = appearance.backgroundTexture ? "on" : "off";
    root.dataset.motion = appearance.reduceMotion ? "reduced" : "full";
    root.dataset.stickyHeader = appearance.stickyHeader ? "on" : "off";
    root.style.setProperty("--ui-zoom", String(scale));
    root.style.setProperty("--scaled-viewport-height", `${100 / scale}vh`);
    root.style.setProperty("--corner-radius", `${appearance.cornerRadius}px`);
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
    window.localStorage.setItem("spectreasy-cytometer", cytometer);
    setSettings((current) => ({
      ...current,
      control: { ...current.control, cytometer },
    }));
    const labels: Record<string, string> = {
      auto: "Auto",
      aurora: "Cytek Aurora",
      aurora_5l: "Cytek Aurora 5L",
      aurora_4l: "Cytek Aurora 4L",
      discover: "Cytek Aurora Discover",
      id7000: "Sony ID7000",
      xenith: "BD FACSDiscover Xenith",
    };
    setProject((current) => ({
      ...current,
      cytometer: labels[cytometer] ?? cytometer,
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

  async function runAction(
    action: "control" | "sample" | "report" | "af",
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
            scc_dir: control.sccDir,
            control_file: control.controlFile,
            output_dir: control.outputDir,
            method: control.method,
            cytometer: control.cytometer,
            auto_create_mapping: control.autoCreateMapping,
            auto_unknown_fluor_policy: control.autoUnknownFluorPolicy,
            gate_file: control.gateFile,
            gating_mode: control.gateFile.trim().length > 0 ? "reuse" : "automatic",
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
            refine: control.refine,
          }
        : action === "sample"
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
                af_min_cluster_events: af.afMinClusterEvents,
                af_min_cluster_proportion: af.afMinClusterProportion,
                seed: af.seed,
              }
            : {
                projectPath: settings.projectPath,
                report_type: "control",
                report_format: control.outputFormat,
                overwrite: window.prompt(
                  "Existing report behavior: version (recommended), overwrite, or cancel",
                  "version",
                ) ?? "cancel",
              };
    if (
      action === "report" &&
      (payload as Record<string, unknown>).overwrite === "cancel"
    ) {
      setJob(emptyJob);
      setToast("Report generation cancelled. No file was changed.");
      return;
    }
    const response = await attemptWorkflowAction(action, payload);
    setToast(response.message);
    if (response.connected) {
      window.setTimeout(() => void refreshProject(), 500);
    }
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

  async function saveProject() {
    const mapping = await persistControlMapping(project.mapping);
    const saved = await persistGuiState(project, settings);
    if (mapping.success)
      setProject((current) => ({ ...current, mappingDirty: false }));
    setToast(
      saved && mapping.success
        ? "Project preferences and control mapping saved to the local R config."
        : "Project preferences saved for this browser session. Connect R to persist them.",
    );
  }

  async function openProject(path: string) {
    const response = await setProjectContext(path);
    setToast(response.message);
    if (response.success) {
      setSettings((current) => ({ ...current, projectPath: path }));
      await refreshProject();
    }
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
      Reports: "reports",
      "QC Metrics": "reports",
      "AF Profiles": "af",
      "Panel Builder": "panel",
      Logs: "settings",
    };
    const section = sectionMap[artifact.group];
    if (section) setActiveSection(section);
  }

  const activeTitle = useMemo(
    () =>
      ({
        overview: "Project setup",
        controls: "Controls",
        samples: "Samples",
        reports: "Reports",
        matrix: "Matrix review",
        panel: "Panel builder",
        af: "AF library",
        comparison: "Method comparison",
        simulator: "Synthetic SCC",
        settings: "Settings & logs",
      })[activeSection],
    [activeSection],
  );

  return (
    <div
      className={`cockpit-app ${activeSection === "controls" && mappingTab === "gating" ? "is-gating-editor" : ""}`}
    >
      <TopBar
        project={project}
        backend={backend}
        cytometer={settings.control.cytometer}
        method={settings.control.method}
        darkMode={darkMode}
        settingsActive={activeSection === "settings"}
        onCytometerChange={changeCytometer}
        onMethodChange={changeMethod}
        onToggleTheme={() =>
          updateSettings("appearance", { theme: darkMode ? "light" : "dark" })
        }
        onRefresh={() => void refreshProject()}
        onSettings={toggleSettings}
      />
      <div className="app-body">
        <main className="main-area">
          <WorkflowRail
            activeSection={activeSection}
            project={project}
            showCounts={settings.appearance.showSectionCounts}
            onChange={setActiveSection}
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
              onSave={() => void saveProject()}
              onSectionChange={setActiveSection}
              settings={settings}
              onSettingsChange={updateSettings}
              onOpenProject={openProject}
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
    </div>
  );
}
