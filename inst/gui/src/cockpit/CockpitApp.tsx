import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { LoaderCircle, Menu } from "lucide-react";
import {
  attemptWorkflowAction,
  createControlMapping,
  createProjectFolder,
  downloadTextFile,
  initialBackendStatus,
  initializeProject,
  loadProjectInputDirectories,
  loadPanelPayload,
  loadProjectSnapshot,
  persistControlMapping,
  persistGuiState,
  selectProjectFolder,
  setProjectContext,
  updateProjectInputDirectory,
} from "./api";
import { emptyProject } from "./projectState";
import { createRefreshSequence, mergeProjectRefresh } from "./projectRefresh";
import { cytometerLabels, normalizeCockpitCytometer } from "./cytometerConfig";
import {
  cockpitExecutionCode,
  createWorkflowPayload,
  normalizeCockpitOutputRoot,
  workflowHasExistingResults,
} from "./workflowExecution";
import { WorkflowRail } from "./components/WorkflowRail";
import { CockpitApplet } from "./components/CockpitApplet";
import { TopBar } from "./components/CockpitTopBar";
import { ProjectFilesDialog } from "./components/ProjectFilesDialog";
import { ProjectInitializationDialog } from "./components/ProjectInitializationDialog";
import { VersionRunDialog } from "./components/VersionRunDialog";
import { GuideDialog } from "./components/GuideDialog";
import { WorkflowWorkspace } from "./workspaces/WorkflowWorkspace";
import { defaultWorkflowSettings, normalizeInterfaceScale, normalizeWorkflowUnmixingMethod } from "./types";
import type {
  Artifact,
  BackendStatus,
  CockpitAppletId,
  ExecutionLogEntry,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  SectionId,
  WorkflowSettings,
} from "./types";
import "./cockpit.css";

const emptyJob: Job = { label: "", state: "idle", progress: 0, subtask: "" };
const persistentDataApplets: CockpitAppletId[] = ["control-gating", "panel-builder", "matrix-adjustment"];


export default function CockpitApp() {
  const [project, setProject] = useState<ProjectState>(() =>
    structuredClone(emptyProject),
  );
  const [backend, setBackend] = useState<BackendStatus>(initialBackendStatus);
  const [activeSection, setActiveSection] = useState<SectionId>("controls");
  const [activeApplet, setActiveApplet] = useState<CockpitAppletId | null>(null);
  const [selectedReportPath, setSelectedReportPath] = useState("");
  const [initializedApplets, setInitializedApplets] = useState<CockpitAppletId[]>([]);
  const [sectionBeforeApplet, setSectionBeforeApplet] =
    useState<SectionId>("controls");
  const [mappingTab, setMappingTab] = useState<
    "mapping" | "gating" | "build" | "qc"
  >("mapping");
  const [sampleTab, setSampleTab] = useState<"unmixing" | "qc">("unmixing");
  const [job, setJob] = useState<Job>(emptyJob);
  const [panelPayload, setPanelPayload] = useState<PanelPayload | null>(null);
  const [terminalOpen, setTerminalOpen] = useState(false);
  const [executionLogs, setExecutionLogs] = useState<ExecutionLogEntry[]>([]);
  const [filesOpen, setFilesOpen] = useState(false);
  const [guideOpen, setGuideOpen] = useState(false);
  const [initializationDismissedFor, setInitializationDismissedFor] = useState("");
  const [initializationBusy, setInitializationBusy] = useState(false);
  const [initializationMessage, setInitializationMessage] = useState("");
  const [projectLoading, setProjectLoading] = useState(false);
  const [projectPickerOpen, setProjectPickerOpen] = useState(false);
  const [pendingVersionRun, setPendingVersionRun] = useState<{ action: "control" | "sample"; label: string } | null>(null);
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
  const refreshSequence = useRef(createRefreshSequence());
  const projectRef = useRef<ProjectState>(structuredClone(emptyProject));
  const darkMode = settings.appearance.theme === "dark";
  const projectIsLoading = projectLoading || !settingsReady;

  function appendExecutionLogs(entries: Array<Pick<ExecutionLogEntry, "kind" | "text">>) {
    setExecutionLogs((current) => {
      const next = [...current];
      entries.forEach((entry, index) => {
        const previous = next.at(-1);
        if (entry.kind === "info" && previous?.kind === "info") {
          next[next.length - 1] = { ...previous, text: `${previous.text}\n${entry.text}` };
        } else {
          next.push({
            ...entry,
            id: `${Date.now()}-${index}-${Math.random().toString(16).slice(2)}`,
          });
        }
      });
      return next.slice(-250);
    });
  }

  const refreshProject = useCallback(async (initial = false, requestedProjectPath = "") => {
    const requestSequence = refreshSequence.current.begin();
    const tabProjectPath = requestedProjectPath || window.sessionStorage.getItem("spectreasy-project-path") || "";
    const snapshot = await loadProjectSnapshot(tabProjectPath, initial ? undefined : projectRef.current);
    if (!refreshSequence.current.isCurrent(requestSequence)) return;
    if (!snapshot.backend.connected) {
      setBackend((current) => ({
        ...current,
        connected: false,
        packageReady: false,
        message: "The local R backend did not answer. Retrying…",
      }));
      return;
    }
    setProject((current) => mergeProjectRefresh(current, snapshot.project));
    if (snapshot.project.projectPath) window.sessionStorage.setItem("spectreasy-project-path", snapshot.project.projectPath);
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
        const savedControl = (saved.control ?? {}) as Partial<WorkflowSettings["control"]>;
        const savedSample: Partial<WorkflowSettings["sample"]> =
          saved.sample ?? {};
        const explicitCytometer = window.localStorage.getItem(
          "spectreasy-cytometer",
        );
        const afFiles = snapshot.project.mapping
          .filter((row) => row.marker.trim().toLowerCase() === "autofluorescence" || /^af(?:$|_|\b)/i.test(row.fluorophore.trim()))
          .map((row) => `${snapshot.project.controlInputDir}/${row.file}`);
        const savedAf: Partial<WorkflowSettings["af"]> = saved.af ?? {};
        const requestedAfFile = savedAf.fcsFile ?? current.af.fcsFile;
        return {
          ...current,
          ...saved,
          projectPath: snapshot.project.projectPath || current.projectPath,
          control: {
            ...current.control,
            ...savedControl,
            sccDir: snapshot.project.controlInputDir,
            controlFile: savedControl.controlFile ?? "fcs_mapping.csv",
            outputDir: normalizeCockpitOutputRoot(savedControl.outputDir),
            // A cockpit project with confirmed gates always reuses its gate CSV.
            manualGateFile: snapshot.project.scan.gates > 0 ? "ssc_gate_config.csv" : "",
            method: normalizeWorkflowUnmixingMethod(
              savedControl.method ??
              snapshot.project.method ??
              current.control.method,
            ),
            cytometer: normalizeCockpitCytometer(
              explicitCytometer ?? savedControl.cytometer ?? snapshot.project.cytometer,
            ),
          },
          sample: {
            ...current.sample,
            ...savedSample,
            sampleDir: snapshot.project.sampleInputDir,
            matrixFile: savedSample.matrixFile ?? "spectreasy_outputs/unmix_controls/scc_reference_matrix.csv",
            detectorNoiseFile: savedSample.detectorNoiseFile ?? "spectreasy_outputs/unmix_controls/scc_detector_noise.csv",
            outputDir: normalizeCockpitOutputRoot(savedSample.outputDir),
            method: normalizeWorkflowUnmixingMethod(
              savedSample.method ??
              savedControl.method ??
              snapshot.project.method ??
              current.sample.method,
            ),
          },
          af: {
            ...current.af,
            ...savedAf,
            fcsFile: afFiles.includes(requestedAfFile) ? requestedAfFile : (afFiles[0] ?? ""),
          },
          appearance: {
            ...current.appearance,
            theme: saved.appearance?.theme === "dark" ? "dark" : (saved.appearance?.theme === "light" ? "light" : current.appearance.theme),
            density: saved.appearance?.density === "compact" || saved.appearance?.density === "spacious" || saved.appearance?.density === "comfortable"
              ? saved.appearance.density
              : current.appearance.density,
            fontFamily:
              String(saved.appearance?.fontFamily ?? "") === "source-sans"
                ? "atkinson"
                : (saved.appearance?.fontFamily ?? current.appearance.fontFamily),
            fontScale: normalizeInterfaceScale(saved.appearance?.fontScale ?? current.appearance.fontScale),
            sidebarWidth:
              saved.appearance?.sidebarWidth === 220
                ? 242
                : (saved.appearance?.sidebarWidth ?? current.appearance.sidebarWidth),
            cornerRadius: typeof saved.appearance?.cornerRadius === "number" ? saved.appearance.cornerRadius : current.appearance.cornerRadius,
            shadows: saved.appearance?.shadows === "none" || saved.appearance?.shadows === "raised" || saved.appearance?.shadows === "subtle"
              ? saved.appearance.shadows
              : current.appearance.shadows,
            backgroundTexture: typeof saved.appearance?.backgroundTexture === "boolean" ? saved.appearance.backgroundTexture : current.appearance.backgroundTexture,
          },
        };
      });
    }
    return snapshot;
  }, []);

  useEffect(() => {
    projectRef.current = project;
  }, [project]);

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
    if (!settingsReady || backend.connected) return;
    const timer = window.setInterval(() => {
      void refreshProject(false);
    }, 1500);
    return () => window.clearInterval(timer);
  }, [backend.connected, refreshProject, settingsReady]);

  useEffect(() => {
    if (!backend.connected || !project.projectPath || projectPickerOpen || projectLoading) return;
    let checking = false;
    const checkLayout = async () => {
      if (checking) return;
      checking = true;
      const layout = await loadProjectInputDirectories(project.projectPath);
      checking = false;
      if (!layout) return;
      if (layout.controlInputDir !== project.controlInputDir || layout.sampleInputDir !== project.sampleInputDir) {
        await refreshProject(false);
      }
    };
    const timer = window.setInterval(() => void checkLayout(), 2000);
    return () => window.clearInterval(timer);
  }, [backend.connected, project.controlInputDir, project.projectPath, project.sampleInputDir, projectLoading, projectPickerOpen, refreshProject]);

  useEffect(() => {
    if (!backend.connected || !project.projectPath || projectPickerOpen || projectLoading) return;
    let refreshing = false;
    const timer = window.setInterval(() => {
      if (refreshing) return;
      refreshing = true;
      void refreshProject(false).finally(() => { refreshing = false; });
    }, 10000);
    return () => window.clearInterval(timer);
  }, [backend.connected, project.projectPath, projectLoading, projectPickerOpen, refreshProject]);

  useEffect(() => {
    if (!settingsReady || projectPickerOpen || projectLoading) return;
    const projectState = {
      projectPath: project.projectPath,
      projectName: project.projectName,
      cytometer: project.cytometer,
      method: project.method,
    };
    const timer = window.setTimeout(() => void persistGuiState(projectState, settings), 450);
    return () => window.clearTimeout(timer);
  }, [project.cytometer, project.method, project.projectName, project.projectPath, projectLoading, projectPickerOpen, settings, settingsReady]);

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
    const root = document.documentElement;
    const appearance = settings.appearance;
    const scale = normalizeInterfaceScale(appearance.fontScale) * 1.4 / 100;
    root.dataset.theme = appearance.theme;
    root.dataset.density = appearance.density;
    root.dataset.shadows = appearance.shadows;
    root.dataset.texture = appearance.backgroundTexture ? "on" : "off";
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

  async function updateInputDirectory(role: "controls" | "samples", path: string) {
    const currentPath = role === "controls" ? project.controlInputDir : project.sampleInputDir;
    if (!settings.projectPath || path.trim() === currentPath) return true;
    const result = await updateProjectInputDirectory(role, path.trim(), settings.projectPath);
    appendExecutionLogs([{ kind: result.success ? "success" : "error", text: result.message }]);
    if (!result.success) {
      setSettings((current) => ({
        ...current,
        control: role === "controls" ? { ...current.control, sccDir: project.controlInputDir } : current.control,
        sample: role === "samples" ? { ...current.sample, sampleDir: project.sampleInputDir } : current.sample,
      }));
      return false;
    }
    await refreshProject(false);
    return true;
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

  function openApplet(applet: CockpitAppletId, section = activeSection, reportPath = "") {
    setSectionBeforeApplet(activeSection);
    setActiveSection(section);
    setSelectedReportPath(reportPath);
    if (persistentDataApplets.includes(applet)) {
      setInitializedApplets((current) => current.includes(applet) ? current : [...current, applet]);
    }
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
    setActiveSection(section);
  }

  async function exitApplet(reason: "exit" | "confirmed" = "exit") {
    const closingApplet = activeApplet;
    setActiveApplet(null);
    setSelectedReportPath("");
    setActiveSection(sectionBeforeApplet);
    const snapshot = await refreshProject();
    if (closingApplet === "control-gating") {
      if (reason === "confirmed" || snapshot?.project.scan.gates) {
        setSettings((current) => ({
          ...current,
          control: { ...current.control, manualGateFile: "ssc_gate_config.csv" },
        }));
      }
      if (reason === "confirmed") {
        setActiveSection("controls");
        setMappingTab("build");
      }
    }
  }

  async function runAction(
    action: "control" | "sample" | "af",
    label: string,
  ) {
    if (action === "control" && !settings.control.manualGateFile) {
      appendExecutionLogs([{ kind: "warning", text: "Review and confirm the control gates before running the control workflow." }]);
      openApplet("control-gating", "controls");
      return false;
    }
    setJob({
      label,
      state: "running",
      progress: 0,
      subtask: "",
      startedAt: "now",
    });
    const payload = createWorkflowPayload(action, settings);
    appendExecutionLogs([
      { kind: "command", text: cockpitExecutionCode(action, payload as Record<string, unknown>) },
      { kind: "info", text: `${label} started.` },
    ]);
    const response = await attemptWorkflowAction(action, payload);
    appendExecutionLogs([
      ...response.logs.map((text) => ({
        kind: /^warning\b/i.test(text) ? "warning" as const : "info" as const,
        text,
      })),
      {
        kind: response.success ? "success" as const : "error" as const,
        text: response.success ? `${label} completed.` : response.message,
      },
    ]);
    if (response.success) {
      setJob({
        label,
        state: "complete",
        progress: 100,
        subtask: "",
        startedAt: "now",
        finishedAt: "now",
        output: response.outputCount ? `${response.outputCount} output${response.outputCount === 1 ? "" : "s"}` : undefined,
      });
      await refreshProject();
      if (action === "control") {
        setActiveSection("controls");
        setMappingTab("qc");
      } else if (action === "sample" && settings.sample.saveReport) {
        setActiveSection("samples");
        setSampleTab("qc");
      }
    } else {
      setJob({
        label,
        state: "failed",
        progress: 0,
        subtask: response.message,
        startedAt: "now",
        finishedAt: "now",
      });
      if (!response.backendReachable) {
        setBackend((current) => ({
          ...current,
          connected: false,
          packageReady: false,
          message: "The local R backend did not answer. Retrying…",
        }));
      }
    }
    return response.success;
  }

  async function requestRunAction(action: "control" | "sample" | "af", label: string) {
    if (action !== "af" && workflowHasExistingResults(action, project)) {
      setPendingVersionRun({ action, label });
      return false;
    }
    return runAction(action, label);
  }

  async function loadPanel() {
    const payload = await loadPanelPayload("aurora");
    if (payload) {
      setPanelPayload(payload);
      appendExecutionLogs([{ kind: "success", text: "Panel library refreshed from the R backend." }]);
    } else {
      appendExecutionLogs([{ kind: "warning", text: "Showing the bundled panel library preview. The R backend is not connected." }]);
    }
  }

  async function confirmMapping() {
    const mapping = await persistControlMapping(project.mapping, project.projectPath);
    if (mapping.success) {
      setProject((current) => ({ ...current, mappingDirty: false }));
      setMappingTab("gating");
      openApplet("control-gating", "controls");
    }
    appendExecutionLogs([{ kind: mapping.success ? "success" : "error", text: mapping.message }]);
  }

  async function chooseProject(mode: "create" | "open") {
    if (projectLoading || projectPickerOpen) return;
    if (project.mappingDirty && !window.confirm("Discard the unsaved control mapping and switch projects?")) return;
    setProjectPickerOpen(true);
    let response: Awaited<ReturnType<typeof selectProjectFolder>>;
    try {
      response = mode === "create"
        ? await createProjectFolder()
        : await selectProjectFolder();
    } finally {
      setProjectPickerOpen(false);
    }
    if (response.cancelled) return;
    if (!response.success || !response.projectPath) {
      appendExecutionLogs([{ kind: "error", text: response.message }]);
      return;
    }

    const selectedProjectPath = response.projectPath;
    setProjectLoading(true);
    try {
      const activated = await setProjectContext(selectedProjectPath);
      if (!activated.success) {
        appendExecutionLogs([{ kind: "error", text: activated.message }]);
        return;
      }
      window.sessionStorage.setItem("spectreasy-project-path", selectedProjectPath);
      await refreshProject(true, selectedProjectPath);
      setInitializationMessage("");
      setInitializationDismissedFor("");
    } finally {
      setProjectLoading(false);
    }
  }

  async function confirmProjectInitialization() {
    if (!project.projectPath || initializationBusy) return;
    setInitializationBusy(true);
    setInitializationMessage("");
    const result = await initializeProject(project.projectPath);
    setInitializationMessage(result.message);
    if (result.success) {
      await refreshProject(true);
      setInitializationDismissedFor("");
    }
    setInitializationBusy(false);
  }

  async function createMapping() {
    const response = await createControlMapping(project.projectPath);
    appendExecutionLogs([{ kind: response.success ? "success" : "error", text: response.message }]);
    if (response.success) await refreshProject();
  }

  function downloadArtifact(artifact: Artifact) {
    const manifest = `Spectreasy artifact\n\nName: ${artifact.name}\nType: ${artifact.type}\nStatus: ${artifact.status}\nPath: ${artifact.path}\nUpdated: ${artifact.updated}\nRun: ${artifact.run ?? "user supplied"}\n`;
    downloadTextFile(`${artifact.name.replaceAll("/", "_")}.txt`, manifest);
  }

  function selectArtifact(artifact: Artifact) {
    if (artifact.group === "Reports" || artifact.group === "QC Metrics") {
      if (artifact.type === "Sample QC report") {
        setActiveSection("samples");
        setSampleTab("qc");
      } else {
        setActiveSection("controls");
        setMappingTab("qc");
      }
      return;
    }
    const sectionMap: Record<string, SectionId> = {
      Controls: "controls",
      Gates: "controls",
      Samples: "samples",
      Matrices: "matrix",
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
        matrix: "Matrix review",
        panel: "Panel builder",
        af: "AF library",
        settings: "Settings & logs",
      })[activeSection],
    [activeSection],
  );

  return (
    <div className="cockpit-app">
      <div className="cockpit-app-content" inert={projectIsLoading} aria-hidden={projectIsLoading}>
        <a className="skip-link" href="#cockpit-main">Skip to workflow</a>
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
        onGuide={() => setGuideOpen(true)}
        onSettings={toggleSettings}
        onTerminal={() => setTerminalOpen((value) => !value)}
        executionLogs={executionLogs}
      />
      <div className="app-body">
        <main className="main-area" id="cockpit-main" tabIndex={-1}>
          <WorkflowRail
            activeSection={activeSection}
            project={project}
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
              sampleTab={sampleTab}
              setSampleTab={setSampleTab}
              onUpdateMapping={updateMapping}
              onRun={requestRunAction}
              onRefresh={() => void refreshProject()}
              onDownload={downloadArtifact}
              onSelectArtifact={selectArtifact}
              onLoadPanel={() => void loadPanel()}
              panelPayload={panelPayload}
              onSave={() => void confirmMapping()}
              onCreateMapping={() => void createMapping()}
              onSectionChange={navigateToSection}
              settings={settings}
              onSettingsChange={updateSettings}
              onInputDirectoryChange={updateInputDirectory}
              onOpenApplet={(applet, reportPath) => openApplet(applet, activeSection, reportPath)}
            />
          </div>
        </main>
      </div>
      {initializedApplets.map((applet) => (
        <CockpitApplet
          key={applet}
          applet={applet}
          theme={settings.appearance.theme}
          project={project}
          projectPath={settings.projectPath}
          outputRoot={normalizeCockpitOutputRoot(settings.control.outputDir)}
          reportPath={selectedReportPath}
          active={activeApplet === applet}
          onExit={exitApplet}
        />
      ))}
      {activeApplet && !persistentDataApplets.includes(activeApplet) && (
        <CockpitApplet
          applet={activeApplet}
          theme={settings.appearance.theme}
          project={project}
          projectPath={settings.projectPath}
          outputRoot={normalizeCockpitOutputRoot(activeApplet === "sample-qc-report" ? settings.sample.outputDir : settings.control.outputDir)}
          reportPath={selectedReportPath}
          onExit={exitApplet}
        />
      )}
      {filesOpen && (
        <ProjectFilesDialog
          projectName={project.projectName}
          projectPath={project.projectPath}
          controlInputDir={project.controlInputDir}
          sampleInputDir={project.sampleInputDir}
          onClose={() => setFilesOpen(false)}
          onChanged={async () => {
            await refreshProject(false);
          }}
        />
      )}
      {guideOpen && <GuideDialog onClose={() => setGuideOpen(false)} />}
      {project.projectPath && project.missingInputDirs.length > 0 && initializationDismissedFor !== project.projectPath && (
        <ProjectInitializationDialog
          projectName={project.projectName}
          missingFolders={project.missingInputDirs}
          busy={initializationBusy}
          message={initializationMessage}
          onConfirm={() => void confirmProjectInitialization()}
          onCancel={() => {
            setInitializationMessage("");
            setInitializationDismissedFor(project.projectPath);
          }}
        />
      )}
      {pendingVersionRun && (
        <VersionRunDialog
          workflow={pendingVersionRun.action}
          onCancel={() => setPendingVersionRun(null)}
          onConfirm={() => {
            const pending = pendingVersionRun;
            setPendingVersionRun(null);
            void runAction(pending.action, pending.label);
          }}
        />
      )}
      </div>
      {projectIsLoading && (
        <div className="project-loading-overlay" role="status" aria-live="polite" aria-label="Loading project files">
          <div className="project-loading-card">
            <LoaderCircle className="is-spinning" size={36} aria-hidden="true" />
            <strong>Loading project files...</strong>
            <span>Preparing the selected project</span>
          </div>
        </div>
      )}
    </div>
  );
}
