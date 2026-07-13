import { lazy, Suspense, useEffect, useMemo, useRef, useState } from "react";
import {
  AlertCircle,
  ArrowRight,
  BarChart3,
  Beaker,
  Check,
  ChevronDown,
  CircleCheckBig,
  Download,
  FileText,
  FlaskConical,
  FolderOpen,
  GitCompare,
  Info,
  Layers3,
  LockKeyhole,
  Play,
  Plus,
  RefreshCcw,
  Save,
  Search,
  Settings2,
  ShieldCheck,
  Sparkles,
  WandSparkles,
  X,
} from "lucide-react";
import type {
  AfSettings,
  Artifact,
  BackendStatus,
  ControlSettings,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  Report,
  SectionId,
  SampleSettings,
  WorkflowSettings,
} from "../types";
import {
  applyAfProfile,
  attemptDiagnosticAction,
  deleteAfProfile,
  downloadBase64File,
  exportPanelOverview,
  importMatrixContent,
  importSampleContent,
  listAfProfiles,
  listMatrixFiles,
  listSampleFiles,
  loadMatrixFile,
  loadPanelMetrics,
  loadProjectReports,
  projectFileUrl,
  rowsFromBackend,
  saveMatrixFile,
} from "../api";
import { demoReports } from "../mockData";
import { StatusPill } from "../components/StatusPill";
import { AppearanceSettings } from "../components/AppearanceSettings";

const GateEditor = lazy(() => import("../../GatingGui.jsx"));

export type WorkflowWorkspaceProps = {
  project: ProjectState;
  backend: BackendStatus;
  job: Job;
  activeSection: SectionId;
  mappingTab: "mapping" | "gating" | "build" | "qc";
  setMappingTab: (tab: "mapping" | "gating" | "build" | "qc") => void;
  onUpdateMapping: (id: string, patch: Partial<MappingRow>) => void;
  onRun: (
    action: "control" | "sample" | "report" | "af",
    label: string,
  ) => void;
  onRefresh: () => void;
  onDownload: (artifact: Artifact) => void;
  onSelectArtifact: (artifact: Artifact) => void;
  onLoadPanel: () => void;
  panelPayload: PanelPayload | null;
  onSave: () => void;
  settings: WorkflowSettings;
  onSettingsChange: (
    section: "projectPath" | "control" | "sample" | "af" | "appearance",
    patch:
      | Partial<WorkflowSettings["control"]>
      | Partial<WorkflowSettings["sample"]>
      | Partial<WorkflowSettings["af"]>
      | Partial<WorkflowSettings["appearance"]>
      | { projectPath: string },
  ) => void;
  onOpenProject: (path: string) => void;
  onSaveMapping?: () => void;
};

function WorkspaceHeader({
  kicker,
  title,
  description,
  action,
  onAction,
}: {
  kicker: string;
  title: string;
  description: string;
  action?: string;
  onAction?: () => void;
}) {
  return (
    <div className="workspace-header">
      <div>
        <span className="eyebrow">{kicker}</span>
        <h1>{title}</h1>
        <p>{description}</p>
      </div>
      {action && onAction && (
        <button className="button button-ghost" onClick={onAction}>
          <RefreshCcw size={15} /> {action}
        </button>
      )}
    </div>
  );
}

function SectionTitle({
  eyebrow,
  title,
  note,
  action,
}: {
  eyebrow: string;
  title: string;
  note?: string;
  action?: React.ReactNode;
}) {
  return (
    <div className="section-title">
      <div>
        <span className="eyebrow">{eyebrow}</span>
        <h2>{title}</h2>
        {note && <p>{note}</p>}
      </div>
      {action}
    </div>
  );
}

function Metric({
  value,
  label,
  accent = "",
}: {
  value: string;
  label: string;
  accent?: string;
}) {
  return (
    <div className={`metric ${accent}`}>
      <strong>{value}</strong>
      <span>{label}</span>
    </div>
  );
}

function JobStrip({ job }: { job: Job }) {
  if (job.state === "idle") return null;
  const state =
    job.state === "complete"
      ? "complete"
      : job.state === "failed"
        ? "warning"
        : "ready";
  return (
    <div className={`job-strip job-${job.state}`}>
      <div className="job-icon">
        {job.state === "running" ? (
          <RefreshCcw size={16} className="spin" />
        ) : job.state === "complete" ? (
          <Check size={16} />
        ) : (
          <AlertCircle size={16} />
        )}
      </div>
      <div className="job-copy">
        <strong>{job.label}</strong>
        <span>{job.subtask}</span>
      </div>
      <div className="job-progress">
        <div className="progress-track">
          <span style={{ width: `${job.progress}%` }} />
        </div>
        <span>{job.progress}%</span>
      </div>
      <StatusPill
        state={state}
        label={
          job.state === "running"
            ? "Running"
            : job.state === "complete"
              ? "Complete"
              : "Failed"
        }
        compact
      />
    </div>
  );
}

function OverviewWorkspace({
  project,
  onSectionChange,
  onRefresh,
}: {
  project: ProjectState;
  onSectionChange: (section: SectionId) => void;
  onRefresh: () => void;
}) {
  const nextSteps = [
    {
      title: "Review 2 new samples",
      detail: "Detector set matches selected matrix",
      state: "warning",
      section: "samples" as SectionId,
    },
    {
      title: "Run sample unmixing",
      detail: "Ready · 28 files · Spectreasy method",
      state: "ready",
      section: "samples" as SectionId,
    },
    {
      title: "Regenerate sample QC",
      detail: "Previous report is stale",
      state: "stale",
      section: "reports" as SectionId,
    },
  ];
  return (
    <>
      <WorkspaceHeader
        kicker="Project overview"
        title={project.projectName}
        description="Review project inputs, outputs, and workflow status."
        action="Rescan project"
        onAction={onRefresh}
      />
      <div className="hero-grid">
        <div className="hero-card">
          <div className="hero-card-top">
            <span className="status-orb">
              <span />
            </span>
            <span className="eyebrow">Project scan · just now</span>
            <StatusPill state="complete" label="Ready for samples" compact />
          </div>
          <h2>Controls are ready for sample unmixing.</h2>
          <p>
            Reference matrix, gates, detector noise, and spectral variants are
            current. Sample unmixing can run with the selected method.
          </p>
          <div className="hero-actions">
            <button
              className="button button-primary"
              onClick={() => onSectionChange("samples")}
            >
              <Play size={14} fill="currentColor" /> Continue to samples
            </button>
            <button
              className="button button-ghost"
              onClick={() => onSectionChange("controls")}
            >
              Review controls <ArrowRight size={14} />
            </button>
          </div>
          <div className="hero-foot">
            <span>
              <ShieldCheck size={14} /> All calculations stay in R
            </span>
            <span>
              <LockKeyhole size={13} /> Local project
            </span>
          </div>
        </div>
        <div className="scan-card">
          <div className="scan-orbit">
            <div className="orbit-ring ring-one" />
            <div className="orbit-ring ring-two" />
            <div className="orbit-ring ring-three" />
            <div className="orbit-core">
              <span>94%</span>
              <small>ready</small>
            </div>
          </div>
          <div className="scan-copy">
            <span className="eyebrow">Artifact coverage</span>
            <strong>14 / 15</strong>
            <span>expected outputs present</span>
            <div className="scan-legend">
              <span>
                <i className="legend-dot current" /> Current
              </span>
              <span>
                <i className="legend-dot stale" /> Stale
              </span>
              <span>
                <i className="legend-dot user" /> Input
              </span>
            </div>
          </div>
        </div>
      </div>
      <div className="metrics-row">
        <Metric value="14" label="control files" accent="accent-teal" />
        <Metric value="28" label="sample files" />
        <Metric value="3" label="matrices" accent="accent-amber" />
        <Metric value="5" label="reports" accent="accent-coral" />
        <Metric value="1" label="attention" accent="accent-red" />
      </div>
      <div className="two-column-grid">
        <section className="surface-card next-card">
          <SectionTitle
            eyebrow="Suggested next"
            title="Recommended next steps"
            note="Actions are ordered by project dependencies."
          />
          <div className="next-list">
            {nextSteps.map((step, index) => (
              <button
                className="next-row"
                key={step.title}
                onClick={() => onSectionChange(step.section)}
              >
                <span className="next-index">0{index + 1}</span>
                <span className="next-row-copy">
                  <strong>{step.title}</strong>
                  <small>{step.detail}</small>
                </span>
                <StatusPill
                  state={step.state as "warning" | "ready" | "stale"}
                  compact
                />
                <ArrowRight size={15} />
              </button>
            ))}
          </div>
        </section>
        <section className="surface-card manifest-card">
          <SectionTitle
            eyebrow="Run manifest"
            title="run-014"
            note="Control workflow · today at 09:48"
          />
          <div className="manifest-list">
            <div>
              <span>Method</span>
              <strong>{project.method}</strong>
            </div>
            <div>
              <span>Cytometer</span>
              <strong>{project.cytometer}</strong>
            </div>
            <div>
              <span>Outputs</span>
              <strong>9 artifacts</strong>
            </div>
            <div>
              <span>Duration</span>
              <strong>04:12</strong>
            </div>
          </div>
          <button className="text-action">
            <FileText size={14} /> View run log <ArrowRight size={13} />
          </button>
        </section>
      </div>
    </>
  );
}

function MappingWorkspace({
  project,
  onUpdateMapping,
  onRun,
  onRefresh,
  mappingTab,
  setMappingTab,
  onSaveMapping,
  settings,
  onSettingsChange,
  onViewReports,
}: Pick<
  WorkflowWorkspaceProps,
  | "project"
  | "onUpdateMapping"
  | "onRun"
  | "onRefresh"
  | "mappingTab"
  | "setMappingTab"
  | "onSaveMapping"
  | "settings"
  | "onSettingsChange"
> & { onViewReports: () => void }) {
  const [selectedRows, setSelectedRows] = useState<string[]>([]);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const toggleRow = (id: string) =>
    setSelectedRows((rows) =>
      rows.includes(id) ? rows.filter((row) => row !== id) : [...rows, id],
    );
  return (
    <>
      <WorkspaceHeader
        kicker="Controls / control stage"
        title="Make the reference trustworthy"
        description="Map, gate, and unmix controls in one pass. Every downstream artifact is tied back to these choices."
      />
      <div className="subnav">
        <button
          className={mappingTab === "mapping" ? "is-active" : ""}
          onClick={() => setMappingTab("mapping")}
        >
          01 Mapping <StatusPill state="complete" compact />
        </button>
        <button
          className={mappingTab === "gating" ? "is-active" : ""}
          onClick={() => setMappingTab("gating")}
        >
          02 Gating <StatusPill state="complete" compact />
        </button>
        <button
          className={mappingTab === "build" ? "is-active" : ""}
          onClick={() => setMappingTab("build")}
        >
          03 Build reference <StatusPill state="complete" compact />
        </button>
        <button
          className={mappingTab === "qc" ? "is-active" : ""}
          onClick={() => setMappingTab("qc")}
        >
          04 Control QC <StatusPill state="complete" compact />
        </button>
      </div>
      {mappingTab === "mapping" && (
        <section className="surface-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Control mapping</span>
              <h2>
                fcs_mapping.csv{" "}
                <span className="dirty-marker">
                  {project.mappingDirty ? "· unsaved" : "· saved"}
                </span>
              </h2>
            </div>
            <div className="toolbar-actions">
              <button className="button button-ghost" onClick={onRefresh}>
                <FolderOpen size={14} /> Use existing
              </button>
              <button className="button button-ghost" onClick={onRefresh}>
                <Plus size={14} /> Create from folder
              </button>
              <button className="button button-primary" onClick={onSaveMapping}>
                <Save size={14} /> Save mapping
              </button>
            </div>
          </div>
          <div className="validation-banner">
            <div className="banner-icon">
              <CircleCheckBig size={17} />
            </div>
            <div>
              <strong>Mapping validated</strong>
              <span>
                {project.mapping.length} rows · editable here
              </span>
            </div>
            <button className="text-action" onClick={onRefresh}>
              Run validation <ArrowRight size={13} />
            </button>
          </div>
          <div className="mapping-tools">
            <label className="search-field compact-search">
              <Search size={14} />
              <input placeholder="Filter controls" />
              <kbd>⌘ F</kbd>
            </label>
            <div className="bulk-actions">
              {selectedRows.length > 0 && (
                <>
                  <span>{selectedRows.length} selected</span>
                  <button
                    className="text-action"
                    onClick={() =>
                      selectedRows.forEach((id) =>
                        onUpdateMapping(id, { controlType: "cell" }),
                      )
                    }
                  >
                    Set as cell
                  </button>
                  <button
                    className="text-action"
                    onClick={() =>
                      selectedRows.forEach((id) =>
                        onUpdateMapping(id, { universalNegative: true }),
                      )
                    }
                  >
                    Set universal negative
                  </button>
                </>
              )}
            </div>
            <span className="table-note">
              <Info size={13} /> Channel choices come from FCS metadata
            </span>
          </div>
          <div className="mapping-table-wrap">
            <table className="mapping-table">
              <thead>
                <tr>
                  <th className="check-col">
                    <input
                      type="checkbox"
                      aria-label="Select all mapping rows"
                      checked={selectedRows.length === project.mapping.length}
                      onChange={() =>
                        setSelectedRows(
                          selectedRows.length === project.mapping.length
                            ? []
                            : project.mapping.map((row) => row.id),
                        )
                      }
                    />
                  </th>
                  <th>File</th>
                  <th>Fluorophore</th>
                  <th>Marker</th>
                  <th>Channel</th>
                  <th>Type</th>
                  <th>Negative</th>
                  <th />
                </tr>
              </thead>
              <tbody>
                {project.mapping.map((row) => (
                  <tr key={row.id} className={row.warning ? "has-warning" : ""}>
                    <td>
                      <input
                        type="checkbox"
                        checked={selectedRows.includes(row.id)}
                        onChange={() => toggleRow(row.id)}
                        aria-label={`Select ${row.file}`}
                      />
                    </td>
                    <td>
                      <span className="file-cell">
                        <span className="file-mini-dot" />
                        {row.file}
                      </span>
                    </td>
                    <td>
                      <input
                        className="table-input"
                        value={row.fluorophore}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            fluorophore: event.target.value,
                          })
                        }
                      />
                    </td>
                    <td>
                      <input
                        className="table-input"
                        value={row.marker}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            marker: event.target.value,
                          })
                        }
                      />
                    </td>
                    <td>
                      <select
                        className="table-select"
                        value={row.channel}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            channel: event.target.value,
                          })
                        }
                      >
                        <option>{row.channel}</option>
                        <option>Violet 450 / B2</option>
                        <option>Blue 530 / B3</option>
                        <option>YellowGreen 586 / YG2</option>
                        <option>Red 670 / R4</option>
                      </select>
                    </td>
                    <td>
                      <select
                        className="type-select"
                        value={row.controlType}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            controlType: event.target
                              .value as MappingRow["controlType"],
                          })
                        }
                      >
                        <option value="cell">Cell SCC</option>
                        <option value="bead">Bead</option>
                        <option value="unstained">Unstained</option>
                        <option value="viability">Viability</option>
                      </select>
                    </td>
                    <td>
                      <input
                        type="checkbox"
                        checked={row.universalNegative}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            universalNegative: event.target.checked,
                          })
                        }
                        aria-label={`${row.file} universal negative`}
                      />
                    </td>
                    <td>
                      {row.warning ? (
                        <span className="row-warning" title={row.warning}>
                          <AlertCircle size={15} />
                        </span>
                      ) : (
                        <Check size={14} className="row-ok" />
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="card-foot">
            <span>
              <AlertCircle size={14} /> Mapping edits are saved directly to the
              project CSV.
            </span>
            <button
              className="text-action"
              onClick={() => setShowAdvanced(!showAdvanced)}
            >
              <Settings2 size={14} /> {showAdvanced ? "Hide" : "Show"}{" "}
              validation details <ChevronDown size={13} />
            </button>
          </div>
          {showAdvanced && (
            <div className="advanced-drawer">
              <div>
                <span className="eyebrow">Blocking checks</span>
                <strong>0 errors</strong>
                <small>
                  Every control has a readable FCS file and detector channel.
                </small>
              </div>
              <div>
                <span className="eyebrow">Advisories</span>
                <strong>Review before run</strong>
                <small>
                  R will validate the mapping against the active SCC files.
                </small>
              </div>
              <div>
                <span className="eyebrow">Detector set</span>
                <strong>Backend checked</strong>
                <small>The chosen channel is sent to the workflow.</small>
              </div>
            </div>
          )}
        </section>
      )}
      {mappingTab === "gating" && <GatingEmbed />}
      {mappingTab === "build" && (
        <BuildReferencePanel
          settings={settings.control}
          onSettingsChange={(patch) => onSettingsChange("control", patch)}
          onRun={onRun}
        />
      )}
      {mappingTab === "qc" && (
        <ControlReportPanel
          onRun={onRun}
          onView={onViewReports}
        />
      )}
    </>
  );
}

function GatingEmbed() {
  return (
    <section className="surface-card gating-embed">
      <div className="gating-embed-note">
        <span className="eyebrow">Live gate editor</span>
        <strong>
          Gate controls with the same CSV used by unmix_controls()
        </strong>
        <small>
          The editor reads real FCS events from the active project and saves the
          gate configuration directly.
        </small>
      </div>
      <Suspense
        fallback={
          <div className="gating-embed-loading">
            Loading the live gate editor…
          </div>
        }
      >
        <GateEditor />
      </Suspense>
    </section>
  );
}

function GatingPanel({ onRun }: { onRun: WorkflowWorkspaceProps["onRun"] }) {
  const [gateMode, setGateMode] = useState("Cell gate");
  return (
    <section className="surface-card gating-workspace">
      <div className="card-toolbar">
        <div>
          <span className="eyebrow">Control gating</span>
          <h2>Gate one control at a time</h2>
        </div>
        <div className="toolbar-actions">
          <button className="button button-ghost">
            <FolderOpen size={14} /> Load gates
          </button>
          <button
            className="button button-primary"
            onClick={() => onRun("control", "Save control gates")}
          >
            <Save size={14} /> Save gates
          </button>
        </div>
      </div>
      <div className="gating-layout">
        <div className="gating-file-column">
          <span className="eyebrow">Selected control</span>
          <button className="file-picker">
            <span>
              <CircleCheckBig size={15} /> unstained_cells.fcs
            </span>
            <ChevronDown size={15} />
          </button>
          <div className="gate-file-meta">
            <div>
              <span>Events</span>
              <strong>148,220</strong>
            </div>
            <div>
              <span>Returned</span>
              <strong>50,000</strong>
            </div>
            <div>
              <span>Type</span>
              <strong>Unstained</strong>
            </div>
          </div>
          <div className="gate-file-list">
            <button className="is-selected">
              <span className="artifact-dot dot-current" />
              unstained_cells.fcs
              <Check size={13} />
            </button>
            <button>
              <span className="artifact-dot dot-current" />
              bead_neg.fcs
            </button>
            <button>
              <span className="artifact-dot dot-current" />
              BV421.fcs
            </button>
            <button>
              <span className="artifact-dot dot-current" />
              FITC.fcs
            </button>
            <button>
              <span className="artifact-dot dot-warning" />
              AF-viability.fcs
              <AlertCircle size={13} />
            </button>
          </div>
        </div>
        <div className="gating-plot-card">
          <div className="plot-toolbar">
            <div className="plot-tabs">
              {["Cell gate", "Singlets", "Positive"].map((mode) => (
                <button
                  className={gateMode === mode ? "is-active" : ""}
                  key={mode}
                  onClick={() => setGateMode(mode)}
                >
                  {mode}
                </button>
              ))}
            </div>
            <span className="plot-caption">
              FSC-A × SSC-A · {gateMode.toLowerCase()}
            </span>
          </div>
          <svg
            className="gating-plot"
            viewBox="0 0 560 340"
            role="img"
            aria-label="Scatter plot preview with cell gate"
          >
            <defs>
              <linearGradient id="plotFade" x1="0" x2="1">
                <stop offset="0" stopColor="#7dc4bf" stopOpacity=".12" />
                <stop offset="1" stopColor="#db745a" stopOpacity=".1" />
              </linearGradient>
            </defs>
            <rect width="560" height="340" fill="url(#plotFade)" />
            {Array.from({ length: 84 }, (_, index) => {
              const x = 74 + ((index * 47) % 420);
              const y = 252 - ((index * 31) % 172);
              return (
                <circle
                  key={index}
                  cx={x}
                  cy={y}
                  r={index % 7 === 0 ? 3 : 2}
                  fill={index % 5 === 0 ? "#d98161" : "#78b9b2"}
                  opacity={0.35 + (index % 4) * 0.12}
                />
              );
            })}
            <polygon
              points="150,245 192,112 332,70 448,143 414,268 240,286"
              fill="#a2d6cf"
              fillOpacity=".13"
              stroke="#d58268"
              strokeWidth="2"
              strokeDasharray="6 5"
            />
            <line
              x1="64"
              y1="300"
              x2="520"
              y2="300"
              stroke="#8c938e"
              strokeWidth="1"
            />
            <line
              x1="64"
              y1="300"
              x2="64"
              y2="22"
              stroke="#8c938e"
              strokeWidth="1"
            />
            <text x="245" y="326" fill="#727b77" fontSize="12">
              FSC-A
            </text>
            <text
              x="14"
              y="168"
              fill="#727b77"
              fontSize="12"
              transform="rotate(-90 14,168)"
            >
              SSC-A
            </text>
          </svg>
          <div className="plot-foot">
            <span>
              <span className="plot-key key-in" /> In gate 82.4%
            </span>
            <span>
              <span className="plot-key key-out" /> Outliers 17.6%
            </span>
            <button className="text-action">
              Refresh preview <RefreshCcw size={13} />
            </button>
          </div>
        </div>
      </div>
      <div className="settings-strip">
        <label>
          Point size <input type="range" min="1" max="4" defaultValue="2" />
        </label>
        <label>
          Max points{" "}
          <select defaultValue="50000">
            <option>10,000</option>
            <option>50,000</option>
            <option>100,000</option>
          </select>
        </label>
        <label>
          Histogram transform{" "}
          <select defaultValue="asinh">
            <option>Asinh</option>
            <option>Linear</option>
            <option>Log10</option>
            <option>Biexponential</option>
          </select>
        </label>
        <span className="settings-note">
          <Info size={14} /> Gate edits are tracked and saved as CSV.
        </span>
      </div>
    </section>
  );
}

function LegacyBuildReferencePanel({
  settings,
  onSettingsChange,
  onRun,
}: {
  settings: ControlSettings;
  onSettingsChange: (patch: Partial<ControlSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const method = settings.method;
  const setMethod = (value: string) => {
    void value;
    return undefined;
  };
  void onSettingsChange;
  const [advanced, setAdvanced] = useState(false);
  return (
    <section className="surface-card run-card">
      <div className="card-toolbar">
        <div>
          <span className="eyebrow">Reference build / control unmixing</span>
          <h2>Turn controls into a trusted reference</h2>
          <p>
            One background job creates the matrix, noise floors, variants,
            unmixed controls, and QC report.
          </p>
        </div>
        <StatusPill state="complete" label="Ready to run" />
      </div>
      <div className="run-summary-grid">
        <div>
          <span className="eyebrow">Method</span>
          <strong>{method}</strong>
          <small>Practical default for spectral cytometry</small>
        </div>
        <div>
          <span className="eyebrow">Input</span>
          <strong>14 mapped controls</strong>
          <small>7 gates · 47 detectors</small>
        </div>
        <div>
          <span className="eyebrow">AF bank</span>
          <strong>Broad · 12 bands</strong>
          <small>Extract from unstained cells</small>
        </div>
        <div>
          <span className="eyebrow">Report</span>
          <strong>HTML</strong>
          <small>Shared QC data object</small>
        </div>
      </div>
      <div className="run-controls">
        <label>
          <span>Unmixing method</span>
          <select
            value={method}
            onChange={(event) => setMethod(event.target.value)}
          >
            <option>Spectreasy</option>
            <option>AutoSpectral</option>
            <option>OLS</option>
            <option>WLS</option>
            <option>RWLS</option>
            <option>NNLS</option>
          </select>
        </label>
        <label className="toggle-label">
          <input type="checkbox" defaultChecked />
          <span className="toggle-ui" />
          <span>Use manual gates</span>
        </label>
        <label className="toggle-label">
          <input type="checkbox" defaultChecked />
          <span className="toggle-ui" />
          <span>Generate reports</span>
        </label>
        <label className="toggle-label">
          <input type="checkbox" />
          <span className="toggle-ui" />
          <span>Overwrite existing outputs</span>
        </label>
      </div>
      <button
        className="advanced-toggle"
        onClick={() => setAdvanced(!advanced)}
      >
        <Settings2 size={15} /> Advanced settings{" "}
        <span>{advanced ? "−" : "+"}</span>
      </button>
      {advanced && (
        <div className="advanced-grid">
          <label>
            Weight quantile
            <input type="number" defaultValue="0.9" step="0.05" />
          </label>
          <label>
            SCC cleanup
            <select defaultValue="median">
              <option>Median background</option>
              <option>None</option>
              <option>Robust spline</option>
            </select>
          </label>
          <label>
            Variant top-k
            <input type="number" defaultValue="5" />
          </label>
          <label>
            Minimum cluster events
            <input type="number" defaultValue="120" />
          </label>
        </div>
      )}
      <div className="run-footer">
        <div className="run-output">
          <span className="output-dot" />
          <span>
            New versioned run folder will be created under{" "}
            <strong>Spectreasy/runs/</strong>
          </span>
        </div>
        <button
          className="button button-primary large-button"
          onClick={() => onRun("control", "Build reference & unmix controls")}
        >
          <Play size={15} fill="currentColor" /> Run control workflow{" "}
          <ArrowRight size={15} />
        </button>
      </div>
    </section>
  );
}

function BuildReferencePanel({
  settings,
  onSettingsChange,
  onRun,
}: {
  settings: ControlSettings;
  onSettingsChange: (patch: Partial<ControlSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const [advanced, setAdvanced] = useState(false);
  return (
    <section className="surface-card run-card">
      <div className="card-toolbar">
        <div>
          <span className="eyebrow">Reference build / control unmixing</span>
          <h2>Turn controls into a trusted reference</h2>
          <p>
            One background job creates the matrix, noise floors, variants,
            unmixed controls, and QC report.
          </p>
        </div>
        <StatusPill state="complete" label="Ready to run" />
      </div>
      <div className="run-summary-grid">
        <div>
          <span className="eyebrow">Method</span>
          <strong>{settings.method}</strong>
          <small>Selected in the control settings</small>
        </div>
        <div>
          <span className="eyebrow">Input</span>
          <strong>{settings.sccDir}</strong>
          <small>
            {settings.controlFile} · auto mapping{" "}
            {settings.autoCreateMapping ? "on" : "off"}
          </small>
        </div>
        <div>
          <span className="eyebrow">AF bank</span>
          <strong>{settings.afNBands} bands</strong>
          <small>{settings.afMinClusterEvents} minimum events</small>
        </div>
        <div>
          <span className="eyebrow">Report</span>
          <strong>
            {settings.saveReport
              ? `${settings.outputFormat.toUpperCase()} enabled`
              : "Report off"}
          </strong>
          <small>
            {settings.saveQcPlots ? "QC plots included" : "QC plots off"}
          </small>
        </div>
      </div>
      <div className="run-controls">
        <label>
          <span>Unmixing method</span>
          <select
            value={settings.method}
            onChange={(event) =>
              onSettingsChange({ method: event.target.value })
            }
          >
            <option>Spectreasy</option>
            <option>AutoSpectral</option>
            <option>OLS</option>
            <option>WLS</option>
            <option>RWLS</option>
            <option>NNLS</option>
          </select>
        </label>
        <label className="toggle-label">
          <input
            type="checkbox"
            checked={settings.gateFile.length > 0}
            onChange={(event) =>
              onSettingsChange({
                gateFile: event.target.checked ? "ssc_gate_config.csv" : "",
              })
            }
          />
          <span className="toggle-ui" />
          <span>Use saved gates</span>
        </label>
        <label className="toggle-label">
          <input
            type="checkbox"
            checked={settings.saveReport}
            onChange={(event) =>
              onSettingsChange({ saveReport: event.target.checked })
            }
          />
          <span className="toggle-ui" />
          <span>Generate report</span>
        </label>
        <label>
          <span>Report format</span>
          <select
            value={settings.outputFormat}
            onChange={(event) =>
              onSettingsChange({
                outputFormat: event.target.value as ControlSettings["outputFormat"],
              })
            }
          >
            <option value="html">HTML</option>
            <option value="pdf">PDF</option>
          </select>
        </label>
        <label className="toggle-label">
          <input
            type="checkbox"
            checked={settings.saveQcPlots}
            onChange={(event) =>
              onSettingsChange({ saveQcPlots: event.target.checked })
            }
          />
          <span className="toggle-ui" />
          <span>Save QC plots</span>
        </label>
      </div>
      <button
        className="advanced-toggle"
        onClick={() => setAdvanced(!advanced)}
      >
        <Settings2 size={15} /> Advanced settings{" "}
        <span>{advanced ? "−" : "+"}</span>
      </button>
      {advanced && (
        <div className="advanced-grid">
          <label>
            AF bands
            <input
              type="number"
              min="1"
              value={settings.afNBands}
              onChange={(event) =>
                onSettingsChange({ afNBands: Number(event.target.value) })
              }
            />
          </label>
          <label>
            AF min events
            <input
              type="number"
              min="1"
              value={settings.afMinClusterEvents}
              onChange={(event) =>
                onSettingsChange({
                  afMinClusterEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Background
            <select
              value={settings.sccBackgroundMethod}
              onChange={(event) =>
                onSettingsChange({
                  sccBackgroundMethod: event.target
                    .value as ControlSettings["sccBackgroundMethod"],
                })
              }
            >
              <option value="scatter_knn">Scatter KNN</option>
              <option value="none">None</option>
            </select>
          </label>
          <label>
            Seed
            <input
              type="number"
              min="1"
              value={settings.seed}
              onChange={(event) =>
                onSettingsChange({ seed: Number(event.target.value) })
              }
            />
          </label>
          <label>
            Threads
            <input
              type="number"
              min="1"
              value={settings.unmixThreads}
              onChange={(event) =>
                onSettingsChange({ unmixThreads: Number(event.target.value) })
              }
            />
          </label>
          <label>
            Variant top-k
            <input
              type="number"
              min="1"
              value={settings.spectralVariantTopK}
              onChange={(event) =>
                onSettingsChange({
                  spectralVariantTopK: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Weight quantile
            <input
              type="number"
              min="0"
              max="1"
              step="0.01"
              value={settings.spectreasyWeightQuantile}
              onChange={(event) =>
                onSettingsChange({
                  spectreasyWeightQuantile: Number(event.target.value),
                })
              }
            />
          </label>
        </div>
      )}
      <div className="run-footer">
        <div className="run-output">
          <span className="output-dot" />
          <span>
            Outputs will be written under <strong>{settings.outputDir}</strong>
          </span>
        </div>
        <button
          className="button button-primary large-button"
          onClick={() => onRun("control", "Build reference & unmix controls")}
        >
          <Play size={15} fill="currentColor" /> Run control workflow{" "}
          <ArrowRight size={15} />
        </button>
      </div>
    </section>
  );
}

function ControlReportPanel({
  onRun,
  onView,
}: {
  onRun: WorkflowWorkspaceProps["onRun"];
  onView: () => void;
}) {
  return (
    <section className="surface-card report-snapshot">
      <div className="card-toolbar">
        <div>
          <span className="eyebrow">Control QC / run-014</span>
          <h2>Controls look healthy</h2>
          <p>
            Generated from the same cached QC data used for the HTML and PDF
            report.
          </p>
        </div>
        <div className="toolbar-actions">
          <button className="button button-ghost" onClick={onView}>
            <ArrowRight size={14} /> View Control Report
          </button>
          <button
            className="button button-primary"
            onClick={() => onRun("report", "Render control QC report")}
          >
            <RefreshCcw size={14} /> Regenerate
          </button>
        </div>
      </div>
      <div className="qc-score-row">
        <div className="qc-score">
          <span>QC score</span>
          <strong>92</strong>
          <small>/ 100 · good</small>
        </div>
        <div className="qc-metric">
          <strong>47 / 47</strong>
          <span>detectors matched</span>
        </div>
        <div className="qc-metric">
          <strong>0.018</strong>
          <span>median residual</span>
        </div>
        <div className="qc-metric">
          <strong>96</strong>
          <span>variants retained</span>
        </div>
        <div className="qc-metric alert">
          <strong>1</strong>
          <span>advisory</span>
        </div>
      </div>
      <div className="report-mini-grid">
        <div>
          <span className="eyebrow">Sections</span>
          <div className="report-section-list">
            <span>
              <Check size={14} /> Mapping status
            </span>
            <span>
              <Check size={14} /> Gating plots
            </span>
            <span>
              <Check size={14} /> Detector noise
            </span>
            <span>
              <Check size={14} /> SCC comparisons
            </span>
          </div>
        </div>
        <div className="mini-spectrum">
          <span className="eyebrow">Reference signature</span>
          <svg viewBox="0 0 280 72">
            <path
              d="M0 54 C18 56 20 21 38 35 S60 58 75 28 S98 12 116 43 S145 55 158 18 S182 38 198 29 S218 42 232 15 S255 40 280 24"
              fill="none"
              stroke="#d27b63"
              strokeWidth="2.5"
            />
            <path
              d="M0 64 C35 58 47 59 71 52 S105 60 134 48 S169 64 197 50 S232 57 280 46"
              fill="none"
              stroke="#79b9b2"
              strokeWidth="1.5"
              opacity=".75"
            />
          </svg>
          <span>32 markers · 47 detectors</span>
        </div>
      </div>
    </section>
  );
}

function ControlsWorkspace(
  props: WorkflowWorkspaceProps & {
    onSectionChange: (section: SectionId) => void;
  },
) {
  return (
    <MappingWorkspace
      {...props}
      onViewReports={() => props.onSectionChange("reports")}
      onSaveMapping={props.onSaveMapping ?? props.onSave}
    />
  );
}

function SamplesWorkspace({
  project,
  onRun,
  onSectionChange,
}: {
  project: ProjectState;
  onRun: WorkflowWorkspaceProps["onRun"];
  onSectionChange: (section: SectionId) => void;
}) {
  const [sampleFilter, setSampleFilter] = useState("");
  const [selectedMatrix, setSelectedMatrix] = useState("reference_matrix.csv");
  const samples = [
    "PBMC_donor_01.fcs",
    "PBMC_donor_02.fcs",
    "PBMC_donor_03.fcs",
    "PBMC_donor_04.fcs",
    "PBMC_donor_05.fcs",
    "PBMC_donor_06.fcs",
  ];
  const filtered = samples.filter((sample) =>
    sample.toLowerCase().includes(sampleFilter.toLowerCase()),
  );
  return (
    <>
      <WorkspaceHeader
        kicker="Samples / sample stage"
        title="Unmix measured samples"
        description="Choose a current matrix, verify detector compatibility, and run the sample batch with the same backend methods."
      />
      <div className="sample-top-grid">
        <section className="surface-card sample-import-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Sample import</span>
              <h2>28 files detected</h2>
            </div>
            <button className="button button-ghost">
              <FolderOpen size={14} /> Add FCS files
            </button>
          </div>
          <div className="sample-filter-row">
            <label className="search-field compact-search">
              <Search size={14} />
              <input
                value={sampleFilter}
                onChange={(event) => setSampleFilter(event.target.value)}
                placeholder="Filter sample files"
              />
            </label>
            <span className="match-note">
              <Check size={13} /> 26 match detector set
            </span>
          </div>
          <div className="sample-list">
            {filtered.map((sample, index) => (
              <div className="sample-row" key={sample}>
                <span className="sample-index">
                  {String(index + 1).padStart(2, "0")}
                </span>
                <span className="file-mini-dot" />
                <strong>{sample}</strong>
                <span className="sample-events">
                  {(112 + index * 18).toFixed(0)}k events
                </span>
                <StatusPill
                  state={index > 3 ? "warning" : "complete"}
                  label={index > 3 ? "New" : "Matched"}
                  compact
                />
              </div>
            ))}
          </div>
          <button className="text-action">
            Show all 28 samples <ArrowRight size={13} />
          </button>
        </section>
        <section className="surface-card matrix-select-card">
          <SectionTitle
            eyebrow="Select matrix"
            title="Use the current reference"
            note="Reference matrices and adjusted matrices stay visibly distinct."
          />
          <label className="matrix-choice">
            <select
              value={selectedMatrix}
              onChange={(event) => setSelectedMatrix(event.target.value)}
            >
              <option>reference_matrix.csv</option>
              <option>adjusted_matrix_v2.csv</option>
              <option>unmixing_matrix.csv</option>
            </select>
            <ChevronDown size={14} />
          </label>
          <div className="matrix-detail-list">
            <div>
              <span>Dimensions</span>
              <strong>32 × 47</strong>
            </div>
            <div>
              <span>Noise file</span>
              <strong>detector_noise.csv</strong>
            </div>
            <div>
              <span>AF rows</span>
              <strong>12 bands</strong>
            </div>
            <div>
              <span>Source run</span>
              <strong>run-014</strong>
            </div>
          </div>
          <div className="matrix-current">
            <span className="artifact-dot dot-current" />
            <div>
              <strong>Current and trusted</strong>
              <small>Inputs have not changed since build.</small>
            </div>
            <button
              className="text-action"
              onClick={() => onSectionChange("matrix")}
            >
              Tune <ArrowRight size={13} />
            </button>
          </div>
        </section>
      </div>
      <section className="surface-card sample-run-card">
        <div className="card-toolbar">
          <div>
            <span className="eyebrow">Sample unmixing</span>
            <h2>Ready for a 28-file run</h2>
            <p>
              Writes to a new versioned folder; existing sample outputs remain
              safe.
            </p>
          </div>
          <StatusPill state="ready" label="Matrix selected" />
        </div>
        <div className="sample-run-options">
          <div>
            <span className="eyebrow">Method</span>
            <strong>{project.method}</strong>
            <small>Inherits control-stage method</small>
          </div>
          <div>
            <span className="eyebrow">Input</span>
            <strong>28 FCS files</strong>
            <small>26 matching · 2 new</small>
          </div>
          <div>
            <span className="eyebrow">QC</span>
            <strong>HTML</strong>
            <small>Report after unmix</small>
          </div>
          <div>
            <span className="eyebrow">Outputs</span>
            <strong>New run folder</strong>
            <small>Overwrite protection on</small>
          </div>
        </div>
        <div className="run-footer">
          <div className="warning-inline">
            <AlertCircle size={15} />
            <span>Sample QC report will be refreshed after this run.</span>
          </div>
          <button
            className="button button-ghost"
            onClick={() => onSectionChange("reports")}
          >
            View Sample Report <ArrowRight size={14} />
          </button>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("sample", "Unmix 28 samples")}
          >
            <Play size={15} fill="currentColor" /> Run sample workflow{" "}
            <ArrowRight size={15} />
          </button>
        </div>
      </section>
    </>
  );
}

function LegacyConfigurableSamplesWorkspace({
  project,
  settings,
  onSettingsChange,
  onRun,
  onSectionChange,
}: {
  project: ProjectState;
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
  onSectionChange: (section: SectionId) => void;
}) {
  const [sampleFilter, setSampleFilter] = useState("");
  const samples = [
    "PBMC_donor_01.fcs",
    "PBMC_donor_02.fcs",
    "PBMC_donor_03.fcs",
    "PBMC_donor_04.fcs",
    "PBMC_donor_05.fcs",
    "PBMC_donor_06.fcs",
  ];
  const filtered = samples.filter((sample) =>
    sample.toLowerCase().includes(sampleFilter.toLowerCase()),
  );
  return (
    <>
      <WorkspaceHeader
        kicker="Samples / sample stage"
        title="Unmix measured samples"
        description="Choose a current matrix, adjust the run settings, and send the batch to Spectreasy R."
      />
      <div className="sample-top-grid">
        <section className="surface-card sample-import-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Sample import</span>
              <h2>{project.scan.samples || 0} files detected</h2>
            </div>
            <button className="button button-ghost">
              <FolderOpen size={14} /> Add FCS files
            </button>
          </div>
          <div className="sample-filter-row">
            <label className="search-field compact-search">
              <Search size={14} />
              <input
                value={sampleFilter}
                onChange={(event) => setSampleFilter(event.target.value)}
                placeholder="Filter sample files"
              />
            </label>
            <span className="match-note">
              <Check size={13} /> Detector compatibility checked in R
            </span>
          </div>
          <div className="sample-list">
            {filtered.map((sample, index) => (
              <div className="sample-row" key={sample}>
                <span className="sample-index">
                  {String(index + 1).padStart(2, "0")}
                </span>
                <span className="file-mini-dot" />
                <strong>{sample}</strong>
                <span className="sample-events">
                  {(112 + index * 18).toFixed(0)}k events
                </span>
                <StatusPill
                  state={index > 3 ? "warning" : "complete"}
                  label={index > 3 ? "New" : "Matched"}
                  compact
                />
              </div>
            ))}
          </div>
          <button className="text-action">
            Show all detected samples <ArrowRight size={13} />
          </button>
        </section>
        <section className="surface-card matrix-select-card">
          <SectionTitle
            eyebrow="Run inputs"
            title="Use a reference matrix"
            note="Paths are relative to the active project folder."
          />
          <div className="settings-form">
            <label>
              Sample folder
              <input
                value={settings.sampleDir}
                onChange={(event) =>
                  onSettingsChange({ sampleDir: event.target.value })
                }
              />
            </label>
            <label>
              Reference matrix
              <input
                value={settings.matrixFile}
                onChange={(event) =>
                  onSettingsChange({ matrixFile: event.target.value })
                }
              />
            </label>
            <label>
              Detector noise file
              <input
                value={settings.detectorNoiseFile}
                onChange={(event) =>
                  onSettingsChange({ detectorNoiseFile: event.target.value })
                }
                placeholder="Optional"
              />
            </label>
          </div>
          <div className="matrix-current">
            <span className="artifact-dot dot-current" />
            <div>
              <strong>{settings.matrixFile}</strong>
              <small>Matrix path sent to the R backend</small>
            </div>
            <button
              className="text-action"
              onClick={() => onSectionChange("matrix")}
            >
              Tune <ArrowRight size={13} />
            </button>
          </div>
        </section>
      </div>
      <section className="surface-card sample-run-card">
        <div className="card-toolbar">
          <div>
            <span className="eyebrow">Sample unmixing</span>
            <h2>Ready for the next batch</h2>
            <p>
              Writes FCS outputs and the optional QC report to the configured
              output folder.
            </p>
          </div>
          <StatusPill state="ready" label={`${settings.method} · configured`} />
        </div>
        <div className="sample-run-options">
          <div>
            <span className="eyebrow">Method</span>
            <strong>{settings.method}</strong>
            <small>Selected in Settings & logs</small>
          </div>
          <div>
            <span className="eyebrow">Plots</span>
            <strong>{settings.plotNEvents.toLocaleString()} events</strong>
            <small>
              {settings.saveQcPlots ? "QC plots on" : "QC plots off"}
            </small>
          </div>
          <div>
            <span className="eyebrow">Chunk</span>
            <strong>{settings.chunkSize.toLocaleString()}</strong>
            <small>events per batch</small>
          </div>
          <div>
            <span className="eyebrow">Outputs</span>
            <strong>
              {settings.writeFcs
                ? `FCS + ${settings.outputFormat.toUpperCase()}`
                : `${settings.outputFormat.toUpperCase()} report only`}
            </strong>
            <small>{settings.outputDir}</small>
          </div>
        </div>
        <details className="settings-section">
          <summary>
            <Settings2 size={14} /> Advanced sample parameters
          </summary>
          <div className="settings-form-grid">
            <label>
              Unmixing method
              <select
                value={settings.method}
                onChange={(event) =>
                  onSettingsChange({ method: event.target.value })
                }
              >
                <option>Spectreasy</option>
                <option>AutoSpectral</option>
                <option>OLS</option>
                <option>WLS</option>
                <option>RWLS</option>
                <option>NNLS</option>
              </select>
            </label>
            <label>
              Threads
              <input
                type="number"
                min="1"
                value={settings.nThreads}
                onChange={(event) =>
                  onSettingsChange({ nThreads: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Plot events
              <input
                type="number"
                min="1"
                value={settings.plotNEvents}
                onChange={(event) =>
                  onSettingsChange({ plotNEvents: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Chunk size
              <input
                type="number"
                min="1"
                value={settings.chunkSize}
                onChange={(event) =>
                  onSettingsChange({ chunkSize: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Variant top-k
              <input
                type="number"
                min="1"
                value={settings.spectralVariantTopK}
                onChange={(event) =>
                  onSettingsChange({
                    spectralVariantTopK: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Weight quantile
              <input
                type="number"
                min="0"
                max="1"
                step="0.01"
                value={settings.spectreasyWeightQuantile}
                onChange={(event) =>
                  onSettingsChange({
                    spectreasyWeightQuantile: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Seed
              <input
                type="number"
                min="1"
                value={settings.seed}
                onChange={(event) =>
                  onSettingsChange({ seed: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Output folder
              <input
                value={settings.outputDir}
                onChange={(event) =>
                  onSettingsChange({ outputDir: event.target.value })
                }
              />
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.estimateAf}
                onChange={(event) =>
                  onSettingsChange({ estimateAf: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Estimate AF from samples</span>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.writeFcs}
                onChange={(event) =>
                  onSettingsChange({ writeFcs: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Write FCS outputs</span>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveReport}
                onChange={(event) =>
                  onSettingsChange({ saveReport: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Generate report</span>
            </label>
            <label>
              Report format
              <select
                value={settings.outputFormat}
                onChange={(event) =>
                  onSettingsChange({
                    outputFormat: event.target.value as SampleSettings["outputFormat"],
                  })
                }
              >
                <option value="html">HTML</option>
                <option value="pdf">PDF</option>
              </select>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveQcPlots}
                onChange={(event) =>
                  onSettingsChange({ saveQcPlots: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Save QC plots</span>
            </label>
          </div>
        </details>
        <div className="run-footer">
          <div className="warning-inline">
            <AlertCircle size={15} />
            <span>
              R will validate the matrix, detector noise, and sample folder
              before writing.
            </span>
          </div>
          <button
            className="button button-ghost"
            onClick={() => onSectionChange("reports")}
          >
            View Sample Report <ArrowRight size={14} />
          </button>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("sample", "Unmix sample workflow")}
          >
            <Play size={15} fill="currentColor" /> Run sample workflow{" "}
            <ArrowRight size={15} />
          </button>
        </div>
      </section>
    </>
  );
}

function ConfigurableSamplesWorkspace({
  project,
  settings,
  onSettingsChange,
  onRun,
  onSectionChange,
}: {
  project: ProjectState;
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
  onSectionChange: (section: SectionId) => void;
}) {
  const [sampleFiles, setSampleFiles] = useState<string[]>([]);
  const [sampleFilter, setSampleFilter] = useState("");
  const [importStatus, setImportStatus] = useState(
    "Loading sample files from the R backend…",
  );
  const uploadRef = useRef<HTMLInputElement>(null);

  const refreshSamples = async () => {
    const files = await listSampleFiles();
    setSampleFiles(files);
    setImportStatus(
      files.length
        ? `${files.length} FCS file${files.length === 1 ? "" : "s"} available in R.`
        : "No FCS files found in the sample folder.",
    );
  };

  useEffect(() => {
    const timer = window.setTimeout(() => {
      void refreshSamples();
    }, 0);
    return () => window.clearTimeout(timer);
  }, []);

  const importFile = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    const bytes = new Uint8Array(await file.arrayBuffer());
    let binary = "";
    bytes.forEach((byte) => {
      binary += String.fromCharCode(byte);
    });
    const imported = await importSampleContent(file.name, window.btoa(binary));
    setImportStatus(
      imported
        ? `${file.name} imported into the active sample folder.`
        : `Could not import ${file.name}.`,
    );
    if (imported) void refreshSamples();
  };

  const filtered = sampleFiles.filter((sample) =>
    sample.toLowerCase().includes(sampleFilter.toLowerCase()),
  );
  return (
    <>
      <WorkspaceHeader
        kicker="Samples / sample stage"
        title="Unmix measured samples"
        description="Choose a current matrix, import FCS files into the active project, and send the batch to Spectreasy R."
      />
      <input
        ref={uploadRef}
        type="file"
        accept=".fcs,application/octet-stream"
        hidden
        onChange={importFile}
      />
      <div className="sample-top-grid">
        <section className="surface-card sample-import-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Sample import</span>
              <h2>
                {sampleFiles.length || project.scan.samples || 0} files detected
              </h2>
            </div>
            <div className="toolbar-actions">
              <button
                className="button button-ghost"
                onClick={() => uploadRef.current?.click()}
              >
                <FolderOpen size={14} /> Add FCS files
              </button>
              <button
                className="button button-ghost"
                onClick={() => void refreshSamples()}
              >
                <RefreshCcw size={14} /> Refresh
              </button>
            </div>
          </div>
          <div className="sample-filter-row">
            <label className="search-field compact-search">
              <Search size={14} />
              <input
                value={sampleFilter}
                onChange={(event) => setSampleFilter(event.target.value)}
                placeholder="Filter sample files"
              />
            </label>
            <span className="match-note">
              <Check size={13} /> {importStatus}
            </span>
          </div>
          <div className="sample-list">
            {filtered.map((sample, index) => (
              <div className="sample-row" key={sample}>
                <span className="sample-index">
                  {String(index + 1).padStart(2, "0")}
                </span>
                <span className="file-mini-dot" />
                <strong>{sample}</strong>
                <span className="sample-events">FCS</span>
                <StatusPill state="complete" label="Available" compact />
              </div>
            ))}
          </div>
          {sampleFiles.length > filtered.length && (
            <span className="table-note">
              Showing {filtered.length} of {sampleFiles.length} files.
            </span>
          )}
        </section>
        <section className="surface-card matrix-select-card">
          <SectionTitle
            eyebrow="Run inputs"
            title="Use a reference matrix"
            note="Paths are relative to the active project folder."
          />
          <div className="settings-form">
            <label>
              Sample folder
              <input
                value={settings.sampleDir}
                onChange={(event) =>
                  onSettingsChange({ sampleDir: event.target.value })
                }
              />
            </label>
            <label>
              Reference matrix
              <input
                value={settings.matrixFile}
                onChange={(event) =>
                  onSettingsChange({ matrixFile: event.target.value })
                }
              />
            </label>
            <label>
              Detector noise file
              <input
                value={settings.detectorNoiseFile}
                onChange={(event) =>
                  onSettingsChange({ detectorNoiseFile: event.target.value })
                }
                placeholder="Optional"
              />
            </label>
          </div>
          <div className="matrix-current">
            <span className="artifact-dot dot-current" />
            <div>
              <strong>{settings.matrixFile}</strong>
              <small>Matrix path sent to the R backend</small>
            </div>
            <button
              className="text-action"
              onClick={() => onSectionChange("matrix")}
            >
              Tune <ArrowRight size={13} />
            </button>
          </div>
        </section>
      </div>
      <section className="surface-card sample-run-card">
        <div className="card-toolbar">
          <div>
            <span className="eyebrow">Sample unmixing</span>
            <h2>Ready for the next batch</h2>
            <p>
              Writes FCS outputs and the optional QC report to the configured
              output folder.
            </p>
          </div>
          <StatusPill state="ready" label={`${settings.method} · configured`} />
        </div>
        <div className="sample-run-options">
          <div>
            <span className="eyebrow">Method</span>
            <strong>{settings.method}</strong>
            <small>Selected in Settings & logs</small>
          </div>
          <div>
            <span className="eyebrow">Plots</span>
            <strong>{settings.plotNEvents.toLocaleString()} events</strong>
            <small>
              {settings.saveQcPlots ? "QC plots on" : "QC plots off"}
            </small>
          </div>
          <div>
            <span className="eyebrow">Chunk</span>
            <strong>{settings.chunkSize.toLocaleString()}</strong>
            <small>events per batch</small>
          </div>
          <div>
            <span className="eyebrow">Outputs</span>
            <strong>
              {settings.writeFcs
                ? `FCS + ${settings.outputFormat.toUpperCase()}`
                : `${settings.outputFormat.toUpperCase()} report only`}
            </strong>
            <small>{settings.outputDir}</small>
          </div>
        </div>
        <details className="settings-section">
          <summary>
            <Settings2 size={14} /> Advanced sample parameters
          </summary>
          <div className="settings-form-grid">
            <label>
              Unmixing method
              <select
                value={settings.method}
                onChange={(event) =>
                  onSettingsChange({ method: event.target.value })
                }
              >
                <option>Spectreasy</option>
                <option>AutoSpectral</option>
                <option>OLS</option>
                <option>WLS</option>
                <option>RWLS</option>
                <option>NNLS</option>
              </select>
            </label>
            <label>
              Threads
              <input
                type="number"
                min="1"
                value={settings.nThreads}
                onChange={(event) =>
                  onSettingsChange({ nThreads: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Plot events
              <input
                type="number"
                min="1"
                value={settings.plotNEvents}
                onChange={(event) =>
                  onSettingsChange({ plotNEvents: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Chunk size
              <input
                type="number"
                min="1"
                value={settings.chunkSize}
                onChange={(event) =>
                  onSettingsChange({ chunkSize: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Variant top-k
              <input
                type="number"
                min="1"
                value={settings.spectralVariantTopK}
                onChange={(event) =>
                  onSettingsChange({
                    spectralVariantTopK: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Weight quantile
              <input
                type="number"
                min="0"
                max="1"
                step="0.01"
                value={settings.spectreasyWeightQuantile}
                onChange={(event) =>
                  onSettingsChange({
                    spectreasyWeightQuantile: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Seed
              <input
                type="number"
                min="1"
                value={settings.seed}
                onChange={(event) =>
                  onSettingsChange({ seed: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Output folder
              <input
                value={settings.outputDir}
                onChange={(event) =>
                  onSettingsChange({ outputDir: event.target.value })
                }
              />
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.estimateAf}
                onChange={(event) =>
                  onSettingsChange({ estimateAf: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Estimate AF from samples</span>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.writeFcs}
                onChange={(event) =>
                  onSettingsChange({ writeFcs: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Write FCS outputs</span>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveReport}
                onChange={(event) =>
                  onSettingsChange({ saveReport: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Generate report</span>
            </label>
            <label>
              Report format
              <select
                value={settings.outputFormat}
                onChange={(event) =>
                  onSettingsChange({
                    outputFormat: event.target.value as SampleSettings["outputFormat"],
                  })
                }
              >
                <option value="html">HTML</option>
                <option value="pdf">PDF</option>
              </select>
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveQcPlots}
                onChange={(event) =>
                  onSettingsChange({ saveQcPlots: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Save QC plots</span>
            </label>
          </div>
        </details>
        <div className="run-footer">
          <div className="warning-inline">
            <AlertCircle size={15} />
            <span>
              R will validate the matrix, detector noise, and sample folder
              before writing.
            </span>
          </div>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("sample", "Unmix sample workflow")}
          >
            <Play size={15} fill="currentColor" /> Run sample workflow{" "}
            <ArrowRight size={15} />
          </button>
        </div>
      </section>
    </>
  );
}

function LegacyReportsWorkspace({
  onDownload,
  onRun,
}: {
  onDownload: WorkflowWorkspaceProps["onDownload"];
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const [selectedReport, setSelectedReport] = useState<Report>(demoReports[0]);
  const grouped = useMemo(
    () => ({
      current: demoReports.filter((report) => report.status === "current"),
      stale: demoReports.filter((report) => report.status === "stale"),
    }),
    [],
  );
  return (
    <>
      <WorkspaceHeader
        kicker="Reports / shared QC data"
        title="Read the run, not the folder tree"
        description="HTML is the in-app view. PDF stays available for archival review, with both formats rendered from one cached report object."
      />
      <div className="reports-layout">
        <section className="surface-card report-library">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Report library</span>
              <h2>5 generated reports</h2>
            </div>
            <button
              className="button button-ghost"
              onClick={() => onRun("report", "Render reports")}
            >
              <RefreshCcw size={14} /> Render reports
            </button>
          </div>
          <div className="report-group">
            <span className="eyebrow">Current</span>
            {grouped.current.map((report) => (
              <button
                className={`report-row ${selectedReport.id === report.id ? "is-selected" : ""}`}
                key={report.id}
                onClick={() => setSelectedReport(report)}
              >
                <span className="report-format">{report.format}</span>
                <span>
                  <strong>{report.title}</strong>
                  <small>
                    {report.type} · {report.created}
                  </small>
                </span>
                <ArrowRight size={14} />
              </button>
            ))}
          </div>
          <div className="report-group stale-group">
            <span className="eyebrow">Needs attention</span>
            {grouped.stale.map((report) => (
              <button
                className={`report-row ${selectedReport.id === report.id ? "is-selected" : ""}`}
                key={report.id}
                onClick={() => setSelectedReport(report)}
              >
                <span className="report-format report-format-stale">
                  {report.format}
                </span>
                <span>
                  <strong>{report.title}</strong>
                  <small>{report.type} · source changed</small>
                </span>
                <StatusPill state="stale" compact />
              </button>
            ))}
          </div>
        </section>
        <section className="surface-card report-viewer">
          <div className="report-viewer-bar">
            <div>
              <span className="eyebrow">Embedded viewer</span>
              <h2>{selectedReport.title}</h2>
            </div>
            <div className="toolbar-actions">
              <button
                className="button button-ghost"
                onClick={() =>
                  onDownload({
                    id: selectedReport.id,
                    name: `${selectedReport.title.replaceAll(" ", "_")}.${selectedReport.format.toLowerCase()}`,
                    type: "QC report",
                    group: "Reports",
                    status: selectedReport.status,
                    detail: selectedReport.format,
                    path: "/reports",
                    updated: selectedReport.created,
                  })
                }
              >
                <Download size={14} /> Download
              </button>
              <button className="button button-ghost">
                Open in browser <ArrowRight size={14} />
              </button>
            </div>
          </div>
          <div className="report-page">
            <div className="report-page-head">
              <span className="report-kicker">SPECTREASY · CONTROL QC</span>
              <h3>Control quality report</h3>
              <p>Lymphocyte panel · run-014 · Cytek Aurora 5L</p>
              <div className="report-badges">
                <StatusPill
                  state={
                    selectedReport.status === "stale" ? "stale" : "complete"
                  }
                />
                <span>Method · Spectreasy</span>
                <span>Generated today, 09:47</span>
              </div>
            </div>
            <div className="report-page-summary">
              <div className="report-score">
                <small>Overall</small>
                <strong>
                  {selectedReport.status === "stale" ? "—" : "92"}
                </strong>
                <span>
                  {selectedReport.status === "stale" ? "Regenerate" : "good"}
                </span>
              </div>
              <div className="report-bars">
                <div>
                  <span>Mapping</span>
                  <i>
                    <b style={{ width: "100%" }} />
                  </i>
                  <strong>100</strong>
                </div>
                <div>
                  <span>Gating</span>
                  <i>
                    <b style={{ width: "88%" }} />
                  </i>
                  <strong>88</strong>
                </div>
                <div>
                  <span>Unmixing</span>
                  <i>
                    <b style={{ width: "94%" }} />
                  </i>
                  <strong>94</strong>
                </div>
                <div>
                  <span>Detector noise</span>
                  <i>
                    <b style={{ width: "96%" }} />
                  </i>
                  <strong>96</strong>
                </div>
              </div>
            </div>
            <div className="report-section-preview">
              <div>
                <span className="eyebrow">01 / Mapping status</span>
                <p>
                  All 14 control rows map to an available file and detector
                  channel.
                </p>
              </div>
              <div>
                <span className="eyebrow">02 / Spectra & residuals</span>
                <p>
                  Reference spectra remain within expected residual bounds for
                  this panel.
                </p>
              </div>
              <div>
                <span className="eyebrow">03 / Advisory</span>
                <p className="advisory-copy">
                  <AlertCircle size={14} /> Viability control has no matched
                  negative.
                </p>
              </div>
            </div>
          </div>
        </section>
      </div>
    </>
  );
}

function LegacyMatrixWorkspace({
  onRun,
}: {
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const [selectedRow, setSelectedRow] = useState("CD45");
  const rows = [
    ["CD45", "0.82", "0.14", "0.06", "0.02"],
    ["CD3", "0.08", "0.74", "0.12", "0.06"],
    ["CD19", "0.04", "0.11", "0.79", "0.06"],
    ["CD8", "0.03", "0.12", "0.18", "0.67"],
    ["AF_01", "0.21", "0.17", "0.11", "0.09"],
  ];
  return (
    <>
      <WorkspaceHeader
        kicker="Matrix review / adjustment"
        title="Tune with context"
        description="Inspect detector signatures and save an adjusted matrix as a new artifact. The R backend remains the source of every value."
      />
      <div className="matrix-header-card surface-card">
        <div>
          <span className="eyebrow">Selected matrix</span>
          <h2>reference_matrix.csv</h2>
          <p>32 markers × 47 detectors · run-014 · current</p>
        </div>
        <div className="toolbar-actions">
          <button className="button button-ghost">Import matrix</button>
          <button
            className="button button-primary"
            onClick={() => onRun("control", "Save adjusted matrix")}
          >
            <Save size={14} /> Save adjusted copy
          </button>
        </div>
      </div>
      <div className="matrix-layout">
        <section className="surface-card matrix-table-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Marker signatures</span>
              <h2>Detector response</h2>
            </div>
            <span className="table-note">
              Values are normalized by the backend
            </span>
          </div>
          <table className="matrix-preview">
            <thead>
              <tr>
                <th>Marker</th>
                <th>Violet</th>
                <th>Blue</th>
                <th>YG</th>
                <th>Red</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr
                  className={selectedRow === row[0] ? "is-selected" : ""}
                  key={row[0]}
                  onClick={() => setSelectedRow(row[0])}
                >
                  <td>
                    <strong>{row[0]}</strong>
                  </td>
                  {row.slice(1).map((value, index) => (
                    <td key={`${row[0]}-${index}`}>
                      <span
                        className="matrix-cell"
                        style={{ opacity: 0.35 + Number(value) * 0.8 }}
                      >
                        {value}
                      </span>
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
          <button className="text-action">
            Show all 32 markers <ArrowRight size={13} />
          </button>
        </section>
        <section className="surface-card tuner-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Adjustment preview</span>
              <h2>{selectedRow}</h2>
            </div>
            <StatusPill state="ready" label="Unsaved edits" />
          </div>
          <div className="signature-chart">
            <svg viewBox="0 0 420 160">
              <path
                d="M0 134 C18 131 24 98 39 109 S62 135 78 72 S97 42 116 96 S145 121 160 80 S180 52 198 76 S216 118 231 49 S248 79 265 65 S286 82 301 32 S321 89 336 64 S366 56 379 24 S400 61 420 45"
                fill="none"
                stroke="#d27860"
                strokeWidth="3"
              />
              <path
                d="M0 145 C55 138 83 132 121 140 S192 125 232 139 S300 122 344 134 S383 127 420 126"
                fill="none"
                stroke="#79bdb5"
                strokeWidth="2"
              />
              <line x1="0" y1="145" x2="420" y2="145" stroke="#d8d2c5" />
            </svg>
            <div className="chart-axis">
              <span>Violet</span>
              <span>47 detector channels</span>
              <span>Red</span>
            </div>
          </div>
          <div className="tuner-controls">
            <label>
              Blend toward source
              <input type="range" min="0" max="100" defaultValue="18" />
            </label>
            <label>
              Outlier threshold
              <input type="range" min="0" max="100" defaultValue="72" />
            </label>
          </div>
          <div className="tuner-note">
            <Info size={14} />
            <span>
              Adjustments are tracked separately and will mark sample outputs
              stale after save.
            </span>
          </div>
        </section>
      </div>
    </>
  );
}

function ReportsWorkspace({
  onRun,
}: {
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const [reports, setReports] = useState<Report[]>(demoReports);
  const [selectedId, setSelectedId] = useState(demoReports[0].id);

  useEffect(() => {
    const timer = window.setTimeout(() => {
      void loadProjectReports().then((liveReports) => {
        if (liveReports.length > 0) {
          setReports(liveReports);
          setSelectedId(liveReports[0].id);
        }
      });
    }, 0);
    return () => window.clearTimeout(timer);
  }, []);

  const selectedReport =
    reports.find((report) => report.id === selectedId) ??
    reports[0] ??
    demoReports[0];
  const grouped = useMemo(
    () => ({
      current: reports.filter((report) => report.status === "current"),
      stale: reports.filter((report) => report.status === "stale"),
    }),
    [reports],
  );
  const openReport = () => {
    if (selectedReport.path)
      window.open(
        projectFileUrl(selectedReport.path),
        "_blank",
        "noopener,noreferrer",
      );
  };

  return (
    <>
      <WorkspaceHeader
        kicker="Reports / shared QC data"
        title="Read the run, not the folder tree"
        description="Reports are discovered from the active project. Open the original HTML or PDF through the local R backend, or render a fresh control report."
      />
      <div className="reports-layout">
        <section className="surface-card report-library">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Report library</span>
              <h2>{reports.length} generated reports</h2>
            </div>
            <button
              className="button button-ghost"
              onClick={() => onRun("report", "Render reports")}
            >
              <RefreshCcw size={14} /> Render reports
            </button>
          </div>
          <div className="report-group">
            <span className="eyebrow">Current</span>
            {grouped.current.length ? (
              grouped.current.map((report) => (
                <button
                  className={`report-row ${selectedReport.id === report.id ? "is-selected" : ""}`}
                  key={report.id}
                  onClick={() => setSelectedId(report.id)}
                >
                  <span className="report-format">{report.format}</span>
                  <span>
                    <strong>{report.title}</strong>
                    <small>
                      {report.type} · {report.created}
                    </small>
                  </span>
                  <ArrowRight size={14} />
                </button>
              ))
            ) : (
              <p className="empty-report-note">
                No current HTML or PDF artifacts found yet.
              </p>
            )}
          </div>
          <div className="report-group stale-group">
            <span className="eyebrow">Needs attention</span>
            {grouped.stale.length ? (
              grouped.stale.map((report) => (
                <button
                  className={`report-row ${selectedReport.id === report.id ? "is-selected" : ""}`}
                  key={report.id}
                  onClick={() => setSelectedId(report.id)}
                >
                  <span className="report-format report-format-stale">
                    {report.format}
                  </span>
                  <span>
                    <strong>{report.title}</strong>
                    <small>{report.type} · source changed</small>
                  </span>
                  <StatusPill state="stale" compact />
                </button>
              ))
            ) : (
              <p className="empty-report-note">No stale reports detected.</p>
            )}
          </div>
        </section>
        <section className="surface-card report-viewer">
          <div className="report-viewer-bar">
            <div>
              <span className="eyebrow">Artifact viewer</span>
              <h2>{selectedReport.title}</h2>
            </div>
            <div className="toolbar-actions">
              <a
                className="button button-ghost"
                href={selectedReport.path ? projectFileUrl(selectedReport.path) : undefined}
                download={selectedReport.path?.split("/").pop()}
                aria-disabled={!selectedReport.path}
                onClick={(event) => {
                  if (!selectedReport.path) event.preventDefault();
                }}
              >
                <Download size={14} /> Download
              </a>
              <button
                className="button button-ghost"
                onClick={openReport}
                disabled={!selectedReport.path}
              >
                Open in browser <ArrowRight size={14} />
              </button>
            </div>
          </div>
          {selectedReport.format === "HTML" && selectedReport.path ? (
            <iframe
              className="report-viewer-frame"
              src={projectFileUrl(selectedReport.path)}
              title={selectedReport.title}
              sandbox="allow-scripts allow-same-origin"
            />
          ) : (
          <div className="report-page">
            <div className="report-page-head">
              <span className="report-kicker">
                SPECTREASY · PROJECT ARTIFACT
              </span>
              <h3>{selectedReport.title}</h3>
              <p>
                {selectedReport.type} · {selectedReport.format} ·{" "}
                {selectedReport.path ?? "Preview artifact"}
              </p>
              <div className="report-badges">
                <StatusPill
                  state={
                    selectedReport.status === "stale" ? "stale" : "complete"
                  }
                />
                <span>Source · R project</span>
                <span>{selectedReport.created}</span>
              </div>
            </div>
            <div className="report-page-summary">
              <div className="report-score">
                <small>Format</small>
                <strong>{selectedReport.format}</strong>
                <span>
                  {selectedReport.status === "stale"
                    ? "needs refresh"
                    : "available"}
                </span>
              </div>
              <div className="report-bars">
                <div>
                  <span>Artifact exists</span>
                  <i>
                    <b style={{ width: selectedReport.path ? "100%" : "0%" }} />
                  </i>
                  <strong>{selectedReport.path ? "yes" : "no"}</strong>
                </div>
                <div>
                  <span>Project-linked</span>
                  <i>
                    <b style={{ width: selectedReport.path ? "100%" : "0%" }} />
                  </i>
                  <strong>{selectedReport.path ? "yes" : "no"}</strong>
                </div>
                <div>
                  <span>QC state</span>
                  <i>
                    <b
                      style={{
                        width:
                          selectedReport.status === "stale" ? "55%" : "94%",
                      }}
                    />
                  </i>
                  <strong>{selectedReport.status}</strong>
                </div>
              </div>
            </div>
            <div className="report-section-preview">
              <div>
                <span className="eyebrow">Selected path</span>
                <p>
                  {selectedReport.path ??
                    "The preview will be replaced when the R backend returns a report artifact."}
                </p>
              </div>
              <div>
                <span className="eyebrow">Next action</span>
                <p>
                  {selectedReport.status === "stale"
                    ? "Render reports after the upstream workflow changes."
                    : "Open the original artifact in a browser tab for full-resolution review."}
                </p>
              </div>
              <div>
                <span className="eyebrow">Data ownership</span>
                <p className="advisory-copy">
                  <ShieldCheck size={14} /> The file remains local to the active
                  Spectreasy project.
                </p>
              </div>
            </div>
          </div>
          )}
        </section>
      </div>
    </>
  );
}

function MatrixWorkspace() {
  const [files, setFiles] = useState<string[]>([]);
  const [filename, setFilename] = useState("");
  const [sourceFilename, setSourceFilename] = useState("");
  const [rows, setRows] = useState<Array<Record<string, unknown>>>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [status, setStatus] = useState("Loading matrices from the R backend…");
  const uploadRef = useRef<HTMLInputElement>(null);

  const load = async (nextFilename: string) => {
    setStatus(`Loading ${nextFilename}…`);
    const result = await loadMatrixFile(nextFilename);
    if (!result) {
      setStatus(`Could not load ${nextFilename}.`);
      return;
    }
    setFilename(result.filename);
    setSourceFilename(result.filename);
    setRows(result.rows);
    setSelectedIndex(0);
    setStatus(
      `${result.rows.length.toLocaleString()} markers loaded from the R backend.`,
    );
  };

  useEffect(() => {
    void listMatrixFiles().then((nextFiles) => {
      setFiles(nextFiles);
      if (nextFiles.length > 0) void load(nextFiles[0]);
      else setStatus("No matrix CSV files found. Import one to begin.");
    });
  }, []);

  const columns = rows.length > 0 ? Object.keys(rows[0]) : [];
  const markerColumn = columns[0] || "Marker";
  const detectorColumns = columns.slice(1, 5);
  const selectedRow = rows[selectedIndex] || {};
  const selectedMarker = String(
    selectedRow[markerColumn] ?? "No marker selected",
  );
  const visibleRows = rows.slice(0, 32);

  const updateCell = (rowIndex: number, column: string, value: string) => {
    setRows((current) =>
      current.map((row, index) => {
        if (index !== rowIndex) return row;
        if (column === markerColumn) return { ...row, [column]: value };
        const numeric = Number(value);
        return {
          ...row,
          [column]:
            value.trim() === "" || Number.isNaN(numeric) ? value : numeric,
        };
      }),
    );
    setStatus("Unsaved matrix edits");
  };

  const save = async () => {
    if (!filename || rows.length === 0) return;
    const saved = await saveMatrixFile(filename, rows, sourceFilename);
    setStatus(
      saved
        ? `${filename} saved through the R backend.`
        : "The R backend rejected the matrix save.",
    );
  };

  const importCsv = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    event.target.value = "";
    if (!file) return;
    const imported = await importMatrixContent(file.name, await file.text());
    if (!imported) {
      setStatus(
        "The matrix import failed. Check that the file is a readable CSV.",
      );
      return;
    }
    const nextFiles = await listMatrixFiles();
    setFiles(nextFiles);
    await load(file.name);
    setStatus(`${file.name} imported and loaded.`);
  };

  const displayFilename =
    filename.replace(/\\/g, "/").split("/").pop() || filename;
  return (
    <div className="matrix-workspace">
      <WorkspaceHeader
        kicker="Matrix review / adjustment"
        title="Tune with context"
        description="Inspect detector signatures and save an adjusted matrix as a new artifact. Every value shown here comes from the active R project."
      />
      <input
        ref={uploadRef}
        type="file"
        accept=".csv,text/csv"
        hidden
        onChange={importCsv}
      />
      <div className="matrix-header-card surface-card">
        <div>
          <span className="eyebrow">Selected matrix</span>
          <h2 title={filename}>{displayFilename || "No matrix selected"}</h2>
          <p>
            {rows.length
              ? `${rows.length.toLocaleString()} markers × ${Math.max(0, columns.length - 1).toLocaleString()} detectors · ${status}`
              : status}
          </p>
        </div>
        <div className="toolbar-actions">
          <select
            className="matrix-file-select"
            value={filename}
            onChange={(event) => void load(event.target.value)}
            disabled={files.length === 0}
          >
            <option value="">Choose matrix…</option>
            {files.map((file) => (
              <option key={file} value={file}>
                {file}
              </option>
            ))}
          </select>
          <button
            className="button button-ghost"
            onClick={() => uploadRef.current?.click()}
          >
            <FolderOpen size={14} /> Import matrix
          </button>
          <button
            className="button button-primary"
            onClick={() => void save()}
            disabled={!rows.length}
          >
            <Save size={14} /> Save adjusted copy
          </button>
        </div>
      </div>
      {rows.length === 0 ? (
        <section className="surface-card empty-workspace">
          <Layers3 size={25} />
          <h2>Bring a matrix into the project</h2>
          <p>
            Import a CSV or build a reference matrix from the Controls
            workspace. The editor will keep the original file and save edits
            through R.
          </p>
          <button
            className="button button-primary"
            onClick={() => uploadRef.current?.click()}
          >
            <FolderOpen size={14} /> Import CSV
          </button>
        </section>
      ) : (
        <div className="matrix-layout">
          <section className="surface-card matrix-table-card">
            <div className="card-toolbar">
              <div>
                <span className="eyebrow">Marker signatures</span>
                <h2>Detector response</h2>
              </div>
              <span className="table-note">
                {rows.length.toLocaleString()} rows · editable values
              </span>
            </div>
            <div className="matrix-scroll">
              <table className="matrix-preview">
                <thead>
                  <tr>
                    <th>{markerColumn}</th>
                    {detectorColumns.map((column) => (
                      <th key={column}>{column}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {visibleRows.map((row, rowIndex) => (
                    <tr
                      className={
                        selectedIndex === rowIndex ? "is-selected" : ""
                      }
                      key={`${String(row[markerColumn])}-${rowIndex}`}
                      onClick={() => setSelectedIndex(rowIndex)}
                    >
                      <td>
                        <strong>{String(row[markerColumn] ?? "")}</strong>
                      </td>
                      {detectorColumns.map((column) => {
                        const value = String(row[column] ?? "");
                        return (
                          <td key={`${rowIndex}-${column}`}>
                            <input
                              className="matrix-edit-cell"
                              value={value}
                              onChange={(event) =>
                                updateCell(rowIndex, column, event.target.value)
                              }
                              aria-label={`${String(row[markerColumn] ?? "")} ${column}`}
                            />
                          </td>
                        );
                      })}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            {rows.length > visibleRows.length && (
              <span className="table-note">
                Showing the first {visibleRows.length} rows. Save preserves all{" "}
                {rows.length.toLocaleString()} rows.
              </span>
            )}
          </section>
          <section className="surface-card tuner-card">
            <div className="card-toolbar">
              <div>
                <span className="eyebrow">Adjustment preview</span>
                <h2>{selectedMarker}</h2>
              </div>
              <StatusPill
                state={status === "Unsaved matrix edits" ? "warning" : "ready"}
                label={
                  status === "Unsaved matrix edits"
                    ? "Unsaved edits"
                    : "Backend values"
                }
              />
            </div>
            <div className="signature-chart">
              <div className="matrix-bars">
                {detectorColumns.map((column) => {
                  const numeric = Number(selectedRow[column] ?? 0);
                  return (
                    <div className="matrix-bar-row" key={column}>
                      <span>{column}</span>
                      <i>
                        <b
                          style={{
                            width: `${Math.min(100, Math.max(0, numeric * 100))}%`,
                          }}
                        />
                      </i>
                      <strong>
                        {Number.isFinite(numeric) ? numeric.toFixed(3) : "—"}
                      </strong>
                    </div>
                  );
                })}
              </div>
              <div className="chart-axis">
                <span>{detectorColumns[0] || "First detector"}</span>
                <span>{Math.max(0, columns.length - 1)} detector channels</span>
                <span>
                  {detectorColumns[detectorColumns.length - 1] ||
                    "Last detector"}
                </span>
              </div>
            </div>
            <div className="tuner-note">
              <Info size={14} />
              <span>
                Edits are sent to R when you save. The source matrix stays
                recoverable through the copy field.
              </span>
            </div>
          </section>
        </div>
      )}
    </div>
  );
}

function LegacyPanelWorkspace({
  panelPayload,
  onLoadPanel,
}: {
  panelPayload: PanelPayload | null;
  onLoadPanel: () => void;
}) {
  const [selected, setSelected] = useState([
    "BV421",
    "FITC",
    "PE",
    "APC",
    "AF700",
    "BUV395",
  ]);
  const fluorophores = panelPayload?.fluorophores?.map(
    (item) => item.fluorophore,
  ) ?? [
    "BV421",
    "FITC",
    "PE",
    "APC",
    "AF700",
    "BUV395",
    "BV510",
    "PerCP-Cy5.5",
    "PE-Cy7",
    "APC-Cy7",
  ];
  return (
    <>
      <WorkspaceHeader
        kicker="Panel builder / spectral design"
        title="See the panel before you run it"
        description="Explore packaged theoretical spectra for Aurora, Discover, ID7000, and Xenith. No controls or samples required."
        action="Refresh library"
        onAction={onLoadPanel}
      />
      <div className="panel-toolbar surface-card">
        <div>
          <span className="eyebrow">Cytometer library</span>
          <h2>Cytek Aurora 5L</h2>
          <p>Violet · Blue · YellowGreen · Red · UV</p>
        </div>
        <div className="panel-selects">
          <select defaultValue="aurora">
            <option value="aurora">Aurora</option>
            <option>Discover</option>
            <option>ID7000</option>
            <option>Xenith</option>
          </select>
          <select defaultValue="5l">
            <option value="5l">5 laser configuration</option>
            <option>4 laser configuration</option>
          </select>
          <button className="button button-primary">
            <Download size={14} /> Export overview
          </button>
        </div>
      </div>
      <div className="panel-builder-grid">
        <section className="surface-card fluor-list">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Fluorophores</span>
              <h2>{selected.length} selected</h2>
            </div>
            <button className="icon-button" aria-label="Add fluorophore">
              <Plus size={16} />
            </button>
          </div>
          <label className="search-field compact-search">
            <Search size={14} />
            <input placeholder="Search library" />
          </label>
          <div className="fluor-options">
            {fluorophores.slice(0, 12).map((fluor) => (
              <button
                className={selected.includes(fluor) ? "is-selected" : ""}
                key={fluor}
                onClick={() =>
                  setSelected((current) =>
                    current.includes(fluor)
                      ? current.filter((item) => item !== fluor)
                      : [...current, fluor],
                  )
                }
              >
                <span className="fluor-swatch" />
                {fluor}
                {selected.includes(fluor) && <Check size={13} />}
              </button>
            ))}
          </div>
        </section>
        <section className="surface-card spectrum-board">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Signatures</span>
              <h2>Overlap at a glance</h2>
            </div>
            <StatusPill state="ready" label={`${selected.length} selected`} />
          </div>
          <svg className="panel-spectrum" viewBox="0 0 650 260">
            <defs>
              <linearGradient id="spectrumArea" x1="0" x2="1">
                <stop offset="0" stopColor="#7cb8ba" stopOpacity=".08" />
                <stop offset="1" stopColor="#d47760" stopOpacity=".18" />
              </linearGradient>
            </defs>
            <rect width="650" height="260" fill="url(#spectrumArea)" />
            {selected.slice(0, 6).map((fluor, index) => (
              <path
                key={fluor}
                d={`M0 ${210 - index * 4} C${60 + index * 10} ${210 - (index % 3) * 90}, ${120 + index * 19} ${48 + index * 17}, ${190 + index * 21} ${146 - index * 10} S${320 + index * 14} ${210 - index * 12}, ${390 + index * 16} ${100 + index * 18} S${520 + index * 6} ${150 - index * 8}, 650 ${82 + index * 18}`}
                fill="none"
                stroke={
                  [
                    "#d37761",
                    "#78b8b3",
                    "#c59a5d",
                    "#7b91ba",
                    "#b47a9c",
                    "#8ca77f",
                  ][index]
                }
                strokeWidth="2.5"
                opacity=".86"
              />
            ))}
            <line x1="0" y1="224" x2="650" y2="224" stroke="#cfc8ba" />
            <text x="0" y="246" fill="#858981" fontSize="12">
              DUV / UV
            </text>
            <text x="297" y="246" fill="#858981" fontSize="12">
              Violet / Blue
            </text>
            <text x="566" y="246" fill="#858981" fontSize="12">
              Red / IR
            </text>
          </svg>
          <div className="panel-metrics">
            <Metric
              value={(panelPayload?.complexity_index ?? 0.68).toFixed(2)}
              label="complexity index"
              accent="accent-amber"
            />
            <Metric value="0.14" label="median overlap" />
            <Metric value="6" label="peak detectors" accent="accent-teal" />
          </div>
        </section>
      </div>
    </>
  );
}

function PanelWorkspace({
  panelPayload,
  onLoadPanel,
}: {
  panelPayload: PanelPayload | null;
  onLoadPanel: () => void;
}) {
  const [cytometer, setCytometer] = useState("aurora");
  const [configuration, setConfiguration] = useState("");
  const [selected, setSelected] = useState<string[]>(
    panelPayload?.selected ?? ["BV421", "FITC", "PE", "APC", "AF700", "BUV395"],
  );
  const [search, setSearch] = useState("");
  const [livePayload, setLivePayload] = useState<PanelPayload | null>(
    panelPayload,
  );
  const payloadForView = livePayload ?? panelPayload;
  const fluorophores = payloadForView?.fluorophores?.map(
    (item) => item.fluorophore,
  ) ?? [
    "BV421",
    "FITC",
    "PE",
    "APC",
    "AF700",
    "BUV395",
    "BV510",
    "PerCP-Cy5.5",
    "PE-Cy7",
    "APC-Cy7",
  ];
  const filteredFluorophores = fluorophores
    .filter((fluor) => fluor.toLowerCase().includes(search.toLowerCase()))
    .slice(0, 18);

  const refresh = async (nextSelected = selected) => {
    const payload = await loadPanelMetrics(
      cytometer,
      configuration,
      nextSelected,
    );
    if (payload) setLivePayload(payload);
  };

  const toggleFluorophore = (fluorophore: string) => {
    const nextSelected = selected.includes(fluorophore)
      ? selected.filter((item) => item !== fluorophore)
      : [...selected, fluorophore];
    setSelected(nextSelected);
    void refresh(nextSelected);
  };

  const exportOverview = async () => {
    const result = await exportPanelOverview(
      cytometer,
      configuration,
      selected,
    );
    if (result)
      downloadBase64File(
        result.filename,
        result.contentBase64,
        result.contentType,
      );
  };

  return (
    <>
      <WorkspaceHeader
        kicker="Panel builder / spectral design"
        title="See the panel before you run it"
        description="Explore packaged theoretical spectra for the active cytometer and update overlap metrics as you select fluorophores."
        action="Refresh library"
        onAction={() => {
          onLoadPanel();
          void refresh();
        }}
      />
      <div className="panel-toolbar surface-card">
        <div>
          <span className="eyebrow">Cytometer library</span>
          <h2>{cytometer === "aurora" ? "Cytek Aurora" : cytometer}</h2>
          <p>Library and metrics are provided by the R backend.</p>
        </div>
        <div className="panel-selects">
          <select
            value={cytometer}
            onChange={(event) => {
              setCytometer(event.target.value);
              void refresh();
            }}
          >
            <option value="aurora">Aurora</option>
            <option value="discover">Discover</option>
            <option value="id7000">ID7000</option>
            <option value="xenith">Xenith</option>
          </select>
          <select
            value={configuration}
            onChange={(event) => {
              setConfiguration(event.target.value);
              void refresh();
            }}
          >
            <option value="">Default configuration</option>
            <option value="5l">5 laser configuration</option>
            <option value="4l">4 laser configuration</option>
          </select>
          <button
            className="button button-primary"
            onClick={() => void exportOverview()}
          >
            <Download size={14} /> Export overview
          </button>
        </div>
      </div>
      <div className="panel-builder-grid">
        <section className="surface-card fluor-list">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Fluorophores</span>
              <h2>{selected.length} selected</h2>
            </div>
            <button
              className="icon-button"
              aria-label="Clear fluorophore selection"
              onClick={() => {
                setSelected([]);
                void refresh([]);
              }}
            >
              <X size={16} />
            </button>
          </div>
          <label className="search-field compact-search">
            <Search size={14} />
            <input
              value={search}
              onChange={(event) => setSearch(event.target.value)}
              placeholder="Search library"
            />
          </label>
          <div className="fluor-options">
            {filteredFluorophores.map((fluor) => (
              <button
                className={selected.includes(fluor) ? "is-selected" : ""}
                key={fluor}
                onClick={() => toggleFluorophore(fluor)}
              >
                <span className="fluor-swatch" />
                {fluor}
                {selected.includes(fluor) && <Check size={13} />}
              </button>
            ))}
          </div>
        </section>
        <section className="surface-card spectrum-board">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Signatures</span>
              <h2>Overlap at a glance</h2>
            </div>
            <StatusPill state="ready" label={`${selected.length} selected`} />
          </div>
          <svg className="panel-spectrum" viewBox="0 0 650 260">
            <defs>
              <linearGradient id="spectrumAreaLive" x1="0" x2="1">
                <stop offset="0" stopColor="#7cb8ba" stopOpacity=".08" />
                <stop offset="1" stopColor="#d47760" stopOpacity=".18" />
              </linearGradient>
            </defs>
            <rect width="650" height="260" fill="url(#spectrumAreaLive)" />
            {selected.slice(0, 8).map((fluor, index) => (
              <path
                key={fluor}
                d={`M0 ${210 - index * 4} C${60 + index * 10} ${210 - (index % 3) * 90}, ${120 + index * 19} ${48 + index * 17}, ${190 + index * 21} ${146 - index * 10} S${320 + index * 14} ${210 - index * 12}, ${390 + index * 16} ${100 + index * 18} S${520 + index * 6} ${150 - index * 8}, 650 ${82 + index * 18}`}
                fill="none"
                stroke={
                  [
                    "#d37761",
                    "#78b8b3",
                    "#c59a5d",
                    "#7b91ba",
                    "#b47a9c",
                    "#8ca77f",
                    "#8a7bb5",
                    "#6c9cb0",
                  ][index]
                }
                strokeWidth="2.5"
                opacity=".86"
              />
            ))}
            <line x1="0" y1="224" x2="650" y2="224" stroke="#cfc8ba" />
            <text x="0" y="246" fill="#858981" fontSize="12">
              DUV / UV
            </text>
            <text x="297" y="246" fill="#858981" fontSize="12">
              Violet / Blue
            </text>
            <text x="566" y="246" fill="#858981" fontSize="12">
              Red / IR
            </text>
          </svg>
          <div className="panel-metrics">
            <Metric
              value={Number(payloadForView?.complexity_index ?? 0).toFixed(2)}
              label="complexity index"
              accent="accent-amber"
            />
            <Metric
              value={
                selected.length ? String(Math.max(0, selected.length - 1)) : "0"
              }
              label="selected overlap pairs"
            />
            <Metric
              value={String(payloadForView?.peak_detectors?.length ?? 0)}
              label="peak detectors"
              accent="accent-teal"
            />
          </div>
        </section>
      </div>
    </>
  );
}

function AfWorkspace({ onRun }: { onRun: WorkflowWorkspaceProps["onRun"] }) {
  return (
    <>
      <WorkspaceHeader
        kicker="AF profile library / global background"
        title="Keep autofluorescence reusable"
        description="Extract once from an unstained file, inspect the bands, then apply a saved profile to the next reference matrix."
      />
      <div className="af-layout">
        <section className="surface-card af-extract-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Extract profile</span>
              <h2>From unstained_cells.fcs</h2>
            </div>
            <StatusPill state="ready" label="12 bands" />
          </div>
          <div className="af-form">
            <label>
              Source FCS
              <select defaultValue="unstained_cells.fcs">
                <option>unstained_cells.fcs</option>
                <option>PBMC_unstained_2026-06-28.fcs</option>
              </select>
            </label>
            <label>
              AF band count
              <input type="number" defaultValue="12" />
            </label>
            <label>
              Maximum cells
              <input type="number" defaultValue="50000" />
            </label>
            <label>
              Seed
              <input type="number" defaultValue="42" />
            </label>
          </div>
          <div className="af-preview">
            <svg viewBox="0 0 500 120">
              <path
                d="M0 99 C18 96 33 82 47 90 S73 108 89 79 S111 61 128 89 S155 101 170 70 S190 37 210 78 S237 95 253 63 S276 56 291 81 S315 105 331 73 S351 52 368 66 S395 91 412 48 S449 73 500 35"
                fill="none"
                stroke="#7cbdb6"
                strokeWidth="3"
              />
              {Array.from({ length: 12 }, (_, index) => (
                <line
                  key={index}
                  x1={30 + index * 39}
                  y1="22"
                  x2={30 + index * 39}
                  y2="103"
                  stroke="#d68a70"
                  strokeOpacity=".25"
                  strokeDasharray="3 4"
                />
              ))}
            </svg>
            <span>Detector-wise AF spectrum · 47 channels</span>
          </div>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("af", "Extract AF profile")}
          >
            <WandSparkles size={15} /> Extract profile
          </button>
        </section>
        <section className="surface-card af-library-card">
          <SectionTitle
            eyebrow="Saved profiles"
            title="1 reusable profile"
            action={
              <button className="button button-ghost">
                <Plus size={14} /> Add profile
              </button>
            }
          />
          <div className="saved-profile">
            <div className="profile-avatar">AF</div>
            <div className="profile-copy">
              <strong>PBMC broad AF</strong>
              <span>12 bands · 47 detectors</span>
              <small>Created Jun 28, 2026 · Spectreasy 0.9.4</small>
            </div>
            <button className="icon-button" aria-label="Profile actions">
              <MoreHorizontalIcon />
            </button>
          </div>
          <div className="profile-actions">
            <button className="text-action">
              <BarChart3 size={14} /> Plot profile
            </button>
            <button className="text-action">
              <Layers3 size={14} /> Apply to matrix
            </button>
            <button className="text-action danger">
              <X size={14} /> Delete
            </button>
          </div>
          <div className="profile-warning">
            <Info size={14} />
            <span>
              Applying this profile will replace existing AF rows after
              confirmation.
            </span>
          </div>
        </section>
      </div>
    </>
  );
}

function LegacyConfigurableAfWorkspace({
  settings,
  onSettingsChange,
  onRun,
}: {
  settings: AfSettings;
  onSettingsChange: (patch: Partial<AfSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  return (
    <>
      <WorkspaceHeader
        kicker="AF profile library / global background"
        title="Keep autofluorescence reusable"
        description="Extract an AF profile from an unstained FCS file with parameters that stay visible and reproducible."
      />
      <div className="af-layout">
        <section className="surface-card af-extract-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Extract profile</span>
              <h2>{settings.fcsFile}</h2>
            </div>
            <StatusPill state="ready" label={`${settings.afNBands} bands`} />
          </div>
          <div className="af-form">
            <label>
              Source FCS
              <input
                value={settings.fcsFile}
                onChange={(event) =>
                  onSettingsChange({ fcsFile: event.target.value })
                }
              />
            </label>
            <label>
              Save as profile
              <input
                value={settings.saveName}
                onChange={(event) =>
                  onSettingsChange({ saveName: event.target.value })
                }
                placeholder="Optional profile name"
              />
            </label>
            <label>
              AF band count
              <input
                type="number"
                min="1"
                value={settings.afNBands}
                onChange={(event) =>
                  onSettingsChange({ afNBands: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Maximum cells
              <input
                type="number"
                min="1"
                value={settings.afMaxCells}
                onChange={(event) =>
                  onSettingsChange({ afMaxCells: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Minimum cluster events
              <input
                type="number"
                min="1"
                value={settings.afMinClusterEvents}
                onChange={(event) =>
                  onSettingsChange({
                    afMinClusterEvents: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Minimum cluster proportion
              <input
                type="number"
                min="0"
                max="1"
                step="0.001"
                value={settings.afMinClusterProportion}
                onChange={(event) =>
                  onSettingsChange({
                    afMinClusterProportion: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Seed
              <input
                type="number"
                min="1"
                value={settings.seed}
                onChange={(event) =>
                  onSettingsChange({ seed: Number(event.target.value) })
                }
              />
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveOverwrite}
                onChange={(event) =>
                  onSettingsChange({ saveOverwrite: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Overwrite saved profile</span>
            </label>
          </div>
          <div className="af-preview">
            <svg viewBox="0 0 500 120">
              <path
                d="M0 99 C18 96 33 82 47 90 S73 108 89 79 S111 61 128 89 S155 101 170 70 S190 37 210 78 S237 95 253 63 S276 56 291 81 S315 105 331 73 S351 52 368 66 S395 91 412 48 S449 73 500 35"
                fill="none"
                stroke="#7cbdb6"
                strokeWidth="3"
              />
              {Array.from(
                { length: Math.min(settings.afNBands, 16) },
                (_, index) => (
                  <line
                    key={index}
                    x1={30 + index * 29}
                    y1="22"
                    x2={30 + index * 29}
                    y2="103"
                    stroke="#d68a70"
                    strokeOpacity=".25"
                    strokeDasharray="3 4"
                  />
                ),
              )}
            </svg>
            <span>Detector-wise AF spectrum · settings are sent to R</span>
          </div>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("af", "Extract AF profile")}
          >
            <WandSparkles size={15} /> Extract profile{" "}
            {settings.saveName ? `& save “${settings.saveName}”` : ""}
          </button>
        </section>
        <section className="surface-card af-library-card">
          <SectionTitle
            eyebrow="Saved profiles"
            title="Reusable background profiles"
            action={
              <button
                className="button button-ghost"
                onClick={() => onRun("af", "Refresh AF profile library")}
              >
                <RefreshCcw size={14} /> Refresh profiles
              </button>
            }
          />
          <div className="saved-profile">
            <div className="profile-avatar">AF</div>
            <div className="profile-copy">
              <strong>Profiles are stored in R</strong>
              <span>
                Extracted profiles are available to the next reference build.
              </span>
              <small>
                Choose a profile name above to save the current extraction.
              </small>
            </div>
            <button className="icon-button" aria-label="Profile actions">
              <MoreHorizontalIcon />
            </button>
          </div>
          <div className="profile-warning">
            <Info size={14} />
            <span>
              Applying a saved profile writes a new matrix artifact so the
              original stays recoverable.
            </span>
          </div>
        </section>
      </div>
    </>
  );
}

function ConfigurableAfWorkspace({
  settings,
  matrixFile,
  onSettingsChange,
  onRun,
}: {
  settings: AfSettings;
  matrixFile: string;
  onSettingsChange: (patch: Partial<AfSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
}) {
  const [profiles, setProfiles] = useState<
    Array<{ name: string; bands: number; detectors: number; created: string }>
  >([]);
  const [profileStatus, setProfileStatus] = useState("Loading saved profiles…");

  const refreshProfiles = async () => {
    const nextProfiles = await listAfProfiles();
    setProfiles(nextProfiles);
    setProfileStatus(
      nextProfiles.length
        ? `${nextProfiles.length} saved profile${nextProfiles.length === 1 ? "" : "s"} found in R.`
        : "No saved profiles yet.",
    );
  };

  useEffect(() => {
    const timer = window.setTimeout(() => {
      void refreshProfiles();
    }, 0);
    return () => window.clearTimeout(timer);
  }, []);

  const removeProfile = async (name: string) => {
    const removed = await deleteAfProfile(name);
    setProfileStatus(
      removed ? `${name} deleted.` : `Could not delete ${name}.`,
    );
    if (removed) void refreshProfiles();
  };

  const applyProfile = async (name: string) => {
    const output = matrixFile.replace(/\.csv$/i, `_${name}.csv`);
    const applied = await applyAfProfile(matrixFile, name, output);
    setProfileStatus(
      applied
        ? `${name} applied to ${output}.`
        : `Could not apply ${name} to the selected matrix.`,
    );
  };

  return (
    <>
      <WorkspaceHeader
        kicker="AF profile library / global background"
        title="Keep autofluorescence reusable"
        description="Extract an AF profile from an unstained FCS file, save it in R, and apply it to a new matrix copy when needed."
      />
      <div className="af-layout">
        <section className="surface-card af-extract-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Extract profile</span>
              <h2>{settings.fcsFile}</h2>
            </div>
            <StatusPill state="ready" label={`${settings.afNBands} bands`} />
          </div>
          <div className="af-form">
            <label>
              Source FCS
              <input
                value={settings.fcsFile}
                onChange={(event) =>
                  onSettingsChange({ fcsFile: event.target.value })
                }
              />
            </label>
            <label>
              Save as profile
              <input
                value={settings.saveName}
                onChange={(event) =>
                  onSettingsChange({ saveName: event.target.value })
                }
                placeholder="Optional profile name"
              />
            </label>
            <label>
              AF band count
              <input
                type="number"
                min="1"
                value={settings.afNBands}
                onChange={(event) =>
                  onSettingsChange({ afNBands: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Maximum cells
              <input
                type="number"
                min="1"
                value={settings.afMaxCells}
                onChange={(event) =>
                  onSettingsChange({ afMaxCells: Number(event.target.value) })
                }
              />
            </label>
            <label>
              Minimum cluster events
              <input
                type="number"
                min="1"
                value={settings.afMinClusterEvents}
                onChange={(event) =>
                  onSettingsChange({
                    afMinClusterEvents: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Minimum cluster proportion
              <input
                type="number"
                min="0"
                max="1"
                step="0.001"
                value={settings.afMinClusterProportion}
                onChange={(event) =>
                  onSettingsChange({
                    afMinClusterProportion: Number(event.target.value),
                  })
                }
              />
            </label>
            <label>
              Seed
              <input
                type="number"
                min="1"
                value={settings.seed}
                onChange={(event) =>
                  onSettingsChange({ seed: Number(event.target.value) })
                }
              />
            </label>
            <label className="toggle-label">
              <input
                type="checkbox"
                checked={settings.saveOverwrite}
                onChange={(event) =>
                  onSettingsChange({ saveOverwrite: event.target.checked })
                }
              />
              <span className="toggle-ui" />
              <span>Overwrite saved profile</span>
            </label>
          </div>
          <div className="af-preview">
            <svg viewBox="0 0 500 120">
              <path
                d="M0 99 C18 96 33 82 47 90 S73 108 89 79 S111 61 128 89 S155 101 170 70 S190 37 210 78 S237 95 253 63 S276 56 291 81 S315 105 331 73 S351 52 368 66 S395 91 412 48 S449 73 500 35"
                fill="none"
                stroke="#7cbdb6"
                strokeWidth="3"
              />
              {Array.from(
                { length: Math.min(settings.afNBands, 16) },
                (_, index) => (
                  <line
                    key={index}
                    x1={30 + index * 29}
                    y1="22"
                    x2={30 + index * 29}
                    y2="103"
                    stroke="#d68a70"
                    strokeOpacity=".25"
                    strokeDasharray="3 4"
                  />
                ),
              )}
            </svg>
            <span>Detector-wise AF spectrum · settings are sent to R</span>
          </div>
          <button
            className="button button-primary large-button"
            onClick={() => onRun("af", "Extract AF profile")}
          >
            <WandSparkles size={15} /> Extract profile{" "}
            {settings.saveName ? `& save “${settings.saveName}”` : ""}
          </button>
        </section>
        <section className="surface-card af-library-card">
          <SectionTitle
            eyebrow="Saved profiles"
            title="Reusable background profiles"
            action={
              <button
                className="button button-ghost"
                onClick={() => void refreshProfiles()}
              >
                <RefreshCcw size={14} /> Refresh profiles
              </button>
            }
          />
          <p className="profile-status">{profileStatus}</p>
          {profiles.length === 0 ? (
            <div className="saved-profile">
              <div className="profile-avatar">AF</div>
              <div className="profile-copy">
                <strong>Profiles are stored in R</strong>
                <span>Extract a profile and give it a name to save it.</span>
                <small>
                  Saved profiles can be applied to a new matrix copy.
                </small>
              </div>
            </div>
          ) : (
            <div className="profile-list">
              {profiles.map((profile) => (
                <div className="saved-profile" key={profile.name}>
                  <div className="profile-avatar">AF</div>
                  <div className="profile-copy">
                    <strong>{profile.name}</strong>
                    <span>
                      {profile.bands} bands · {profile.detectors} detectors
                    </span>
                    <small>
                      {profile.created ||
                        "Saved in the local R profile library"}
                    </small>
                  </div>
                  <div className="profile-actions">
                    <button
                      className="text-action"
                      onClick={() => void applyProfile(profile.name)}
                    >
                      <Layers3 size={14} /> Apply
                    </button>
                    <button
                      className="text-action danger"
                      onClick={() => void removeProfile(profile.name)}
                    >
                      <X size={14} /> Delete
                    </button>
                  </div>
                </div>
              ))}
            </div>
          )}
          <div className="profile-warning">
            <Info size={14} />
            <span>
              Applying a saved profile writes{" "}
              {matrixFile.replace(/\.csv$/i, "_profile.csv")} and leaves the
              source matrix unchanged.
            </span>
          </div>
        </section>
      </div>
    </>
  );
}

function MoreHorizontalIcon() {
  return <span className="more-dots">•••</span>;
}

function LegacyExperimentalWorkspace({
  kind,
}: {
  kind: "comparison" | "simulator";
}) {
  const comparison = kind === "comparison";
  return (
    <>
      <WorkspaceHeader
        kicker={
          comparison
            ? "Method comparison / experimental"
            : "Synthetic SCC simulator / experimental"
        }
        title={
          comparison
            ? "Compare methods without losing the baseline"
            : "Benchmark synthetic controls carefully"
        }
        description={
          comparison
            ? "A reserved diagnostic surface for per-marker spread, residuals, NPS, and gate stability."
            : "Paint trusted spectra onto real unstained events in R, with truth metadata and an explicit synthetic label."
        }
      />
      <section className="surface-card experimental-card">
        <div className="experimental-mark">
          {comparison ? <GitCompare size={23} /> : <Sparkles size={23} />}
        </div>
        <StatusPill
          state="blocked"
          label={comparison ? "Experimental module" : "Needs trusted inputs"}
        />
        <h2>
          {comparison
            ? "The comparison bench is ready for a backend action."
            : "Add a trusted matrix and unstained source to begin."}
        </h2>
        <p>
          {comparison
            ? "Select Spectreasy, AutoSpectral, OLS, WLS, RWLS, or NNLS, then run the same inputs through each R method. This panel will compare outputs without changing your active workflow."
            : "Synthetic SCCs are generated from a prior spectrum and real unstained events. They can benchmark a workflow, but they cannot discover an unknown true spectrum."}
        </p>
        <div className="experimental-chips">
          <span>{comparison ? "6 methods reserved" : "Truth metadata"}</span>
          <span>{comparison ? "NPS + residuals" : "Benchmark mode"}</span>
          <span>
            {comparison ? "Diagnostic only" : "Explicit confirmation"}
          </span>
        </div>
        <button className="button button-ghost" disabled>
          <LockKeyhole size={14} />{" "}
          {comparison
            ? "Backend action not available yet"
            : "Select trusted inputs first"}
        </button>
      </section>
    </>
  );
}

function ExperimentalWorkspace({
  kind,
  project,
  settings,
}: {
  kind: "comparison" | "simulator";
  project: ProjectState;
  settings: WorkflowSettings;
}) {
  const comparison = kind === "comparison";
  const [methods, setMethods] = useState(["Spectreasy", "OLS", "WLS", "NNLS"]);
  const [comparisonRows, setComparisonRows] = useState<
    Array<Record<string, unknown>>
  >([]);
  const [status, setStatus] = useState("Ready for a backend diagnostic run.");
  const [marker, setMarker] = useState("");
  const [events, setEvents] = useState(1000);
  const [noise, setNoise] = useState(0.02);

  const runComparison = async () => {
    if (methods.length === 0) {
      setStatus("Select at least one method.");
      return;
    }
    setStatus("Running selected methods in R…");
    const response = await attemptDiagnosticAction("compare", {
      projectPath: project.projectPath,
      matrix_file: settings.sample.matrixFile,
      sample_dir: settings.sample.sampleDir,
      methods,
      seed: settings.sample.seed,
    });
    setStatus(response.message);
    const result = response.result as Record<string, unknown> | undefined;
    setComparisonRows(rowsFromBackend(result?.rows));
  };

  const runSynthetic = async () => {
    setStatus("Generating synthetic SCC FCS in R…");
    const response = await attemptDiagnosticAction("synthetic", {
      projectPath: project.projectPath,
      matrix_file: settings.sample.matrixFile,
      marker,
      events,
      noise,
      seed: settings.sample.seed,
    });
    setStatus(response.message);
    const result = response.result as Record<string, unknown> | undefined;
    if (result?.output_file)
      setComparisonRows([{ ...result, status: "complete" }]);
  };

  const toggleMethod = (method: string) =>
    setMethods((current) =>
      current.includes(method)
        ? current.filter((item) => item !== method)
        : [...current, method],
    );
  return (
    <>
      <WorkspaceHeader
        kicker={
          comparison
            ? "Method comparison / diagnostic"
            : "Synthetic SCC / benchmark"
        }
        title={
          comparison
            ? "Compare methods on the same inputs"
            : "Generate a trusted synthetic SCC"
        }
        description={
          comparison
            ? "Run selected unmixing methods against one active matrix and sample, then save a compact comparison table."
            : "Generate a real FCS benchmark from one selected matrix signature, with detector noise and truth metadata saved beside it."
        }
      />
      <section className="surface-card diagnostic-card">
        <div className="diagnostic-topline">
          <div className="experimental-mark">
            {comparison ? <GitCompare size={23} /> : <Sparkles size={23} />}
          </div>
          <StatusPill state="ready" label="Backend action available" />
        </div>
        {comparison ? (
          <>
            <div className="diagnostic-grid">
              <label>
                Reference matrix
                <input value={settings.sample.matrixFile} readOnly />
              </label>
              <label>
                Sample folder
                <input value={settings.sample.sampleDir} readOnly />
              </label>
            </div>
            <div className="method-picker">
              <span className="eyebrow">Methods to compare</span>
              <div>
                {[
                  "Spectreasy",
                  "AutoSpectral",
                  "OLS",
                  "WLS",
                  "RWLS",
                  "NNLS",
                ].map((method) => (
                  <label className="diagnostic-check" key={method}>
                    <input
                      type="checkbox"
                      checked={methods.includes(method)}
                      onChange={() => toggleMethod(method)}
                    />
                    <span className="toggle-ui" />
                    <span>{method}</span>
                  </label>
                ))}
              </div>
            </div>
            <button
              className="button button-primary"
              onClick={() => void runComparison()}
            >
              <Play size={15} fill="currentColor" /> Run comparison
            </button>
            {comparisonRows.length > 0 && (
              <div className="diagnostic-results">
                <div className="card-toolbar">
                  <div>
                    <span className="eyebrow">R results</span>
                    <h2>Method residuals</h2>
                  </div>
                  <span className="table-note">
                    Saved under spectreasy_outputs/method_comparison
                  </span>
                </div>
                <table className="diagnostic-table">
                  <thead>
                    <tr>
                      <th>Method</th>
                      <th>Status</th>
                      <th>Residual RMS</th>
                      <th>Message</th>
                    </tr>
                  </thead>
                  <tbody>
                    {comparisonRows.map((row, index) => (
                      <tr key={`${String(row.method)}-${index}`}>
                        <td>
                          <strong>{String(row.method ?? "")}</strong>
                        </td>
                        <td>{String(row.status ?? "")}</td>
                        <td>
                          {row.residual_rms == null || row.residual_rms === "NA"
                            ? "—"
                            : Number(row.residual_rms).toFixed(5)}
                        </td>
                        <td>{String(row.message ?? "")}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </>
        ) : (
          <>
            <div className="diagnostic-grid">
              <label>
                Reference matrix
                <input value={settings.sample.matrixFile} readOnly />
              </label>
              <label>
                Matrix marker
                <input
                  value={marker}
                  onChange={(event) => setMarker(event.target.value)}
                  placeholder="Blank = first marker"
                />
              </label>
              <label>
                Event count
                <input
                  type="number"
                  min="1"
                  max="100000"
                  value={events}
                  onChange={(event) => setEvents(Number(event.target.value))}
                />
              </label>
              <label>
                Multiplicative noise
                <input
                  type="number"
                  min="0"
                  step="0.005"
                  value={noise}
                  onChange={(event) => setNoise(Number(event.target.value))}
                />
              </label>
            </div>
            <button
              className="button button-primary"
              onClick={() => void runSynthetic()}
            >
              <Sparkles size={15} /> Generate synthetic FCS
            </button>
            {comparisonRows.length > 0 && (
              <div className="diagnostic-result-note">
                <Check size={14} />
                <span>
                  {String(
                    comparisonRows[0].output_file ?? "Synthetic output written",
                  )}{" "}
                  · truth metadata saved alongside it.
                </span>
              </div>
            )}
          </>
        )}
        <p className="diagnostic-status">
          <Info size={14} /> {status}
        </p>
      </section>
    </>
  );
}

function ConfigurableSettingsWorkspace({
  backend,
  job,
  settings,
  onSettingsChange,
  onOpenProject,
  onSave,
}: {
  backend: BackendStatus;
  job: Job;
  settings: WorkflowSettings;
  onSettingsChange: WorkflowWorkspaceProps["onSettingsChange"];
  onOpenProject: WorkflowWorkspaceProps["onOpenProject"];
  onSave: () => void;
}) {
  const control = settings.control;
  const sample = settings.sample;
  const af = settings.af;
  return (
    <>
      <WorkspaceHeader
        kicker="Settings, logs & advanced tools"
        title="Make every run reproducible"
        description="The browser owns the project context and user-facing parameters. Spectreasy R remains responsible for validation and numerical work."
      />
      <div className="settings-grid">
        <section className="surface-card settings-card">
          <SectionTitle
            eyebrow="Active project"
            title="Choose the folder once"
            note="All relative paths below are resolved from this folder."
          />
          <div className="settings-form">
            <label>
              Project path
              <input
                value={settings.projectPath}
                onChange={(event) =>
                  onSettingsChange(
                    "projectPath" as never,
                    { projectPath: event.target.value } as never,
                  )
                }
              />
            </label>
            <button
              className="button button-primary"
              onClick={() => onOpenProject(settings.projectPath)}
            >
              <FolderOpen size={14} /> Open project in R
            </button>
          </div>
          <div className="health-status">
            <span
              className={`health-orb ${backend.connected ? "is-connected" : ""}`}
            />
            <div>
              <strong>
                {backend.connected
                  ? "Connected to Spectreasy R"
                  : "Not connected"}
              </strong>
              <span>{backend.message}</span>
            </div>
            <StatusPill
              state={backend.connected ? "connected" : "offline"}
              compact
            />
          </div>
          <div className="settings-list">
            <div>
              <span>Host</span>
              <strong>127.0.0.1</strong>
            </div>
            <div>
              <span>API port</span>
              <strong>{backend.apiPort}</strong>
            </div>
            <div>
              <span>Package</span>
              <strong>{backend.version}</strong>
            </div>
            <div>
              <span>Assets</span>
              <strong>bundled locally</strong>
            </div>
          </div>
        </section>
        <section className="surface-card settings-card">
          <SectionTitle
            eyebrow="Defaults"
            title="Workflow defaults"
            note="These preferences persist in the local Spectreasy config."
          />
          <div className="settings-list">
            <div>
              <span>Controls</span>
              <strong>
                {control.method} · {control.cytometer}
              </strong>
            </div>
            <div>
              <span>Samples</span>
              <strong>
                {sample.method} · {sample.nThreads} thread(s)
              </strong>
            </div>
            <div>
              <span>AF extraction</span>
              <strong>
                {af.afNBands} bands · {af.afMaxCells.toLocaleString()} cells
              </strong>
            </div>
            <div>
              <span>Reports</span>
              <strong>
                {control.saveReport && sample.saveReport
                  ? "Control + sample QC on"
                  : "Selective reports"}
              </strong>
            </div>
          </div>
          <button
            className="button button-primary settings-save"
            onClick={onSave}
          >
            <Save size={14} /> Save all settings
          </button>
        </section>
      </div>
      <AppearanceSettings
        value={settings.appearance}
        onChange={(patch) => onSettingsChange("appearance", patch)}
      />
      <details className="surface-card settings-section" open>
        <summary>
          <CircleCheckBig size={15} /> Control-stage parameters{" "}
          <span>unmix_controls · build_reference_matrix</span>
        </summary>
        <div className="settings-form-grid">
          <label>
            SCC folder
            <input
              value={control.sccDir}
              onChange={(event) =>
                onSettingsChange("control", { sccDir: event.target.value })
              }
            />
          </label>
          <label>
            Mapping CSV
            <input
              value={control.controlFile}
              onChange={(event) =>
                onSettingsChange("control", { controlFile: event.target.value })
              }
            />
          </label>
          <label>
            Output folder
            <input
              value={control.outputDir}
              onChange={(event) =>
                onSettingsChange("control", { outputDir: event.target.value })
              }
            />
          </label>
          <label>
            Cytometer
            <select
              value={control.cytometer}
              onChange={(event) =>
                onSettingsChange("control", { cytometer: event.target.value })
              }
            >
              <option value="auto">Auto-detect</option>
              <option value="aurora">Cytek Aurora</option>
              <option value="aurora_5l">Cytek Aurora 5L</option>
              <option value="aurora_4l">Cytek Aurora 4L</option>
              <option value="discover">Cytek Aurora Discover</option>
              <option value="id7000">Sony ID7000</option>
            </select>
          </label>
          <label>
            Unmixing method
            <select
              value={control.method}
              onChange={(event) =>
                onSettingsChange("control", { method: event.target.value })
              }
            >
              <option>Spectreasy</option>
              <option>AutoSpectral</option>
              <option>OLS</option>
              <option>WLS</option>
              <option>RWLS</option>
              <option>NNLS</option>
            </select>
          </label>
          <label>
            Unknown fluor policy
            <select
              value={control.autoUnknownFluorPolicy}
              onChange={(event) =>
                onSettingsChange("control", {
                  autoUnknownFluorPolicy: event.target
                    .value as ControlSettings["autoUnknownFluorPolicy"],
                })
              }
            >
              <option value="by_channel">Infer by channel</option>
              <option value="empty">Leave empty</option>
              <option value="filename">Infer by filename</option>
            </select>
          </label>
          <label>
            Gate CSV
            <input
              value={control.gateFile}
              onChange={(event) =>
                onSettingsChange("control", { gateFile: event.target.value })
              }
              placeholder="Optional"
            />
          </label>
          <label>
            AF bands
            <input
              type="number"
              min="1"
              value={control.afNBands}
              onChange={(event) =>
                onSettingsChange("control", {
                  afNBands: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            AF max cluster events
            <input
              type="number"
              min="1"
              value={control.afMinClusterEvents}
              onChange={(event) =>
                onSettingsChange("control", {
                  afMinClusterEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            AF min proportion
            <input
              type="number"
              min="0"
              max="1"
              step="0.001"
              value={control.afMinClusterProportion}
              onChange={(event) =>
                onSettingsChange("control", {
                  afMinClusterProportion: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Scatter panel width (mm)
            <input
              type="number"
              min="1"
              value={control.unmixScatterPanelSizeMm}
              onChange={(event) =>
                onSettingsChange("control", {
                  unmixScatterPanelSizeMm: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            RWLS iterations
            <input
              type="number"
              min="1"
              value={control.rwlsMaxIter}
              onChange={(event) =>
                onSettingsChange("control", {
                  rwlsMaxIter: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Unmix threads
            <input
              type="number"
              min="1"
              value={control.unmixThreads}
              onChange={(event) =>
                onSettingsChange("control", {
                  unmixThreads: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Seed
            <input
              type="number"
              min="1"
              value={control.seed}
              onChange={(event) =>
                onSettingsChange("control", {
                  seed: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Background method
            <select
              value={control.sccBackgroundMethod}
              onChange={(event) =>
                onSettingsChange("control", {
                  sccBackgroundMethod: event.target
                    .value as ControlSettings["sccBackgroundMethod"],
                })
              }
            >
              <option value="scatter_knn">Scatter KNN</option>
              <option value="none">None</option>
            </select>
          </label>
          <label>
            Background K
            <input
              type="number"
              min="1"
              value={control.sccBackgroundK}
              onChange={(event) =>
                onSettingsChange("control", {
                  sccBackgroundK: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Variant SOM nodes
            <input
              type="number"
              min="1"
              value={control.spectralVariantSomNodes}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectralVariantSomNodes: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Variant top-k
            <input
              type="number"
              min="1"
              value={control.spectralVariantTopK}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectralVariantTopK: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Cosine threshold
            <input
              type="number"
              min="0"
              max="1"
              step="0.001"
              value={control.spectralVariantCosineThreshold}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectralVariantCosineThreshold: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Maximum variants
            <input
              type="number"
              min="1"
              value={control.spectralVariantMaxVariants}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectralVariantMaxVariants: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Variant min events
            <input
              type="number"
              min="1"
              value={control.spectralVariantMinEvents}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectralVariantMinEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Weight quantile
            <input
              type="number"
              min="0"
              max="1"
              step="0.01"
              value={control.spectreasyWeightQuantile}
              onChange={(event) =>
                onSettingsChange("control", {
                  spectreasyWeightQuantile: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            AutoSpectral candidates
            <input
              type="number"
              min="1"
              value={control.autospectralNCandidates}
              onChange={(event) =>
                onSettingsChange("control", {
                  autospectralNCandidates: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            AutoSpectral spectra
            <input
              type="number"
              min="1"
              value={control.autospectralNSpectral}
              onChange={(event) =>
                onSettingsChange("control", {
                  autospectralNSpectral: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            AutoSpectral min events
            <input
              type="number"
              min="1"
              value={control.autospectralMinEvents}
              onChange={(event) =>
                onSettingsChange("control", {
                  autospectralMinEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={control.autoCreateMapping}
              onChange={(event) =>
                onSettingsChange("control", {
                  autoCreateMapping: event.target.checked,
                })
              }
            />
            <span className="toggle-ui" />
            <span>Auto-create mapping</span>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={control.saveQcPlots}
              onChange={(event) =>
                onSettingsChange("control", {
                  saveQcPlots: event.target.checked,
                })
              }
            />
            <span className="toggle-ui" />
            <span>Save QC plots</span>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={control.saveReport}
              onChange={(event) =>
                onSettingsChange("control", {
                  saveReport: event.target.checked,
                })
              }
            />
            <span className="toggle-ui" />
            <span>Generate QC report</span>
          </label>
          <label>
            Control report format
            <select
              value={control.outputFormat}
              onChange={(event) =>
                onSettingsChange("control", {
                  outputFormat: event.target.value as ControlSettings["outputFormat"],
                })
              }
            >
              <option value="html">HTML</option>
              <option value="pdf">PDF</option>
            </select>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={control.refine}
              onChange={(event) =>
                onSettingsChange("control", { refine: event.target.checked })
              }
            />
            <span className="toggle-ui" />
            <span>Refine AutoSpectral matrix</span>
          </label>
        </div>
      </details>
      <details className="surface-card settings-section">
        <summary>
          <Beaker size={15} /> Sample-stage parameters{" "}
          <span>unmix_samples</span>
        </summary>
        <div className="settings-form-grid">
          <label>
            Sample folder
            <input
              value={sample.sampleDir}
              onChange={(event) =>
                onSettingsChange("sample", { sampleDir: event.target.value })
              }
            />
          </label>
          <label>
            Reference matrix
            <input
              value={sample.matrixFile}
              onChange={(event) =>
                onSettingsChange("sample", { matrixFile: event.target.value })
              }
            />
          </label>
          <label>
            Detector noise file
            <input
              value={sample.detectorNoiseFile}
              onChange={(event) =>
                onSettingsChange("sample", {
                  detectorNoiseFile: event.target.value,
                })
              }
            />
          </label>
          <label>
            Output folder
            <input
              value={sample.outputDir}
              onChange={(event) =>
                onSettingsChange("sample", { outputDir: event.target.value })
              }
            />
          </label>
          <label>
            Method
            <select
              value={sample.method}
              onChange={(event) =>
                onSettingsChange("sample", { method: event.target.value })
              }
            >
              <option>Spectreasy</option>
              <option>AutoSpectral</option>
              <option>OLS</option>
              <option>WLS</option>
              <option>RWLS</option>
              <option>NNLS</option>
            </select>
          </label>
          <label>
            RWLS iterations
            <input
              type="number"
              min="1"
              value={sample.rwlsMaxIter}
              onChange={(event) =>
                onSettingsChange("sample", {
                  rwlsMaxIter: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Threads
            <input
              type="number"
              min="1"
              value={sample.nThreads}
              onChange={(event) =>
                onSettingsChange("sample", {
                  nThreads: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Variant top-k
            <input
              type="number"
              min="1"
              value={sample.spectralVariantTopK}
              onChange={(event) =>
                onSettingsChange("sample", {
                  spectralVariantTopK: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Variant min abundance
            <input
              type="number"
              min="0"
              value={sample.spectralVariantMinAbundance}
              onChange={(event) =>
                onSettingsChange("sample", {
                  spectralVariantMinAbundance: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Positive fraction
            <input
              type="number"
              min="0"
              max="1"
              step="0.001"
              value={sample.spectralVariantPositiveFraction}
              onChange={(event) =>
                onSettingsChange("sample", {
                  spectralVariantPositiveFraction: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Min improvement
            <input
              type="number"
              min="0"
              step="0.001"
              value={sample.spectralVariantMinImprovement}
              onChange={(event) =>
                onSettingsChange("sample", {
                  spectralVariantMinImprovement: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Weight quantile
            <input
              type="number"
              min="0"
              max="1"
              step="0.01"
              value={sample.spectreasyWeightQuantile}
              onChange={(event) =>
                onSettingsChange("sample", {
                  spectreasyWeightQuantile: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Plot events
            <input
              type="number"
              min="1"
              value={sample.plotNEvents}
              onChange={(event) =>
                onSettingsChange("sample", {
                  plotNEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Chunk size
            <input
              type="number"
              min="1"
              value={sample.chunkSize}
              onChange={(event) =>
                onSettingsChange("sample", {
                  chunkSize: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Seed
            <input
              type="number"
              min="1"
              value={sample.seed}
              onChange={(event) =>
                onSettingsChange("sample", { seed: Number(event.target.value) })
              }
            />
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={sample.estimateAf}
              onChange={(event) =>
                onSettingsChange("sample", { estimateAf: event.target.checked })
              }
            />
            <span className="toggle-ui" />
            <span>Estimate AF</span>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={sample.writeFcs}
              onChange={(event) =>
                onSettingsChange("sample", { writeFcs: event.target.checked })
              }
            />
            <span className="toggle-ui" />
            <span>Write FCS</span>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={sample.saveReport}
              onChange={(event) =>
                onSettingsChange("sample", { saveReport: event.target.checked })
              }
            />
            <span className="toggle-ui" />
            <span>Generate report</span>
          </label>
          <label>
            Sample report format
            <select
              value={sample.outputFormat}
              onChange={(event) =>
                onSettingsChange("sample", {
                  outputFormat: event.target.value as SampleSettings["outputFormat"],
                })
              }
            >
              <option value="html">HTML</option>
              <option value="pdf">PDF</option>
            </select>
          </label>
          <label className="toggle-label">
            <input
              type="checkbox"
              checked={sample.saveQcPlots}
              onChange={(event) =>
                onSettingsChange("sample", {
                  saveQcPlots: event.target.checked,
                })
              }
            />
            <span className="toggle-ui" />
            <span>Save QC plots</span>
          </label>
        </div>
      </details>
      <details className="surface-card settings-section">
        <summary>
          <FlaskConical size={15} /> AF profile parameters{" "}
          <span>extract_af_profile</span>
        </summary>
        <div className="settings-form-grid">
          <label>
            Source FCS
            <input
              value={af.fcsFile}
              onChange={(event) =>
                onSettingsChange("af", { fcsFile: event.target.value })
              }
            />
          </label>
          <label>
            AF bands
            <input
              type="number"
              min="1"
              value={af.afNBands}
              onChange={(event) =>
                onSettingsChange("af", { afNBands: Number(event.target.value) })
              }
            />
          </label>
          <label>
            Maximum cells
            <input
              type="number"
              min="1"
              value={af.afMaxCells}
              onChange={(event) =>
                onSettingsChange("af", {
                  afMaxCells: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Minimum cluster events
            <input
              type="number"
              min="1"
              value={af.afMinClusterEvents}
              onChange={(event) =>
                onSettingsChange("af", {
                  afMinClusterEvents: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Minimum cluster proportion
            <input
              type="number"
              min="0"
              max="1"
              step="0.001"
              value={af.afMinClusterProportion}
              onChange={(event) =>
                onSettingsChange("af", {
                  afMinClusterProportion: Number(event.target.value),
                })
              }
            />
          </label>
          <label>
            Seed
            <input
              type="number"
              min="1"
              value={af.seed}
              onChange={(event) =>
                onSettingsChange("af", { seed: Number(event.target.value) })
              }
            />
          </label>
        </div>
      </details>
      <section className="surface-card job-history">
        <div className="card-toolbar">
          <div>
            <span className="eyebrow">Job history</span>
            <h2>Recent actions</h2>
          </div>
          <button className="button button-ghost">
            <Download size={14} /> Download logs
          </button>
        </div>
        <div className="job-history-row">
          <span className="job-history-status complete">
            <Check size={13} />
          </span>
          <div>
            <strong>Control QC report generated</strong>
            <small>run-014 · today, 09:48 · 4m 12s</small>
          </div>
          <span className="job-history-output">9 outputs</span>
          <ArrowRight size={14} />
        </div>
        <div className="job-history-row">
          <span className="job-history-status complete">
            <Check size={13} />
          </span>
          <div>
            <strong>Reference matrix built</strong>
            <small>run-014 · today, 09:44 · 3m 38s</small>
          </div>
          <span className="job-history-output">6 outputs</span>
          <ArrowRight size={14} />
        </div>
        {job.state === "running" && (
          <div className="job-history-live">
            <RefreshCcw size={14} className="spin" />
            <span>
              {job.label} · {job.subtask}
            </span>
            <strong>{job.progress}%</strong>
          </div>
        )}
      </section>
    </>
  );
}

function ControlReferenceTuning({
  settings,
  onSettingsChange,
}: {
  settings: ControlSettings;
  onSettingsChange: (patch: Partial<ControlSettings>) => void;
}) {
  return (
    <section className="surface-card settings-section">
      <div className="settings-section-title">
        <div>
          <span className="eyebrow">Reference builder</span>
          <h2>Gating & clustering parameters</h2>
        </div>
        <span>build_reference_matrix</span>
      </div>
      <div className="settings-form-grid">
        <label>
          AF max cells
          <input
            type="number"
            min="1"
            value={settings.afMaxCells}
            onChange={(event) =>
              onSettingsChange({ afMaxCells: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Default sample type
          <select
            value={settings.defaultSampleType}
            onChange={(event) =>
              onSettingsChange({
                defaultSampleType: event.target
                  .value as ControlSettings["defaultSampleType"],
              })
            }
          >
            <option value="beads">Beads</option>
            <option value="cells">Cells</option>
          </select>
        </label>
        <label>
          Bead histogram width
          <input
            type="number"
            min="0"
            max="1"
            step="0.01"
            value={settings.histogramPctBeads}
            onChange={(event) =>
              onSettingsChange({
                histogramPctBeads: Number(event.target.value),
              })
            }
          />
        </label>
        <label>
          Bead histogram direction
          <select
            value={settings.histogramDirectionBeads}
            onChange={(event) =>
              onSettingsChange({
                histogramDirectionBeads: event.target
                  .value as ControlSettings["histogramDirectionBeads"],
              })
            }
          >
            <option value="right">Right</option>
            <option value="both">Both</option>
            <option value="left">Left</option>
          </select>
        </label>
        <label>
          Cell histogram width
          <input
            type="number"
            min="0"
            max="1"
            step="0.01"
            value={settings.histogramPctCells}
            onChange={(event) =>
              onSettingsChange({
                histogramPctCells: Number(event.target.value),
              })
            }
          />
        </label>
        <label>
          Cell histogram direction
          <select
            value={settings.histogramDirectionCells}
            onChange={(event) =>
              onSettingsChange({
                histogramDirectionCells: event.target
                  .value as ControlSettings["histogramDirectionCells"],
              })
            }
          >
            <option value="right">Right</option>
            <option value="both">Both</option>
            <option value="left">Left</option>
          </select>
        </label>
        <label>
          Outlier percentile
          <input
            type="number"
            min="0"
            max="1"
            step="0.001"
            value={settings.outlierPercentile}
            onChange={(event) =>
              onSettingsChange({
                outlierPercentile: Number(event.target.value),
              })
            }
          />
        </label>
        <label>
          Debris percentile
          <input
            type="number"
            min="0"
            max="1"
            step="0.001"
            value={settings.debrisPercentile}
            onChange={(event) =>
              onSettingsChange({ debrisPercentile: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Bead gate scale
          <input
            type="number"
            min="0"
            step="0.1"
            value={settings.beadGateScale}
            onChange={(event) =>
              onSettingsChange({ beadGateScale: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Maximum clusters
          <input
            type="number"
            min="1"
            value={settings.maxClusters}
            onChange={(event) =>
              onSettingsChange({ maxClusters: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Minimum cluster proportion
          <input
            type="number"
            min="0"
            max="1"
            step="0.001"
            value={settings.minClusterProportion}
            onChange={(event) =>
              onSettingsChange({
                minClusterProportion: Number(event.target.value),
              })
            }
          />
        </label>
        <label>
          Bead gate contour
          <input
            type="number"
            min="0"
            max="1"
            step="0.01"
            value={settings.gateContourBeads}
            onChange={(event) =>
              onSettingsChange({ gateContourBeads: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Cell gate contour
          <input
            type="number"
            min="0"
            max="1"
            step="0.01"
            value={settings.gateContourCells}
            onChange={(event) =>
              onSettingsChange({ gateContourCells: Number(event.target.value) })
            }
          />
        </label>
        <label>
          Reference subsample
          <input
            type="number"
            min="1"
            value={settings.subsampleN}
            onChange={(event) =>
              onSettingsChange({ subsampleN: Number(event.target.value) })
            }
          />
        </label>
      </div>
    </section>
  );
}

function SampleOutputTuning({
  settings,
  onSettingsChange,
}: {
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
}) {
  return (
    <section className="surface-card settings-section">
      <div className="settings-section-title">
        <div>
          <span className="eyebrow">Sample outputs</span>
          <h2>Library & return type</h2>
        </div>
        <span>unmix_samples</span>
      </div>
      <div className="settings-form-grid">
        <label>
          Spectral variant library
          <input
            value={settings.spectralVariantLibraryFile}
            onChange={(event) =>
              onSettingsChange({
                spectralVariantLibraryFile: event.target.value,
              })
            }
            placeholder="Optional .rds file"
          />
        </label>
        <label>
          Return type
          <select
            value={settings.returnType}
            onChange={(event) =>
              onSettingsChange({
                returnType: event.target.value as SampleSettings["returnType"],
              })
            }
          >
            <option value="list">List</option>
            <option value="flowSet">flowSet</option>
            <option value="SingleCellExperiment">SingleCellExperiment</option>
          </select>
        </label>
      </div>
    </section>
  );
}

function SettingsWorkspace({
  backend,
  job,
  onSave,
}: {
  backend: BackendStatus;
  job: Job;
  onSave: () => void;
}) {
  return (
    <>
      <WorkspaceHeader
        kicker="Settings, logs & advanced tools"
        title="Keep every run explainable"
        description="Project preferences, backend health, and every long-running action live here."
      />
      <div className="settings-grid">
        <section className="surface-card settings-card">
          <SectionTitle
            eyebrow="Backend health"
            title="Local R connection"
            note="The browser only orchestrates. Numerical work stays in Spectreasy."
          />
          <div className="health-status">
            <span
              className={`health-orb ${backend.connected ? "is-connected" : ""}`}
            />
            <div>
              <strong>
                {backend.connected
                  ? "Connected to Spectreasy R"
                  : "Not connected"}
              </strong>
              <span>{backend.message}</span>
            </div>
            <StatusPill
              state={backend.connected ? "connected" : "offline"}
              compact
            />
          </div>
          <div className="settings-list">
            <div>
              <span>Host</span>
              <strong>localhost</strong>
            </div>
            <div>
              <span>API port</span>
              <strong>{backend.apiPort}</strong>
            </div>
            <div>
              <span>Package</span>
              <strong>{backend.version}</strong>
            </div>
            <div>
              <span>GUI assets</span>
              <strong>bundled</strong>
            </div>
          </div>
        </section>
        <section className="surface-card settings-card">
          <SectionTitle eyebrow="Project defaults" title="Lymphocyte panel" />
          <div className="settings-form">
            <label>
              Project path
              <input defaultValue="/projects/lymphocyte-panel" />
            </label>
            <label>
              Default cytometer
              <select defaultValue="Cytek Aurora 5L">
                <option>Cytek Aurora 5L</option>
                <option>Cytek Aurora 4L</option>
                <option>Cytek Northern Lights</option>
              </select>
            </label>
            <label>
              Default report format
              <select defaultValue="HTML">
                <option>HTML</option>
                <option>PDF</option>
              </select>
            </label>
            <label>
              Existing output behavior
              <select defaultValue="New run folder">
                <option>New run folder</option>
                <option>Ask every time</option>
                <option>Overwrite with confirmation</option>
              </select>
            </label>
          </div>
          <button className="button button-primary" onClick={onSave}>
            <Save size={14} /> Save project settings
          </button>
        </section>
      </div>
      <section className="surface-card job-history">
        <div className="card-toolbar">
          <div>
            <span className="eyebrow">Job history</span>
            <h2>Recent actions</h2>
          </div>
          <button className="button button-ghost">
            <Download size={14} /> Download logs
          </button>
        </div>
        <div className="job-history-row">
          <span className="job-history-status complete">
            <Check size={13} />
          </span>
          <div>
            <strong>Control QC report generated</strong>
            <small>run-014 · today, 09:48 · 4m 12s</small>
          </div>
          <span className="job-history-output">9 outputs</span>
          <ArrowRight size={14} />
        </div>
        <div className="job-history-row">
          <span className="job-history-status complete">
            <Check size={13} />
          </span>
          <div>
            <strong>Reference matrix built</strong>
            <small>run-014 · today, 09:44 · 3m 38s</small>
          </div>
          <span className="job-history-output">6 outputs</span>
          <ArrowRight size={14} />
        </div>
        <div className="job-history-row">
          <span className="job-history-status warning">
            <AlertCircle size={13} />
          </span>
          <div>
            <strong>Sample QC report</strong>
            <small>run-013 · yesterday, 16:06 · stale</small>
          </div>
          <span className="job-history-output">1 report</span>
          <ArrowRight size={14} />
        </div>
        {job.state === "running" && (
          <div className="job-history-live">
            <RefreshCcw size={14} className="spin" />
            <span>
              {job.label} · {job.subtask}
            </span>
            <strong>{job.progress}%</strong>
          </div>
        )}
      </section>
    </>
  );
}

export function WorkflowWorkspace(
  props: WorkflowWorkspaceProps & {
    onSectionChange: (section: SectionId) => void;
  },
) {
  const { activeSection, project, backend, job, onSectionChange } = props;
  return (
    <div className="workspace-content">
      <JobStrip job={job} />
      {activeSection === "overview" && (
        <OverviewWorkspace
          project={project}
          onSectionChange={onSectionChange}
          onRefresh={props.onRefresh}
        />
      )}
      {activeSection === "controls" && <ControlsWorkspace {...props} />}
      {activeSection === "samples" && (
        <ConfigurableSamplesWorkspace
          project={project}
          settings={props.settings.sample}
          onSettingsChange={(patch) => props.onSettingsChange("sample", patch)}
          onRun={props.onRun}
          onSectionChange={onSectionChange}
        />
      )}
      {activeSection === "reports" && (
        <ReportsWorkspace onRun={props.onRun} />
      )}
      {activeSection === "matrix" && <MatrixWorkspace />}
      {activeSection === "panel" && (
        <PanelWorkspace
          panelPayload={props.panelPayload}
          onLoadPanel={props.onLoadPanel}
        />
      )}
      {activeSection === "af" && (
        <ConfigurableAfWorkspace
          settings={props.settings.af}
          matrixFile={props.settings.sample.matrixFile}
          onSettingsChange={(patch) => props.onSettingsChange("af", patch)}
          onRun={props.onRun}
        />
      )}
      {(activeSection === "comparison" || activeSection === "simulator") && (
        <ExperimentalWorkspace
          kind={activeSection}
          project={project}
          settings={props.settings}
        />
      )}
      {activeSection === "settings" && (
        <>
          <ConfigurableSettingsWorkspace
            backend={backend}
            job={job}
            settings={props.settings}
            onSettingsChange={props.onSettingsChange}
            onOpenProject={props.onOpenProject}
            onSave={props.onSave}
          />
          <ControlReferenceTuning
            settings={props.settings.control}
            onSettingsChange={(patch) =>
              props.onSettingsChange("control", patch)
            }
          />
          <SampleOutputTuning
            settings={props.settings.sample}
            onSettingsChange={(patch) =>
              props.onSettingsChange("sample", patch)
            }
          />
        </>
      )}
    </div>
  );
}

void LegacyBuildReferencePanel;
void LegacyMatrixWorkspace;
void LegacyPanelWorkspace;
void LegacyConfigurableAfWorkspace;
void LegacyConfigurableSamplesWorkspace;
void LegacyReportsWorkspace;
void LegacyExperimentalWorkspace;
void GatingPanel;
void SamplesWorkspace;
void AfWorkspace;
void SettingsWorkspace;
