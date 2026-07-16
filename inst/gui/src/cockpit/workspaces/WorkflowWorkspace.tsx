import { useEffect, useMemo, useState } from "react";
import { createPortal } from "react-dom";
import {
  AlertCircle,
  ArrowRight,
  Beaker,
  CircleCheckBig,
  FlaskConical,
  FolderOpen,
  Layers3,
  Play,
  Plus,
  RefreshCcw,
  Settings2,
  WandSparkles,
  X,
} from "lucide-react";
import type {
  AfSettings,
  Artifact,
  BackendStatus,
  CockpitAppletId,
  ControlSettings,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  SectionId,
  SampleSettings,
  WorkflowSettings,
} from "../types";
import {
  activateAfProfile,
  deactivateAfProfile,
  deleteAfProfile,
  listAfProfiles,
  loadAfProfileData,
  selectAfSourceFile,
} from "../api";
import { StatusPill } from "../components/StatusPill";
import { AppearanceSettings } from "../components/AppearanceSettings";
import { GuiSelect } from "../components/GuiSelect";
import { InlineProjectFiles } from "../components/ProjectFileList";
import { JobStrip } from "../components/JobStrip";
import { ControlReportPanel } from "../components/ControlReportPanel";
import { ResetSettingsButton, SettingsCardSummary } from "../components/SettingsCardSummary";
import { defaultWorkflowSettings } from "../types";

export type WorkflowWorkspaceProps = {
  project: ProjectState;
  backend: BackendStatus;
  job: Job;
  activeSection: SectionId;
  mappingTab: "mapping" | "gating" | "build" | "qc";
  setMappingTab: (tab: "mapping" | "gating" | "build" | "qc") => void;
  onUpdateMapping: (id: string, patch: Partial<MappingRow>) => void;
  onRun: (
    action: "control" | "sample" | "af",
    label: string,
  ) => Promise<boolean>;
  onRefresh: () => void;
  onDownload: (artifact: Artifact) => void;
  onSelectArtifact: (artifact: Artifact) => void;
  onLoadPanel: () => void;
  panelPayload: PanelPayload | null;
  onSave: () => void;
  onCreateMapping: () => void;
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
  onOpenApplet: (applet: CockpitAppletId) => void;
  onSaveMapping?: () => void;
};

function WorkspaceHeader({
  kicker,
  action,
  onAction,
}: {
  kicker: string;
  action?: string;
  onAction?: () => void;
}) {
  return (
    <div className="workspace-header">
      <div>
        <span className="eyebrow">{kicker}</span>
      </div>
      {action && onAction && (
        <button className="button button-ghost" onClick={onAction}>
          <RefreshCcw size={15} /> {action}
        </button>
      )}
    </div>
  );
}

function MappingWorkspace({
  project,
  onUpdateMapping,
  onRun,
  mappingTab,
  setMappingTab,
  onSaveMapping,
  onCreateMapping,
  onRefresh,
  settings,
  onSettingsChange,
  onViewReports,
  onOpenApplet,
}: Pick<
  WorkflowWorkspaceProps,
  | "project"
  | "onUpdateMapping"
  | "onRun"
  | "onRefresh"
  | "mappingTab"
  | "setMappingTab"
  | "onSaveMapping"
  | "onCreateMapping"
  | "settings"
  | "onSettingsChange"
  | "onOpenApplet"
> & { onViewReports: () => void }) {
  const [selectedRows, setSelectedRows] = useState<string[]>([]);
  const [bulkFillColumn, setBulkFillColumn] = useState<"controlType" | "universalNegative" | null>(null);
  const negativeCandidates = useMemo(
    () => project.mapping.filter((row) =>
      row.marker.trim().toLowerCase() === "autofluorescence" || /^af(?:$|_|\b)/i.test(row.fluorophore.trim()),
    ),
    [project.mapping],
  );
  const toggleRow = (id: string) =>
    setSelectedRows((rows) =>
      rows.includes(id) ? rows.filter((row) => row !== id) : [...rows, id],
    );
  const firstMappingRow = project.mapping[0];
  const confirmBulkFill = () => {
    if (!bulkFillColumn || !firstMappingRow) return;
    const value = firstMappingRow[bulkFillColumn];
    project.mapping.forEach((row) => onUpdateMapping(row.id, { [bulkFillColumn]: value }));
    setBulkFillColumn(null);
  };
  return (
    <>
      <div className="subnav">
        <button
          className={mappingTab === "mapping" ? "is-active" : ""}
          onClick={() => setMappingTab("mapping")}
        >
          01 Mapping <StatusPill state={project.mapping.length ? "complete" : "idle"} compact />
        </button>
        <button
          className={mappingTab === "gating" ? "is-active" : ""}
          onClick={() => onOpenApplet("control-gating")}
          disabled={project.mapping.length === 0}
          title={project.mapping.length === 0 ? "Create the control mapping first" : undefined}
        >
          02 Gating <StatusPill state={project.scan.gates > 0 ? "complete" : "ready"} compact />
        </button>
        <button
          className={mappingTab === "build" ? "is-active" : ""}
          onClick={() => setMappingTab("build")}
          disabled={project.scan.gates === 0}
          title={project.scan.gates === 0 ? "Confirm control gates first" : undefined}
        >
          03 Unmixing <StatusPill state={project.scan.matrices > 0 ? "complete" : "ready"} compact />
        </button>
        <button
          className={mappingTab === "qc" ? "is-active" : ""}
          onClick={() => setMappingTab("qc")}
          disabled={project.mapping.length === 0}
        >
          04 Quality Control <StatusPill state={project.scan.reports > 0 ? "complete" : "idle"} compact />
        </button>
      </div>
      {mappingTab === "mapping" && project.mapping.length === 0 && (
        <section className="surface-card mapping-empty-state">
          <div>
            <h2>No control mapping created</h2>
            <p>Create <code>fcs_mapping.csv</code> from the FCS files in the project's <code>scc</code> folder.</p>
          </div>
          <button className="button button-primary create-control-file-button" onClick={onCreateMapping}>
            <Plus size={15} /> Create control file
          </button>
        </section>
      )}
      {mappingTab === "mapping" && project.mapping.length > 0 && (
        <section className="surface-card">
          <div className="card-toolbar">
            <div>
              <h2>
                fcs_mapping.csv
                {project.mappingDirty && <span className="dirty-marker"> · unsaved</span>}
              </h2>
            </div>
            <div className="toolbar-actions">
              <button className="button button-primary" onClick={onSaveMapping}>
                <CircleCheckBig size={14} /> Confirm
              </button>
            </div>
          </div>
          {selectedRows.length > 0 && <div className="mapping-tools">
            <div className="bulk-actions">
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
            </div>
          </div>}
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
                  <th>
                    <span className="mapping-column-head">
                      Type
                      <button type="button" onClick={() => setBulkFillColumn("controlType")} aria-label="Set all type values from the first row" title="Apply first value to all rows">
                        <WandSparkles size={13} />
                      </button>
                    </span>
                  </th>
                  <th>
                    <span className="mapping-column-head">
                      Univ. neg.
                      <button type="button" onClick={() => setBulkFillColumn("universalNegative")} aria-label="Set all universal negative values from the first row" title="Apply first value to all rows">
                        <WandSparkles size={13} />
                      </button>
                    </span>
                  </th>
                </tr>
              </thead>
              <tbody>
                {project.mapping.map((row) => (
                  <tr key={row.id} className={`${row.warning ? "has-warning" : ""} ${row.ignored ? "is-ignored" : ""}`} title={row.ignoredReason || undefined}>
                    <td>
                      <input
                        type="checkbox"
                        checked={selectedRows.includes(row.id)}
                        onChange={() => toggleRow(row.id)}
                        aria-label={`Select ${row.file}`}
                      />
                    </td>
                    <td>
                      <span className="file-cell" title={row.ignoredReason || undefined}>
                        <span className="file-mini-dot" />
                        {row.file}
                        {row.warning && <AlertCircle size={14} className="row-warning" aria-label={row.warning} />}
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
                      <input
                        className="table-input"
                        value={row.channel}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            channel: event.target.value,
                          })
                        }
                      />
                    </td>
                    <td>
                      <GuiSelect
                        className="type-select"
                        value={row.controlType}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            controlType: event.target
                              .value as MappingRow["controlType"],
                          })
                        }
                      >
                        <option value="cell">Cell</option>
                        <option value="bead">Bead</option>
                      </GuiSelect>
                    </td>
                    <td>
                      <GuiSelect
                        className="negative-select"
                        aria-label={`${row.file} universal negative`}
                        value={row.universalNegative}
                        onChange={(event) =>
                          onUpdateMapping(row.id, {
                            universalNegative: event.target.value,
                          })
                        }
                      >
                        <option value="">None</option>
                        {negativeCandidates.map((candidate) => (
                          <option value={candidate.file} key={candidate.file}>{candidate.file}</option>
                        ))}
                      </GuiSelect>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}
      {bulkFillColumn && firstMappingRow && createPortal(
        <div className="cockpit-confirm-overlay" role="presentation" onMouseDown={() => setBulkFillColumn(null)}>
          <div className="cockpit-confirm" role="dialog" aria-modal="true" aria-labelledby="mapping-fill-title" onMouseDown={(event) => event.stopPropagation()}>
            <WandSparkles size={20} />
            <h2 id="mapping-fill-title">Apply the first value to all rows?</h2>
            <p>
              Set every {bulkFillColumn === "controlType" ? "Type" : "Univ. neg."} entry to <strong>{bulkFillColumn === "controlType" ? (firstMappingRow.controlType === "cell" ? "Cell" : "Bead") : (firstMappingRow.universalNegative || "None")}</strong>?
            </p>
            <div>
              <button className="button button-ghost" onClick={() => setBulkFillColumn(null)}>Cancel</button>
              <button className="button button-primary" onClick={confirmBulkFill}>Apply to all</button>
            </div>
          </div>
        </div>,
        document.body,
      )}
      {mappingTab === "mapping" && (
        <InlineProjectFiles
          kind="controls"
          refreshKey={`${project.projectPath}:${project.scan.controls}`}
          onChanged={onRefresh}
        />
      )}
      {mappingTab === "build" && (
        <BuildReferencePanel
          settings={settings.control}
          onSettingsChange={(patch) => onSettingsChange("control", patch)}
          onRun={onRun}
        />
      )}
      {mappingTab === "qc" && (
        <ControlReportPanel
          project={project}
          onView={onViewReports}
        />
      )}
    </>
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
  const method = settings.method;
  const useSpectralPipeline = method === "Spectreasy" || method === "AutoSpectral";
  const [advanced, setAdvanced] = useState(false);
  const defaults = defaultWorkflowSettings("").control;
  return (
    <section className="surface-card run-card streamlined-run-card">
      <div className="settings-card-plain-header">
        <strong>Settings</strong>
        <ResetSettingsButton label="control settings" onReset={() => onSettingsChange({ ...defaults, manualGateFile: "ssc_gate_config.csv" })} />
      </div>
      <div className="run-controls">
        <label>
          <span>Unmixing method</span>
          <GuiSelect
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
          </GuiSelect>
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
          <GuiSelect
            value={settings.outputFormat}
            onChange={(event) =>
              onSettingsChange({
                outputFormat: event.target.value as ControlSettings["outputFormat"],
              })
            }
          >
            <option value="html">HTML</option>
            <option value="pdf">PDF</option>
          </GuiSelect>
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
          <span>Save standalone QC plots</span>
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
          {useSpectralPipeline && <label>
            Background
            <GuiSelect
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
            </GuiSelect>
          </label>}
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
          {useSpectralPipeline && <label>
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
          </label>}
          {method === "Spectreasy" && <label>
            Spectreasy weight quantile
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
          </label>}
        </div>
      )}
      <div className="run-footer run-footer-actions-only">
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

function ControlsWorkspace(
  props: WorkflowWorkspaceProps & {
    onSectionChange: (section: SectionId) => void;
  },
) {
  return (
    <MappingWorkspace
      {...props}
      onViewReports={() => props.onSectionChange("control-reports")}
      onSaveMapping={props.onSaveMapping ?? props.onSave}
    />
  );
}

function ConfigurableSamplesWorkspace({
  project,
  settings,
  onSettingsChange,
  onRun,
  onRefresh,
}: {
  project: ProjectState;
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
  onRefresh: WorkflowWorkspaceProps["onRefresh"];
}) {
  const defaults = defaultWorkflowSettings("").sample;
  return (
    <>
      <section className="surface-card sample-run-card streamlined-run-card">
        <details className="settings-section" open>
          <summary>
            <SettingsCardSummary icon={<Settings2 size={14} />} title="Settings" onReset={() => onSettingsChange(defaults)} />
          </summary>
          <div className="settings-form-grid">
            <label>
              Unmixing method
              <GuiSelect
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
              </GuiSelect>
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
              Spectreasy weight quantile
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
              <GuiSelect
                value={settings.outputFormat}
                onChange={(event) =>
                  onSettingsChange({
                    outputFormat: event.target.value as SampleSettings["outputFormat"],
                  })
                }
              >
                <option value="html">HTML</option>
                <option value="pdf">PDF</option>
              </GuiSelect>
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
              <span>Save standalone QC plots</span>
            </label>
          </div>
        </details>
        <div className="run-footer run-footer-actions-only">
          <button
            className="button button-primary large-button"
            onClick={() => onRun("sample", "Unmix sample workflow")}
          >
            <Play size={15} fill="currentColor" /> Unmix samples{" "}
            <ArrowRight size={15} />
          </button>
        </div>
      </section>
      <InlineProjectFiles
        kind="samples"
        refreshKey={`${project.projectPath}:${project.scan.samples}`}
        onChanged={onRefresh}
      />
    </>
  );
}

function ConfigurableAfWorkspace({
  settings,
  onSettingsChange,
  onRun,
  onRefresh,
}: {
  settings: AfSettings;
  onSettingsChange: (patch: Partial<AfSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
  onRefresh: () => void;
}) {
  const defaults = defaultWorkflowSettings("").af;
  const [profiles, setProfiles] = useState<
    Array<{ name: string; bands: number; detectors: number; created: string; active: boolean }>
  >([]);
  const [preview, setPreview] = useState<Awaited<ReturnType<typeof loadAfProfileData>>>(null);
  const [confirmAction, setConfirmAction] = useState<{ type: "link" | "unlink" | "delete"; name: string } | null>(null);

  const refreshProfiles = async () => {
    const nextProfiles = await listAfProfiles();
    setProfiles(nextProfiles);
    const previewName = nextProfiles.find((profile) => profile.active)?.name ?? nextProfiles[0]?.name;
    setPreview(previewName ? await loadAfProfileData(previewName) : null);
  };

  useEffect(() => {
    const timer = window.setTimeout(() => {
      void refreshProfiles();
    }, 0);
    return () => window.clearTimeout(timer);
  }, []);

  const removeProfile = async (name: string) => {
    const removed = await deleteAfProfile(name);
    if (removed) await refreshProfiles();
  };

  const linkProfile = async (name: string) => {
    const result = await activateAfProfile(name);
    if (result.success) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const unlinkProfile = async (name: string) => {
    const removed = await deactivateAfProfile(name);
    if (removed) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const chooseSource = async () => {
    const result = await selectAfSourceFile();
    if (result.success && result.path) onSettingsChange({ fcsFile: result.path });
  };

  const extractProfile = async () => {
    const saved = await onRun("af", "Extract AF profile");
    if (saved) await refreshProfiles();
  };

  const chartWidth = Math.max(760, (preview?.detectors.length ?? 0) * 22 + 92);
  const chartHeight = 355;
  const chartLeft = 58;
  const chartRight = 18;
  const chartTop = 22;
  const chartBottom = 92;
  const plotWidth = chartWidth - chartLeft - chartRight;
  const plotHeight = chartHeight - chartTop - chartBottom;
  const detectorX = (index: number, length: number) => chartLeft + (index / Math.max(1, length - 1)) * plotWidth;
  const spectrumPath = (values: number[]) => values.map((value, index) => {
    const x = detectorX(index, values.length);
    const y = chartTop + (1 - Math.max(0, Math.min(1, value))) * plotHeight;
    return `${index === 0 ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)}`;
  }).join(" ");

  return (
    <>
      <WorkspaceHeader
        kicker="AF profile library"
      />
      <div className="af-layout">
        <section className="surface-card af-extract-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Extract profile</span>
              <h2>{settings.fcsFile}</h2>
            </div>
            <div className="toolbar-actions">
              <ResetSettingsButton label="AF settings" onReset={() => onSettingsChange(defaults)} />
              <StatusPill state="ready" label={`${settings.afNBands} bands`} />
            </div>
          </div>
          <div className="af-form">
            <label>
              Source FCS
              <span className="file-picker-field">
                <input value={settings.fcsFile} readOnly placeholder="Choose an unstained FCS file" />
                <button type="button" className="button button-ghost" onClick={() => void chooseSource()}><FolderOpen size={14} /> Browse</button>
              </span>
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
              Downsampled events
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
          <button
            className="button button-primary large-button"
            onClick={() => void extractProfile()}
          >
            <WandSparkles size={15} /> Extract profile{" "}
            {settings.saveName ? `& save “${settings.saveName}”` : ""}
          </button>
        </section>
        <section className="surface-card af-library-card">
          {profiles.length === 0 ? (
            <div className="profile-list" aria-label="No saved AF profiles" />
          ) : (
            <div className="profile-list">
              {profiles.map((profile) => (
                <div className={`saved-profile ${preview?.name === profile.name ? "is-selected" : ""}`} key={profile.name}>
                  <button
                    className="profile-copy profile-preview-button"
                    aria-pressed={preview?.name === profile.name}
                    onClick={() => void loadAfProfileData(profile.name).then(setPreview)}
                  >
                    <strong>{profile.name}</strong>
                    <span>
                      {profile.bands} bands · {profile.detectors} detectors
                    </span>
                    <small>
                      {profile.created ||
                        "Saved in the local R profile library"}
                    </small>
                  </button>
                  <div className="profile-actions profile-actions-stacked">
                    <button
                      className="text-action"
                      onClick={() => setConfirmAction({ type: profile.active ? "unlink" : "link", name: profile.name })}
                    >
                      <Layers3 size={14} /> {profile.active ? "Unlink from dataset" : "Link to dataset"}
                    </button>
                    <button
                      className="text-action danger"
                      onClick={() => setConfirmAction({ type: "delete", name: profile.name })}
                    >
                      <X size={14} /> Delete
                    </button>
                  </div>
                </div>
              ))}
            </div>
          )}
        </section>
        {preview && preview.spectra.length > 0 && <section className="surface-card af-spectrum-card">
          <div className="card-toolbar">
            <div>
              <span className="eyebrow">Spectrum preview</span>
              <h2>{preview.name}</h2>
              <p>{preview.spectra.length} AF spectra across {preview.detectors.length} detectors</p>
            </div>
          </div>
          <div className="af-spectrum-scroll">
            <svg
              className="af-spectrum-chart"
              viewBox={`0 0 ${chartWidth} ${chartHeight}`}
              role="img"
              aria-label={`${preview.name}: ${preview.spectra.length} AF spectra across ${preview.detectors.length} detectors`}
            >
              {[0, .2, .4, .6, .8, 1].map((tick) => {
                const y = chartTop + (1 - tick) * plotHeight;
                return <g key={tick}>
                  <line x1={chartLeft} x2={chartWidth - chartRight} y1={y} y2={y} className="af-grid-line" />
                  <text x={chartLeft - 9} y={y + 3} textAnchor="end" className="af-tick-label">{tick.toFixed(1)}</text>
                </g>;
              })}
              {preview.detectors.map((detector, index) => {
                const x = detectorX(index, preview.detectors.length);
                return <g key={detector}>
                  <line x1={x} x2={x} y1={chartTop} y2={chartTop + plotHeight} className="af-grid-line af-grid-line-vertical" />
                  <text transform={`translate(${x + 3} ${chartTop + plotHeight + 12}) rotate(68)`} textAnchor="start" className="af-detector-label">{detector}</text>
                </g>;
              })}
              {preview.spectra.map((spectrum) => <path
                key={spectrum.name}
                d={spectrumPath(spectrum.values)}
                fill="none"
                stroke="var(--cockpit-accent)"
                strokeWidth="1.35"
                strokeLinecap="round"
                strokeLinejoin="round"
                opacity=".2"
                vectorEffect="non-scaling-stroke"
              />)}
              <line x1={chartLeft} x2={chartWidth - chartRight} y1={chartTop + plotHeight} y2={chartTop + plotHeight} className="af-axis-line" />
              <line x1={chartLeft} x2={chartLeft} y1={chartTop} y2={chartTop + plotHeight} className="af-axis-line" />
              <text x={chartLeft + plotWidth / 2} y={chartHeight - 8} textAnchor="middle" className="af-axis-title">Detector</text>
              <text transform={`translate(15 ${chartTop + plotHeight / 2}) rotate(-90)`} textAnchor="middle" className="af-axis-title">Normalized intensity</text>
            </svg>
          </div>
        </section>}
      </div>
      {confirmAction && createPortal(<div className="cockpit-confirm-overlay" role="presentation" onMouseDown={() => setConfirmAction(null)}>
        <div className="cockpit-confirm" role="dialog" aria-modal="true" aria-labelledby="af-confirm-title" onMouseDown={(event) => event.stopPropagation()}>
          <h2 id="af-confirm-title">{confirmAction.type === "link" ? "Link this AF profile?" : confirmAction.type === "unlink" ? "Unlink this AF profile?" : "Delete this AF profile?"}</h2>
          <p>{confirmAction.type === "link" ? `Use ${confirmAction.name} as this dataset's unstained cell control? The mapped unstained cell file will be ignored.` : confirmAction.type === "unlink" ? `${confirmAction.name} will no longer replace this dataset's mapped unstained cell control.` : `${confirmAction.name} will be permanently removed from the local profile library.`}</p>
          <div>
            <button className="button button-ghost" onClick={() => setConfirmAction(null)}>Cancel</button>
            <button className={`button ${confirmAction.type === "delete" ? "button-danger" : "button-primary"}`} onClick={() => {
              const action = confirmAction;
              setConfirmAction(null);
              if (action.type === "link") void linkProfile(action.name);
              else if (action.type === "unlink") void unlinkProfile(action.name);
              else void removeProfile(action.name);
            }}>{confirmAction.type === "link" ? "Link to dataset" : confirmAction.type === "unlink" ? "Unlink from dataset" : "Delete"}</button>
          </div>
        </div>
      </div>, document.body)}
    </>
  );
}

function ConfigurableSettingsWorkspace({
  settings,
  onSettingsChange,
}: {
  settings: WorkflowSettings;
  onSettingsChange: WorkflowWorkspaceProps["onSettingsChange"];
}) {
  const defaults = defaultWorkflowSettings(settings.projectPath);
  const control = settings.control;
  const sample = settings.sample;
  const af = settings.af;
  const useSpectralControlPipeline = control.method === "Spectreasy" || control.method === "AutoSpectral";
  return (
    <>
      <WorkspaceHeader
        kicker="Settings"
      />
      <AppearanceSettings
        value={settings.appearance}
        onChange={(patch) => onSettingsChange("appearance", patch)}
        onReset={() => onSettingsChange("appearance", defaults.appearance)}
      />
      <details className="surface-card settings-section" open>
        <summary>
          <SettingsCardSummary icon={<CircleCheckBig size={15} />} title="Control-stage parameters" onReset={() => onSettingsChange("control", defaults.control)} />
        </summary>
        <div className="settings-form-grid">
          <label>
            Cytometer
            <GuiSelect
              value={control.cytometer}
              onChange={(event) =>
                onSettingsChange("control", { cytometer: event.target.value })
              }
            >
              <option value="auto">Auto-detect</option>
              <option value="aurora">Cytek Aurora</option>
              <option value="northern_lights">Cytek Northern Lights</option>
              <option value="id7000">Sony ID7000</option>
              <option value="discover_s8">BD FACSDiscover S8</option>
              <option value="discover_a8">BD FACSDiscover A8</option>
              <option value="a5se">BD FACSymphony A5 SE</option>
              <option value="opteon">Agilent NovoCyte Opteon</option>
              <option value="mosaic">Beckman Coulter CytoFLEX Mosaic</option>
              <option value="xenith">Thermo Fisher Attune Xenith</option>
            </GuiSelect>
          </label>
          <label>
            Unmixing method
            <GuiSelect
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
            </GuiSelect>
          </label>
          <label>
            Unknown fluor policy
            <GuiSelect
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
            </GuiSelect>
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
          {useSpectralControlPipeline && <>
          <label>
            Background method
            <GuiSelect
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
            </GuiSelect>
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
            Spectreasy weight quantile
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
          </>}
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
            <span>Save standalone QC plots</span>
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
            <GuiSelect
              value={control.outputFormat}
              onChange={(event) =>
                onSettingsChange("control", {
                  outputFormat: event.target.value as ControlSettings["outputFormat"],
                })
              }
            >
              <option value="html">HTML</option>
              <option value="pdf">PDF</option>
            </GuiSelect>
          </label>
          {control.method === "AutoSpectral" && <label className="toggle-label">
            <input
              type="checkbox"
              checked={control.refine}
              onChange={(event) =>
                onSettingsChange("control", { refine: event.target.checked })
              }
            />
            <span className="toggle-ui" />
            <span>Refine AutoSpectral matrix</span>
          </label>}
        </div>
      </details>
      <details className="surface-card settings-section">
        <summary>
          <SettingsCardSummary icon={<Beaker size={15} />} title="Sample-stage parameters" onReset={() => onSettingsChange("sample", defaults.sample)} />
        </summary>
        <div className="settings-form-grid">
          <label>
            Method
            <GuiSelect
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
            </GuiSelect>
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
            Spectreasy weight quantile
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
            <GuiSelect
              value={sample.outputFormat}
              onChange={(event) =>
                onSettingsChange("sample", {
                  outputFormat: event.target.value as SampleSettings["outputFormat"],
                })
              }
            >
              <option value="html">HTML</option>
              <option value="pdf">PDF</option>
            </GuiSelect>
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
            <span>Save standalone QC plots</span>
          </label>
        </div>
      </details>
      <details className="surface-card settings-section">
        <summary>
          <SettingsCardSummary icon={<FlaskConical size={15} />} title="AF profile parameters" detail="extract_af_profile" onReset={() => onSettingsChange("af", defaults.af)} />
        </summary>
        <div className="settings-form-grid">
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
            Downsampled events
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
    </>
  );
}

function ControlReferenceTuning({
  settings,
  onSettingsChange,
  onReset,
}: {
  settings: ControlSettings;
  onSettingsChange: (patch: Partial<ControlSettings>) => void;
  onReset: () => void;
}) {
  return (
    <section className="surface-card settings-section">
      <div className="settings-section-title">
        <div>
          <span className="eyebrow">Reference builder</span>
          <h2>Gating & clustering parameters</h2>
        </div>
        <ResetSettingsButton label="gating and clustering parameters" onReset={onReset} />
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
          <GuiSelect
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
          </GuiSelect>
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
          <GuiSelect
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
          </GuiSelect>
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
          <GuiSelect
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
          </GuiSelect>
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
  onReset,
}: {
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
  onReset: () => void;
}) {
  return (
    <section className="surface-card settings-section">
      <div className="settings-section-title">
        <div>
          <span className="eyebrow">Sample outputs</span>
          <h2>Library & return type</h2>
        </div>
        <ResetSettingsButton label="sample output parameters" onReset={onReset} />
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
          <GuiSelect
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
          </GuiSelect>
        </label>
      </div>
    </section>
  );
}

export function WorkflowWorkspace(
  props: WorkflowWorkspaceProps & {
    onSectionChange: (section: SectionId) => void;
  },
) {
  const { activeSection, job } = props;
  return (
    <div className="workspace-content">
      <JobStrip job={job} />
      {activeSection === "controls" && <ControlsWorkspace {...props} />}
      {activeSection === "samples" && (
        <ConfigurableSamplesWorkspace
          project={props.project}
          settings={props.settings.sample}
          onSettingsChange={(patch) => props.onSettingsChange("sample", patch)}
          onRun={props.onRun}
          onRefresh={props.onRefresh}
        />
      )}
      {activeSection === "af" && (
        <ConfigurableAfWorkspace
          settings={props.settings.af}
          onSettingsChange={(patch) => props.onSettingsChange("af", patch)}
          onRun={props.onRun}
          onRefresh={props.onRefresh}
        />
      )}
      {activeSection === "settings" && (
        <>
          <ConfigurableSettingsWorkspace
            settings={props.settings}
            onSettingsChange={props.onSettingsChange}
          />
          <ControlReferenceTuning
            settings={props.settings.control}
            onSettingsChange={(patch) =>
              props.onSettingsChange("control", patch)
            }
            onReset={() => props.onSettingsChange("control", defaultWorkflowSettings(props.settings.projectPath).control)}
          />
          <SampleOutputTuning
            settings={props.settings.sample}
            onSettingsChange={(patch) =>
              props.onSettingsChange("sample", patch)
            }
            onReset={() => props.onSettingsChange("sample", defaultWorkflowSettings(props.settings.projectPath).sample)}
          />
        </>
      )}
    </div>
  );
}
