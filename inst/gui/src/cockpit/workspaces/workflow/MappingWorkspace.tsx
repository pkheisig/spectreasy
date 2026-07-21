import { useMemo, useState } from "react";
import { createPortal } from "react-dom";
import { AlertCircle, CircleCheckBig, Plus, WandSparkles } from "lucide-react";
import type { MappingRow } from "../../types";
import { GuiSelect } from "../../components/GuiSelect";
import { InlineProjectFiles } from "../../components/ProjectFileList";
import { QcReportPanel } from "../../components/ControlReportPanel";
import { StatusPill } from "../../components/StatusPill";
import { BuildReferencePanel } from "./BuildReferencePanel";
import type { WorkflowWorkspaceProps } from "./WorkflowShared";

export function MappingWorkspace({
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
  onInputDirectoryChange,
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
  | "onInputDirectoryChange"
  | "onOpenApplet"
> & { onViewReports: (reportPath: string) => void }) {
  const [selectedRows, setSelectedRows] = useState<string[]>([]);
  const [bulkFillColumn, setBulkFillColumn] = useState<"controlType" | "universalNegative" | null>(null);
  const negativeCandidates = useMemo(
    () => project.mapping.filter((row) =>
      row.marker.trim().toLowerCase() === "autofluorescence" || /^af(?:$|_|\b)/i.test(row.fluorophore.trim()),
    ),
    [project.mapping],
  );
  const controlReportCount = project.artifacts.filter((artifact) => artifact.type === "Control QC report").length;
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
        >
          04 Quality Control <StatusPill state={controlReportCount > 0 ? "complete" : "idle"} compact />
        </button>
      </div>
      {mappingTab === "mapping" && project.mapping.length === 0 && (
        <section className="surface-card mapping-empty-state">
          <div>
            <h2>No control mapping created</h2>
            <p>Create <code>fcs_mapping.csv</code> from the FCS files in <code>{settings.control.sccDir}</code>.</p>
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
          projectPath={project.projectPath}
          directory={project.controlInputDir}
          refreshKey={`${project.projectPath}:${project.scan.controls}`}
          onChanged={onRefresh}
        />
      )}
      {mappingTab === "build" && (
        <>
          <BuildReferencePanel
            settings={settings.control}
            onSettingsChange={(patch) => onSettingsChange("control", patch)}
            onInputDirectoryChange={onInputDirectoryChange}
            onRun={onRun}
          />
          <InlineProjectFiles
            kind="controls"
            projectPath={project.projectPath}
            directory={project.controlInputDir}
            refreshKey={`${project.projectPath}:${project.scan.controls}`}
            onChanged={onRefresh}
          />
        </>
      )}
      {mappingTab === "qc" && (
        <QcReportPanel
          project={project}
          kind="control"
          onView={onViewReports}
          onAiQc={() => onOpenApplet("ai-ready-qc")}
        />
      )}
    </>
  );
}
