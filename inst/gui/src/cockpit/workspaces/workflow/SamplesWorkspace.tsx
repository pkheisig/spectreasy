import { useState } from "react";
import { ArrowRight, Play, Settings2 } from "lucide-react";
import type { ProjectState, SampleSettings, SectionId } from "../../types";
import { defaultWorkflowSettings } from "../../types";
import { GuiSelect } from "../../components/GuiSelect";
import { InlineProjectFiles } from "../../components/ProjectFileList";
import { ProjectDirectoryField } from "../../components/ProjectDirectoryField";
import { QcReportPanel } from "../../components/ControlReportPanel";
import { ResetSettingsButton } from "../../components/SettingsCardSummary";
import { StatusPill } from "../../components/StatusPill";
import { MappingWorkspace } from "./MappingWorkspace";
import type { WorkflowWorkspaceProps } from "./WorkflowShared";

export function ControlsWorkspace(
  props: WorkflowWorkspaceProps & {
    onSectionChange: (section: SectionId) => void;
  },
) {
  return (
    <MappingWorkspace
      {...props}
      onViewReports={(reportPath) => props.onOpenApplet("control-qc-report", reportPath)}
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
  onInputDirectoryChange,
}: {
  project: ProjectState;
  settings: SampleSettings;
  onSettingsChange: (patch: Partial<SampleSettings>) => void;
  onRun: WorkflowWorkspaceProps["onRun"];
  onRefresh: WorkflowWorkspaceProps["onRefresh"];
  onInputDirectoryChange: WorkflowWorkspaceProps["onInputDirectoryChange"];
}) {
  const defaults = defaultWorkflowSettings("").sample;
  const [advanced, setAdvanced] = useState(false);
  return (
    <>
      <section className="surface-card run-card streamlined-run-card">
        <div className="settings-card-plain-header">
          <strong>Settings</strong>
          <ResetSettingsButton label="sample settings" onReset={() => onSettingsChange({ ...defaults, sampleDir: settings.sampleDir })} />
        </div>
        <div className="run-controls">
          <label>
            <span>Unmixing method</span>
            <GuiSelect value={settings.method} onChange={(event) => onSettingsChange({ method: event.target.value })}>
              <option>Spectreasy</option>
              <option>AutoSpectral</option>
              <option>OLS</option>
              <option>WLS</option>
              <option>RWLS</option>
              <option>NNLS</option>
            </GuiSelect>
          </label>
          <label className="toggle-label">
            <input type="checkbox" checked={settings.saveReport} onChange={(event) => onSettingsChange({ saveReport: event.target.checked })} />
            <span className="toggle-ui" />
            <span>Generate report</span>
          </label>
          <label>
            <span>Report format</span>
            <GuiSelect value={settings.outputFormat} onChange={(event) => onSettingsChange({ outputFormat: event.target.value as SampleSettings["outputFormat"] })}>
              <option value="html">HTML</option>
              <option value="pdf">PDF</option>
            </GuiSelect>
          </label>
          <label className="toggle-label">
            <input type="checkbox" checked={settings.saveQcPlots} onChange={(event) => onSettingsChange({ saveQcPlots: event.target.checked })} />
            <span className="toggle-ui" />
            <span>Save standalone QC plots</span>
          </label>
        </div>
        <button className="advanced-toggle" onClick={() => setAdvanced(!advanced)}>
          <Settings2 size={15} /> Advanced settings <span>{advanced ? "−" : "+"}</span>
        </button>
        {advanced && (
          <div className="advanced-grid">
            <ProjectDirectoryField
              key={`samples-${settings.sampleDir}`}
              label="Sample FCS folder"
              role="samples"
              value={settings.sampleDir}
              onCommit={onInputDirectoryChange}
            />
            <label>Unmixing matrix<input value={settings.matrixFile} onChange={(event) => onSettingsChange({ matrixFile: event.target.value })} /></label>
            <label>Detector noise file<input value={settings.detectorNoiseFile} onChange={(event) => onSettingsChange({ detectorNoiseFile: event.target.value })} placeholder="Optional" /></label>
            <label>Output folder<input value={settings.outputDir} onChange={(event) => onSettingsChange({ outputDir: event.target.value })} /></label>
            <label>Threads<input type="number" min="1" value={settings.nThreads} onChange={(event) => onSettingsChange({ nThreads: Number(event.target.value) })} /></label>
            <label>Chunk size<input type="number" min="1" value={settings.chunkSize} onChange={(event) => onSettingsChange({ chunkSize: Number(event.target.value) })} /></label>
            <label>Plot events<input type="number" min="1" value={settings.plotNEvents} onChange={(event) => onSettingsChange({ plotNEvents: Number(event.target.value) })} /></label>
            {settings.method === "RWLS" && <label>RWLS max iterations<input type="number" min="1" value={settings.rwlsMaxIter} onChange={(event) => onSettingsChange({ rwlsMaxIter: Number(event.target.value) })} /></label>}
            {(settings.method === "Spectreasy" || settings.method === "AutoSpectral") && <label>Variant top-k<input type="number" min="1" value={settings.spectralVariantTopK} onChange={(event) => onSettingsChange({ spectralVariantTopK: Number(event.target.value) })} /></label>}
            {(settings.method === "Spectreasy" || settings.method === "AutoSpectral") && <label>Variant min abundance<input type="number" min="0" step="0.01" value={settings.spectralVariantMinAbundance} onChange={(event) => onSettingsChange({ spectralVariantMinAbundance: Number(event.target.value) })} /></label>}
            {(settings.method === "Spectreasy" || settings.method === "AutoSpectral") && <label>Variant positive fraction<input type="number" min="0" max="1" step="0.01" value={settings.spectralVariantPositiveFraction} onChange={(event) => onSettingsChange({ spectralVariantPositiveFraction: Number(event.target.value) })} /></label>}
            {(settings.method === "Spectreasy" || settings.method === "AutoSpectral") && <label>Variant min improvement<input type="number" min="0" step="0.01" value={settings.spectralVariantMinImprovement} onChange={(event) => onSettingsChange({ spectralVariantMinImprovement: Number(event.target.value) })} /></label>}
            {(settings.method === "Spectreasy" || settings.method === "AutoSpectral") && <label>Spectral variant library<input value={settings.spectralVariantLibraryFile} onChange={(event) => onSettingsChange({ spectralVariantLibraryFile: event.target.value })} placeholder="Optional .rds file" /></label>}
            {settings.method === "Spectreasy" && <label>Spectreasy weight quantile<input type="number" min="0" max="1" step="0.01" value={settings.spectreasyWeightQuantile} onChange={(event) => onSettingsChange({ spectreasyWeightQuantile: Number(event.target.value) })} /></label>}
            <label>Seed<input type="number" min="1" value={settings.seed} onChange={(event) => onSettingsChange({ seed: Number(event.target.value) })} /></label>
            <label>
              Return type
              <GuiSelect value={settings.returnType} onChange={(event) => onSettingsChange({ returnType: event.target.value as SampleSettings["returnType"] })}>
                <option value="list">List</option>
                <option value="flowSet">flowSet</option>
                <option value="SingleCellExperiment">SingleCellExperiment</option>
              </GuiSelect>
            </label>
            <label className="toggle-label">
              <input type="checkbox" checked={settings.reportPerSample} onChange={(event) => onSettingsChange({ reportPerSample: event.target.checked })} />
              <span className="toggle-ui" /><span>One QC report per sample</span>
            </label>
            <label className="toggle-label">
              <input type="checkbox" checked={settings.estimateAf} onChange={(event) => onSettingsChange({ estimateAf: event.target.checked })} />
              <span className="toggle-ui" /><span>Estimate AF from samples</span>
            </label>
            <label className="toggle-label">
              <input type="checkbox" checked={settings.writeFcs} onChange={(event) => onSettingsChange({ writeFcs: event.target.checked })} />
              <span className="toggle-ui" /><span>Write FCS outputs</span>
            </label>
          </div>
        )}
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
        projectPath={project.projectPath}
        directory={project.sampleInputDir}
        refreshKey={`${project.projectPath}:${project.scan.samples}`}
        onChanged={onRefresh}
      />
    </>
  );
}

export function SamplesWorkspace(props: WorkflowWorkspaceProps) {
  const sampleReports = props.project.artifacts.filter((artifact) => artifact.type === "Sample QC report");
  return (
    <>
      <div className="subnav">
        <button
          className={props.sampleTab === "unmixing" ? "is-active" : ""}
          onClick={() => props.setSampleTab("unmixing")}
        >
          01 Unmixing <StatusPill state={props.project.scan.samples > 0 ? "ready" : "idle"} compact />
        </button>
        <button
          className={props.sampleTab === "qc" ? "is-active" : ""}
          onClick={() => props.setSampleTab("qc")}
        >
          02 Quality Control <StatusPill state={sampleReports.length > 0 ? "complete" : "idle"} compact />
        </button>
      </div>
      {props.sampleTab === "unmixing" ? (
        <ConfigurableSamplesWorkspace
          project={props.project}
          settings={props.settings.sample}
          onSettingsChange={(patch) => props.onSettingsChange("sample", patch)}
          onRun={props.onRun}
          onRefresh={props.onRefresh}
          onInputDirectoryChange={props.onInputDirectoryChange}
        />
      ) : (
        <QcReportPanel
          project={props.project}
          kind="sample"
          onView={(reportPath) => props.onOpenApplet("sample-qc-report", reportPath)}
        />
      )}
    </>
  );
}
