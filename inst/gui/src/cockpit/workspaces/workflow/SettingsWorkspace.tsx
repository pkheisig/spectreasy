import { Beaker, CircleCheckBig, FlaskConical } from "lucide-react";
import type { ControlSettings, SampleSettings, WorkflowSettings } from "../../types";
import { defaultWorkflowSettings } from "../../types";
import { AppearanceSettings } from "../../components/AppearanceSettings";
import { GuiSelect } from "../../components/GuiSelect";
import { SettingsCardSummary } from "../../components/SettingsCardSummary";
import { WorkspaceHeader } from "./WorkflowShared";
import type { WorkflowWorkspaceProps } from "./WorkflowShared";

export function ConfigurableSettingsWorkspace({
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
  const useSpectralControlPipeline = control.method === "AutoSpectral";
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
          <SettingsCardSummary icon={<CircleCheckBig size={15} />} title="Control-stage parameters" onReset={() => onSettingsChange("control", { ...defaults.control, sccDir: control.sccDir })} />
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
      <details className="surface-card settings-section" open>
        <summary>
          <SettingsCardSummary icon={<Beaker size={15} />} title="Sample-stage parameters" onReset={() => onSettingsChange("sample", { ...defaults.sample, sampleDir: sample.sampleDir })} />
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
