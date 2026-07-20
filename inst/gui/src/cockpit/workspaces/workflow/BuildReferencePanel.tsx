import { useState } from "react";
import { ArrowRight, Play, Settings2 } from "lucide-react";
import type { ControlSettings } from "../../types";
import { defaultWorkflowSettings } from "../../types";
import { GuiSelect } from "../../components/GuiSelect";
import { ProjectDirectoryField } from "../../components/ProjectDirectoryField";
import { ResetSettingsButton } from "../../components/SettingsCardSummary";
import type { WorkflowWorkspaceProps } from "./WorkflowShared";

export function BuildReferencePanel({
  settings,
  onSettingsChange,
  onInputDirectoryChange,
  onRun,
}: {
  settings: ControlSettings;
  onSettingsChange: (patch: Partial<ControlSettings>) => void;
  onInputDirectoryChange: WorkflowWorkspaceProps["onInputDirectoryChange"];
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
        <ResetSettingsButton label="control settings" onReset={() => onSettingsChange({ ...defaults, sccDir: settings.sccDir, manualGateFile: "ssc_gate_config.csv" })} />
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
          <ProjectDirectoryField
            key={`controls-${settings.sccDir}`}
            label="Control FCS folder"
            role="controls"
            value={settings.sccDir}
            onCommit={onInputDirectoryChange}
          />
          <label>
            Output folder
            <input
              value={settings.outputDir}
              onChange={(event) => onSettingsChange({ outputDir: event.target.value })}
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
          {method === "RWLS" && <label>
            RWLS max iterations
            <input
              type="number"
              min="1"
              value={settings.rwlsMaxIter}
              onChange={(event) =>
                onSettingsChange({ rwlsMaxIter: Number(event.target.value) })
              }
            />
          </label>}
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
            AF maximum cells
            <input
              type="number"
              min="1"
              value={settings.afMaxCells}
              onChange={(event) =>
                onSettingsChange({ afMaxCells: Number(event.target.value) })
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
          {useSpectralPipeline && settings.sccBackgroundMethod === "scatter_knn" && <label>
            Background neighbours
            <input
              type="number"
              min="1"
              value={settings.sccBackgroundK}
              onChange={(event) =>
                onSettingsChange({ sccBackgroundK: Number(event.target.value) })
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
          {useSpectralPipeline && <label>
            Variant SOM nodes
            <input
              type="number"
              min="1"
              value={settings.spectralVariantSomNodes}
              onChange={(event) =>
                onSettingsChange({ spectralVariantSomNodes: Number(event.target.value) })
              }
            />
          </label>}
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
          {useSpectralPipeline && <label>
            Variant cosine threshold
            <input
              type="number"
              min="0"
              max="1"
              step="0.01"
              value={settings.spectralVariantCosineThreshold}
              onChange={(event) =>
                onSettingsChange({
                  spectralVariantCosineThreshold: Number(event.target.value),
                })
              }
            />
          </label>}
          {useSpectralPipeline && <label>
            Maximum variants
            <input
              type="number"
              min="1"
              value={settings.spectralVariantMaxVariants}
              onChange={(event) =>
                onSettingsChange({ spectralVariantMaxVariants: Number(event.target.value) })
              }
            />
          </label>}
          {useSpectralPipeline && <label>
            Minimum variant events
            <input
              type="number"
              min="1"
              value={settings.spectralVariantMinEvents}
              onChange={(event) =>
                onSettingsChange({ spectralVariantMinEvents: Number(event.target.value) })
              }
            />
          </label>}
          {useSpectralPipeline && <label>
            Candidate events
            <input
              type="number"
              min="1"
              value={settings.autospectralNCandidates}
              onChange={(event) =>
                onSettingsChange({ autospectralNCandidates: Number(event.target.value) })
              }
            />
          </label>}
          {useSpectralPipeline && <label>
            Spectrum events
            <input
              type="number"
              min="1"
              value={settings.autospectralNSpectral}
              onChange={(event) =>
                onSettingsChange({ autospectralNSpectral: Number(event.target.value) })
              }
            />
          </label>}
          {useSpectralPipeline && <label>
            Minimum selector events
            <input
              type="number"
              min="1"
              value={settings.autospectralMinEvents}
              onChange={(event) =>
                onSettingsChange({ autospectralMinEvents: Number(event.target.value) })
              }
            />
          </label>}
          <label>
            Scatter panel size (mm)
            <input
              type="number"
              min="1"
              value={settings.unmixScatterPanelSizeMm}
              onChange={(event) =>
                onSettingsChange({ unmixScatterPanelSizeMm: Number(event.target.value) })
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
          {method === "AutoSpectral" && <label className="toggle-label">
            <input
              type="checkbox"
              checked={settings.refine}
              onChange={(event) => onSettingsChange({ refine: event.target.checked })}
            />
            <span className="toggle-ui" />
            <span>Refine AF bank</span>
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
