import type { ControlSettings, SampleSettings } from "../../types";
import { GuiSelect } from "../../components/GuiSelect";
import { ResetSettingsButton } from "../../components/SettingsCardSummary";

export function ControlReferenceTuning({
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

export function SampleOutputTuning({
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
