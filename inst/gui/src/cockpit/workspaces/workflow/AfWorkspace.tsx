import { useCallback, useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { FolderOpen, Layers3, Pencil, WandSparkles, X } from "lucide-react";
import type { AfSettings } from "../../types";
import { defaultWorkflowSettings } from "../../types";
import {
  activateAfProfile,
  deactivateAfProfile,
  deleteAfProfile,
  listAfProfiles,
  loadAfProfileData,
  renameAfProfile,
  selectAfSourceFile,
} from "../../api";
import { ResetSettingsButton } from "../../components/SettingsCardSummary";
import { StatusPill } from "../../components/StatusPill";
import { WorkspaceHeader } from "./WorkflowShared";
import type { WorkflowWorkspaceProps } from "./WorkflowShared";

export function ConfigurableAfWorkspace({
  projectPath,
  settings,
  onSettingsChange,
  onRun,
  onRefresh,
}: {
  projectPath: string;
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
  const [confirmAction, setConfirmAction] = useState<{ type: "link" | "unlink" | "delete" | "overwrite"; name: string } | null>(null);
  const [renameAction, setRenameAction] = useState<{ original: string; value: string; error: string } | null>(null);
  const refreshGeneration = useRef(0);

  const refreshProfiles = useCallback(async () => {
    const generation = ++refreshGeneration.current;
    const nextProfiles = await listAfProfiles(projectPath);
    if (generation !== refreshGeneration.current) return;
    setProfiles(nextProfiles);
    const previewName = nextProfiles.find((profile) => profile.active)?.name ?? nextProfiles[0]?.name;
    const nextPreview = previewName ? await loadAfProfileData(previewName) : null;
    if (generation !== refreshGeneration.current) return;
    setPreview(nextPreview);
  }, [projectPath]);

  useEffect(() => {
    void refreshProfiles();
    return () => {
      refreshGeneration.current += 1;
    };
  }, [refreshProfiles]);

  const removeProfile = async (name: string) => {
    const removed = await deleteAfProfile(name, projectPath);
    if (removed) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const linkProfile = async (name: string) => {
    const result = await activateAfProfile(name, projectPath);
    if (result.success) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const unlinkProfile = async (name: string) => {
    const removed = await deactivateAfProfile(name, projectPath);
    if (removed) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const renameProfile = async () => {
    if (!renameAction) return;
    const newName = renameAction.value.trim();
    if (!newName) {
      setRenameAction({ ...renameAction, error: "Enter a profile name." });
      return;
    }
    if (newName === renameAction.original) {
      setRenameAction({ ...renameAction, error: "Enter a different profile name." });
      return;
    }
    const result = await renameAfProfile(renameAction.original, newName, projectPath);
    if (!result.success) {
      setRenameAction({ ...renameAction, error: result.message });
      return;
    }
    setRenameAction(null);
    await refreshProfiles();
    setPreview(await loadAfProfileData(result.name ?? newName));
    onRefresh();
  };

  const chooseSource = async () => {
    const result = await selectAfSourceFile(projectPath);
    if (result.success && result.path) onSettingsChange({ fcsFile: result.path });
  };

  const runExtractProfile = async () => {
    const saved = await onRun("af", "Extract AF profile");
    if (saved) {
      await refreshProfiles();
      onRefresh();
    }
  };

  const extractProfile = () => {
    if (settings.saveOverwrite) {
      setConfirmAction({ type: "overwrite", name: settings.saveName || settings.fcsFile });
      return;
    }
    void runExtractProfile();
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
                placeholder="Defaults to source filename"
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
            onClick={extractProfile}
          >
            <WandSparkles size={15} /> Extract & save profile
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
                      className="text-action"
                      onClick={() => setRenameAction({ original: profile.name, value: profile.name, error: "" })}
                    >
                      <Pencil size={14} /> Rename
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
          <h2 id="af-confirm-title">{confirmAction.type === "link" ? "Link this AF profile?" : confirmAction.type === "unlink" ? "Unlink this AF profile?" : confirmAction.type === "overwrite" ? "Overwrite this AF profile?" : "Delete this AF profile?"}</h2>
          <p>{confirmAction.type === "link" ? `Use ${confirmAction.name} as this dataset's unstained cell control? The mapped unstained cell file will be ignored.` : confirmAction.type === "unlink" ? `${confirmAction.name} will no longer replace this dataset's mapped unstained cell control.` : confirmAction.type === "overwrite" ? `The existing ${confirmAction.name} profile will be replaced with the newly extracted profile.` : `${confirmAction.name} will be permanently removed from the package root.`}</p>
          <div>
            <button className="button button-ghost" onClick={() => setConfirmAction(null)}>Cancel</button>
            <button className={`button ${confirmAction.type === "delete" ? "button-danger" : "button-primary"}`} onClick={() => {
              const action = confirmAction;
              setConfirmAction(null);
              if (action.type === "link") void linkProfile(action.name);
              else if (action.type === "unlink") void unlinkProfile(action.name);
              else if (action.type === "overwrite") void runExtractProfile();
              else void removeProfile(action.name);
            }}>{confirmAction.type === "link" ? "Link to dataset" : confirmAction.type === "unlink" ? "Unlink from dataset" : confirmAction.type === "overwrite" ? "Extract & overwrite" : "Delete"}</button>
          </div>
        </div>
      </div>, document.body)}
      {renameAction && createPortal(<div className="cockpit-confirm-overlay" role="presentation" onMouseDown={() => setRenameAction(null)}>
        <div className="cockpit-confirm cockpit-rename-profile" role="dialog" aria-modal="true" aria-labelledby="af-rename-title" onMouseDown={(event) => event.stopPropagation()}>
          <h2 id="af-rename-title">Rename AF profile</h2>
          <label>
            Profile name
            <input
              autoFocus
              value={renameAction.value}
              onFocus={(event) => event.currentTarget.select()}
              onChange={(event) => setRenameAction({ ...renameAction, value: event.target.value, error: "" })}
              onKeyDown={(event) => {
                if (event.key === "Enter") {
                  event.preventDefault();
                  void renameProfile();
                }
              }}
            />
          </label>
          {renameAction.error && <p className="cockpit-confirm-error" role="alert">{renameAction.error}</p>}
          <div>
            <button className="button button-ghost" onClick={() => setRenameAction(null)}>Cancel</button>
            <button className="button button-primary" onClick={() => void renameProfile()}>Rename</button>
          </div>
        </div>
      </div>, document.body)}
    </>
  );
}
