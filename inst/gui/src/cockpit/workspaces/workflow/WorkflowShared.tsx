import { RefreshCcw } from "lucide-react";
import type {
  Artifact,
  BackendStatus,
  CockpitAppletId,
  Job,
  MappingRow,
  PanelPayload,
  ProjectState,
  SectionId,
  WorkflowSettings,
} from "../../types";

export type WorkflowWorkspaceProps = {
  project: ProjectState;
  backend: BackendStatus;
  job: Job;
  activeSection: SectionId;
  mappingTab: "mapping" | "gating" | "build" | "qc";
  setMappingTab: (tab: "mapping" | "gating" | "build" | "qc") => void;
  sampleTab: "unmixing" | "qc";
  setSampleTab: (tab: "unmixing" | "qc") => void;
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
  onInputDirectoryChange: (role: "controls" | "samples", path: string) => Promise<boolean>;
  onOpenApplet: (applet: CockpitAppletId, reportPath?: string) => void;
  onSaveMapping?: () => void;
};

export function WorkspaceHeader({
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
