import { JobStrip } from "../components/JobStrip";
import { defaultWorkflowSettings } from "../types";
import type { SectionId } from "../types";
import { ConfigurableAfWorkspace } from "./workflow/AfWorkspace";
import { ControlsWorkspace, SamplesWorkspace } from "./workflow/SamplesWorkspace";
import { ConfigurableSettingsWorkspace } from "./workflow/SettingsWorkspace";
import { ControlReferenceTuning, SampleOutputTuning } from "./workflow/TuningPanels";
import type { WorkflowWorkspaceProps } from "./workflow/WorkflowShared";

export type { WorkflowWorkspaceProps } from "./workflow/WorkflowShared";

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
      {activeSection === "samples" && <SamplesWorkspace {...props} />}
      {activeSection === "af" && (
        <ConfigurableAfWorkspace
          projectPath={props.project.projectPath}
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
            onReset={() => props.onSettingsChange("control", { ...defaultWorkflowSettings(props.settings.projectPath).control, sccDir: props.settings.control.sccDir })}
          />
          <SampleOutputTuning
            settings={props.settings.sample}
            onSettingsChange={(patch) =>
              props.onSettingsChange("sample", patch)
            }
            onReset={() => props.onSettingsChange("sample", { ...defaultWorkflowSettings(props.settings.projectPath).sample, sampleDir: props.settings.sample.sampleDir })}
          />
        </>
      )}
    </div>
  );
}
