import type { ProjectState, SectionId } from "../types";

type RailProps = {
  activeSection: SectionId;
  project: ProjectState;
  onChange: (section: SectionId) => void;
};

const sections: Array<{ id: SectionId; number: string; title: string }> = [
  { id: "overview", number: "01", title: "Overview" },
  { id: "controls", number: "02", title: "Controls" },
  { id: "samples", number: "03", title: "Samples" },
  { id: "reports", number: "04", title: "Reports" },
  { id: "matrix", number: "05", title: "Matrix review" },
  { id: "panel", number: "06", title: "Panel builder" },
  { id: "af", number: "07", title: "AF library" },
  { id: "comparison", number: "08", title: "Compare methods" },
  { id: "simulator", number: "09", title: "Synthetic SCC" },
  { id: "settings", number: "10", title: "Settings & logs" },
];

function sectionCount(section: SectionId, project: ProjectState) {
  const counts: Partial<Record<SectionId, number>> = {
    controls: project.scan.controls,
    samples: project.scan.samples,
    reports: project.scan.reports,
    matrix: project.scan.matrices,
    af: project.artifacts.filter((artifact) => artifact.group === "AF Profiles")
      .length,
  };
  return counts[section];
}

export function WorkflowRail({ activeSection, project, onChange }: RailProps) {
  return (
    <nav className="workflow-rail" aria-label="Workflow sections">
      <div className="rail-intro">
        <span className="eyebrow">Workflow</span>
        <strong>{project.projectName}</strong>
      </div>
      <div className="rail-list">
        {sections.map(({ id, number, title }) => {
          const count = sectionCount(id, project);
          return (
            <button
              key={id}
              className={`rail-item ${activeSection === id ? "is-active" : ""}`}
              onClick={() => onChange(id)}
              aria-current={activeSection === id ? "page" : undefined}
            >
              <span className="rail-number">{number}</span>
              <span className="rail-copy">
                <strong>{title}</strong>
              </span>
              {count !== undefined && (
                <span className="rail-count">{count}</span>
              )}
            </button>
          );
        })}
      </div>
    </nav>
  );
}
