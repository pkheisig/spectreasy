import { useState } from "react";

type ProjectDirectoryFieldProps = {
  label: string;
  role: "controls" | "samples";
  value: string;
  onCommit: (role: "controls" | "samples", path: string) => Promise<boolean>;
};

export function ProjectDirectoryField({
  label,
  role,
  value,
  onCommit,
}: ProjectDirectoryFieldProps) {
  const [draft, setDraft] = useState(value);
  const [busy, setBusy] = useState(false);

  const normalizedDraft = draft.trim();
  const changed = normalizedDraft.length > 0 && normalizedDraft !== value;

  async function commitRename() {
    if (!changed || busy) return;
    const confirmed = window.confirm(
      `Rename the project folder “${value}” to “${normalizedDraft}”? Its contents will be moved.`,
    );
    if (!confirmed) return;

    setBusy(true);
    const success = await onCommit(role, normalizedDraft);
    setBusy(false);
    if (!success) setDraft(value);
  }

  return (
    <div className="project-directory-field">
      <label>
        {label}
        <input
          value={draft}
          onChange={(event) => setDraft(event.target.value)}
          disabled={busy}
        />
      </label>
      <button
        type="button"
        className="button button-ghost directory-rename-button"
        disabled={!changed || busy}
        onClick={() => void commitRename()}
      >
        {busy ? "Renaming…" : "Rename folder"}
      </button>
    </div>
  );
}
