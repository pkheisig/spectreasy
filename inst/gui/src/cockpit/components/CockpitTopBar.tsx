type TopBarProps = {
  project: ProjectState;
  cytometer: string;
  method: string;
  darkMode: boolean;
  settingsActive: boolean;
  terminalActive: boolean;
  backendConnected: boolean;
  onCytometerChange: (value: string) => void;
  onMethodChange: (value: string) => void;
  onToggleTheme: () => void;
  onFiles: () => void;
  onCreateProject: () => void;
  onOpenProject: () => void;
  onGuide: () => void;
  onSettings: () => void;
  onTerminal: () => void;
  executionLogs: ExecutionLogEntry[];
};

export function TopBar({
  project,
  cytometer,
  method,
  darkMode,
  settingsActive,
  terminalActive,
  backendConnected,
  onCytometerChange,
  onMethodChange,
  onToggleTheme,
  onFiles,
  onCreateProject,
  onOpenProject,
  onGuide,
  onSettings,
  onTerminal,
  executionLogs,
}: TopBarProps) {
  const [projectMenuOpen, setProjectMenuOpen] = useState(false);
  const projectMenuRef = useRef<HTMLDivElement>(null);
  const logsMenuRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!projectMenuOpen) return;
    const closeMenu = (event: MouseEvent | KeyboardEvent) => {
      if (event instanceof KeyboardEvent && event.key !== "Escape") return;
      if (event instanceof MouseEvent && projectMenuRef.current?.contains(event.target as Node)) return;
      setProjectMenuOpen(false);
    };
    document.addEventListener("mousedown", closeMenu);
    document.addEventListener("keydown", closeMenu);
    return () => {
      document.removeEventListener("mousedown", closeMenu);
      document.removeEventListener("keydown", closeMenu);
    };
  }, [projectMenuOpen]);

  useEffect(() => {
    if (!terminalActive) return;
    const closeLogs = (event: MouseEvent | KeyboardEvent) => {
      if (event instanceof KeyboardEvent && event.key !== "Escape") return;
      if (event instanceof MouseEvent && logsMenuRef.current?.contains(event.target as Node)) return;
      onTerminal();
    };
    document.addEventListener("mousedown", closeLogs);
    document.addEventListener("keydown", closeLogs);
    return () => {
      document.removeEventListener("mousedown", closeLogs);
      document.removeEventListener("keydown", closeLogs);
    };
  }, [onTerminal, terminalActive]);

  const chooseProjectAction = (action: () => void) => {
    setProjectMenuOpen(false);
    action();
  };

  return (
    <header className="topbar">
      <div className="brand-lockup">
        <div className="brand-mark">
          <span />
          <span />
          <span />
        </div>
        <div>
          <strong>spectreasy</strong>
        </div>
      </div>
      <div className="project-actions">
        <div className="project-menu-shell" ref={projectMenuRef}>
          <button
            className="project-switcher"
            type="button"
            aria-expanded={projectMenuOpen}
            aria-haspopup="menu"
            onClick={() => setProjectMenuOpen((open) => !open)}
          >
            <FolderOpen size={16} />
            <strong>Project</strong>
          </button>
          <div
            className={`project-action-flyout ${projectMenuOpen ? "is-open" : ""}`}
            role="menu"
            aria-hidden={!projectMenuOpen}
          >
            <button type="button" role="menuitem" tabIndex={projectMenuOpen ? 0 : -1} onClick={() => chooseProjectAction(onCreateProject)}>
              Create new project
            </button>
            <button type="button" role="menuitem" tabIndex={projectMenuOpen ? 0 : -1} onClick={() => chooseProjectAction(onOpenProject)}>
              Open project
            </button>
          </div>
        </div>
        <button className="project-files-button" type="button" onClick={onFiles} disabled={!project.projectPath}>
          <FilesIcon size={15} />
          <span>Files</span>
        </button>
      </div>
      <button className="guide-trigger" type="button" onClick={onGuide}>
        <CircleHelp size={16} />
        <span>Guide</span>
      </button>
      <div className="topbar-spacer" />
      <label className="context-select cytometer-select">
        <span className="chip-label">Cytometer</span>
        <GuiSelect
          aria-label="Cytometer"
          className="cytometer-picker"
          value={cytometer}
          onChange={(event) => onCytometerChange(event.target.value)}
        >
          {cytometerOptions.map(([value, label]) => (
            <option value={value} key={value}>{label}</option>
          ))}
        </GuiSelect>
      </label>
      <label className="context-select method-select">
        <span className="chip-label">Method</span>
        <GuiSelect
          aria-label="Method"
          value={method}
          onChange={(event) => onMethodChange(event.target.value)}
        >
          <option>Spectreasy</option>
          <option>AutoSpectral</option>
          <option>OLS</option>
          <option>WLS</option>
          <option>RWLS</option>
          <option>NNLS</option>
        </GuiSelect>
      </label>
      <button
        className="topbar-icon"
        onClick={onToggleTheme}
        aria-label={darkMode ? "Use light mode" : "Use dark mode"}
        title={darkMode ? "Use light mode" : "Use dark mode"}
      >
        {darkMode ? <Sun size={17} /> : <Moon size={17} />}
      </button>
      <div className="logs-menu-shell" ref={logsMenuRef}>
        <button
          className={`topbar-icon terminal-trigger ${terminalActive ? "is-active" : ""}`}
          onClick={onTerminal}
          aria-label={terminalActive ? "Close execution logs" : "Open execution logs"}
          aria-expanded={terminalActive}
          title={backendConnected ? "Execution logs" : "Execution logs · R backend offline"}
        >
          <TerminalSquare size={18} />
          <span className={`terminal-status-dot ${backendConnected ? "is-connected" : ""}`} />
        </button>
        {terminalActive && <TerminalPanel connected={backendConnected} entries={executionLogs} />}
      </div>
      <button
        className={`topbar-icon ${settingsActive ? "is-active" : ""}`}
        onClick={onSettings}
        aria-label={settingsActive ? "Return to workflow" : "Settings"}
        aria-pressed={settingsActive}
        title={settingsActive ? "Return to workflow" : "Settings"}
      >
        <Settings size={18} />
      </button>
    </header>
  );
}
import { useEffect, useRef, useState } from "react";
import {
  CircleHelp,
  Files as FilesIcon,
  FolderOpen,
  Moon,
  Settings,
  Sun,
  TerminalSquare,
} from "lucide-react";
import { cytometerOptions } from "../cytometerConfig";
import type { ExecutionLogEntry, ProjectState } from "../types";
import { GuiSelect } from "./GuiSelect";
import { TerminalPanel } from "./TerminalPanel";
