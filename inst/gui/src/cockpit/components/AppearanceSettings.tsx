import {
  Contrast,
  LayoutPanelLeft,
  Moon,
  Sun,
  Type,
} from 'lucide-react'
import type { AppearanceSettings as AppearanceSettingsValue } from '../types'

type Props = {
  value: AppearanceSettingsValue
  onChange: (patch: Partial<AppearanceSettingsValue>) => void
}

export function AppearanceSettings({ value, onChange }: Props) {
  return (
    <details className="surface-card settings-section appearance-settings" open>
      <summary>
        <Type size={16} /> Appearance
        <span>Shared cockpit display preferences</span>
      </summary>
      <div className="appearance-settings-body">
        <section className="appearance-group">
          <div className="appearance-group-heading">
            {value.theme === 'dark' ? <Moon size={16} /> : <Sun size={16} />}
            <div>
              <strong>Theme</strong>
              <small>Uses the same light and dark modes as the analysis GUIs.</small>
            </div>
          </div>
          <div className="choice-row" role="group" aria-label="Color theme">
            {(['light', 'dark'] as const).map((theme) => (
              <button
                key={theme}
                className={`choice-button ${value.theme === theme ? 'is-selected' : ''}`}
                onClick={() => onChange({ theme })}
                type="button"
              >
                {theme === 'light' ? <Sun size={15} /> : <Moon size={15} />}
                {theme === 'light' ? 'Light' : 'Dark'}
              </button>
            ))}
          </div>
        </section>

        <section className="appearance-group">
          <div className="appearance-group-heading">
            <Type size={16} />
            <div>
              <strong>Scale and spacing</strong>
              <small>Avenir Next and the Spectreasy accent colors stay fixed.</small>
            </div>
          </div>
          <div className="appearance-field-grid">
            <label>
              Interface size <strong>{value.fontScale}%</strong>
              <input
                aria-label="Interface size"
                type="range"
                min="90"
                max="140"
                step="5"
                value={value.fontScale}
                onChange={(event) => onChange({ fontScale: Number(event.target.value) })}
              />
            </label>
            <label>
              Density
              <select
                aria-label="Interface density"
                value={value.density}
                onChange={(event) => onChange({ density: event.target.value as AppearanceSettingsValue['density'] })}
              >
                <option value="compact">Compact</option>
                <option value="comfortable">Comfortable</option>
                <option value="spacious">Spacious</option>
              </select>
            </label>
            <label>
              Corner radius <strong>{value.cornerRadius}px</strong>
              <input
                aria-label="Corner radius"
                type="range"
                min="0"
                max="16"
                step="2"
                value={value.cornerRadius}
                onChange={(event) => onChange({ cornerRadius: Number(event.target.value) })}
              />
            </label>
            <label>
              Card depth
              <select
                aria-label="Card depth"
                value={value.shadows}
                onChange={(event) => onChange({ shadows: event.target.value as AppearanceSettingsValue['shadows'] })}
              >
                <option value="none">Flat</option>
                <option value="subtle">Subtle</option>
                <option value="raised">Raised</option>
              </select>
            </label>
          </div>
        </section>

        <section className="appearance-group">
          <div className="appearance-group-heading">
            <LayoutPanelLeft size={16} />
            <div>
              <strong>Navigation and behavior</strong>
              <small>Control how much supporting information stays visible.</small>
            </div>
          </div>
          <div className="appearance-field-grid">
            <label>
              Sidebar width
              <select
                aria-label="Sidebar width"
                value={value.sidebarWidth}
                onChange={(event) => onChange({ sidebarWidth: event.target.value as AppearanceSettingsValue['sidebarWidth'] })}
              >
                <option value="narrow">Narrow</option>
                <option value="standard">Standard</option>
                <option value="wide">Wide</option>
              </select>
            </label>
          </div>
          <div className="preference-toggles">
            <label><input type="checkbox" checked={value.showSectionCounts} onChange={(event) => onChange({ showSectionCounts: event.target.checked })} /><span>Show file counts</span></label>
            <label><input type="checkbox" checked={value.stickyHeader} onChange={(event) => onChange({ stickyHeader: event.target.checked })} /><span>Keep header visible</span></label>
            <label><input type="checkbox" checked={value.backgroundTexture} onChange={(event) => onChange({ backgroundTexture: event.target.checked })} /><span>Subtle background texture</span></label>
            <label><input type="checkbox" checked={value.reduceMotion} onChange={(event) => onChange({ reduceMotion: event.target.checked })} /><span>Reduce motion</span></label>
            <label><input type="checkbox" checked={value.highContrast} onChange={(event) => onChange({ highContrast: event.target.checked })} /><Contrast size={14} /><span>Higher contrast</span></label>
          </div>
        </section>
      </div>
    </details>
  )
}
