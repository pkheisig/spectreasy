import { LoaderCircle } from 'lucide-react'
import './ModuleLoadingState.css'

type ModuleLoadingStateProps = {
  label: string
  theme?: 'light' | 'dark' | null
}

export function ModuleLoadingState({ label, theme = 'light' }: ModuleLoadingStateProps) {
  return (
    <div
      className={`module-loading-state theme-${theme === 'dark' ? 'dark' : 'light'}`}
      role="status"
      aria-live="polite"
      aria-label={label}
    >
      <div className="module-loading-card">
        <LoaderCircle className="module-loading-spinner" size={34} aria-hidden="true" />
        <strong>{label}</strong>
      </div>
    </div>
  )
}
