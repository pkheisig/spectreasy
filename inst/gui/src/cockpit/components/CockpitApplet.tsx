import { Component, lazy, Suspense, useEffect } from 'react'
import type { ComponentType, ErrorInfo, ReactNode } from 'react'
import { createPortal } from 'react-dom'
import { ModuleLoadingState } from '../../ModuleLoadingState'
import { isStaleChunkLoadError } from '../../staleChunkRecovery'
import type { CockpitAppletId, ProjectState } from '../types'

const GatingGui = lazy(() => import('../../GatingGui.jsx')) as ComponentType<{
  embedded?: boolean
  cockpitTheme?: 'light' | 'dark'
  projectPath?: string
  projectRevision?: string
  initialFiles?: Array<Record<string, unknown>> | null
  initialMetadata?: Record<string, unknown>
  onRequestExit?: () => void
  onRequestClose?: () => void
}>
const PanelBuilder = lazy(() => import('../../PanelBuilder.tsx')) as ComponentType<{ embedded?: boolean; cockpitTheme?: 'light' | 'dark'; projectPath?: string; projectRevision?: string; onRequestExit?: () => void }>
const MatrixAdjustment = lazy(() => import('../../MatrixAdjustment.tsx')) as ComponentType<{ embedded?: boolean; cockpitTheme?: 'light' | 'dark'; projectPath?: string; projectRevision?: string; initialMatrixFiles?: string[]; initialSampleFiles?: string[]; initialUnmixingMethod?: string; onRequestExit?: () => void }>
const QcReportApplet = lazy(() => import('./QcReportApplet.tsx'))
const AnalysisWorkspace = lazy(() => import('../../analysis/AnalysisWorkspace.tsx'))

const appletLabels: Record<CockpitAppletId, string> = {
  'control-gating': 'control gating',
  'panel-builder': 'panel builder',
  'matrix-adjustment': 'matrix adjustment',
  'sample-analysis': 'sample analysis',
  'control-qc-report': 'controls QC report',
  'sample-qc-report': 'samples QC report',
}

const appletLoadingLabels: Record<CockpitAppletId, string> = {
  'control-gating': 'Loading Control Gating',
  'panel-builder': 'Loading Spectral Panel Builder',
  'matrix-adjustment': 'Loading Matrix Adjustment',
  'sample-analysis': 'Loading Sample Analysis',
  'control-qc-report': 'Loading Controls QC Report',
  'sample-qc-report': 'Loading Samples QC Report',
}

type AppletErrorBoundaryProps = {
  label: string
  theme: 'light' | 'dark'
  onClose: () => void
  children: ReactNode
}

class AppletErrorBoundary extends Component<AppletErrorBoundaryProps, { error: Error | null }> {
  state = { error: null as Error | null }

  static getDerivedStateFromError(error: Error) {
    return { error }
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    console.error(`Embedded ${this.props.label} failed to render.`, error, info)
  }

  render() {
    if (!this.state.error) return this.props.children
    const staleChunk = isStaleChunkLoadError(this.state.error)
    return (
      <div className={`cockpit-applet-error theme-${this.props.theme}`} role="alert">
        <div className="cockpit-applet-error-card">
          <strong>{this.props.label} could not open</strong>
          <span>{this.state.error.message || 'The module encountered an unexpected error.'}</span>
          <div className="cockpit-applet-error-actions">
            {staleChunk && <button type="button" onClick={() => window.location.reload()}>Reload app</button>}
            <button type="button" className={staleChunk ? 'is-secondary' : ''} onClick={this.props.onClose}>Return to cockpit</button>
          </div>
        </div>
      </div>
    )
  }
}

type CockpitAppletProps = {
  applet: CockpitAppletId
  theme: 'light' | 'dark'
  projectPath?: string
  project?: ProjectState
  outputRoot?: string
  reportPath?: string
  onExit: (reason?: 'exit' | 'confirmed') => void
  active?: boolean
}

export function CockpitApplet({ applet, theme, projectPath = '', project, outputRoot = 'spectreasy_outputs', reportPath = '', onExit, active = true }: CockpitAppletProps) {
  useEffect(() => {
    if (!active) return
    const root = document.getElementById('root')
    const rootWasInert = root?.hasAttribute('inert') ?? false
    const bodyOverflow = document.body.style.overflow
    const bodyClass = document.body.className
    const bodyStyle = document.body.getAttribute('style')
    const htmlStyle = document.documentElement.getAttribute('style')
    const htmlTheme = document.documentElement.dataset.theme

    root?.setAttribute('inert', '')
    document.body.style.overflow = 'hidden'

    return () => {
      if (!rootWasInert) root?.removeAttribute('inert')
      document.body.className = bodyClass
      if (bodyStyle === null) document.body.removeAttribute('style')
      else document.body.setAttribute('style', bodyStyle)
      document.body.style.overflow = bodyOverflow
      if (htmlStyle === null) document.documentElement.removeAttribute('style')
      else document.documentElement.setAttribute('style', htmlStyle)
      if (htmlTheme === undefined) delete document.documentElement.dataset.theme
      else document.documentElement.dataset.theme = htmlTheme
    }
  }, [active])

  const label = appletLabels[applet]
  const loadingLabel = appletLoadingLabels[applet]
  const projectRevision = project?.dataRevision || 'empty'
  const appletDataKey = applet === 'panel-builder' ? projectPath : `${projectPath}:${projectRevision}`
  return createPortal(
    <div className={`cockpit-applet theme-${theme} ${active ? '' : 'is-hidden'}`} role={active ? 'dialog' : undefined} aria-modal={active ? 'true' : undefined} aria-hidden={!active} aria-label={active ? `Embedded ${label}` : undefined}>
      <div className="cockpit-applet-content">
        <AppletErrorBoundary key={`${applet}:${appletDataKey}`} label={label} theme={theme} onClose={() => onExit('exit')}>
          <Suspense fallback={<ModuleLoadingState label={loadingLabel} theme={theme} />}>
            {applet === 'control-gating' && <GatingGui embedded cockpitTheme={theme} projectPath={projectPath} projectRevision={projectRevision} initialFiles={project?.gatingFiles ?? null} initialMetadata={project?.gatingMetadata ?? {}} onRequestExit={() => onExit('confirmed')} onRequestClose={() => onExit('exit')} />}
            {applet === 'panel-builder' && <PanelBuilder embedded cockpitTheme={theme} projectPath={projectPath} projectRevision={projectRevision} onRequestExit={() => onExit('exit')} />}
            {applet === 'matrix-adjustment' && <MatrixAdjustment embedded cockpitTheme={theme} projectPath={projectPath} projectRevision={projectRevision} initialMatrixFiles={project?.matrixFiles ?? undefined} initialSampleFiles={project?.sampleFiles ?? undefined} initialUnmixingMethod={project?.method} onRequestExit={() => onExit('exit')} />}
            {applet === 'sample-analysis' && <AnalysisWorkspace cockpitTheme={theme} projectPath={projectPath} onRequestExit={() => onExit('exit')} />}
            {applet === 'control-qc-report' && <QcReportApplet kind="control" theme={theme} projectPath={projectPath} outputRoot={outputRoot} initialReportPath={reportPath} onRequestExit={() => onExit('exit')} />}
            {applet === 'sample-qc-report' && <QcReportApplet kind="sample" theme={theme} projectPath={projectPath} outputRoot={outputRoot} initialReportPath={reportPath} onRequestExit={() => onExit('exit')} />}
          </Suspense>
        </AppletErrorBoundary>
      </div>
    </div>,
    document.body,
  )
}
