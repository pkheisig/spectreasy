import { lazy, Suspense, useEffect, useRef } from 'react'
import type { ComponentType } from 'react'
import { createPortal } from 'react-dom'
import { X } from 'lucide-react'
import type { CockpitAppletId } from '../types'

const GatingGui = lazy(() => import('../../GatingGui.jsx')) as ComponentType<{
  embedded?: boolean
  cockpitTheme?: 'light' | 'dark'
  onRequestExit?: () => void
}>
const PanelBuilder = lazy(() => import('../../PanelBuilder.tsx')) as ComponentType<{ embedded?: boolean; cockpitTheme?: 'light' | 'dark' }>
const MatrixAdjustment = lazy(() => import('../../MatrixAdjustment.tsx')) as ComponentType<{ embedded?: boolean; cockpitTheme?: 'light' | 'dark' }>
const QcReportApplet = lazy(() => import('./QcReportApplet.tsx'))
const AiReadyQcApplet = lazy(() => import('./AiReadyQcApplet.tsx'))

const appletLabels: Record<CockpitAppletId, string> = {
  'control-gating': 'control gating',
  'panel-builder': 'panel builder',
  'matrix-adjustment': 'matrix adjustment',
  'control-qc-report': 'controls QC report',
  'sample-qc-report': 'samples QC report',
  'ai-ready-qc': 'AI-ready QC',
}

type CockpitAppletProps = {
  applet: CockpitAppletId
  theme: 'light' | 'dark'
  projectPath?: string
  outputRoot?: string
  reportPath?: string
  onExit: (reason?: 'exit' | 'confirmed') => void
  active?: boolean
}

export function CockpitApplet({ applet, theme, projectPath = '', outputRoot = 'spectreasy_outputs', reportPath = '', onExit, active = true }: CockpitAppletProps) {
  const exitButtonRef = useRef<HTMLButtonElement>(null)

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
    exitButtonRef.current?.focus()

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
  return createPortal(
    <div className={`cockpit-applet theme-${theme} ${active ? '' : 'is-hidden'}`} role={active ? 'dialog' : undefined} aria-modal={active ? 'true' : undefined} aria-hidden={!active} aria-label={active ? `Embedded ${label}` : undefined}>
      <button
        ref={exitButtonRef}
        type="button"
        className="cockpit-applet-exit"
        onClick={() => onExit('exit')}
        aria-label={`Exit ${label} and return to cockpit`}
      >
        <X size={15} /> Exit
      </button>
      <div className="cockpit-applet-content">
        <Suspense fallback={<div className="cockpit-applet-loading">Loading {label}…</div>}>
          {applet === 'control-gating' && <GatingGui embedded cockpitTheme={theme} onRequestExit={() => onExit('confirmed')} />}
          {applet === 'panel-builder' && <PanelBuilder embedded cockpitTheme={theme} />}
          {applet === 'matrix-adjustment' && <MatrixAdjustment embedded cockpitTheme={theme} />}
          {applet === 'control-qc-report' && <QcReportApplet kind="control" theme={theme} projectPath={projectPath} outputRoot={outputRoot} initialReportPath={reportPath} />}
          {applet === 'sample-qc-report' && <QcReportApplet kind="sample" theme={theme} projectPath={projectPath} outputRoot={outputRoot} initialReportPath={reportPath} />}
          {applet === 'ai-ready-qc' && <AiReadyQcApplet theme={theme} projectPath={projectPath} outputRoot={outputRoot} />}
        </Suspense>
      </div>
    </div>,
    document.body,
  )
}
