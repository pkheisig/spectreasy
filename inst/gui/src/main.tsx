import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'

const params = new URLSearchParams(window.location.search)
const mode = params.get('mode')
const root = createRoot(document.getElementById('root')!)

async function renderRoot() {
  if (mode === 'panel-builder') {
    const { default: PanelBuilder } = await import('./PanelBuilder.tsx')
    root.render(<StrictMode><PanelBuilder /></StrictMode>)
    return
  }
  if (mode === 'control-gating') {
    const { default: GatingGui } = await import('./GatingGui.jsx')
    root.render(<StrictMode><GatingGui /></StrictMode>)
    return
  }
  root.render(<StrictMode><App /></StrictMode>)
}

void renderRoot()
