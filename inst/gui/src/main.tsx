import { StrictMode } from 'react'
import axios from 'axios'
import { createRoot } from 'react-dom/client'
import { resolveApiToken } from './apiBase.ts'
import './index.css'
import App from './App.tsx'

const params = new URLSearchParams(window.location.search)
const mode = params.get('mode')
const apiToken = resolveApiToken()
if (apiToken) axios.defaults.headers.common['X-Spectreasy-Token'] = apiToken
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
  if (mode === 'tuner') {
    const { default: MatrixAdjustment } = await import('./MatrixAdjustment.tsx')
    root.render(<StrictMode><MatrixAdjustment /></StrictMode>)
    return
  }
  root.render(<StrictMode><App /></StrictMode>)
}

void renderRoot()
