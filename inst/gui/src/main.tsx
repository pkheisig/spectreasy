import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
import App from './App.tsx'
import PanelBuilder from './PanelBuilder.tsx'

const params = new URLSearchParams(window.location.search)
const mode = params.get('mode')
const RootComponent = mode === 'panel-builder' ? PanelBuilder : App

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <RootComponent />
  </StrictMode>,
)
