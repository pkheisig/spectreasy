import axios from 'axios'
import { demoProject } from './mockData'
import type { Artifact, BackendStatus, PanelPayload, ProjectState, Report, WorkflowSettings } from './types'

const API_BASE = (() => {
  const envBase = (import.meta.env.VITE_API_BASE as string | undefined)?.trim()
  if (envBase) return envBase.replace(/\/$/, '')
  if (typeof window !== 'undefined' && window.location.port === '5174') return 'http://localhost:8000'
  if (typeof window !== 'undefined') return window.location.origin.replace(/\/$/, '')
  return 'http://localhost:8000'
})()

const client = axios.create({ baseURL: API_BASE, timeout: 900 })

function scalarValue(value: unknown, fallback = ''): string {
  if (Array.isArray(value)) return scalarValue(value[0], fallback)
  if (value == null) return fallback
  return String(value)
}

function unboxGuiState(value: unknown): unknown {
  if (Array.isArray(value)) {
    if (value.length === 1) return unboxGuiState(value[0])
    return value.map((item) => unboxGuiState(item))
  }
  if (value && typeof value === 'object') {
    return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, unboxGuiState(item)]))
  }
  return value
}

function projectNameFromPath(path: string): string {
  const normalized = path.replace(/\\/g, '/').replace(/\/+$/, '')
  return normalized.split('/').pop()?.trim() || 'Spectreasy project'
}

export function rowsFromBackend(value: unknown): Array<Record<string, unknown>> {
  if (Array.isArray(value)) return value.filter((row): row is Record<string, unknown> => Boolean(row && typeof row === 'object'))
  if (!value || typeof value !== 'object') return []
  const columns = value as Record<string, unknown>
  const lengths = Object.values(columns).filter(Array.isArray).map((column) => column.length)
  const length = Math.max(0, ...lengths)
  return Array.from({ length }, (_, index) => Object.fromEntries(Object.entries(columns).map(([key, column]) => [key, Array.isArray(column) ? column[index] : column])))
}

function normalizeMappingRows(value: unknown): ProjectState['mapping'] {
  return rowsFromBackend(value).map((row, index) => {
    const controlType = String(row['control.type'] ?? row.controlType ?? 'cells').toLowerCase()
    const normalizedType: ProjectState['mapping'][number]['controlType'] = controlType.includes('viab') ? 'viability' : controlType.includes('bead') ? 'bead' : controlType.includes('unstain') ? 'unstained' : 'cell'
    const universal = String(row['universal.negative'] ?? row.universalNegative ?? '').toLowerCase()
    return {
      id: String(row.id ?? row.filename ?? `mapping-${index + 1}`),
      file: String(row.filename ?? row.file ?? ''),
      fluorophore: String(row.fluorophore ?? ''),
      marker: String(row.marker ?? ''),
      channel: String(row.channel ?? ''),
      controlType: normalizedType,
      universalNegative: ['true', 't', '1', 'yes'].includes(universal),
    }
  }).filter((row) => row.file.length > 0)
}

function liveArtifacts(projectPath: string, value: unknown): Artifact[] {
  const files = Array.isArray(value) ? value.map((file) => String(file).replace(/\\/g, '/')).filter(Boolean) : []
  const absolute = (relative: string) => `${projectPath.replace(/\\/g, '/').replace(/\/$/, '')}/${relative}`
  const artifact = (file: string, group: string, type: string, detail: string, status: Artifact['status'] = 'current'): Artifact => ({
    id: `live-${group}-${file}`,
    name: file.split('/').pop() || file,
    type,
    group,
    status,
    detail,
    path: absolute(file),
    updated: 'Present in active project',
  })
  const out: Artifact[] = []
  const controls = files.filter((file) => /(^|\/)(scc|controls?)\/.*\.fcs$/i.test(file))
  const samples = files.filter((file) => /(^|\/)samples?\/.*\.fcs$/i.test(file))
  if (controls.length) out.push(artifact('scc/', 'Controls', 'Folder', `${controls.length} FCS file${controls.length === 1 ? '' : 's'}`))
  files.filter((file) => /^fcs_mapping\.csv$/i.test(file)).forEach((file) => out.push(artifact(file, 'Controls', 'Mapping', 'Control mapping CSV')))
  files.filter((file) => /gate.*\.csv$/i.test(file)).forEach((file) => out.push(artifact(file, 'Gates', 'Gates', 'Saved gate configuration')))
  if (samples.length) out.push(artifact('samples/', 'Samples', 'Folder', `${samples.length} FCS file${samples.length === 1 ? '' : 's'}`))
  const matrices = files.filter((file) => /(matrix|unmixing|detector_noise|variant).*\.(csv|rds)$/i.test(file)).slice(0, 12)
  matrices.forEach((file) => out.push(artifact(file, 'Matrices', /variant/i.test(file) ? 'Spectral library' : 'Matrix', 'Project matrix artifact')))
  const reports = files.filter((file) => /(^|\/)(reports|spectreasy_outputs)\//i.test(file) && /\.(html?|pdf)$/i.test(file)).slice(-8)
  reports.forEach((file) => out.push(artifact(file, 'Reports', /sample|qc_samples/i.test(file) ? 'Sample QC report' : 'QC report', 'Project report')))
  files.filter((file) => /metric.*\.csv$/i.test(file)).slice(0, 6).forEach((file) => out.push(artifact(file, 'QC Metrics', 'Metrics', 'QC metric table', /sample|013/i.test(file) ? 'stale' : 'current')))
  return out
}

export const initialBackendStatus: BackendStatus = {
  connected: false,
  version: '—',
  mode: 'Local preview',
  message: 'Using a local preview state until the R backend connects.',
  apiPort: '8000',
  packageReady: false,
}

export async function loadProjectSnapshot(): Promise<{ project: ProjectState; backend: BackendStatus; savedSettings?: Partial<WorkflowSettings> }> {
    try {
    const [statusResponse, matricesResponse, samplesResponse, gatesResponse, guiStateResponse] = await Promise.all([
      client.get('/status'),
      client.get('/matrices'),
      client.get('/samples'),
      client.get('/gate_files'),
      client.get('/gui_state', { params: { module: 'spectreasy_cockpit' } }).catch(() => null),
    ])
    const projectResponse = await client.get('/project/status').catch(() => null)
    const status = statusResponse.data as Record<string, unknown>
    const matrices = Array.isArray(matricesResponse.data) ? matricesResponse.data : []
    const samples = Array.isArray(samplesResponse.data) ? samplesResponse.data : []
    const gateRows = rowsFromBackend(gatesResponse.data?.files)
    const liveMapping = normalizeMappingRows(gatesResponse.data?.files)
    const statusOk = scalarValue(status.status) === 'ok'
    const rawScan = (projectResponse?.data?.scan ?? {}) as Record<string, unknown>
    const backendScan = Object.fromEntries(Object.entries(rawScan).map(([key, value]) => [key, Number(scalarValue(value, '0'))])) as Partial<ProjectState['scan']>
    const backend: BackendStatus = {
      connected: statusOk,
      version: scalarValue(status.version, 'Spectreasy R package'),
      mode: scalarValue(status.gui_mode, 'Local R backend'),
      message: 'R backend connected. Project artifacts refreshed from disk.',
      apiPort: API_BASE.split(':').pop() ?? '8000',
      packageReady: statusOk,
    }
    const projectPath = scalarValue(projectResponse?.data?.project_path, scalarValue(status.matrix_dir, scalarValue(status.wd, demoProject.projectPath)))
    const liveProjectArtifacts = liveArtifacts(projectPath, projectResponse?.data?.files)
    const project: ProjectState = {
      ...demoProject,
      projectName: scalarValue(status.project_name, projectNameFromPath(projectPath)),
      projectPath,
      method: scalarValue(status.unmixing_method, demoProject.method),
      cytometer: scalarValue(status.panel_cytometer, demoProject.cytometer),
      artifacts: liveProjectArtifacts.length > 0 ? liveProjectArtifacts : demoProject.artifacts,
      mapping: liveMapping.length > 0 ? liveMapping : demoProject.mapping,
      scan: { ...demoProject.scan, ...backendScan, matrices: matrices.length, samples: samples.length, gates: gateRows.length },
    }
    const savedConfig = unboxGuiState(guiStateResponse?.data?.config) as Record<string, unknown> | undefined
    const savedSettings = (savedConfig?.settings ?? savedConfig ?? undefined) as Partial<WorkflowSettings> | undefined
    return { project, backend, savedSettings }
  } catch {
    return { project: structuredClone(demoProject), backend: initialBackendStatus }
  }
}

export async function loadPanelPayload(cytometer: string): Promise<PanelPayload | null> {
  return loadPanelMetrics(cytometer, '', [])
}

export async function loadPanelMetrics(cytometer: string, configuration: string, fluorophores: string[]): Promise<PanelPayload | null> {
  try {
    const response = await client.post('/spectral_panel_metrics', { cytometer, configuration: configuration || null, fluorophores })
    if (response.data?.error) return null
    return response.data as PanelPayload
  } catch {
    return null
  }
}

export async function exportPanelOverview(cytometer: string, configuration: string, fluorophores: string[], markers: string[] = []): Promise<{ filename: string; contentType: string; contentBase64: string } | null> {
  try {
    const response = await client.post('/export_spectral_panel_overview', { cytometer, configuration: configuration || null, fluorophores, markers }, { timeout: 120000 })
    if (response.data?.error || !response.data?.content_base64) return null
    return {
      filename: scalarValue(response.data.filename, 'spectreasy_panel_overview.pdf'),
      contentType: scalarValue(response.data.content_type, 'application/pdf'),
      contentBase64: scalarValue(response.data.content_base64),
    }
  } catch {
    return null
  }
}

export function downloadBase64File(filename: string, contentBase64: string, contentType = 'application/octet-stream') {
  const binary = window.atob(contentBase64)
  const bytes = Uint8Array.from(binary, (character) => character.charCodeAt(0))
  const blob = new Blob([bytes], { type: contentType })
  const href = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = href
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(href)
}

export async function listMatrixFiles(): Promise<string[]> {
  try {
    const response = await client.get('/matrices')
    return Array.isArray(response.data) ? response.data.map((filename) => String(filename)) : []
  } catch {
    return []
  }
}

export async function listSampleFiles(): Promise<string[]> {
  try {
    const response = await client.get('/samples')
    return Array.isArray(response.data) ? response.data.map((filename: unknown) => String(filename)) : []
  } catch {
    return []
  }
}

export async function importSampleContent(filename: string, contentBase64: string): Promise<boolean> {
  try {
    const response = await client.post('/import_sample_content', { filename, content_base64: contentBase64 })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function loadMatrixFile(filename: string): Promise<{ filename: string; rows: Array<Record<string, unknown>> } | null> {
  try {
    const response = await client.get('/load_matrix', { params: { filename } })
    if (response.data?.error) return null
    return { filename, rows: rowsFromBackend(response.data) }
  } catch {
    return null
  }
}

export async function saveMatrixFile(filename: string, rows: Array<Record<string, unknown>>, sourceFilename = ''): Promise<boolean> {
  try {
    const response = await client.post('/save_matrix', { filename, source_filename: sourceFilename, matrix_json: rows })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function importMatrixContent(filename: string, content: string): Promise<boolean> {
  try {
    const response = await client.post('/import_matrix_content', { filename, content })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function listAfProfiles(): Promise<Array<{ name: string; bands: number; detectors: number; created: string; path: string }>> {
  try {
    const response = await client.get('/af_profiles')
    return rowsFromBackend(response.data?.profiles).map((row) => ({
      name: String(row.name ?? ''),
      bands: Number(row.bands ?? 0),
      detectors: Number(row.detectors ?? 0),
      created: String(row.created ?? ''),
      path: String(row.path ?? ''),
    })).filter((profile) => profile.name.length > 0)
  } catch {
    return []
  }
}

export async function loadProjectReports(): Promise<Report[]> {
  try {
    const response = await client.get('/project/reports')
    return rowsFromBackend(response.data?.reports).map((row) => {
      const normalized = String(row.path ?? '').replace(/\\/g, '/')
      const format = String(row.format ?? (normalized.toLowerCase().endsWith('.pdf') ? 'PDF' : 'HTML')) as 'HTML' | 'PDF'
      const basename = normalized.split('/').pop() ?? normalized
      const title = basename.replace(/\.[^.]+$/, '').replace(/[_-]+/g, ' ')
      return {
        id: `live-${normalized}`,
        title,
        type: String(row.report_type ?? (/sample|qc_samples/i.test(normalized) ? 'Sample QC' : 'Control QC')) as Report['type'],
        format,
        run: normalized.split('/').find((part: string) => /^run[-_]/i.test(part)) ?? 'project',
        created: String(row.created ?? 'Present in active project'),
        status: String(row.status ?? 'current') as Report['status'],
        matrix: '—',
        path: normalized,
      } satisfies Report
    }).filter((report) => report.path.length > 0)
  } catch {
    return []
  }
}

export function projectFileUrl(path: string): string {
  return `${API_BASE}/project/file?path=${encodeURIComponent(path)}`
}

export async function deleteAfProfile(name: string): Promise<boolean> {
  try {
    const response = await client.delete('/af_profiles', { params: { name } })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function applyAfProfile(matrixFilename: string, profileName: string, outputFilename = ''): Promise<boolean> {
  try {
    const response = await client.post('/af_profiles/apply', { matrix_filename: matrixFilename, profile_name: profileName, output_filename: outputFilename || matrixFilename })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function persistGuiState(project: ProjectState, settings?: WorkflowSettings): Promise<boolean> {
  try {
    await client.post('/gui_state', {
      module: 'spectreasy_cockpit',
      config_json: {
        projectPath: project.projectPath,
        projectName: project.projectName,
        cytometer: project.cytometer,
        method: project.method,
        settings,
      },
    })
    return true
  } catch {
    return false
  }
}

export async function persistControlMapping(mapping: ProjectState['mapping']): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/control_mapping', { rows: mapping })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    const error = scalarValue(response.data?.error, '')
    return success ? { success: true, message: 'Control mapping saved to the project CSV.' } : { success: false, message: error || 'The R backend rejected the control mapping.' }
  } catch {
    return { success: false, message: 'Control mapping could not be saved. Connect the R backend and retry.' }
  }
}

export async function setProjectContext(projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/project/context', { projectPath })
    if (response.data?.success === false || response.data?.error) {
      return { success: false, message: response.data?.error ?? 'The R backend rejected this project folder.' }
    }
    return { success: true, message: 'Project folder opened. Artifacts refreshed from disk.' }
  } catch {
    return { success: false, message: 'The project folder could not be opened. Check the path and R connection.' }
  }
}

export async function attemptWorkflowAction(action: string, payload: Record<string, unknown>): Promise<{ connected: boolean; message: string }> {
  const endpoint = action === 'control' ? '/workflow/control' : action === 'sample' ? '/workflow/sample' : action === 'af' ? '/workflow/af' : '/workflow/report'
  try {
    const response = await client.post(endpoint, payload, { timeout: 900000 })
    const success = scalarValue(response.data?.success, 'true')
    const error = scalarValue(response.data?.error, '')
    if (success === 'false' || error) {
      return { connected: false, message: error || 'The R backend rejected this workflow action.' }
    }
    return { connected: true, message: 'Job accepted by the R backend.' }
  } catch {
    return { connected: false, message: 'Preview job queued locally. Connect the R backend to run the full Spectreasy workflow.' }
  }
}

export async function attemptDiagnosticAction(action: 'compare' | 'synthetic', payload: Record<string, unknown>): Promise<{ connected: boolean; message: string; result?: unknown }> {
  const endpoint = action === 'compare' ? '/workflow/compare' : '/workflow/synthetic'
  try {
    const response = await client.post(endpoint, payload, { timeout: 900000 })
    const success = scalarValue(response.data?.success, 'true')
    const error = scalarValue(response.data?.error, '')
    if (success === 'false' || error) return { connected: false, message: error || 'The R backend rejected this diagnostic action.' }
    return { connected: true, message: action === 'compare' ? 'Method comparison completed in the R backend.' : 'Synthetic SCC FCS written by the R backend.', result: response.data?.result }
  } catch {
    return { connected: false, message: 'The diagnostic action could not reach the R backend.' }
  }
}

export function downloadTextFile(filename: string, content: string) {
  const blob = new Blob([content], { type: 'text/plain;charset=utf-8' })
  const href = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = href
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(href)
}
