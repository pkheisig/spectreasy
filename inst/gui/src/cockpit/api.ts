import axios from 'axios'
import { emptyProject } from './projectState'
import { resolveApiBase, resolveApiToken } from '../apiBase'
import type { Artifact, BackendStatus, PanelPayload, ProjectState, Report, WorkflowSettings } from './types'

const API_BASE = resolveApiBase()

const client = axios.create({
  baseURL: API_BASE,
  timeout: 30000,
  headers: { 'X-Spectreasy-Token': resolveApiToken() },
})

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
  if (!normalized) return 'Choose a project folder'
  return normalized.split('/').pop()?.trim() || 'Spectreasy project'
}

function normalizeOutputRoot(value: unknown): string {
  const output = String(value ?? 'spectreasy_outputs').trim().replace(/[\\/]+$/, '')
  return output.replace(/[\\/](?:unmix_controls|unmix_samples)$/i, '') || 'spectreasy_outputs'
}

export function rowsFromBackend(value: unknown): Array<Record<string, unknown>> {
  if (Array.isArray(value)) return value.filter((row): row is Record<string, unknown> => Boolean(row && typeof row === 'object'))
  if (!value || typeof value !== 'object') return []
  const columns = value as Record<string, unknown>
  const lengths = Object.values(columns).filter(Array.isArray).map((column) => column.length)
  const length = Math.max(0, ...lengths)
  return Array.from({ length }, (_, index) => Object.fromEntries(Object.entries(columns).map(([key, column]) => [key, Array.isArray(column) ? column[index] : column])))
}

function stringsFromBackend(value: unknown): string[] {
  if (Array.isArray(value)) return value.map((item) => scalarValue(item)).filter(Boolean)
  const scalar = scalarValue(value)
  return scalar ? [scalar] : []
}

function normalizeMappingRows(value: unknown): ProjectState['mapping'] {
  return rowsFromBackend(value).map((row, index) => {
    const controlType = String(row['control.type'] ?? row.controlType ?? 'cells').toLowerCase()
    const normalizedType: ProjectState['mapping'][number]['controlType'] = controlType.includes('bead') ? 'bead' : 'cell'
    return {
      id: String(row.id ?? row.filename ?? `mapping-${index + 1}`),
      file: String(row.filename ?? row.file ?? ''),
      fluorophore: String(row.fluorophore ?? ''),
      marker: String(row.marker ?? ''),
      channel: String(row.channel ?? ''),
      controlType: normalizedType,
      universalNegative: String(row['universal.negative'] ?? row.universalNegative ?? ''),
      isViability: ['true', 't', '1', 'yes'].includes(String(row['is.viability'] ?? row.isViability ?? '').toLowerCase()),
      ignored: ['true', 't', '1', 'yes'].includes(String(row.ignored ?? '').toLowerCase()),
      ignoredReason: String(row.ignored_reason ?? row.ignoredReason ?? ''),
    }
  }).filter((row) => row.file.length > 0)
}

function liveArtifacts(projectPath: string, value: unknown, reportValue: unknown, controlInputDir: string, sampleInputDir: string): Artifact[] {
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
    relativePath: file,
    updated: 'Present in active project',
  })
  const out: Artifact[] = []
  const directoryFiles = (directory: string) => {
    const prefix = `${directory.replace(/\\/g, '/').replace(/^\/+|\/+$/g, '')}/`.toLowerCase()
    return files.filter((file) => file.toLowerCase().startsWith(prefix) && /\.fcs$/i.test(file))
  }
  const controls = directoryFiles(controlInputDir)
  const samples = directoryFiles(sampleInputDir)
  if (controls.length) out.push(artifact(`${controlInputDir}/`, 'Controls', 'Folder', `${controls.length} FCS file${controls.length === 1 ? '' : 's'}`))
  files.filter((file) => /^fcs_mapping\.csv$/i.test(file)).forEach((file) => out.push(artifact(file, 'Controls', 'Mapping', 'Control mapping CSV')))
  files.filter((file) => /gate.*\.csv$/i.test(file)).forEach((file) => out.push(artifact(file, 'Gates', 'Gates', 'Saved gate configuration')))
  if (samples.length) out.push(artifact(`${sampleInputDir}/`, 'Samples', 'Folder', `${samples.length} FCS file${samples.length === 1 ? '' : 's'}`))
  const matrices = files.filter((file) => /(matrix|unmixing|detector_noise|variant).*\.(csv|rds)$/i.test(file)).slice(0, 12)
  matrices.forEach((file) => out.push(artifact(file, 'Matrices', /variant/i.test(file) ? 'Spectral library' : 'Matrix', 'Project matrix artifact')))
  rowsFromBackend(reportValue)
    .map((row) => ({
        path: String(row.path ?? '').replace(/\\/g, '/'),
        type: String(row.report_type ?? ''),
        status: String(row.status ?? 'current'),
        created: String(row.created ?? 'Present in active project'),
        createdEpoch: Number(row.created_epoch),
    }))
    .filter((report) => report.path.length > 0)
    .sort((left, right) => (Number.isFinite(left.createdEpoch) ? left.createdEpoch : 0) - (Number.isFinite(right.createdEpoch) ? right.createdEpoch : 0))
    .forEach((report) => {
        const reportArtifact = artifact(
          report.path,
          'Reports',
          report.type === 'Sample QC' ? 'Sample QC report' : 'Control QC report',
          'Project report',
          report.status === 'stale' ? 'stale' : 'current',
        )
        reportArtifact.updated = report.created
        reportArtifact.updatedEpoch = Number.isFinite(report.createdEpoch) ? report.createdEpoch : undefined
      out.push(reportArtifact)
    })
  files.filter((file) => /metric.*\.csv$/i.test(file)).slice(0, 6).forEach((file) => out.push(artifact(file, 'QC Metrics', 'Metrics', 'QC metric table')))
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

export async function loadProjectSnapshot(requestedProjectPath = ''): Promise<{ project: ProjectState; backend: BackendStatus; savedSettings?: Partial<WorkflowSettings> }> {
  try {
    // Establish connectivity from the lightweight health check first. Artifact
    // endpoints are independent: one unavailable artifact must not make a
    // running R session appear offline.
    const statusResponse = await client.get('/status', { timeout: 5000 })
    const [projectResponse, matricesResponse, samplesResponse, mappingResponse, guiStateResponse] = await Promise.all([
      client.get('/project/status', { params: { project_path: requestedProjectPath }, timeout: 10000 }).catch(() => null),
      client.get('/matrices', { params: { project_path: requestedProjectPath }, timeout: 10000 }).catch(() => null),
      client.get('/samples', { params: { project_path: requestedProjectPath }, timeout: 10000 }).catch(() => null),
      client.get('/control_mapping', { params: { project_path: requestedProjectPath }, timeout: 10000 }).catch(() => null),
      client.get('/gui_state', { params: { module: 'spectreasy_cockpit', project_path: requestedProjectPath }, timeout: 5000 }).catch(() => null),
    ])
    const status = statusResponse.data as Record<string, unknown>
    const matrices = Array.isArray(matricesResponse?.data) ? matricesResponse.data : []
    const samples = Array.isArray(samplesResponse?.data) ? samplesResponse.data : []
    const liveMapping = normalizeMappingRows(mappingResponse?.data?.rows)
    const statusOk = scalarValue(status.status) === 'ok'
    const rawScan = (projectResponse?.data?.scan ?? {}) as Record<string, unknown>
    const rawLayout = (projectResponse?.data?.layout ?? {}) as Record<string, unknown>
    const controlInputDir = scalarValue(rawLayout.control_input_dir, 'scc')
    const sampleInputDir = scalarValue(rawLayout.sample_input_dir, 'samples')
    const scanCount = (snakeCase: string, camelCase = snakeCase) =>
      Number(scalarValue(rawScan[snakeCase] ?? rawScan[camelCase], '0')) || 0
    const backendScan: ProjectState['scan'] = {
      controls: scanCount('controls'),
      samples: scanCount('samples'),
      matrices: scanCount('matrices'),
      reports: scanCount('reports'),
      gates: scanCount('gates'),
      qcMetrics: scanCount('qc_metrics', 'qcMetrics'),
      spectralVariants: scanCount('spectral_variants', 'spectralVariants'),
    }
    const backend: BackendStatus = {
      connected: statusOk,
      version: scalarValue(status.version, 'Spectreasy R package'),
      mode: scalarValue(status.gui_mode, 'Local R backend'),
      message: 'R backend connected. Project artifacts refreshed from disk.',
      apiPort: API_BASE.split(':').pop() ?? '8000',
      packageReady: statusOk,
    }
    const projectPath: string = projectResponse
      ? scalarValue(projectResponse.data?.project_path, '')
      : scalarValue(status.matrix_dir, '')
    const reportedProjectName = scalarValue(status.project_name, '').trim()
    // On a fresh browser session the first request cannot be project-scoped
    // because the project path is not in sessionStorage yet. The cockpit
    // persists settings in the project, so reload that config after the
    // project endpoint gives us the canonical path. Keep the initial global
    // config as a fallback for projects that have no saved state yet.
    let scopedGuiStateResponse = guiStateResponse
    if (projectPath && !requestedProjectPath) {
      scopedGuiStateResponse = await client.get('/gui_state', {
        params: { module: 'spectreasy_cockpit', project_path: projectPath },
        timeout: 5000,
      }).catch(() => guiStateResponse)
    }
    const scopedConfig = scopedGuiStateResponse?.data?.config
    const hasScopedConfig = scopedConfig && typeof scopedConfig === 'object' && Object.keys(scopedConfig).length > 0
    const savedConfig = unboxGuiState(hasScopedConfig ? scopedConfig : guiStateResponse?.data?.config) as Record<string, unknown> | undefined
    const savedSettings = (savedConfig?.settings ?? savedConfig ?? undefined) as Partial<WorkflowSettings> | undefined
    const [controlReportsResponse, sampleReportsResponse] = projectPath
      ? await Promise.all([
          client.get('/project/reports', {
            params: { project_path: projectPath, output_root: normalizeOutputRoot(savedSettings?.control?.outputDir), report_type: 'control' },
            timeout: 10000,
          }).catch(() => null),
          client.get('/project/reports', {
            params: { project_path: projectPath, output_root: normalizeOutputRoot(savedSettings?.sample?.outputDir), report_type: 'sample' },
            timeout: 10000,
          }).catch(() => null),
        ])
      : [null, null]
    const exactReports = [
      ...rowsFromBackend(controlReportsResponse?.data?.reports),
      ...rowsFromBackend(sampleReportsResponse?.data?.reports),
    ]
    const liveProjectArtifacts = liveArtifacts(
      projectPath,
      projectResponse?.data?.files,
      exactReports,
      controlInputDir,
      sampleInputDir,
    )
    const project: ProjectState = {
      ...emptyProject,
      projectName: reportedProjectName || projectNameFromPath(projectPath),
      projectPath,
      method: scalarValue(status.unmixing_method, emptyProject.method),
      cytometer: scalarValue(status.panel_cytometer, emptyProject.cytometer),
      artifacts: liveProjectArtifacts,
      mapping: liveMapping,
      missingInputDirs: stringsFromBackend(projectResponse?.data?.missing_input_dirs),
      controlInputDir,
      sampleInputDir,
      scan: { ...emptyProject.scan, ...backendScan, matrices: matrices.length, samples: samples.length },
    }
    return { project, backend, savedSettings }
  } catch {
    return { project: structuredClone(emptyProject), backend: initialBackendStatus }
  }
}

export async function updateProjectInputDirectory(
  role: 'controls' | 'samples',
  path: string,
  projectPath: string,
): Promise<{ success: boolean; message: string; controlInputDir?: string; sampleInputDir?: string }> {
  try {
    const response = await client.post('/project/layout', { role, path, projectPath }, { timeout: 30000 })
    if (!response.data?.success) return { success: false, message: scalarValue(response.data?.error, 'Could not update the input directory.') }
    return {
      success: true,
      message: `${role === 'controls' ? 'Control' : 'Sample'} input directory updated.`,
      controlInputDir: scalarValue(response.data?.layout?.control_input_dir),
      sampleInputDir: scalarValue(response.data?.layout?.sample_input_dir),
    }
  } catch (error) {
    const detail = axios.isAxiosError(error) ? scalarValue(error.response?.data?.error) : ''
    return { success: false, message: detail || 'Could not update the input directory.' }
  }
}

export async function loadProjectInputDirectories(projectPath: string): Promise<{ controlInputDir: string; sampleInputDir: string } | null> {
  try {
    const response = await client.get('/project/layout', { params: { project_path: projectPath }, timeout: 5000 })
    if (!response.data?.success) return null
    return {
      controlInputDir: scalarValue(response.data?.layout?.control_input_dir, 'scc'),
      sampleInputDir: scalarValue(response.data?.layout?.sample_input_dir, 'samples'),
    }
  } catch {
    return null
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

export async function listMatrixFiles(projectPath = ''): Promise<string[]> {
  try {
    const response = await client.get('/matrices', { params: { project_path: projectPath } })
    return Array.isArray(response.data) ? response.data.map((filename) => String(filename)) : []
  } catch {
    return []
  }
}

export async function listSampleFiles(projectPath = ''): Promise<string[]> {
  try {
    const response = await client.get('/samples', { params: { project_path: projectPath } })
    return Array.isArray(response.data) ? response.data.map((filename: unknown) => String(filename)) : []
  } catch {
    return []
  }
}

export type ProjectFileKind = 'controls' | 'samples'

export type ProjectFileEntry = {
  name: string
  size: number
  modified: string
  modifiedEpoch: number
  kind: ProjectFileKind
}

function blobAsBase64(blob: Blob, label: string): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader()
    reader.onerror = () => reject(reader.error ?? new Error(`Could not read ${label}.`))
    reader.onload = () => {
      const result = String(reader.result ?? '')
      const separator = result.indexOf(',')
      if (separator < 0) reject(new Error(`Could not encode ${label}.`))
      else resolve(result.slice(separator + 1))
    }
    reader.readAsDataURL(blob)
  })
}

export async function listProjectFiles(kind: ProjectFileKind, projectPath: string): Promise<{ success: boolean; files: ProjectFileEntry[]; message?: string }> {
  try {
    const response = await client.get('/project/files', { params: { kind, project_path: projectPath }, timeout: 10000 })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    const files = rowsFromBackend(response.data?.files).map((row) => ({
      name: scalarValue(row.name),
      size: Number(scalarValue(row.size, '0')) || 0,
      modified: scalarValue(row.modified),
      modifiedEpoch: Number(scalarValue(row.modified_epoch, '0')) || 0,
      kind: scalarValue(row.kind, kind) === 'samples' ? 'samples' : 'controls',
    } satisfies ProjectFileEntry)).filter((file) => file.name.length > 0)
    return { success, files, message: success ? undefined : scalarValue(response.data?.error, 'Project files could not be loaded.') }
  } catch {
    return { success: false, files: [], message: 'The local R backend did not answer.' }
  }
}

export async function uploadProjectFile(kind: ProjectFileKind, file: File, projectPath: string): Promise<{ success: boolean; message: string }> {
  if (!/\.fcs$/i.test(file.name)) return { success: false, message: `${file.name} is not an FCS file.` }
  let uploadId = ''
  try {
    const start = await client.post('/project/upload-start', {
      kind,
      filename: file.name,
      size: file.size,
      projectPath,
    }, { timeout: 30000 })
    const started = scalarValue(start.data?.success, 'false') === 'true'
    uploadId = scalarValue(start.data?.upload_id)
    if (!started || !uploadId) {
      return { success: false, message: scalarValue(start.data?.error, `${file.name} could not be added.`) }
    }
    const chunkSize = 4 * 1024 * 1024
    for (let offset = 0; offset < file.size; offset += chunkSize) {
      const chunk = file.slice(offset, Math.min(file.size, offset + chunkSize))
      const response = await client.post('/project/upload-chunk', {
        upload_id: uploadId,
        offset,
        content_base64: await blobAsBase64(chunk, file.name),
      }, { timeout: 120000 })
      if (scalarValue(response.data?.success, 'false') !== 'true') {
        return { success: false, message: scalarValue(response.data?.error, `${file.name} could not be added.`) }
      }
    }
    const response = await client.post('/project/upload-finish', { upload_id: uploadId }, { timeout: 120000 })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    return {
      success,
      message: success
        ? `${file.name} added.`
        : scalarValue(response.data?.error, `${file.name} could not be added.`),
    }
  } catch {
    if (uploadId) void client.post('/project/upload-abort', { upload_id: uploadId }, { timeout: 5000 }).catch(() => undefined)
    return { success: false, message: `${file.name} could not be added. The local R backend did not answer.` }
  }
}

export async function deleteProjectFile(kind: ProjectFileKind, filename: string, projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.delete('/project/files', { params: { kind, filename, project_path: projectPath }, timeout: 10000 })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    return {
      success,
      message: success
        ? `${filename} deleted.`
        : scalarValue(response.data?.error, `${filename} could not be deleted.`),
    }
  } catch {
    return { success: false, message: `${filename} could not be deleted. The local R backend did not answer.` }
  }
}

export async function deleteAllProjectFiles(kind: ProjectFileKind, projectPath: string): Promise<{ success: boolean; deleted: number; message: string }> {
  try {
    const response = await client.delete('/project/files/all', { params: { kind, project_path: projectPath }, timeout: 30000 })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    const deleted = Number(scalarValue(response.data?.deleted, '0')) || 0
    return {
      success,
      deleted,
      message: success
        ? `${deleted} file${deleted === 1 ? '' : 's'} deleted.`
        : scalarValue(response.data?.error, 'The files could not be deleted.'),
    }
  } catch {
    return { success: false, deleted: 0, message: 'The files could not be deleted. The local R backend did not answer.' }
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

export type AfProfileSummary = { name: string; bands: number; detectors: number; created: string; path: string; active: boolean }
export type AfProfileData = { name: string; detectors: string[]; spectra: Array<{ name: string; values: number[] }> }

export async function listAfProfiles(projectPath = ''): Promise<AfProfileSummary[]> {
  try {
    const response = await client.get('/af_profiles', { params: { project_path: projectPath } })
    return rowsFromBackend(response.data?.profiles).map((row) => ({
      name: String(row.name ?? ''),
      bands: Number(row.bands ?? 0),
      detectors: Number(row.detectors ?? 0),
      created: String(row.created ?? ''),
      path: String(row.path ?? ''),
      active: ['true', 't', '1', 'yes'].includes(String(row.active ?? '').toLowerCase()),
    })).filter((profile) => profile.name.length > 0)
  } catch {
    return []
  }
}

export async function loadAfProfileData(name: string): Promise<AfProfileData | null> {
  try {
    const response = await client.get('/af_profiles/data', { params: { name }, timeout: 5000 })
    if (response.data?.error) return null
    const detectors = Array.isArray(response.data?.detectors) ? response.data.detectors.map(String) : []
    const spectra = rowsFromBackend(response.data?.spectra).map((row) => ({
      name: String(row.name ?? ''),
      values: Array.isArray(row.values) ? row.values.map(Number) : [],
    })).filter((row) => row.name && row.values.length === detectors.length)
    return { name: String(response.data?.name ?? name), detectors, spectra }
  } catch {
    return null
  }
}

export async function selectAfSourceFile(projectPath = ''): Promise<{ success: boolean; cancelled: boolean; path?: string; message: string }> {
  try {
    const response = await client.post('/af_profiles/select-source', { projectPath }, { timeout: 0 })
    const cancelled = scalarValue(response.data?.cancelled, 'false') === 'true'
    if (cancelled) return { success: false, cancelled: true, message: 'File selection cancelled.' }
    if (response.data?.error) return { success: false, cancelled: false, message: scalarValue(response.data.error) }
    return { success: true, cancelled: false, path: scalarValue(response.data?.path), message: 'Source FCS selected.' }
  } catch {
    return { success: false, cancelled: false, message: 'Connect the local R backend to choose an FCS file.' }
  }
}

export async function loadProjectReports(projectPath = '', outputRoot = 'spectreasy_outputs', kind?: 'control' | 'sample'): Promise<Report[]> {
  const response = await client.get('/project/reports', {
    params: { project_path: projectPath, output_root: outputRoot, report_type: kind ?? '' },
  })
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
      createdEpoch: Number.isFinite(Number(row.created_epoch)) ? Number(row.created_epoch) : undefined,
      status: String(row.status ?? 'current') as Report['status'],
      matrix: '—',
      path: normalized,
    } satisfies Report
  }).filter((report) => report.path.length > 0)
}

export function projectFileUrl(path: string, projectPath = ''): string {
  const params = new URLSearchParams({ path, project_path: projectPath, token: resolveApiToken() })
  return `${API_BASE}/project/file?${params.toString()}`
}

export async function initializeProject(projectPath: string): Promise<{ success: boolean; created: string[]; message: string }> {
  try {
    const response = await client.post('/project/initialize', { projectPath }, { timeout: 30000 })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    const created = stringsFromBackend(response.data?.created)
    return {
      success,
      created,
      message: success
        ? `${created.length ? created.join(' and ') : 'Project input folders'} created.`
        : scalarValue(response.data?.error, 'Project input folders could not be created.'),
    }
  } catch {
    return { success: false, created: [], message: 'The local R backend did not answer.' }
  }
}

function downloadBlob(filename: string, blob: Blob) {
  const href = URL.createObjectURL(blob)
  const anchor = document.createElement('a')
  anchor.href = href
  anchor.download = filename
  anchor.click()
  URL.revokeObjectURL(href)
}

export async function downloadProjectReport(path: string, projectPath = ''): Promise<boolean> {
  try {
    const response = await client.get('/project/file', { params: { path, project_path: projectPath }, responseType: 'blob', timeout: 120000 })
    downloadBlob(path.split('/').pop() || 'spectreasy_qc_report.html', response.data as Blob)
    return true
  } catch {
    return false
  }
}

export async function exportProjectReportPdf(path: string, projectPath = ''): Promise<boolean> {
  try {
    const response = await client.post('/project/report/export-pdf', { path, projectPath }, { timeout: 120000 })
    if (response.data?.error || !response.data?.content_base64) return false
    downloadBase64File(
      scalarValue(response.data.filename, 'spectreasy_qc_report.pdf'),
      scalarValue(response.data.content_base64),
      scalarValue(response.data.content_type, 'application/pdf'),
    )
    return true
  } catch {
    return false
  }
}

export async function deleteAfProfile(name: string, projectPath = ''): Promise<boolean> {
  try {
    const response = await client.delete('/af_profiles/delete', { params: { name, project_path: projectPath } })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function renameAfProfile(profileName: string, newName: string, projectPath: string): Promise<{ success: boolean; message: string; name?: string }> {
  try {
    const response = await client.post('/af_profiles/rename', { profile_name: profileName, new_name: newName, projectPath })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    return {
      success,
      name: success ? scalarValue(response.data?.new_name, newName) : undefined,
      message: success ? `${profileName} renamed to ${newName}.` : scalarValue(response.data?.error, `Could not rename ${profileName}.`),
    }
  } catch {
    return { success: false, message: `Could not rename ${profileName}.` }
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

export async function activateAfProfile(profileName: string, projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/af_profiles/activate', { profile_name: profileName, use_as_unstained: true, projectPath })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    return { success, message: success ? `${profileName} is linked to this dataset.` : scalarValue(response.data?.error, `Could not link ${profileName}.`) }
  } catch {
    return { success: false, message: `Could not link ${profileName}.` }
  }
}

export async function deactivateAfProfile(profileName: string, projectPath: string): Promise<boolean> {
  try {
    const response = await client.post('/af_profiles/deactivate', { profile_name: profileName, projectPath })
    return scalarValue(response.data?.success, 'false') === 'true'
  } catch {
    return false
  }
}

export async function persistGuiState(project: ProjectState, settings?: WorkflowSettings): Promise<boolean> {
  try {
    await client.post('/gui_state', {
      module: 'spectreasy_cockpit',
      projectPath: project.projectPath,
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

export async function persistControlMapping(mapping: ProjectState['mapping'], projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/control_mapping', { rows: mapping, projectPath })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    const error = scalarValue(response.data?.error, '')
    return success ? { success: true, message: 'Control mapping saved to the project CSV.' } : { success: false, message: error || 'The R backend rejected the control mapping.' }
  } catch {
    return { success: false, message: 'Control mapping could not be saved. Connect the R backend and retry.' }
  }
}

export async function createControlMapping(projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/control_mapping/create', { projectPath })
    const success = scalarValue(response.data?.success, 'false') === 'true'
    return success
      ? { success: true, message: 'Created fcs_mapping.csv from the SCC folder.' }
      : { success: false, message: scalarValue(response.data?.error, 'The control mapping could not be created.') }
  } catch {
    return { success: false, message: 'The control mapping could not be created. Select a project with an SCC folder and retry.' }
  }
}

export async function setProjectContext(projectPath: string): Promise<{ success: boolean; message: string }> {
  try {
    const response = await client.post('/project/context', { projectPath })
    if (response.data?.success === false || response.data?.error) {
      return { success: false, message: response.data?.error ?? 'The R backend rejected this project folder.' }
    }
    return { success: true, message: '' }
  } catch {
    return { success: false, message: 'The project folder could not be opened. Check the path and R connection.' }
  }
}

async function pickProjectFolder(endpoint: '/project/select' | '/project/create'): Promise<{ success: boolean; cancelled: boolean; message: string; projectPath?: string }> {
  try {
    // Native file dialogs intentionally wait for the user. They must not inherit
    // the short timeout used by ordinary API calls.
    const response = await client.post(endpoint, {}, { timeout: 0 })
    const cancelled = scalarValue(response.data?.cancelled, 'false') === 'true'
    if (cancelled) return { success: false, cancelled: true, message: 'Project selection cancelled.' }
    if (response.data?.success === false || response.data?.error) {
      return { success: false, cancelled: false, message: scalarValue(response.data?.error, 'The project folder could not be opened.') }
    }
    return { success: true, cancelled: false, message: '', projectPath: scalarValue(response.data?.project?.project_path) }
  } catch (error) {
    if (axios.isAxiosError(error) && error.response?.status === 403) {
      return {
        success: false,
        cancelled: false,
        message: 'This browser tab belongs to an earlier R session. Use the cockpit tab opened by the current spectreasy_gui() command.',
      }
    }
    return {
      success: false,
      cancelled: false,
      message: 'The local R backend did not answer. Keep spectreasy_gui() running and retry.',
    }
  }
}

export function selectProjectFolder() {
  return pickProjectFolder('/project/select')
}

export function createProjectFolder() {
  return pickProjectFolder('/project/create')
}

export type WorkflowActionResult = {
  success: boolean
  backendReachable: boolean
  message: string
  outputCount: number
  logs: string[]
}

export async function attemptWorkflowAction(action: string, payload: Record<string, unknown>): Promise<WorkflowActionResult> {
  const endpoint = action === 'control' ? '/workflow/control' : action === 'sample' ? '/workflow/sample' : action === 'af' ? '/workflow/af' : '/workflow/report'
  try {
    // Scientific workflows can legitimately take longer than 15 minutes. A
    // client-side timeout abandons only the HTTP request; R keeps computing and
    // can finish successfully, leaving the cockpit with a false failure.
    const response = await client.post(endpoint, payload, { timeout: 0 })
    const success = scalarValue(response.data?.success, 'false')
    const error = scalarValue(response.data?.error, '')
    const logs = stringsFromBackend(response.data?.logs)
    if (success === 'false' || error) {
      return { success: false, backendReachable: true, message: error || 'The R backend rejected this workflow action.', outputCount: 0, logs }
    }
    const result = response.data?.result
    const outputCount = result && typeof result === 'object'
      ? Object.values(result as Record<string, unknown>).filter((value) => value != null && scalarValue(value).length > 0).length
      : 0
    return { success: true, backendReachable: true, message: 'Workflow completed in the local R session.', outputCount, logs }
  } catch (error) {
    const status = axios.isAxiosError(error) ? error.response?.status : undefined
    return {
      success: false,
      backendReachable: status != null,
      message: status === 403
        ? 'This cockpit tab belongs to a different R session. Reopen the URL printed by the active spectreasy_gui() process.'
        : 'The local R backend did not answer. Keep spectreasy_gui() running and retry.',
      outputCount: 0,
      logs: [],
    }
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
