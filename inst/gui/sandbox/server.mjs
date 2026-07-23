import { createServer } from 'node:http'
import { spawn } from 'node:child_process'
import { mkdir, readFile, readdir, rename, rm, stat } from 'node:fs/promises'
import { arch, homedir, platform } from 'node:os'
import { dirname, join, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const here = dirname(fileURLToPath(import.meta.url))
const repositoryRoot = resolve(here, '../../..')
const sandboxRoot = resolve(
  process.env.SPECTREASY_SANDBOX_ROOT || join(homedir(), '.cache', 'spectreasy-onboarding-sandbox'),
)
const environmentPrefix = join(sandboxRoot, 'environment')
const mambaRoot = join(sandboxRoot, 'mamba-root')
const toolsRoot = join(sandboxRoot, 'tools')
const micromamba = join(toolsRoot, 'bin', 'micromamba')
const projectRoot = join(sandboxRoot, 'project')
const userLibrary = join(sandboxRoot, 'user-library')
const packageBuildRoot = join(sandboxRoot, 'package-build')
const port = Number(process.env.SPECTREASY_SANDBOX_PORT || 8787)
const minimumRVersion = '4.5.0'
const currentRVersion = process.env.SPECTREASY_SANDBOX_R_VERSION || '4.5.3'
const oldRVersion = process.env.SPECTREASY_SANDBOX_OLD_R_VERSION || '4.4.3'
const guiDependencies = ['httpuv', 'later', 'plumber', 'jsonlite']
const maxLogLines = 500

const job = {
  busy: false,
  action: '',
  phase: 'idle',
  error: '',
  log: [],
  startedAt: null,
  finishedAt: null,
}

function appendLog(value) {
  const lines = String(value).replace(/\r/g, '').split('\n').filter(Boolean)
  job.log.push(...lines)
  if (job.log.length > maxLogLines) job.log.splice(0, job.log.length - maxLogLines)
}

async function exists(path) {
  try {
    await stat(path)
    return true
  } catch {
    return false
  }
}

function executable(name) {
  return platform() === 'win32' ? join(environmentPrefix, name === 'R' ? 'R.exe' : `${name}.exe`) : join(environmentPrefix, 'bin', name)
}

function isolatedEnvironment(extra = {}) {
  const nullDevice = platform() === 'win32' ? 'NUL' : '/dev/null'
  const pathSeparator = platform() === 'win32' ? ';' : ':'
  return {
    ...process.env,
    MAMBA_ROOT_PREFIX: mambaRoot,
    CONDA_PREFIX: environmentPrefix,
    PATH: `${join(environmentPrefix, 'bin')}${pathSeparator}${process.env.PATH || ''}`,
    PKG_CONFIG_PATH: join(environmentPrefix, 'lib', 'pkgconfig'),
    R_LIBS_USER: userLibrary,
    R_ENVIRON_USER: nullDevice,
    R_PROFILE_USER: nullDevice,
    ...extra,
  }
}

function platformArchive() {
  const key = `${platform()}-${arch()}`
  return {
    'darwin-arm64': 'osx-arm64',
    'darwin-x64': 'osx-64',
    'linux-x64': 'linux-64',
    'linux-arm64': 'linux-aarch64',
  }[key] || null
}

function run(command, args, options = {}) {
  return new Promise((resolvePromise, reject) => {
    appendLog(`$ ${[command, ...args].join(' ')}`)
    const child = spawn(command, args, {
      cwd: options.cwd || repositoryRoot,
      env: isolatedEnvironment(options.env),
      stdio: ['ignore', 'pipe', 'pipe'],
    })
    child.stdout.on('data', appendLog)
    child.stderr.on('data', appendLog)
    child.on('error', reject)
    child.on('close', (code) => {
      if (code === 0) resolvePromise()
      else reject(new Error(`${command} exited with status ${code}`))
    })
  })
}

async function capture(command, args) {
  let output = ''
  await new Promise((resolvePromise, reject) => {
    const child = spawn(command, args, {
      cwd: repositoryRoot,
      env: isolatedEnvironment(),
      stdio: ['ignore', 'pipe', 'pipe'],
    })
    child.stdout.on('data', (chunk) => { output += chunk })
    child.stderr.on('data', (chunk) => { output += chunk })
    child.on('error', reject)
    child.on('close', (code) => {
      if (code === 0) return resolvePromise()
      const error = new Error(output.trim() || `${command} exited with status ${code}`)
      error.exitCode = code
      reject(error)
    })
  })
  return output.trim()
}

function compareVersions(left, right) {
  const a = left.split('.').map(Number)
  const b = right.split('.').map(Number)
  for (let index = 0; index < Math.max(a.length, b.length); index += 1) {
    const difference = (a[index] || 0) - (b[index] || 0)
    if (difference !== 0) return difference
  }
  return 0
}

async function inspectEnvironment() {
  const rscript = executable('Rscript')
  if (!(await exists(rscript))) {
    return { r: 'missing', rVersion: null, spectreasy: 'missing', spectreasyVersion: null, detail: 'No R executable exists in the sandbox prefix.' }
  }

  let rVersion = null
  try {
    rVersion = await capture(rscript, ['--vanilla', '-e', 'cat(as.character(getRversion()))'])
  } catch (error) {
    return { r: 'broken', rVersion: null, spectreasy: 'unknown', spectreasyVersion: null, detail: error.message }
  }

  if (compareVersions(rVersion, minimumRVersion) < 0) {
    return { r: 'outdated', rVersion, spectreasy: 'missing', spectreasyVersion: null, detail: `Spectreasy requires R ${minimumRVersion} or newer.` }
  }

  try {
    const result = await capture(rscript, ['--vanilla', '-e', [
      'package_path <- suppressWarnings(find.package("spectreasy", quiet = TRUE))',
      'if (!length(package_path) || !nzchar(package_path)) quit(status = 42)',
      'load_error <- tryCatch({ invisible(loadNamespace("spectreasy")); NULL }, error = function(error) conditionMessage(error))',
      'if (!is.null(load_error)) { cat(load_error); quit(status = 44) }',
      `gui_dependencies <- c(${guiDependencies.map((name) => JSON.stringify(name)).join(', ')})`,
      'missing_gui_dependencies <- gui_dependencies[!vapply(gui_dependencies, requireNamespace, logical(1), quietly = TRUE)]',
      'if (length(missing_gui_dependencies)) { cat(paste(missing_gui_dependencies, collapse = ", ")); quit(status = 43) }',
      'cat(as.character(utils::packageVersion("spectreasy")))',
    ].join(';')])
    return { r: 'ready', rVersion, spectreasy: 'ready', spectreasyVersion: result, detail: 'The isolated package namespace loads successfully.' }
  } catch (error) {
    const missing = error.exitCode === 42
    const missingGuiDependencies = error.exitCode === 43
    return {
      r: 'ready',
      rVersion,
      spectreasy: missing ? 'missing' : 'broken',
      spectreasyVersion: null,
      detail: missing
        ? 'Spectreasy is not installed in the isolated R library.'
        : missingGuiDependencies
          ? `Spectreasy is installed, but GUI dependencies are missing: ${error.message}.`
          : error.message,
    }
  }
}

async function ensureMicromamba() {
  if (await exists(micromamba)) return
  const archivePlatform = platformArchive()
  if (!archivePlatform) throw new Error(`Real sandbox installs are not supported on ${platform()} ${arch()}. Use preview states instead.`)
  job.phase = 'bootstrap'
  appendLog(`Downloading the repo-local Micromamba runtime for ${archivePlatform}.`)
  await mkdir(toolsRoot, { recursive: true })
  const archive = join(toolsRoot, 'micromamba.tar.bz2')
  const temporary = `${archive}.download`
  await run('curl', ['-fL', `https://micro.mamba.pm/api/micromamba/${archivePlatform}/latest`, '-o', temporary])
  await rename(temporary, archive)
  await run('tar', ['-xjf', archive, '-C', toolsRoot, 'bin/micromamba'])
  await rm(archive, { force: true })
}

async function createREnvironment(version) {
  await ensureMicromamba()
  job.phase = 'install-r'
  await rm(environmentPrefix, { recursive: true, force: true })
  await rm(userLibrary, { recursive: true, force: true })
  await mkdir(userLibrary, { recursive: true })
  await run(micromamba, [
    'create', '--yes', '--override-channels', '--prefix', environmentPrefix, '--channel', 'conda-forge',
    `r-base=${version}`, 'r-remotes', 'libsodium', 'libuv', 'pkg-config',
  ])
}

async function ensureSystemDependencies() {
  await ensureMicromamba()
  job.phase = 'install-system-dependencies'
  await run(micromamba, [
    'install', '--yes', '--freeze-installed', '--override-channels', '--prefix', environmentPrefix,
    '--channel', 'conda-forge', 'libsodium', 'libuv', 'pkg-config',
  ])
}

async function installSpectreasy() {
  const rscript = executable('Rscript')
  if (!(await exists(rscript))) throw new Error('Install R in the sandbox before installing Spectreasy.')
  const inspected = await inspectEnvironment()
  if (inspected.r !== 'ready') throw new Error(`R ${minimumRVersion} or newer is required before installing Spectreasy.`)

  await ensureSystemDependencies()
  job.phase = 'build-package'
  await rm(packageBuildRoot, { recursive: true, force: true })
  await mkdir(packageBuildRoot, { recursive: true })
  await run(executable('R'), ['CMD', 'build', '--no-build-vignettes', '--no-manual', repositoryRoot], { cwd: packageBuildRoot })
  const tarballs = (await readdir(packageBuildRoot)).filter((file) => /^spectreasy_.*\.tar\.gz$/.test(file))
  if (tarballs.length !== 1) throw new Error('The sandbox package build did not produce exactly one source tarball.')

  job.phase = 'install-dependencies'
  const packageTarball = JSON.stringify(join(packageBuildRoot, tarballs[0]))
  const installExpression = [
    'options(repos = c(CRAN = "https://cloud.r-project.org"))',
    `gui_dependencies <- c(${guiDependencies.map((name) => JSON.stringify(name)).join(', ')})`,
    'missing_gui_dependencies <- gui_dependencies[!vapply(gui_dependencies, requireNamespace, logical(1), quietly = TRUE)]',
    'if (length(missing_gui_dependencies)) install.packages(missing_gui_dependencies)',
    `remotes::install_local(${packageTarball}, dependencies = NA, upgrade = "never", force = TRUE)`,
  ].join(';')
  await run(rscript, ['--vanilla', '-e', installExpression])
  job.phase = 'verify-spectreasy'
  const verified = await inspectEnvironment()
  if (verified.spectreasy !== 'ready') throw new Error(`Installation finished, but Spectreasy did not load: ${verified.detail}`)
}

async function resetSandbox() {
  job.phase = 'reset'
  await rm(environmentPrefix, { recursive: true, force: true })
  await rm(userLibrary, { recursive: true, force: true })
  await rm(projectRoot, { recursive: true, force: true })
  await mkdir(projectRoot, { recursive: true })
  appendLog('Removed the isolated R prefix and sandbox project. Download caches were retained.')
}

async function purgeSandbox() {
  job.phase = 'purge'
  await rm(sandboxRoot, { recursive: true, force: true })
  await mkdir(projectRoot, { recursive: true })
  await mkdir(userLibrary, { recursive: true })
  appendLog('Deleted the sandbox R runtime, libraries, tools, project, and download caches.')
}

async function prepareROnly() {
  const inspected = await inspectEnvironment()
  if (inspected.r !== 'ready') {
    await createREnvironment(currentRVersion)
    return
  }
  job.phase = 'remove-spectreasy'
  await rm(join(userLibrary, 'spectreasy'), { recursive: true, force: true })
  await rm(join(environmentPrefix, 'lib', 'R', 'library', 'spectreasy'), { recursive: true, force: true })
  appendLog('Removed Spectreasy only; R and dependency packages were retained.')
}

async function breakSpectreasy() {
  const packagePath = (await exists(join(userLibrary, 'spectreasy')))
    ? join(userLibrary, 'spectreasy')
    : join(environmentPrefix, 'lib', 'R', 'library', 'spectreasy')
  if (!(await exists(packagePath))) throw new Error('Install Spectreasy before creating a broken-package state.')
  job.phase = 'break-spectreasy'
  const libraryDirectory = join(packagePath, 'libs')
  const libraries = (await readdir(libraryDirectory)).filter((file) => /\.(so|dll|dylib)$/.test(file))
  if (!libraries.length) throw new Error('Could not find the installed Spectreasy shared library to disable.')
  await rm(join(libraryDirectory, libraries[0]), { force: true })
  appendLog('Removed the sandbox package shared library to create a real load failure.')
}

async function prepareReady() {
  const inspected = await inspectEnvironment()
  if (inspected.r !== 'ready') await createREnvironment(currentRVersion)
  const afterR = await inspectEnvironment()
  if (afterR.spectreasy !== 'ready') await installSpectreasy()
}

async function prepareBrokenSpectreasy() {
  const inspected = await inspectEnvironment()
  if (inspected.spectreasy !== 'ready') await prepareReady()
  await breakSpectreasy()
}

function startJob(action, operation) {
  if (job.busy) throw new Error(`The sandbox is already running: ${job.action}`)
  Object.assign(job, {
    busy: true,
    action,
    phase: 'starting',
    error: '',
    log: [],
    startedAt: new Date().toISOString(),
    finishedAt: null,
  })
  Promise.resolve()
    .then(operation)
    .then(() => appendLog('Sandbox action completed successfully.'))
    .catch((error) => {
      job.error = error instanceof Error ? error.message : String(error)
      appendLog(`ERROR: ${job.error}`)
    })
    .finally(() => {
      job.busy = false
      job.phase = job.error ? 'failed' : 'complete'
      job.finishedAt = new Date().toISOString()
    })
}

function sendJson(response, status, value) {
  const body = JSON.stringify(value)
  response.writeHead(status, {
    'Content-Type': 'application/json; charset=utf-8',
    'Content-Length': Buffer.byteLength(body),
    'Cache-Control': 'no-store',
    'Access-Control-Allow-Origin': '*',
    'Access-Control-Allow-Headers': 'Content-Type, X-Spectreasy-Token',
    'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
  })
  response.end(body)
}

async function requestBody(request) {
  let body = ''
  for await (const chunk of request) body += chunk
  return body ? JSON.parse(body) : {}
}

async function sandboxStatus() {
  return {
    environment: await inspectEnvironment(),
    job: { ...job, log: job.log.slice(-120) },
    paths: { sandboxRoot, environmentPrefix, repositoryRoot },
    minimumRVersion,
    supported: Boolean(platformArchive()),
  }
}

await mkdir(projectRoot, { recursive: true })
await mkdir(userLibrary, { recursive: true })

const server = createServer(async (request, response) => {
  const url = new URL(request.url || '/', `http://${request.headers.host}`)
  if (request.method === 'OPTIONS') return sendJson(response, 204, {})

  try {
    if (request.method === 'GET' && url.pathname === '/sandbox/status') {
      return sendJson(response, 200, await sandboxStatus())
    }
    if (request.method === 'POST' && url.pathname === '/sandbox/action') {
      const { action } = await requestBody(request)
      const operations = {
        reset: resetSandbox,
        purge: purgeSandbox,
        'install-r': () => createREnvironment(currentRVersion),
        'install-old-r': () => createREnvironment(oldRVersion),
        'prepare-r-only': prepareROnly,
        'break-spectreasy': prepareBrokenSpectreasy,
        'install-spectreasy': installSpectreasy,
        'prepare-ready': prepareReady,
      }
      if (!operations[action]) return sendJson(response, 400, { error: 'Unknown sandbox action.' })
      startJob(action, operations[action])
      return sendJson(response, 202, { accepted: true, action })
    }

    // Minimal, harmless cockpit API used after the sandbox reaches the ready state.
    if (request.method === 'GET' && url.pathname === '/status') {
      const environment = await inspectEnvironment()
      return sendJson(response, 200, {
        status: 'ok',
        version: environment.spectreasyVersion || 'Sandbox preview',
        gui_mode: 'Onboarding sandbox',
        project_name: 'Sandbox project',
        panel_cytometer: 'aurora',
        unmixing_method: 'AutoSpectral',
      })
    }
    if (request.method === 'GET' && ['/matrices', '/samples'].includes(url.pathname)) return sendJson(response, 200, [])
    if (request.method === 'GET' && url.pathname === '/control_mapping') return sendJson(response, 200, { rows: [] })
    if (request.method === 'GET' && url.pathname === '/gui_state') return sendJson(response, 200, { config: {} })
    if (request.method === 'POST' && url.pathname === '/gui_state') return sendJson(response, 200, { success: true })
    if (request.method === 'POST' && url.pathname === '/spectral_panel_metrics') return sendJson(response, 200, {})
    if (request.method === 'GET' && url.pathname === '/project/status') {
      return sendJson(response, 200, { project_path: projectRoot, files: [], scan: { controls: 0, samples: 0, matrices: 0, reports: 0 } })
    }
    if (request.method === 'GET' && url.pathname === '/sandbox/log') {
      const logPath = join(sandboxRoot, 'last.log')
      return sendJson(response, 200, { log: await readFile(logPath, 'utf8').catch(() => job.log.join('\n')) })
    }
    return sendJson(response, 404, { error: 'Not found' })
  } catch (error) {
    return sendJson(response, 500, { error: error instanceof Error ? error.message : String(error) })
  }
})

server.listen(port, '127.0.0.1', () => {
  console.log(`Spectreasy onboarding sandbox listening at http://127.0.0.1:${port}`)
  console.log(`Isolated files: ${sandboxRoot}`)
})
