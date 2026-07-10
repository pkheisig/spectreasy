import type { Artifact, MappingRow, ProjectState, Report } from './types'

export const demoArtifacts: Artifact[] = [
  { id: 'controls', name: 'controls/', type: 'Folder', group: 'Controls', status: 'current', detail: '14 FCS files · 1 unstained', path: '/projects/lymphocyte-panel/controls', updated: 'Today, 09:18' },
  { id: 'mapping', name: 'fcs_mapping.csv', type: 'Mapping', group: 'Controls', status: 'current', detail: '14 rows · validated 09:22', path: '/projects/lymphocyte-panel/fcs_mapping.csv', updated: 'Today, 09:22', run: 'run-014' },
  { id: 'gates', name: 'ssc_gate_config.csv', type: 'Gates', group: 'Gates', status: 'current', detail: '7 gates · manual', path: '/projects/lymphocyte-panel/ssc_gate_config.csv', updated: 'Today, 09:27', run: 'run-014' },
  { id: 'samples', name: 'samples/', type: 'Folder', group: 'Samples', status: 'user', detail: '28 FCS files · 2 new', path: '/projects/lymphocyte-panel/samples', updated: 'Today, 10:04' },
  { id: 'reference', name: 'reference_matrix.csv', type: 'Reference matrix', group: 'Matrices', status: 'current', detail: '32 markers × 47 detectors', path: '/projects/lymphocyte-panel/Spectreasy/reference_matrix.csv', updated: 'Today, 09:41', run: 'run-014' },
  { id: 'unmixing', name: 'unmixing_matrix.csv', type: 'Static unmixing matrix', group: 'Matrices', status: 'current', detail: '47 detectors × 32 markers', path: '/projects/lymphocyte-panel/Spectreasy/unmixing_matrix.csv', updated: 'Today, 09:41', run: 'run-014' },
  { id: 'noise', name: 'detector_noise.csv', type: 'Detector noise', group: 'Matrices', status: 'current', detail: '47 detector floors', path: '/projects/lymphocyte-panel/Spectreasy/detector_noise.csv', updated: 'Today, 09:41', run: 'run-014' },
  { id: 'variants', name: 'spectral_variants.rds', type: 'Spectral library', group: 'Matrices', status: 'current', detail: '96 variants · 18 clusters', path: '/projects/lymphocyte-panel/Spectreasy/spectral_variants.rds', updated: 'Today, 09:43', run: 'run-014' },
  { id: 'control-report', name: 'control_qc_report.html', type: 'QC report', group: 'Reports', status: 'current', detail: 'HTML · 6 sections', path: '/projects/lymphocyte-panel/reports/run-014/control_qc_report.html', updated: 'Today, 09:47', run: 'run-014' },
  { id: 'control-pdf', name: 'control_qc_report.pdf', type: 'QC report', group: 'Reports', status: 'current', detail: 'PDF · 18 pages', path: '/projects/lymphocyte-panel/reports/run-014/control_qc_report.pdf', updated: 'Today, 09:48', run: 'run-014' },
  { id: 'sample-report', name: 'sample_qc_report.html', type: 'QC report', group: 'Reports', status: 'stale', detail: 'Stale · 2 new samples', path: '/projects/lymphocyte-panel/reports/run-013/sample_qc_report.html', updated: 'Yesterday, 16:06', run: 'run-013' },
  { id: 'metrics', name: 'qc_metrics.csv', type: 'Metrics', group: 'QC Metrics', status: 'stale', detail: '1,248 rows · needs refresh', path: '/projects/lymphocyte-panel/reports/run-013/qc_metrics.csv', updated: 'Yesterday, 16:06', run: 'run-013' },
  { id: 'af-profile', name: 'PBMC broad AF', type: 'AF profile', group: 'AF Profiles', status: 'user', detail: '12 bands · 47 detectors', path: '/Users/pkheisig/.spectreasy/af_profiles/PBMC_broad_AF.rds', updated: 'Jun 28, 2026' },
  { id: 'panel-export', name: 'panel_overview.pdf', type: 'Panel export', group: 'Panel Builder', status: 'user', detail: 'Aurora · 18 fluorophores', path: '/projects/lymphocyte-panel/exports/panel_overview.pdf', updated: 'Jun 26, 2026' },
  { id: 'job-log', name: 'run-014.json', type: 'Run manifest', group: 'Logs', status: 'current', detail: 'Control workflow · 4m 12s', path: '/projects/lymphocyte-panel/logs/run-014.json', updated: 'Today, 09:48', run: 'run-014' },
]

export const demoMapping: MappingRow[] = [
  { id: 'r1', file: 'unstained_cells.fcs', fluorophore: 'Unstained', marker: 'AF', channel: 'Violet 405 / B1', controlType: 'unstained', universalNegative: true },
  { id: 'r2', file: 'bead_neg.fcs', fluorophore: 'Bead negative', marker: 'Bead negative', channel: 'Violet 405 / B1', controlType: 'bead', universalNegative: true },
  { id: 'r3', file: 'BV421.fcs', fluorophore: 'BV421', marker: 'CD45', channel: 'Violet 450 / B2', controlType: 'cell', universalNegative: false },
  { id: 'r4', file: 'FITC.fcs', fluorophore: 'FITC', marker: 'CD3', channel: 'Blue 530 / B3', controlType: 'cell', universalNegative: false },
  { id: 'r5', file: 'PE.fcs', fluorophore: 'PE', marker: 'CD19', channel: 'YellowGreen 586 / YG2', controlType: 'cell', universalNegative: false },
  { id: 'r6', file: 'APC.fcs', fluorophore: 'APC', marker: 'CD8', channel: 'Red 670 / R4', controlType: 'cell', universalNegative: false },
  { id: 'r7', file: 'AF-viability.fcs', fluorophore: 'Zombie NIR', marker: 'Viability', channel: 'Red 780 / R7', controlType: 'viability', universalNegative: false, warning: 'No matched viability negative' },
]

export const demoReports: Report[] = [
  { id: 'r-control-html', title: 'Control QC · run-014', type: 'Control QC', format: 'HTML', run: 'run-014', created: 'Today, 09:47', status: 'current', matrix: 'reference_matrix.csv' },
  { id: 'r-control-pdf', title: 'Control QC · run-014', type: 'Control QC', format: 'PDF', run: 'run-014', created: 'Today, 09:48', status: 'current', matrix: 'reference_matrix.csv' },
  { id: 'r-sample-html', title: 'Sample QC · run-013', type: 'Sample QC', format: 'HTML', run: 'run-013', created: 'Yesterday, 16:06', status: 'stale', matrix: 'reference_matrix.csv' },
  { id: 'r-sample-pdf', title: 'Sample QC · run-013', type: 'Sample QC', format: 'PDF', run: 'run-013', created: 'Yesterday, 16:07', status: 'stale', matrix: 'reference_matrix.csv' },
  { id: 'r-panel-pdf', title: 'Aurora panel overview', type: 'Panel overview', format: 'PDF', run: 'panel-006', created: 'Jun 26, 2026', status: 'current', matrix: '—' },
]

export const demoProject: ProjectState = {
  projectName: 'Lymphocyte panel',
  projectPath: '/projects/lymphocyte-panel',
  cytometer: 'Cytek Aurora 5L',
  method: 'Spectreasy',
  lastAction: 'Control QC report generated',
  lastActionAt: 'Today, 09:48',
  artifacts: demoArtifacts,
  mapping: demoMapping,
  mappingDirty: false,
  gatesDirty: false,
  scan: { controls: 14, samples: 28, matrices: 3, reports: 5, qcMetrics: 1, spectralVariants: 1, gates: 1 },
}

export const workflowStatus: Record<string, { state: 'complete' | 'ready' | 'warning' | 'blocked' | 'stale' | 'idle'; note: string }> = {
  setup: { state: 'complete', note: 'Project manifest loaded' },
  controlImport: { state: 'complete', note: '14 controls detected' },
  mapping: { state: 'complete', note: 'Validated 09:22' },
  gating: { state: 'complete', note: '7 saved gates' },
  controlRun: { state: 'complete', note: 'run-014 · 4m 12s' },
  controlReport: { state: 'complete', note: 'HTML report ready' },
  sampleImport: { state: 'warning', note: '2 new samples since last run' },
  sampleRun: { state: 'ready', note: 'Matrix selected · 28 samples' },
  sampleReport: { state: 'stale', note: 'Regenerate after sample run' },
}
