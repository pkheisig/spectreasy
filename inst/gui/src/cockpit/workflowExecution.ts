import type { WorkflowSettings } from "./types";

export type WorkflowAction = "control" | "sample" | "af";

function rValue(value: unknown): string {
  if (value == null) return "NULL";
  if (typeof value === "boolean") return value ? "TRUE" : "FALSE";
  if (typeof value === "number") return String(value);
  return JSON.stringify(String(value));
}

export function cockpitExecutionCode(action: WorkflowAction, payload: Record<string, unknown>): string {
  const projectPath = rValue(payload.projectPath);
  const argumentAliases: Record<string, string> = {
    method: "unmixing_method",
    matrix_file: "unmixing_matrix_file",
  };
  const functionName = action === "control"
    ? "unmix_controls"
    : action === "sample"
      ? "unmix_samples"
      : "extract_af_profile";
  const excluded = new Set(["projectPath", "save_name", "save_overwrite"]);
  const pathArguments = new Set(
    action === "control"
      ? ["scc_dir", "control_file", "output_dir", "manual_gate_file"]
      : action === "sample"
        ? ["sample_dir", "matrix_file", "detector_noise_file", "output_dir", "spectral_variant_library_file"]
        : ["fcs_file"],
  );
  const argumentValue = (key: string, value: unknown): string => {
    if (!pathArguments.has(key) || typeof value !== "string") return rValue(value);
    const path = value.trim();
    if (!path) return "NULL";
    if (/^(?:\/|[A-Za-z]:[\\/])/.test(path)) return rValue(path);
    return `file.path(project_dir, ${rValue(path)})`;
  };
  const args = Object.entries(payload)
    .filter(([key]) => !excluded.has(key))
    .map(([key, value]) => `  ${argumentAliases[key] ?? key} = ${argumentValue(key, value)}`)
    .join(",\n");
  return `# Cockpit-generated R operation\nproject_dir <- ${projectPath}\n${functionName}(\n${args}\n)`;
}

export function normalizeCockpitOutputRoot(value: unknown): string {
  const output = String(value ?? "spectreasy_outputs").trim().replace(/[\\/]+$/, "");
  return output.replace(/[\\/](?:unmix_controls|unmix_samples)$/i, "") || "spectreasy_outputs";
}

export function createWorkflowPayload(action: WorkflowAction, settings: WorkflowSettings): Record<string, unknown> {
  const { control, sample, af } = settings;
  if (action === "control") {
    return {
      projectPath: settings.projectPath,
      scc_dir: control.sccDir,
      control_file: control.controlFile,
      output_dir: normalizeCockpitOutputRoot(control.outputDir),
      method: control.method,
      cytometer: control.cytometer,
      auto_create_mapping: control.autoCreateMapping,
      auto_unknown_fluor_policy: control.autoUnknownFluorPolicy,
      manual_gate_file: "ssc_gate_config.csv",
      gating_mode: "reuse",
      af_n_bands: control.afNBands,
      af_max_cells: control.afMaxCells,
      default_sample_type: control.defaultSampleType,
      histogram_pct_beads: control.histogramPctBeads,
      histogram_direction_beads: control.histogramDirectionBeads,
      histogram_pct_cells: control.histogramPctCells,
      histogram_direction_cells: control.histogramDirectionCells,
      outlier_percentile: control.outlierPercentile,
      debris_percentile: control.debrisPercentile,
      bead_gate_scale: control.beadGateScale,
      max_clusters: control.maxClusters,
      min_cluster_proportion: control.minClusterProportion,
      gate_contour_beads: control.gateContourBeads,
      gate_contour_cells: control.gateContourCells,
      subsample_n: control.subsampleN,
      unmix_scatter_panel_size_mm: control.unmixScatterPanelSizeMm,
      rwls_max_iter: control.rwlsMaxIter,
      n_threads: control.unmixThreads,
      seed: control.seed,
      save_qc_png: control.saveQcPlots,
      save_report: control.saveReport,
      report_format: control.outputFormat,
      scc_background_method: control.sccBackgroundMethod,
      scc_background_k: control.sccBackgroundK,
      spectral_variant_som_nodes: control.spectralVariantSomNodes,
      spectral_variant_top_k: control.spectralVariantTopK,
      spectral_variant_cosine_threshold: control.spectralVariantCosineThreshold,
      spectral_variant_max_variants: control.spectralVariantMaxVariants,
      spectral_variant_min_events: control.spectralVariantMinEvents,
      spectreasy_weight_quantile: control.spectreasyWeightQuantile,
      autospectral_n_candidates: control.autospectralNCandidates,
      autospectral_n_spectral: control.autospectralNSpectral,
      autospectral_min_events: control.autospectralMinEvents,
      autospectral_refine: control.refine,
    };
  }
  if (action === "sample") {
    return {
      projectPath: settings.projectPath,
      sample_dir: sample.sampleDir,
      matrix_file: sample.matrixFile,
      detector_noise_file: sample.detectorNoiseFile,
      output_dir: normalizeCockpitOutputRoot(sample.outputDir),
      method: sample.method,
      rwls_max_iter: sample.rwlsMaxIter,
      n_threads: sample.nThreads,
      spectral_variant_top_k: sample.spectralVariantTopK,
      spectral_variant_min_abundance: sample.spectralVariantMinAbundance,
      spectral_variant_positive_fraction: sample.spectralVariantPositiveFraction,
      spectral_variant_min_improvement: sample.spectralVariantMinImprovement,
      spectral_variant_library_file: sample.spectralVariantLibraryFile,
      spectreasy_weight_quantile: sample.spectreasyWeightQuantile,
      estimate_af: sample.estimateAf,
      write_fcs: sample.writeFcs,
      save_report: sample.saveReport,
      report_format: sample.outputFormat,
      report_per_sample: sample.reportPerSample,
      save_qc_plots: sample.saveQcPlots,
      plot_n_events: sample.plotNEvents,
      chunk_size: sample.chunkSize,
      seed: sample.seed,
      return_type: sample.returnType,
    };
  }
  return {
    projectPath: settings.projectPath,
    fcs_file: af.fcsFile,
    save_name: af.saveName,
    save_overwrite: af.saveOverwrite,
    af_n_bands: af.afNBands,
    af_max_cells: af.afMaxCells,
    seed: af.seed,
  };
}
