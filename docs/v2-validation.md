# v2 population-analysis validation

This file records the maintained validation contract and the latest local
development measurements. Timings are hardware- and runtime-specific and are
not product guarantees.

## Scientific contract

- Every enabled reduction adapter runs on the same seeded 180-cell continuum,
  returns finite coordinates, and preserves continuum distances. HSNE uses a
  separate, lower global-distance threshold because it optimizes a landmark
  hierarchy.
- FlowSOM and PhenoGraph recover known separated synthetic populations with
  adjusted Rand index above 0.8.
- DPT, Slingshot, TSCAN, Palantir, PAGA+DPT, Wanderlust, and Wishbone order a
  known seeded lineage from the selected root with Spearman correlation above
  0.6. The selected root must remain within the first 12% of normalized
  pseudotime.
- Every enabled adapter also satisfies finite-output, dimensionality,
  prerequisite, same-seed reproducibility, artifact-cache, and provenance
  checks on a second small fixture.

These tests live in `tests/testthat/test-analysis-workspace-v2.R`; maintained
Wanderlust and Wishbone implementations also have independent noisy-lineage
tests in `tests/python/test_spectreasy_trajectory.py`.

## Performance harness

Run:

```sh
Rscript tests/performance/analysis_v2_full_benchmark.R /tmp/analysis-v2-full full
```

The full profile runs all 15 methods cold and warm at 500, 2,000, and 5,000
events. It additionally runs PCA, UMAP, diffusion map, FlowSOM, and PhenoGraph
at 20,000 events. Each run verifies row count, finite coordinates, normalized
pseudotime when applicable, and fitted-object reuse. RSS values are diagnostic
before/after deltas rather than peak-memory guarantees.

Latest local measurements on 2026-07-23, using R 4.6.1 on Apple silicon and 12
markers:

| Events | Cold method | Seconds | Events/second |
| ---: | --- | ---: | ---: |
| 20,000 | PCA | 1.783 | 11,217 |
| 20,000 | UMAP | 23.235 | 861 |
| 20,000 | Diffusion map | 11.658 | 1,716 |
| 20,000 | FlowSOM | 4.612 | 4,337 |
| 20,000 | PhenoGraph | 11.125 | 1,798 |

Across 50 configurations and 100 executions, all scientific-output checks
passed and every warm run reused the requested fitted object. Median cold times
were 1.783 seconds for reductions, 2.463 seconds for clustering, and 2.438
seconds for trajectories. At 5,000 events the trajectory cold times were DPT
1.255 s, Slingshot 4.070 s, TSCAN 0.644 s, Palantir 5.152 s, PAGA+DPT 8.755 s,
Wanderlust 3.477 s, and Wishbone 3.594 s.

The first full run identified two scale-specific failures that the small method
contracts did not: covariance-aware Slingshot distances could create an NA
cluster graph, and Scanpy PAGA could fail when the event kNN graph contained no
inter-cluster edges. The maintained defaults now use Slingshot's supported
Euclidean-center distance, and PAGA records a centroid minimum-spanning-tree
fallback while DPT continues to use its event-level diffusion graph. The
100-execution matrix above is the clean rerun after both corrections.

## Real PKH acceptance

Run against a local dataset without writing to it:

```sh
Rscript tests/acceptance/analysis_v2_real_pkh.R DATASET_ROOT /tmp/pkh-acceptance
```

The 2026-07-23 acceptance run used
`Erlangen_P13_test_culture_repeat/samples_full/1.fcs` (301,304 events, 71
parameters), copied it to an isolated scratch project, and selected eight real
marker channels. Rectangle, ellipse, polygon, and histogram-range gates all
returned nonempty, proper subsets; the POS/NEG sibling gates produced a finite
robust-MAD staining index; FCS and CSV export checksums and workspace
round-tripping passed. All 15 available methods returned 600 aligned real
events with finite coordinates and finite pseudotime where applicable.

The corresponding authenticated, read-only OMIQ inspection was limited to
datasets whose names contained `PKH`. The OMIQ pregating dataset's four sample
event counts—301,304; 220,912; 123,328; and 302,288—exactly matched the local
`samples_full` files. No OMIQ task, gate, plot, or dataset was saved or changed.

## Browser and native UI acceptance

`tests/browser/analysis_v2_smoke.py`, `analysis_v2_compact.py`, and
`setup_smoke.py` verify the Other tools entry, fixed square plots, searchable
external axis selectors, plot deletion, undo/redo, complete workspace JSON
round-trip, clustering/map cache reuse, advanced settings, 2D/3D result plots,
marker palettes, cell identities, plot/table export, compact layout, and the
cross-platform Setup instructions. At compact widths the analysis dialog is
full-screen, its file/marker and pipeline sections stack without overlap, and
the result remains a fixed square instead of stretching to fill the viewport.

The macOS export-folder button was also exercised in a real Chrome session. It
opened Finder's native “Choose a Folder” dialog at the active scratch project,
selected `spectreasy_outputs/analysis`, returned that project-relative path to
the GUI, and left the original default export field editable. The server also
rejects picker results outside the active project.

The package workflow builds and installs the same source tarball on Linux,
macOS, and Windows. No platform-specific application bundle is required or
published; Python-backed analysis remains optional and uses the managed runtime
installed from the Setup page or `install_analysis_dependencies()`.
