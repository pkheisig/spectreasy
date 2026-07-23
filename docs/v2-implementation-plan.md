# Spectreasy v2 post-unmixing analysis implementation plan

## Outcome

Spectreasy v2 adds a project-scoped analysis workspace for raw or unmixed FCS
sample files. It follows the interaction model of a conventional FlowJo/OMIQ
workspace: a population hierarchy on the left, a user-configurable plot canvas
in the center, and contextual settings/statistics on the right. It does not
reuse the fixed control-gating plot sequence.

The workspace covers four connected jobs:

1. Hierarchical manual gating and population statistics.
2. Reproducible event downsampling and FCS/CSV export at any population.
3. Dimensional reduction, clustering, and trajectory inference.
4. Publication-ready plots, analysis provenance, and method citations.

All implementation belongs on the `v2` branch. It must not be merged into
`master` as part of this work.

## Local OMIQ evidence

The Erlangen OMIQ statistics export at
`Erlangen_HD_CMV_testrun/omiq.csv` records population ancestry as slash-delimited
paths, for example `cells/singlets/live/CD3_CD8/tet+`. The available
10,000-event OMIQ FCS export in this checkout records
`WRITTEN_BY=OMIQ (www.omiq.ai)` but contains no population-parent or gate-geometry
keywords. Therefore FCS alone cannot reconstruct the OMIQ workspace tree.

Spectreasy consequently keeps editable gate geometry and ancestry in
`.spectreasy/analysis-v2/workspace.json`, exports the selected events to FCS,
and writes the complete slash-delimited population path plus source, seed, and
writer provenance into that FCS. This follows the observed OMIQ export semantics
without pretending that gate geometry is embedded in the file.

## Scientific and interaction contract

### Data sources

- Discover the project's configured sample input folder.
- Discover unmixed FCS outputs under the selected output root.
- Permit an explicit project-contained FCS folder.
- Never infer raw versus unmixed from the number of parameters alone; label the
  selected source and preserve its path relative to the project.
- Keep SCC/control and sample sources separately selectable.

### Event identity and reproducibility

- Assign each event a stable source file ID and one-based source row index.
- Random downsampling is without replacement.
- Store a master seed and derive a stable per-file seed from the master seed and
  normalized file identity. Reordering or adding files must not alter existing
  samples.
- Persist source file size, modification time, event count, channel schema,
  selected markers, transform, target count, seed, method, package version, and
  output checksum.
- Downsampling is optional for trajectory analyses because rare transitional
  states can be lost.

### Gate hierarchy

Each gate is a project-persistent node with:

- stable ID, name, parent ID, and source scope;
- rectangle, polygon, or one-dimensional range geometry;
- x/y channels and per-axis transform;
- global or per-file scope;
- optional semantic role (`positive`, `negative`, `root`, `terminal`);
- creation/update timestamps and software version.

Gate membership is derived from the source data and ancestor chain. Cached
memberships may accelerate display, but geometry plus ancestry remain the source
of truth. Double-clicking a population opens it as the active parent. Deleting a
node requires explicit handling of descendants.

### Plot canvas

- Start with one configurable plot, not a preset scatter sequence.
- Add plots with a visible plus button; remove and reorder plots independently.
- Support scatter, density, contour, hexbin, and histogram presentations over
  the same event selection.
- Support marker/density/population colouring, backgating overlays, adjustable
  point size/opacity, and linear/logicle/asinh display transforms.
- Plot state belongs to the workspace, not to an FCS export.

### Staining index

Staining index uses two ordinary gates assigned the `positive` and `negative`
roles for a selected marker:

`SI = (MFI_positive - MFI_negative) / (2 * robust_SD_negative)`

The default MFI is the median. The default robust standard deviation is
`1.4826 * MAD`, with both choices recorded in the result. The calculation must
fail clearly if either role is missing, the gates come from incompatible parent
populations, or the negative population has insufficient events.

### Population annotation

Clustering, annotation, and gating are distinct:

- clustering discovers numerical groups;
- gating defines explicit geometric populations;
- annotation assigns a biological label from marker evidence.

The first implementation provides marker summaries and user-confirmed labels.
Rule-based suggestions must expose the exact marker rule and may never be
presented as verified biological identity.

## Method registry

Every method entry records its primary paper, implementation URL, package and
version, license, runtime capability, required inputs, seed support, and output
contract. A method is executable only when its backend is installed and has
passed a Spectreasy adapter test.

### Dimensional reduction

| Method | Backend | Initial status | Reference |
| --- | --- | --- | --- |
| PCA | base R | enabled | Pearson 1901 / standard linear decomposition |
| t-SNE | `Rtsne` | enabled | van der Maaten and Hinton 2008, JMLR |
| UMAP | `uwot` | enabled | McInnes et al. 2018, arXiv:1802.03426 |
| Diffusion map | `destiny` | enabled | Angerer et al. 2016, Bioinformatics, doi:10.1093/bioinformatics/btv715 |
| PHATE | `phate` Python | enabled when managed runtime is present | Moon et al. 2019, Nature Biotechnology, doi:10.1038/s41587-019-0336-3 |
| HSNE | `nptsne` Python | enabled when managed runtime is present | Pezzotti et al. 2016, Computer Graphics Forum, doi:10.1111/cgf.12878 |

### Clustering

| Method | Backend | Initial status | Reference |
| --- | --- | --- | --- |
| FlowSOM | `FlowSOM` | enabled | Van Gassen et al. 2015, Cytometry A, doi:10.1002/cyto.a.22625 |
| PhenoGraph | `Rphenograph` | enabled | Levine et al. 2015, Cell, doi:10.1016/j.cell.2015.05.047 |

### Trajectory and pseudotime

| Method | Backend | Initial status | Reference |
| --- | --- | --- | --- |
| DPT | `destiny` | enabled | Haghverdi et al. 2016, Nature Methods, doi:10.1038/nmeth.3971 |
| Slingshot | `slingshot` | enabled when installed | Street et al. 2018, BMC Genomics, doi:10.1186/s12864-018-4772-0 |
| TSCAN | `TSCAN` | enabled when installed | Ji and Ji 2016, Nucleic Acids Research, doi:10.1093/nar/gkw430 |
| Palantir | `palantir` Python | enabled when managed runtime is present | Setty et al. 2019, Nature Biotechnology, doi:10.1038/s41587-019-0068-4 |
| PAGA + DPT | `scanpy` Python | enabled when managed runtime is present | Wolf et al. 2019, Genome Biology, doi:10.1186/s13059-019-1663-x |
| Wanderlust | maintained Spectreasy Python implementation | enabled when managed runtime is present | Bendall et al. 2014, Cell, doi:10.1016/j.cell.2014.04.005 |
| Wishbone | maintained Spectreasy Python implementation | enabled when managed runtime is present | Setty et al. 2016, Nature Biotechnology, doi:10.1038/nbt.3569 |

The maintained Wanderlust and Wishbone module is an independent implementation
of the published algorithms. It does not copy the historical GPL Wishbone
repository into the MIT package. Synthetic linear and bifurcating trajectory
contracts, noisy-seed stability tests, event alignment, determinism, and
20,000-cell performance benchmarks gate the adapters.

## API and storage architecture

### Modules

- `inst/api/modules/analysis_workspace.R`: source discovery, state persistence,
  FCS access, gate evaluation, statistics, export, and method execution.
- `R/population_analysis.R`: exported R interface and the shared request
  boundary used by both code and GUI callers.
- `inst/gui/src/analysis/`: analysis workspace components, plot rendering,
  geometry tools, method registry display, and API client.
- `inst/python/`: maintained Python adapters and the independent
  Wanderlust/Wishbone trajectory implementation.
- `spectreasy_outputs/analysis/`: user exports and method results.
- `.spectreasy/analysis-v2/workspace.json`: editable gate tree, plots,
  annotations, and analysis settings.

### R code interface

The browser is an orchestration layer, not an independent scientific
implementation. The method, run, and annotation endpoints call the same package
backend used in R sessions:

```r
methods <- analysis_methods()

fit <- analyze_population(
  project_path = project,
  file = "samples/sample.fcs",
  markers = c("CD3-A", "CD19-A", "CD56-A"),
  clustering = "flowsom",
  reduction = "umap",
  cluster_settings = list(clusters = 12, xdim = 10, ydim = 10),
  reduction_settings = list(
    neighbors = 20,
    min_dist = 0.05,
    metric = "cosine"
  ),
  max_events = 20000,
  seed = 20260723
)

plot_population_analysis(
  fit,
  x = "UMAP 1",
  y = "UMAP 2",
  color_by = "CD3",
  palette = "sunset"
)

export_population_analysis(fit, "analysis.svg", color_by = "cluster")
export_population_analysis(fit, "analysis.csv")
```

`analysis_methods()` exposes the complete validated parameter schema, including
defaults, ranges, choices, prerequisites, and runtime availability.
`analyze_population()` accepts method-specific named setting lists, applies the
same validation as the GUI, and records resolved values in provenance and cache
keys. `load_population_analysis()`, `annotate_population()`,
`plot_population_analysis()`, `population_analysis_palettes()`, and
`export_population_analysis()` provide programmatic access to saved objects,
marker identities, square 2D plots, conditional Plotly 3D plots, and
CSV/RDS/PNG/SVG/PDF/HTML output.

### Endpoint families

- `GET /analysis/sources`
- `GET /analysis/workspace`
- `POST /analysis/workspace`
- `GET /analysis/events`
- `POST /analysis/statistics`
- `POST /analysis/export`
- `GET /analysis/methods`
- `POST /analysis/run`

Mutating routes reuse the existing local-session token protection. Every path is
resolved beneath the active project root.

### Common method output

```text
analysis_id
method_id
method_version
method_citation
source_fingerprint
population_id
markers
transform
seed
parameters
event_ids
coordinates
cluster_id
pseudotime
lineage_membership
root_event_id
terminal_state
warnings
runtime_seconds
```

Unavailable fields are explicit nulls, not silently invented values.
Each completed method run also produces a coordinate CSV, metadata JSON, PNG,
PDF, HTML, and SVG when the `svglite` device is installed. A separate hierarchy
statistics export contains per-file/per-population counts, percentages, marker
median/mean, and robust SD.

## Performance strategy

- Read FCS metadata without loading event matrices.
- Cache one selected source matrix and derived population masks on the backend.
- Return compact, bounded plot payloads; keep full membership and calculations
  in R.
- Use canvas for event clouds and SVG only for axes/gate overlays.
- Derive display samples independently from analytical samples.
- Lazy-load the analysis applet and optional method backends.
- Start independent API requests concurrently and cancel stale plot requests.

## Validation gates

### Contract tests

- Project containment for every path-taking endpoint.
- Stable per-file sampling under file reordering.
- Gate ancestry and rectangle/polygon/range membership.
- Correct population counts and percentages.
- Staining-index calculation against a hand-calculated fixture.
- FCS/CSV export schema, event count, keywords, and re-readability.
- Method availability reflects the real runtime.

### Scientific adapter tests

- DPT: reproduce ordering on a documented hematopoietic/reference fixture.
- Slingshot/TSCAN: reproduce a seeded bifurcating synthetic trajectory before
  enabling.
- Wanderlust/Wishbone: recover seeded linear and bifurcating trajectories,
  remain stable across noisy seeds, and preserve event alignment.
- Record package versions and tolerate only documented numerical variation.

### UI and system tests

- Add plot, change axes, draw gate, create child gate, navigate hierarchy.
- Assign positive/negative roles and calculate staining index.
- Backgate a descendant onto an ancestor plot.
- Persist/reload the workspace.
- Export selected population and inspect the resulting files.
- Run every available method; unavailable methods are disabled with a concrete
  setup action. Paper citations remain in provenance metadata, not the GUI.
- Tune method-specific advanced settings and verify parameter-dependent cache
  keys, cross-parameter validation, reset behavior, and disabled Run states.
- Add independent result plots, use 2D and 3D coordinate controls, and export
  PNG, SVG, interactive HTML, and the complete plotting table.
- Frontend lint, tests, production build; R focused tests and package check;
  headless Chromium smoke at desktop and compact widths.

## Delivery sequence

1. Persistent workspace schema and method registry.
2. Source discovery, compact event payloads, and gate evaluation.
3. OMIQ-style hierarchy plus configurable multi-plot canvas.
4. Statistics, staining index, reproducible downsampling, and exports.
5. Available reduction/clustering/DPT adapters.
6. Provenance metadata, dynamic method settings, and capability-gated runtimes.
7. Focused, full, and browser validation.
