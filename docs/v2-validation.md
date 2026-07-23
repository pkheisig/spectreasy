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
Rscript tests/performance/analysis_v2_benchmark.R 20000 /tmp/analysis-v2-20k.csv
Rscript tests/performance/analysis_v2_benchmark.R 100000 /tmp/analysis-v2-100k.csv
```

Latest local measurements on 2026-07-23, using 12 markers:

| Events | Pipeline | Seconds | Events/second | Cache note |
| ---: | --- | ---: | ---: | --- |
| 20,000 | PCA | 9.667 | 2,069 | new embedding |
| 20,000 | FlowSOM → PCA | 5.085 | 3,933 | PCA reused; timing is added clustering work |
| 20,000 | UMAP | 26.021 | 769 | new embedding |
| 100,000 | PCA | 15.558 | 6,428 | new embedding |
| 100,000 | FlowSOM → PCA | 15.405 | 6,491 | PCA reused; timing is added clustering work |

The benchmark asserts that every requested event is returned and deliberately
records cache reuse rather than conflating it with a cold run.
