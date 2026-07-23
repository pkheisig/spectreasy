#!/usr/bin/env python3
"""Maintained Python adapters for Spectreasy sample analysis.

The R API owns sampling, transformations, provenance, and cache keys. This
bridge only fits methods whose maintained reference implementation is Python
and writes deterministic, event-aligned numeric results back to the R process.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def _parameters(path: str | None) -> dict:
    if not path:
        return {}
    with open(path, encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("Advanced method settings must be a JSON object.")
    return value


def _matrix(path: str) -> np.ndarray:
    value = np.loadtxt(path, delimiter=",", ndmin=2)
    if value.ndim != 2 or value.shape[0] < 20 or value.shape[1] < 2:
        raise ValueError("The Python analysis adapter requires at least 20 events and two markers.")
    if not np.isfinite(value).all():
        raise ValueError("The Python analysis adapter received non-finite marker values.")
    return np.asarray(value, dtype=np.float64)


def _write_table(path: str, coordinates: np.ndarray, **columns: np.ndarray) -> None:
    coordinates = np.asarray(coordinates, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError("The analysis method did not return at least two coordinates.")
    if not np.isfinite(coordinates).all():
        raise ValueError("The analysis method returned non-finite coordinates.")
    frame = pd.DataFrame(
        coordinates,
        columns=[f"dimension_{index + 1}" for index in range(coordinates.shape[1])],
    )
    for name, values in columns.items():
        frame[name] = np.asarray(values)
    frame.to_csv(path, index=False)


def reduce_phate(args: argparse.Namespace) -> None:
    import joblib
    import phate

    data = _matrix(args.input)
    parameters = _parameters(args.parameters)
    operator = phate.PHATE(
        n_components=args.dimensions,
        knn=min(int(parameters.get("neighbors", args.neighbors)), data.shape[0] - 1),
        decay=int(parameters.get("decay", 40)),
        gamma=float(parameters.get("gamma", 1)),
        n_landmark=min(int(parameters.get("landmarks", 2000)), data.shape[0]),
        n_pca=min(int(parameters.get("pca_components", 100)), data.shape[1]),
        mds=str(parameters.get("mds", "metric")),
        knn_dist=str(parameters.get("distance", "euclidean")),
        random_state=args.seed,
        n_jobs=1,
        verbose=0,
    )
    coordinates = operator.fit_transform(data)
    joblib.dump(operator, args.model, compress=3)
    _write_table(args.output, coordinates)


def reduce_hsne(args: argparse.Namespace) -> None:
    import nptsne

    data = np.asarray(_matrix(args.input), dtype=np.float32)
    parameters = _parameters(args.parameters)
    # HSNE is a hierarchical landmark method. The overview embeds the highest
    # scale; every original event is mapped to its representative landmark so
    # event IDs, marker values, and cached cluster labels remain aligned.
    requested_scales = max(2, min(5, int(parameters.get("scales", 3))))
    np.random.seed(args.seed)
    # nptsne always embeds the highest scale. Small populations can collapse
    # that scale to a single landmark, for which t-SNE is undefined and emits
    # NaN coordinates. Build the deepest requested hierarchy whose overview
    # still contains at least three landmarks.
    hierarchy = None
    scales = requested_scales
    while scales >= 2:
        candidate = nptsne.HSne(False)
        if not candidate.create_hsne(data, scales):
            raise RuntimeError("nptsne could not construct the HSNE hierarchy.")
        top_points = int(candidate.get_scale(scales - 1).num_points)
        if top_points >= 3:
            hierarchy = candidate
            break
        scales -= 1
    if hierarchy is None:
        raise RuntimeError("HSNE requires enough events to create at least three overview landmarks.")
    hierarchy.save(args.model)
    analysis = nptsne.hsne_analysis.Analysis(
        hierarchy,
        nptsne.hsne_analysis.EmbedderType.CPU,
    )
    latest_finite_embedding = None
    for _ in range(int(parameters.get("iterations", 500))):
        analysis.do_iteration()
        candidate = np.asarray(analysis.embedding, dtype=np.float64)
        if candidate.ndim == 2 and candidate.shape[1] >= 2 and np.isfinite(candidate).all():
            latest_finite_embedding = candidate.copy()
    landmark_coordinates = np.asarray(analysis.embedding, dtype=np.float64)
    if not np.isfinite(landmark_coordinates).all():
        if latest_finite_embedding is None:
            raise RuntimeError("nptsne did not produce a finite HSNE landmark embedding.")
        landmark_coordinates = latest_finite_embedding
    coordinates = np.zeros((data.shape[0], 2), dtype=np.float64)
    assigned = np.zeros(data.shape[0], dtype=bool)
    for landmark in range(analysis.number_of_points):
        mask = np.asarray(
            analysis.get_fast_area_of_influence([landmark]),
            dtype=np.float32,
        ) > 0
        coordinates[mask] = landmark_coordinates[landmark]
        assigned |= mask
    if not assigned.all():
        # This should be rare and reflects thresholded areas of influence.
        # Fall back to the closest original landmark for only those events.
        landmark_rows = np.asarray(analysis.landmark_orig_indexes, dtype=int)
        missing = np.flatnonzero(~assigned)
        distances = (
            (data[missing, None, :] - data[landmark_rows][None, :, :]) ** 2
        ).sum(axis=2)
        coordinates[missing] = landmark_coordinates[np.argmin(distances, axis=1)]
    _write_table(args.output, coordinates)


def trajectory_palantir(args: argparse.Namespace) -> None:
    import palantir.core
    import palantir.utils

    data = _matrix(args.input)
    parameters = _parameters(args.parameters)
    index = pd.Index([str(value) for value in range(data.shape[0])], name="event")
    frame = pd.DataFrame(data, index=index)
    components = min(
        max(int(parameters.get("diffusion_components", 10)), args.dimensions + 1),
        data.shape[0] - 2,
    )
    diffusion = palantir.utils.run_diffusion_maps(
        frame,
        n_components=components,
        knn=min(int(parameters.get("neighbors", args.neighbors)), data.shape[0] - 1),
        alpha=float(parameters.get("diffusion_alpha", 0)),
        seed=args.seed,
    )
    multiscale = palantir.utils.determine_multiscale_space(diffusion)
    result = palantir.core.run_palantir(
        multiscale,
        early_cell=str(args.root_index),
        knn=min(int(parameters.get("neighbors", args.neighbors)), data.shape[0] - 1),
        num_waypoints=min(
            int(parameters.get("waypoints", 1200)),
            max(20, data.shape[0]),
        ),
        n_jobs=1,
        scale_components=bool(parameters.get("scale_components", True)),
        use_early_cell_as_start=True,
        max_iterations=int(parameters.get("iterations", 25)),
        seed=args.seed,
    )
    eigenvectors = np.asarray(diffusion["EigenVectors"], dtype=np.float64)
    coordinates = eigenvectors[:, : args.dimensions]
    pseudotime = np.asarray(result.pseudotime.reindex(index), dtype=np.float64)
    with open(args.model, "wb") as handle:
        pickle.dump(
            {"diffusion": diffusion, "multiscale": multiscale, "result": result},
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    _write_table(args.output, coordinates, pseudotime=pseudotime)


def trajectory_paga_dpt(args: argparse.Namespace) -> None:
    import anndata
    import scanpy as sc
    import scipy.sparse
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.spatial.distance import cdist

    data = _matrix(args.input)
    parameters = _parameters(args.parameters)
    clusters = np.loadtxt(args.clusters, delimiter=",", dtype=str, ndmin=1)
    if clusters.shape[0] != data.shape[0]:
        raise ValueError("PAGA+DPT requires one cached cluster label per event.")
    adata = anndata.AnnData(data)
    adata.obs_names = [str(value) for value in range(data.shape[0])]
    adata.obs["spectreasy_cluster"] = pd.Categorical(clusters)
    sc.pp.neighbors(
        adata,
        n_neighbors=min(int(parameters.get("neighbors", args.neighbors)), data.shape[0] - 1),
        metric=str(parameters.get("metric", "euclidean")),
        random_state=args.seed,
    )
    try:
        sc.tl.paga(
            adata,
            groups="spectreasy_cluster",
            model=str(parameters.get("paga_model", "v1.2")),
        )
    except ValueError as error:
        # Scanpy/igraph currently fails when the kNN graph contains no
        # inter-cluster edges. Preserve a valid, explicit PAGA graph by
        # connecting cluster centroids with their minimum spanning tree.
        # DPT itself continues to use the event-level diffusion graph.
        message = str(error).lower()
        if "index arrays" not in message and "adjacency" not in message:
            raise
        categories = list(adata.obs["spectreasy_cluster"].cat.categories)
        codes = np.asarray(adata.obs["spectreasy_cluster"].cat.codes)
        if len(categories) < 2:
            raise ValueError("PAGA+DPT requires clustering with at least two populated clusters.") from error
        centroids = np.vstack([data[codes == index].mean(axis=0) for index in range(len(categories))])
        distances = cdist(centroids, centroids, metric="euclidean")
        positive = distances[np.isfinite(distances) & (distances > 0)]
        scale = float(np.median(positive)) if positive.size else 1.0
        connectivities = np.exp(-distances / max(scale, np.finfo(float).eps))
        np.fill_diagonal(connectivities, 0.0)
        tree_cost = np.where(connectivities > 0, 1.0 / connectivities, 0.0)
        tree = minimum_spanning_tree(scipy.sparse.csr_matrix(tree_cost))
        tree = tree + tree.T
        adata.uns["paga"] = {
            "connectivities": scipy.sparse.csr_matrix(connectivities),
            "connectivities_tree": scipy.sparse.csr_matrix(tree),
            "groups": "spectreasy_cluster",
        }
        adata.uns["spectreasy_paga_fallback"] = "centroid_mst_no_intercluster_edges"
    diffusion_components = min(
        max(int(parameters.get("diffusion_components", 10)), args.dimensions + 1),
        data.shape[0] - 1,
    )
    sc.tl.diffmap(adata, n_comps=diffusion_components)
    adata.uns["iroot"] = int(args.root_index)
    sc.tl.dpt(
        adata,
        n_dcs=diffusion_components,
        n_branchings=int(parameters.get("branchings", 0)),
        min_group_size=float(parameters.get("minimum_group_size", 0.01)),
    )
    diffusion = np.asarray(adata.obsm["X_diffmap"], dtype=np.float64)
    start = 1 if diffusion.shape[1] > args.dimensions else 0
    coordinates = diffusion[:, start : start + args.dimensions]
    if coordinates.shape[1] < args.dimensions:
        coordinates = diffusion[:, : args.dimensions]
    pseudotime = np.asarray(adata.obs["dpt_pseudotime"], dtype=np.float64)
    adata.write_h5ad(args.model)
    _write_table(args.output, coordinates, pseudotime=pseudotime)


def trajectory_wanderlust(args: argparse.Namespace) -> None:
    from spectreasy_trajectory import run_wanderlust

    data = _matrix(args.input)
    parameters = _parameters(args.parameters)
    result = run_wanderlust(
        data,
        args.root_index,
        neighbors=int(parameters.get("neighbors", 15)),
        candidate_neighbors=int(parameters.get("candidate_neighbors", 20)),
        graphs=int(parameters.get("graphs", 5)),
        waypoints=min(int(parameters.get("waypoints", 64)), data.shape[0]),
        iterations=int(parameters.get("iterations", 15)),
        weighting=str(parameters.get("weighting", "exponential")),
        dimensions=args.dimensions,
        seed=args.seed,
    )
    with open(args.model, "wb") as handle:
        pickle.dump(
            {
                "method": "wanderlust",
                "parameters": parameters,
                "pseudotime": result.pseudotime,
                "coordinates": result.coordinates,
                "waypoints": result.waypoints,
            },
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    _write_table(args.output, result.coordinates, pseudotime=result.pseudotime)


def trajectory_wishbone(args: argparse.Namespace) -> None:
    from spectreasy_trajectory import run_wishbone

    data = _matrix(args.input)
    parameters = _parameters(args.parameters)
    result = run_wishbone(
        data,
        args.root_index,
        neighbors=int(parameters.get("neighbors", 15)),
        candidate_neighbors=int(parameters.get("candidate_neighbors", 20)),
        graphs=int(parameters.get("graphs", 5)),
        waypoints=min(int(parameters.get("waypoints", 96)), data.shape[0]),
        iterations=int(parameters.get("iterations", 15)),
        diffusion_components=min(
            int(parameters.get("diffusion_components", 10)),
            data.shape[0] - 2,
        ),
        branch_confidence=float(parameters.get("branch_confidence", 0.65)),
        dimensions=args.dimensions,
        seed=args.seed,
    )
    with open(args.model, "wb") as handle:
        pickle.dump(
            {
                "method": "wishbone",
                "parameters": parameters,
                "pseudotime": result.pseudotime,
                "coordinates": result.coordinates,
                "branches": result.branches,
                "branch_probabilities": result.branch_probabilities,
                "branch_point": result.branch_point,
                "waypoints": result.waypoints,
            },
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    _write_table(
        args.output,
        result.coordinates,
        pseudotime=result.pseudotime,
        trajectory_branch=result.branches,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("operation", choices=("reduce", "trajectory", "probe"))
    parser.add_argument("--method", default="")
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--model")
    parser.add_argument("--clusters")
    parser.add_argument("--dimensions", type=int, default=2)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--root-index", type=int, default=0)
    parser.add_argument("--parameters")
    args = parser.parse_args()

    if args.operation == "probe":
        import importlib.metadata
        import importlib.util

        modules = ("phate", "nptsne", "scanpy", "palantir")
        status = {
            module: {
                "available": importlib.util.find_spec(module) is not None,
                "version": (
                    importlib.metadata.version(module)
                    if importlib.util.find_spec(module) is not None
                    else ""
                ),
            }
            for module in modules
        }
        built_in_dependencies = ("numpy", "pandas", "scipy", "sklearn")
        status["spectreasy_builtin"] = {
            "available": all(
                importlib.util.find_spec(module) is not None
                for module in built_in_dependencies
            ),
            "version": "1.0.0",
        }
        print(
            json.dumps(status)
        )
        return

    if args.operation == "reduce" and args.method == "phate":
        reduce_phate(args)
    elif args.operation == "reduce" and args.method == "hsne":
        reduce_hsne(args)
    elif args.operation == "trajectory" and args.method == "palantir":
        trajectory_palantir(args)
    elif args.operation == "trajectory" and args.method == "paga-dpt":
        trajectory_paga_dpt(args)
    elif args.operation == "trajectory" and args.method == "wanderlust":
        trajectory_wanderlust(args)
    elif args.operation == "trajectory" and args.method == "wishbone":
        trajectory_wishbone(args)
    else:
        raise ValueError(f"Unsupported Python analysis adapter: {args.operation}/{args.method}")


if __name__ == "__main__":
    main()
