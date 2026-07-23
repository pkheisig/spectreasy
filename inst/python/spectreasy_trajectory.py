"""Independent maintained Wanderlust and Wishbone implementations.

This module is an original implementation based on the algorithms described in:

* Bendall et al. (2014), Cell, doi:10.1016/j.cell.2014.04.005
* Setty et al. (2016), Nature Biotechnology, doi:10.1038/nbt.3569

It does not contain code from the historical GPL-2 Wishbone repository.  The
implementation keeps the published core semantics—bootstrapped kNN graphs,
waypoint-perspective consensus, diffusion geometry, and bifurcation inference—
while using current NumPy, SciPy, and scikit-learn APIs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors


def _normalize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(values)
    if not finite.any():
        raise ValueError("Trajectory calculation returned no finite values.")
    if not finite.all():
        values[~finite] = np.nanmedian(values[finite])
    low, high = np.min(values), np.max(values)
    if high <= low:
        return np.zeros(values.shape, dtype=np.float64)
    return (values - low) / (high - low)


def _standardize(data: np.ndarray) -> np.ndarray:
    data = np.asarray(data, dtype=np.float64)
    center = np.median(data, axis=0)
    scale = np.median(np.abs(data - center), axis=0) * 1.4826
    fallback = np.std(data, axis=0)
    scale[(~np.isfinite(scale)) | (scale <= np.finfo(float).eps)] = fallback[
        (~np.isfinite(scale)) | (scale <= np.finfo(float).eps)
    ]
    scale[(~np.isfinite(scale)) | (scale <= np.finfo(float).eps)] = 1.0
    return np.clip((data - center) / scale, -8, 8)


def _connect_components(graph: csr_matrix, data: np.ndarray) -> csr_matrix:
    count, labels = connected_components(graph, directed=False)
    if count <= 1:
        return graph
    graph = graph.tolil(copy=True)
    # Connect adjacent component centroids, then use the closest observed pair
    # for each bridge. This avoids inventing a zero-cost edge.
    components = [np.flatnonzero(labels == value) for value in range(count)]
    centroids = np.vstack([data[index].mean(axis=0) for index in components])
    centroid_order = np.argsort(np.linalg.norm(centroids - centroids[0], axis=1))
    for left_id, right_id in zip(centroid_order[:-1], centroid_order[1:]):
        left = components[int(left_id)]
        right = components[int(right_id)]
        search = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(data[right])
        distances, matches = search.kneighbors(data[left])
        best = int(np.argmin(distances[:, 0]))
        source = int(left[best])
        target = int(right[int(matches[best, 0])])
        weight = max(float(distances[best, 0]), np.finfo(float).eps)
        graph[source, target] = weight
        graph[target, source] = weight
    return graph.tocsr()


def _candidate_neighbors(
    data: np.ndarray,
    candidate_neighbors: int,
) -> tuple[np.ndarray, np.ndarray]:
    n = data.shape[0]
    candidate_neighbors = min(max(3, candidate_neighbors), n - 1)
    search = NearestNeighbors(
        n_neighbors=candidate_neighbors + 1,
        metric="euclidean",
        n_jobs=1,
    ).fit(data)
    distances, indices = search.kneighbors(data)
    return distances[:, 1:], indices[:, 1:]


def _neighbor_graph(
    data: np.ndarray,
    candidate_distances: np.ndarray,
    candidate_indices: np.ndarray,
    retained_neighbors: int,
    rng: np.random.Generator,
    bootstrap: bool,
) -> csr_matrix:
    n, candidate_neighbors = candidate_indices.shape
    retained_neighbors = min(max(2, retained_neighbors), candidate_neighbors)
    rows: list[int] = []
    cols: list[int] = []
    weights: list[float] = []
    for row in range(n):
        if bootstrap and retained_neighbors < candidate_neighbors:
            selected = np.sort(
                rng.choice(candidate_neighbors, retained_neighbors, replace=False)
            )
        else:
            selected = np.arange(retained_neighbors)
        rows.extend([row] * selected.size)
        cols.extend(candidate_indices[row, selected].tolist())
        weights.extend(np.maximum(candidate_distances[row, selected], np.finfo(float).eps).tolist())
    directed = coo_matrix((weights, (rows, cols)), shape=(n, n)).tocsr()
    graph = directed.maximum(directed.T)
    graph.setdiag(0)
    graph.eliminate_zeros()
    return _connect_components(graph, data)


def _waypoints(graph: csr_matrix, root: int, count: int) -> tuple[np.ndarray, np.ndarray]:
    n = graph.shape[0]
    count = min(max(2, count), n)
    selected = [int(root)]
    minimum = np.asarray(dijkstra(graph, directed=False, indices=root), dtype=float)
    minimum[~np.isfinite(minimum)] = np.nanmax(minimum[np.isfinite(minimum)])
    while len(selected) < count:
        candidate = int(np.argmax(minimum))
        if candidate in selected:
            break
        selected.append(candidate)
        distance = np.asarray(
            dijkstra(graph, directed=False, indices=candidate),
            dtype=float,
        )
        finite = np.isfinite(distance)
        distance[~finite] = np.nanmax(distance[finite]) if finite.any() else minimum.max()
        minimum = np.minimum(minimum, distance)
    landmarks = np.asarray(selected, dtype=int)
    distances = np.asarray(
        dijkstra(graph, directed=False, indices=landmarks),
        dtype=np.float64,
    )
    finite = np.isfinite(distances)
    replacement = np.nanmax(distances[finite]) * 1.05 if finite.any() else 1.0
    distances[~finite] = replacement
    return landmarks, distances


def _weights(distances: np.ndarray, scheme: str) -> np.ndarray:
    if scheme == "uniform":
        weights = np.ones_like(distances)
    elif scheme == "linear":
        weights = np.maximum(np.max(distances, axis=0, keepdims=True) - distances, 0)
    else:
        sigma = max(float(np.mean(np.std(distances, axis=0))) * 3, np.finfo(float).eps)
        weights = np.exp(-0.5 * (distances / sigma) ** 2)
    denominator = weights.sum(axis=0, keepdims=True)
    denominator[denominator <= 0] = 1
    return weights / denominator


def _consensus(
    landmarks: np.ndarray,
    distances: np.ndarray,
    root_row: int,
    iterations: int,
    weighting: str,
) -> np.ndarray:
    trajectory = _normalize(distances[root_row])
    weights = _weights(distances, weighting)
    for _ in range(max(2, iterations)):
        perspectives = np.empty_like(distances)
        for row, landmark in enumerate(landmarks):
            position = trajectory[landmark]
            direction = np.where(trajectory < position, -1.0, 1.0)
            perspectives[row] = position + direction * distances[row]
        updated = _normalize(np.sum(perspectives * weights, axis=0))
        correlation = np.corrcoef(trajectory, updated)[0, 1]
        trajectory = updated
        if np.isfinite(correlation) and correlation >= 0.9999:
            break
    return trajectory


def _diffusion_embedding(
    data: np.ndarray,
    neighbors: int,
    components: int,
) -> np.ndarray:
    n = data.shape[0]
    neighbors = min(max(3, neighbors), n - 1)
    search = NearestNeighbors(n_neighbors=neighbors + 1, n_jobs=1).fit(data)
    distances, indices = search.kneighbors(data)
    distances, indices = distances[:, 1:], indices[:, 1:]
    local_scale = np.maximum(distances[:, -1], np.finfo(float).eps)
    rows = np.repeat(np.arange(n), neighbors)
    cols = indices.reshape(-1)
    denominator = local_scale[rows] * local_scale[cols]
    affinity = np.exp(-(distances.reshape(-1) ** 2) / np.maximum(denominator, np.finfo(float).eps))
    matrix = coo_matrix((affinity, (rows, cols)), shape=(n, n)).tocsr()
    matrix = matrix.maximum(matrix.T)
    degree = np.asarray(matrix.sum(axis=1)).ravel()
    inverse = np.zeros_like(degree)
    valid = degree > 0
    inverse[valid] = 1 / np.sqrt(degree[valid])
    normalized = diags(inverse) @ matrix @ diags(inverse)
    requested = min(max(components + 1, 3), n - 1)
    values, vectors = eigsh(normalized, k=requested, which="LA")
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    embedding = vectors[:, 1:] * values[None, 1:]
    if embedding.shape[1] < components:
        embedding = np.pad(embedding, ((0, 0), (0, components - embedding.shape[1])))
    return np.asarray(embedding[:, :components], dtype=np.float64)


@dataclass
class WanderlustResult:
    pseudotime: np.ndarray
    coordinates: np.ndarray
    waypoints: list[np.ndarray]


def run_wanderlust(
    data: np.ndarray,
    root: int,
    *,
    neighbors: int = 15,
    candidate_neighbors: int = 20,
    graphs: int = 5,
    waypoints: int = 100,
    iterations: int = 15,
    weighting: str = "exponential",
    dimensions: int = 3,
    seed: int = 1,
    coordinate_embedding: np.ndarray | None = None,
) -> WanderlustResult:
    data = _standardize(data)
    if not 0 <= root < data.shape[0]:
        raise ValueError("Wanderlust root index is outside the event matrix.")
    rng = np.random.default_rng(seed)
    candidate_neighbors = min(
        max(candidate_neighbors, neighbors + 1),
        data.shape[0] - 1,
    )
    candidate_distances, candidate_indices = _candidate_neighbors(
        data,
        candidate_neighbors,
    )
    trajectories: list[np.ndarray] = []
    waypoint_sets: list[np.ndarray] = []
    for _ in range(max(1, graphs)):
        graph = _neighbor_graph(
            data,
            candidate_distances=candidate_distances,
            candidate_indices=candidate_indices,
            retained_neighbors=neighbors,
            rng=rng,
            bootstrap=True,
        )
        landmarks, distances = _waypoints(graph, root, waypoints)
        root_row = int(np.flatnonzero(landmarks == root)[0])
        trajectories.append(
            _consensus(landmarks, distances, root_row, iterations, weighting)
        )
        waypoint_sets.append(landmarks)
    pseudotime = _normalize(np.mean(np.vstack(trajectories), axis=0))
    coordinates = (
        np.asarray(coordinate_embedding[:, :dimensions], dtype=np.float64)
        if coordinate_embedding is not None
        else _diffusion_embedding(data, candidate_neighbors, dimensions)
    )
    # Orient the displayed diffusion geometry consistently with the selected root.
    for column in range(coordinates.shape[1]):
        if np.corrcoef(coordinates[:, column], pseudotime)[0, 1] < 0:
            coordinates[:, column] *= -1
    return WanderlustResult(pseudotime, coordinates, waypoint_sets)


@dataclass
class WishboneResult:
    pseudotime: np.ndarray
    coordinates: np.ndarray
    branches: np.ndarray
    branch_probabilities: np.ndarray
    branch_point: float
    waypoints: list[np.ndarray]


def run_wishbone(
    data: np.ndarray,
    root: int,
    *,
    neighbors: int = 15,
    candidate_neighbors: int = 20,
    graphs: int = 5,
    waypoints: int = 150,
    iterations: int = 15,
    diffusion_components: int = 10,
    branch_confidence: float = 0.65,
    dimensions: int = 3,
    seed: int = 1,
) -> WishboneResult:
    standardized = _standardize(data)
    diffusion = _diffusion_embedding(
        standardized,
        candidate_neighbors,
        max(diffusion_components, dimensions),
    )
    wanderlust = run_wanderlust(
        diffusion[:, :diffusion_components],
        root,
        neighbors=neighbors,
        candidate_neighbors=candidate_neighbors,
        graphs=graphs,
        waypoints=waypoints,
        iterations=iterations,
        weighting="exponential",
        dimensions=dimensions,
        seed=seed,
        coordinate_embedding=diffusion,
    )
    pseudotime = wanderlust.pseudotime
    # Use a deliberately broad post-root candidate set. Restricting terminal
    # discovery to only the highest pseudotime quartile can omit a slower fate
    # when the two branches progress at different graph distances.
    late = np.flatnonzero(pseudotime >= np.quantile(pseudotime, 0.45))
    if late.size < 6:
        late = np.argsort(pseudotime)[-max(6, data.shape[0] // 4) :]
    terminal = KMeans(n_clusters=2, random_state=seed, n_init=20).fit(
        standardized[late]
    )
    centers = terminal.cluster_centers_
    distances = np.linalg.norm(
        standardized[:, None, :] - centers[None, :, :],
        axis=2,
    )
    scale = max(float(np.median(np.abs(distances[:, 0] - distances[:, 1]))), 0.05)
    logits = -distances / scale
    logits -= logits.max(axis=1, keepdims=True)
    probabilities = np.exp(logits)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    confidence = probabilities.max(axis=1)
    assignment = probabilities.argmax(axis=1)

    # The first sustained interval with confident representation of both fates
    # marks the bifurcation; earlier uncertain cells remain trunk cells.
    order = np.argsort(pseudotime)
    window = min(max(20, data.shape[0] // 10), data.shape[0])
    branch_point = float(np.quantile(pseudotime, 0.45))
    for offset in range(0, max(1, data.shape[0] - window + 1), max(1, window // 5)):
        selected = order[offset : offset + window]
        counts = np.bincount(assignment[selected], minlength=2) / max(1, selected.size)
        if (
            counts.min() >= 0.2
            and np.median(confidence[selected]) >= branch_confidence
            and np.median(pseudotime[selected]) >= 0.15
        ):
            branch_point = float(np.quantile(pseudotime[selected], 0.1))
            break
    branches = assignment.astype(int) + 2
    branches[(pseudotime < branch_point) | (confidence < branch_confidence)] = 1
    return WishboneResult(
        pseudotime=pseudotime,
        coordinates=diffusion[:, :dimensions],
        branches=branches,
        branch_probabilities=probabilities,
        branch_point=branch_point,
        waypoints=wanderlust.waypoints,
    )
