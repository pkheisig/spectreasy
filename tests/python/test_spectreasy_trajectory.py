import unittest

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

from spectreasy_trajectory import run_wanderlust, run_wishbone


class SpectreasyTrajectoryTests(unittest.TestCase):
    def test_wanderlust_recovers_seeded_linear_progression_deterministically(self):
        rng = np.random.default_rng(731)
        truth = np.linspace(0, 1, 240)
        data = np.column_stack(
            [truth, np.sin(2 * truth), truth**2, np.cos(truth)]
        ) + rng.normal(0, 0.025, (truth.size, 4))
        first = run_wanderlust(
            data,
            0,
            graphs=3,
            waypoints=36,
            iterations=10,
            seed=731,
        )
        second = run_wanderlust(
            data,
            0,
            graphs=3,
            waypoints=36,
            iterations=10,
            seed=731,
        )
        self.assertGreater(spearmanr(truth, first.pseudotime).statistic, 0.95)
        self.assertTrue(np.isfinite(first.coordinates).all())
        np.testing.assert_allclose(first.pseudotime, second.pseudotime)

    def test_wishbone_recovers_two_terminal_fates(self):
        rng = np.random.default_rng(917)
        trunk_time = np.linspace(0, 0.45, 100)
        branch_time = np.linspace(0.45, 1, 100)
        trunk = np.column_stack(
            [trunk_time, np.zeros(100), trunk_time * 0.2, np.zeros(100)]
        )
        positive = np.column_stack(
            [
                branch_time,
                (branch_time - 0.45) * 1.4,
                branch_time * 0.2,
                (branch_time - 0.45) * 0.3,
            ]
        )
        negative = np.column_stack(
            [
                branch_time,
                -(branch_time - 0.45) * 1.4,
                branch_time * 0.2,
                -(branch_time - 0.45) * 0.3,
            ]
        )
        data = np.vstack([trunk, positive, negative])
        data += rng.normal(0, 0.018, data.shape)
        result = run_wishbone(
            data,
            0,
            graphs=3,
            waypoints=50,
            iterations=10,
            diffusion_components=6,
            branch_confidence=0.55,
            seed=917,
        )
        truth_time = np.concatenate([trunk_time, branch_time, branch_time])
        terminal_truth = np.concatenate([np.zeros(100), np.ones(100)])
        self.assertGreater(spearmanr(truth_time, result.pseudotime).statistic, 0.85)
        self.assertGreater(
            adjusted_rand_score(terminal_truth, result.branches[100:]),
            0.65,
        )
        self.assertEqual(set(np.unique(result.branches)), {1, 2, 3})

    def test_wishbone_branch_recovery_is_stable_across_noisy_seeds(self):
        trunk_time = np.linspace(0, 0.45, 100)
        branch_time = np.linspace(0.45, 1, 100)
        base = np.vstack(
            [
                np.column_stack([trunk_time, np.zeros(100), trunk_time * 0.2, np.zeros(100)]),
                np.column_stack([branch_time, (branch_time - 0.45) * 1.4, branch_time * 0.2, (branch_time - 0.45) * 0.3]),
                np.column_stack([branch_time, -(branch_time - 0.45) * 1.4, branch_time * 0.2, -(branch_time - 0.45) * 0.3]),
            ]
        )
        terminal_truth = np.concatenate([np.zeros(100), np.ones(100)])
        for seed in (3, 7, 10):
            rng = np.random.default_rng(seed)
            data = base + rng.normal(0, 0.025, base.shape)
            result = run_wishbone(
                data,
                0,
                graphs=3,
                waypoints=50,
                iterations=10,
                diffusion_components=6,
                branch_confidence=0.55,
                seed=seed,
            )
            self.assertGreater(
                adjusted_rand_score(terminal_truth, result.branches[100:]),
                0.65,
                msg=f"seed={seed}",
            )


if __name__ == "__main__":
    unittest.main()
