import unittest
from pathlib import Path

from jacobian_trace import cantera_attempt_trace, latest_stat


class JacobianTraceTests(unittest.TestCase):
    def test_stale_retry_refreshes_are_per_attempt_not_cumulative(self):
        fixture = Path(__file__).parent / "fixtures" / "nakamura_diffusion_jacobian_trace.log"
        attempts = cantera_attempt_trace(fixture)
        self.assertEqual([row["jacobians"] for row in attempts], [1, 1, 1])
        self.assertEqual([row["jacobians_cumulative"] for row in attempts], [1, 2, 3])
        self.assertEqual([row["jacobian_age_start"] for row in attempts], [1, 1, 1])
        self.assertEqual([row["accepted"] for row in attempts], [False, False, False])

    def test_attempt_closes_before_the_next_ten_step_batch(self):
        fixture = Path(__file__).parent / "fixtures" / "steady_batch_jacobian_trace.log"
        attempts = cantera_attempt_trace(fixture)
        self.assertEqual([row["dt"] for row in attempts], [1e-5, 1e-5, 1.5e-5])
        self.assertEqual([row["jacobians"] for row in attempts], [2, 0, 1])
        self.assertEqual([row["jacobians_cumulative"] for row in attempts], [2, 2, 3])
        self.assertEqual([row["accepted"] for row in attempts], [True, True, False])

    def test_cumulative_jacobian_stat_uses_latest_snapshot(self):
        stats = {"jacobian_evals": [3, 4, 4], "steps": [1, 1, 1]}
        self.assertEqual(latest_stat(stats, "jacobian_evals"), 4)
        self.assertEqual(sum(stats["steps"]), 3)


if __name__ == "__main__":
    unittest.main()
