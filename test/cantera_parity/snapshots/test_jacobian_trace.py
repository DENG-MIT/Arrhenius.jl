import unittest
from pathlib import Path

from jacobian_trace import cantera_attempt_trace


class JacobianTraceTests(unittest.TestCase):
    def test_stale_retry_refreshes_are_per_attempt_not_cumulative(self):
        fixture = Path(__file__).parent / "fixtures" / "nakamura_diffusion_jacobian_trace.log"
        attempts = cantera_attempt_trace(fixture)
        self.assertEqual([row["jacobians"] for row in attempts], [1, 1, 1])
        self.assertEqual([row["jacobians_cumulative"] for row in attempts], [1, 2, 3])
        self.assertEqual([row["jacobian_age_start"] for row in attempts], [1, 1, 1])
        self.assertEqual([row["accepted"] for row in attempts], [False, False, False])


if __name__ == "__main__":
    unittest.main()
