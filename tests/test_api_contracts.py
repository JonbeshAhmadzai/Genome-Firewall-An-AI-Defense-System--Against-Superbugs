from __future__ import annotations

import unittest

from api import config, health


class ApiContractTests(unittest.TestCase):
    def test_health_and_config_contracts(self) -> None:
        self.assertEqual(health()["status"], "ok")
        payload = config()
        self.assertIn("Escherichia coli", payload["species"])
        self.assertIn("ciprofloxacin", payload["drugs"]["Escherichia coli"])
        self.assertIn("llm_available", payload)


if __name__ == "__main__":
    unittest.main()
