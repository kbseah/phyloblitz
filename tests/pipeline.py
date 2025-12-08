import unittest

from phyloblitz.pipeline import parse_spoa_r2


class TestReportFunctions(unittest.TestCase):

    def test_parse_spoa_r2(self):
        f = """>S1\nAATAC\n>S2\nAAATC\n>Consensus\nAAAAC"""
        self.assertEqual(
            parse_spoa_r2(f), {"S1": "AATAC", "S2": "AAATC", "Consensus": "AAAAC"}
        )


if __name__ == "__main__":
    unittest.main()
