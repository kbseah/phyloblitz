import unittest

from phyloblitz.pipeline import parse_spoa_r2, count_spoa_aln_vars


class TestReportFunctions(unittest.TestCase):

    def test_parse_spoa_r2(self):
        f = """>S1\nAATAC\n>S2\nAAATC\n>Consensus\nAAAAC"""
        self.assertEqual(
            parse_spoa_r2(f), {"S1": "AATAC", "S2": "AAATC", "Consensus": "AAAAC"}
        )

    def test_count_spoa_aln_vars(self):
        s = {"S1": "---AT-CG-A--", "Consensus": "ATAA-ACG-TGC"}
        o = count_spoa_aln_vars(s)
        self.assertEqual(
            [
                o["S1"]["query_lead_gap"],
                o["S1"]["query_trail_gap"],
                o["S1"]["match"],
                o["S1"]["cons_gap"],
                o["S1"]["query_gap"],
            ],
            [3, 2, 3, 1, 1],
        )


if __name__ == "__main__":
    unittest.main()
