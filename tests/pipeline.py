import unittest

from phyloblitz.pipeline import (
    count_spoa_aln_vars,
    filter_paf_overhang,
    lists_common_prefix,
    parse_spoa_r2,
)


class TestReportFunctions(unittest.TestCase):

    def test_parse_spoa_r2(self):
        f = """>S1\nAATAC\n>S2\nAAATC\n>Consensus\nAAAAC"""
        self.assertEqual(
            parse_spoa_r2(f),
            {"S1": "AATAC", "S2": "AAATC", "Consensus": "AAAAC"},
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

    def test_lists_common_prefix(self):
        lol = [
            ["a", "b", "c", "d"],
            ["a", "b", "c", "e", "p"],
            ["a", "b", "c", "f"],
            ["a", "b", "z", "g", "q", "r"],
        ]
        self.assertEqual(lists_common_prefix(lol), ["a", "b"])

    def test_filter_paf_overhang(self):
        rec_ok = "qname\t1000\t5\t800\t+\ttname\t1000\t200\t990\t750\t850\t0\ttp:A:S\n"
        rec_no = "qname\t1000\t5\t800\t+\ttname\t1000\t5\t800\t750\t850\t0\ttp:A:S\n"
        self.assertEqual(rec_ok, filter_paf_overhang(rec_ok, max_overhang_frac=0.05))
        self.assertIsNone(filter_paf_overhang(rec_no, max_overhang_frac=0.05))


if __name__ == "__main__":
    unittest.main()
