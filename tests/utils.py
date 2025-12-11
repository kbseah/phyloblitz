#!/usr/bin/env python3

import unittest

from phyloblitz.utils import lists_common_prefix, filter_paf_overhang


class TestReportFunctions(unittest.TestCase):

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
