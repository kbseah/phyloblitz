#!/usr/bin/env python3

import unittest

from phyloblitz.report import lists_common_prefix


class TestReportFunctions(unittest.TestCase):

    def test_lists_common_prefix(self):
        lol = [
            ["a", "b", "c", "d"],
            ["a", "b", "c", "e", "p"],
            ["a", "b", "c", "f"],
            ["a", "b", "z", "g", "q", "r"],
        ]
        self.assertEqual(lists_common_prefix(lol), ["a", "b"])


if __name__ == "__main__":
    unittest.main()
