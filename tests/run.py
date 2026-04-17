import unittest

from phyloblitz.run import merge_intervals


class TestReportFunctions(unittest.TestCase):

    def test_merge_intervals(self):
        intervals = [
            (1, 5),
            (20, 25),
            (24, 30),
            (80, 85),
            (85, 90),
            (92, 95),
            (96, 100),
        ]
        result = [(1, 5), (20, 30), (80, 90), (92, 95), (96, 100)]
        out = merge_intervals(intervals)
        self.assertEqual(
            sorted(result),
            sorted(out),
        )


if __name__ == "__main__":
    unittest.main()
