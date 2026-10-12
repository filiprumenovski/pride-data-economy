"""Regression checks for the release's numerical comparison gate."""
import tempfile
import unittest
from pathlib import Path
from compare_reconstruction import compare


class ComparisonTests(unittest.TestCase):
    def test_required_missing_output_cannot_pass(self):
        with tempfile.TemporaryDirectory() as d:
            reference, rebuilt = Path(d) / 'reference', Path(d) / 'rebuilt'
            reference.mkdir(); rebuilt.mkdir()
            (reference / 'table.csv').write_text('year,value\n2025,1\n')
            result = compare(reference, rebuilt, ['table.csv'])
            self.assertEqual(result[0]['status'], 'missing')

    def test_order_is_ignored_but_changed_count_is_detected(self):
        with tempfile.TemporaryDirectory() as d:
            reference, rebuilt = Path(d) / 'reference', Path(d) / 'rebuilt'
            reference.mkdir(); rebuilt.mkdir()
            (reference / 'table.csv').write_text('year,value\n2024,10\n2025,20\n')
            (rebuilt / 'table.csv').write_text('year,value\n2025,20\n2024,10\n')
            self.assertEqual(compare(reference, rebuilt, ['table.csv'])[0]['status'], 'match')
            (rebuilt / 'table.csv').write_text('year,value\n2025,21\n2024,10\n')
            self.assertEqual(compare(reference, rebuilt, ['table.csv'])[0]['status'], 'changed')


if __name__ == '__main__':
    unittest.main()
