"""Equal-sized FASTA files must not change date selection across row orders."""
import sys
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'work/stack-and-reuse'))
from database_rules import main_fasta


class DatabaseSelectionTests(unittest.TestCase):
    def test_size_then_name_determine_selection(self):
        candidates = pd.DataFrame([
            ('PXD1', 100, '2024-01-01.fasta', 'top'),
            ('PXD1', 100, '2023-01-01.fasta', 'top'),
            ('PXD1', 99, '2022-01-01.fasta', 'top'),
            ('PXD2', 10, 'x.fasta', 'archive'),
        ], columns=['accession', 'sz', 'file_name', 'where'])
        for seed in range(10):
            chosen = main_fasta(candidates.sample(frac=1, random_state=seed)).set_index('accession')
            self.assertEqual(chosen.loc['PXD1', 'file_name'], '2023-01-01.fasta')
            self.assertEqual(len(chosen), 2)


if __name__ == '__main__':
    unittest.main()
