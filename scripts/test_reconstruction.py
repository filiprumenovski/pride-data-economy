"""Checks for the inference boundaries in the public reconstruction runner."""
import unittest
import pandas as pd
from release_geography import reconstruct


class GeographyEvidenceTests(unittest.TestCase):
    @staticmethod
    def contact(acc, repo='MassIVE', country='USA', affinity='A', email='E', role='pi'):
        return dict(repository=repo, accession=acc, role=role, affinity=affinity, email=email,
                    person=acc, raw_country=country, person_raw_country=country,
                    project_country=None, affiliation_suffix=None)

    def target(self, seeds, **changes):
        row = self.contact('PXD1', repo='PRIDE', country=None); row.update(changes)
        _, result = reconstruct(pd.DataFrame([*seeds, row]))
        return result.set_index('accession').loc['PXD1']

    def test_one_reference_is_not_enough(self):
        self.assertTrue(pd.isna(self.target([self.contact('MSV1')]).country))

    def test_two_distinct_sources_support_inference(self):
        result = self.target([self.contact('MSV1'), self.contact('MSV2')])
        self.assertEqual(result.country, 'United States')
        self.assertTrue(result.pi_geo_inferred_only)

    def test_duplicate_contact_does_not_count_as_second_source(self):
        self.assertTrue(pd.isna(self.target([self.contact('MSV1'), self.contact('MSV1')]).country))

    def test_conflicting_sources_remain_unresolved(self):
        self.assertTrue(pd.isna(self.target([self.contact('MSV1'), self.contact('MSV2'), self.contact('MSV3', country='Canada')]).country))

    def test_reported_country_takes_priority(self):
        result = self.target([self.contact('MSV1'), self.contact('MSV2')], raw_country='Canada')
        self.assertEqual(result.country, 'Canada')
        self.assertFalse(result.pi_geo_inferred_only)

    def test_email_can_resolve_unmatched_affiliation(self):
        self.assertEqual(self.target([self.contact('MSV1'), self.contact('MSV2')], affinity='Different').country, 'United States')


if __name__ == '__main__':
    unittest.main()
