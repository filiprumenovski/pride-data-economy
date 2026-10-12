"""Deterministic selection of a dataset's main dated sequence database."""


def main_fasta(candidates):
    # Equal-sized files can have different date-bearing names. Break ties by
    # filename and location, not the platform-dependent order of an unstable sort.
    return candidates.sort_values(
        ['sz', 'accession', 'file_name', 'where'],
        ascending=[False, True, True, True], kind='stable',
    ).drop_duplicates('accession')
