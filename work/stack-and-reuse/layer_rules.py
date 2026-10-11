"""File-layer classification rules shared by layers.py and audit_sheet.py (SQL conditions on
n = lower-cased file name, b = n without one compression/archive suffix, cat = PRIDE file category)."""
RULES = [  # (layer, SQL condition) — order matters
    ('standard_results', r"regexp_matches(n, '\.(mzid|mzidentml|mztab)(\.txt)?(\.(gz|zip|bz2))?$')"),
    # SDRF tables only; quantms writes derived files whose names start with 'sdrf_openms'
    ('sample_metadata', r"n like '%sdrf%' and n not like '%openms%' and n not like '%msstats%' and regexp_matches(b, '\.(tsv|txt|csv|xlsx)$')"),
    ('open_spectra', r"regexp_matches(b, '\.(mzml|mzxml|mgf|ms2|mzdata|mz5|imzml|ibd|dta|pkl)$')"),
    ('spectral_library', r"regexp_matches(b, '\.(speclib|blib|msp|sptxt|dlib|elib|kit|tsv\.speclib)$') or n like '%spectral_lib%' or n like '%speclib%'"),
    ('tool_results', r"""regexp_matches(b, '(^|/)(proteingroups|evidence|msms|msmsscans|peptides|allpeptides|modificationspecificpeptides|summary|parameters)\.txt$')
        or n like '%mqpar%' or regexp_matches(b, '\.(msf|pdresult|pdstudy|dat|sne|sf3|group|sky|skyd|cpsx|psdb|idxml|featurexml|consensusxml|pepxml|protxml|byrslt|pin)$')
        or regexp_matches(b, '\.(pep|prot)\.xml$')
        or regexp_matches(b, '\.(sqt|msgfdb_out|inspect_out|xcorr_hit_list_best)$')
        or (regexp_matches(b, '\.(out|xml)$') and cat in ('SEARCH', 'RESULT'))
        or regexp_matches(b, '(^|[/_.-])report(\.pr_matrix|\.pg_matrix|\.gg_matrix|\.unique_genes_matrix|\.stats)?\.(tsv|parquet)$')
        or regexp_matches(b, '(^|[/_.-])report\.(pr|pg|gg)_matrix') or b like '%.apl'
        or (regexp_matches(n, '\.(zip|tgz|tar|rar|7z|gz)$') and regexp_matches(n, 'maxquant|txt\.(zip|tgz|tar|rar|7z|gz)$|combined|proteingroups'))
        or regexp_matches(b, '(^|/)(psm|peptide|protein|ion|combined_protein|combined_peptide|combined_ion|combined_modified_peptide)\.tsv$')
        or (n like '%txt%' and n like '%combined%' and cat <> 'RAW')"""),
    # other archives filed by the submitter as search/result files: results of unknown format, packed
    ('archived_results', r"regexp_matches(n, '\.(zip|tgz|tar|rar|7z|gz)$') and b not like '%.%' and cat in ('SEARCH', 'RESULT')"),
    ('vendor_raw', r"""cat = 'RAW' or regexp_matches(b, '\.(raw|wiff|wiff2|scan|d|tdf|tdf_bin|baf|lcd|yep|t2d|fid|ser|uimf)$')"""),
    ('sequence_db', r"regexp_matches(b, '\.(fasta|fa|faa)$')"),
    # submission-tool checksums and run logs are bookkeeping, not results
    ('bookkeeping', r"n like '%checksum%' or regexp_matches(b, '\.(log|md5)$') or regexp_matches(b, '\.log\.txt$')"),
    ('unidentified_tables', r"regexp_matches(b, '\.(tsv|csv|xlsx|xls|txt|parquet)$')"),
    ('documents', r"regexp_matches(b, '\.(pdf|docx?|pptx?|md|html?|png|jpe?g|tiff?|svg)$') or n like '%readme%'"),
]
