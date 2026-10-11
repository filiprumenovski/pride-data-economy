"""Registry of analysis software: how a tool is named, what native result files it writes,
and whether it can export PSI standard formats (mzIdentML / mzTab).

`native` patterns: regular expressions on the lower-cased file name with one archive/compression
suffix removed. `text`: regular expression on the dataset's software field + processing protocol.

Export capability (`export`):
  'native'   — the tool's own result file IS a PSI standard (written without an extra step)
  'built-in' — the tool can write mzIdentML or mzTab itself, as an optional export
  'none'     — no standard export documented in the tool's own docs/source
  None       — unknown
`since`: first year the export is documented (None = not dated). `default_on`: export is
enabled by default. Evidence grades: V read in primary source (docs/source/rules file);
C shown by real tool-written files or wrappers; S search-snippet only; U unverified.
Sources: PRIDE-2026 = reference/pride-submission-validation.yaml v2.2.2; agent research notes of
2 October 2026 summarised in outputs/Stack-and-reuse/Software-registry.md.
"""

TOOLS = {
    # --- commercial / DDA ----------------------------------------------------------------
    'Mascot': dict(text=r'mascot', native=[r'\.dat$', r'mascot'], archive=r'mascot',
                   export='built-in', formats='mzIdentML, mzTab', since=2013, default_on=False, grade='C/S'),
    'Proteome Discoverer': dict(text=r'proteome discoverer', native=[r'\.(msf|pdresult|pdstudy|pdresultview|msfview|pdresultdetails)$'],
                   archive=r'pdresult|proteome.?discoverer', export='built-in', formats='mzIdentML (mzTab S)', since=2017,
                   default_on=False, grade='V (PRIDE recommends) / S (version)'),
    'MaxQuant': dict(text=r'maxquant', native=[r'(^|/)[^/]*(evidence|peptides|proteingroups|msms|allpeptides|modificationspecificpeptides|msmsscans|summary|parameters)\.txt$',
                   r'mqpar[^/]*\.xml$', r'\.apl$'], archive=r'maxquant|(^|[/_.-])txt$|combined',
                   export='built-in', formats='mzTab (no mzIdentML found)', since=2018, default_on=False, grade='C/S'),
    'PEAKS': dict(text=r'peaks ?(studio|x\b|xpro|x\+|db|online|software|ab\b|au\b|\d)|bioinformatics solutions|"name":"peaks|\(peaks\)',
                  native=[r'peaks', r'db search psm', r'protein-peptides'], archive=r'peaks',
                  export='built-in', formats='mzIdentML', since=None, default_on=False, grade='S'),
    'Scaffold': dict(text=r'scaffold', native=[r'\.(sf3|sfd)$', r'scaffold'], archive=r'scaffold',
                     export='built-in', formats='mzIdentML', since=2014, default_on=False, grade='C'),
    'ProteinPilot': dict(text=r'proteinpilot|protein pilot|paragon', native=[r'\.group$', r'proteinpilot'], archive=r'proteinpilot|\.group$',
                         export='built-in', formats='mzIdentML (reported to fail PRIDE validation)', since=None, default_on=False, grade='C/S'),
    'Progenesis': dict(text=r'progenesis', native=[r'progenesis'], archive=r'progenesis',
                       export='built-in', formats='mzIdentML (Mascot/Ion Accounting searches)', since=None, default_on=False, grade='S'),
    'Byonic': dict(text=r'byonic', native=[r'\.byrslt$', r'byonic'], archive=r'byonic',
                   export='built-in', formats='mzIdentML', since=None, default_on=False, grade='S'),
    'PLGS': dict(text=r'plgs|proteinlynx', native=[r'final_peptide', r'plgs'], archive=r'plgs',
                 export='built-in', formats='mzIdentML (v3.0.3+)', since=None, default_on=False, grade='S'),
    # --- DIA / targeted ---------------------------------------------------------------------
    'DIA-NN': dict(text=r'dia-?nn', native=[r'(^|[/_.-])report\.(tsv|parquet)$', r'\.(pr|pg|gg|unique_genes)_matrix\.tsv$', r'report[^/]*_matrix',
                   r'report\.stats\.tsv$', r'\.speclib$', r'site_report\.(tsv|parquet)$', r'diann|dia-nn'], archive=r'dia-?nn|diann',
                   export='none', formats='-', since=None, default_on=False, grade='V (absent from docs)'),
    'Spectronaut': dict(text=r'spectronaut', native=[r'\.(sne|kit|bgsfasta)$', r'spectronaut'], archive=r'spectronaut',
                        export=None, formats='unknown', since=None, default_on=False, grade='U'),
    'FragPipe': dict(text=r'msfragger|fragpipe', native=[r'(^|[/_.-])(psm|peptide|protein|ion)\.tsv$', r'combined_(protein|peptide|ion|modified_peptide)[^/]*\.tsv$',
                     r'fragpipe\.workflow$|fp-manifest$|fragger\.params$', r'fragpipe', r'\.pepxml$'], archive=r'fragpipe|msfragger',
                     export='none', formats='mzIdentML only via Philosopher CLI; not in current GUI', since=None, default_on=False, grade='V'),
    'Skyline': dict(text=r'skyline', native=[r'\.(sky|skyd|blib|skyl|skyr)$'], archive=r'skyline|\.sky$',
                    export='none', formats='-', since=None, default_on=False, grade='V'),
    'OpenSWATH': dict(text=r'openswath|pyprophet', native=[r'\.(osw|oswpq|pqp|sqmass)$'], archive=r'openswath|pyprophet',
                      export='none', formats='-', since=None, default_on=False, grade='V'),
    # EncyclopeDIA text pattern excludes KEGG and the Cancer Cell Line Encyclopedia
    'EncyclopeDIA': dict(text=r'(?<!kyoto )(?<!line )encyclopedia(?! of)', native=[r'\.(elib|dlib)$', r'encyclopedia'], archive=r'encyclopedia',
                         export=None, formats='unknown', since=None, default_on=False, grade='U'),
    # --- open-source DDA ----------------------------------------------------------------------
    'MS-GF+': dict(text=r'ms-?gf', native=[r'\.mzid$'], archive=r'msgf',
                   export='native', formats='mzIdentML', since=None, default_on=True, grade='V'),
    'MetaMorpheus': dict(text=r'metamorpheus', native=[r'\.psmtsv$', r'allquantified'], archive=r'metamorpheus',
                         export='native', formats='mzIdentML (on by default)', since=None, default_on=True, grade='V'),
    'Comet': dict(text=r'\bcomet\b', native=[r'\.(sqt|pin)$', r'comet'], archive=r'comet',
                  export='built-in', formats='mzIdentML (from 2020)', since=2020, default_on=False, grade='V'),
    'X!Tandem': dict(text=r'x!? ?tandem|xtandem', native=[r'\.t\.xml$', r'tandem'], archive=r'tandem',
                     export='built-in', formats='mzIdentML (snippet)', since=2015, default_on=False, grade='S'),
    'PeptideShaker': dict(text=r'peptideshaker|searchgui', native=[r'\.(cpsx|psdb)$', r'peptideshaker|searchgui'], archive=r'peptideshaker|searchgui',
                          export='built-in', formats='mzIdentML', since=None, default_on=False, grade='V'),
    'OpenMS': dict(text=r'openms', native=[r'\.(idxml|consensusxml|featurexml)$'], archive=r'openms',
                   export='built-in', formats='mzIdentML, mzTab', since=None, default_on=False, grade='V'),
    'TPP': dict(text=r'trans-proteomic|peptideprophet|proteinprophet|\btpp\b', native=[r'\.(pep|prot)\.?xml$', r'interact'], archive=r'interact|tpp',
                export='built-in', formats='mzIdentML via tpp2mzid', since=None, default_on=False, grade='S'),
    'Percolator': dict(text=r'percolator', native=[r'\.pout$', r'percolator'], archive=r'percolator',
                       export='none', formats='-', since=None, default_on=False, grade='V'),
    # Sage text pattern excludes 'Sage-N Research' (SORCERER vendor)
    'Sage': dict(text=r'\bsage\b(?! group)(?!-n)', native=[r'results\.sage\.(tsv|parquet)$', r'\.sage\.'], archive=r'sage',
                 export='none', formats='-', since=None, default_on=False, grade='V'),
    'pFind': dict(text=r'pfind|plink', native=[r'\.spectra$', r'pfind|plink'], archive=r'pfind|plink',
                  export='none', formats='-', since=None, default_on=False, grade='V/U'),
}


def compile_registry():
    import re
    out = {}
    for k, t in TOOLS.items():
        out[k] = dict(t, text_re=re.compile(t['text']), native_re=re.compile('|'.join(f'(?:{p})' for p in t['native'])),
                      archive_re=re.compile(t['archive']) if t.get('archive') else None)
    return out
