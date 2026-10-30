"""Arsenal H5, H6 and H8 (October 2026): quantification/PTM strategy, declared vs evidenced
software, and documentation richness. Descriptive only; offline (reads the release databases,
the tracked cohort and the outputs of archive_contents.py / header_sniff.py).

Evidence sources (each dataset's denominator for a source is recorded):
  top     top-level, non-generated file names (pride_files, '/generated/' paths excluded)
  zip     member paths inside listed results zips (Archive-contents/Members.parquet; rar, 7z and
          tar archives are not listed, so absence there is "not detected")
  header  header rows of top-level tables that the file-name rules could not identify
          (Header-sniff/Headers.jsonl, fixed sample of <= 20 tables per dataset; tables inside zips
          and canonically named tool tables such as proteinGroups.txt were NOT read)
  text    metadata: title, description and keywords ("design text"); sample-processing and
          data-processing protocols ("protocol text"); PRIDE's structured fields (quantification
          methods, software, identified PTMs)

H5  Quantification and PTM strategy: modification site tables (MaxQuant '<mod>Sites.txt',
    FragPipe 'combined_site_*', DIA-NN '*.unimod_<id>_sites_*' and TMT-Integrator single/multi-site
    tables), MaxQuant processing steps (proc/ folder), and quantification columns in table
    headers (isobaric reporter channels, LFQ intensity, iBAQ, SILAC ratios), compared with the
    metadata flags (tmt_expanded, lfq_expanded, silac_expanded) and text mentions.
H6  Software named in metadata (software_registry TOOLS text regexes on software field + data
    processing protocol, as in common.software_flags) vs software evidenced by conservative,
    tool-specific file signatures (SIG below) and table-header signatures (HEADER_SIG below).
    '.peaks'/'.peaks_t2' files (MaxQuant intermediates in '<run>/p0/' folders) are NOT PEAKS evidence;
    generic names (peptides.txt, psm.tsv, protein.tsv, report.tsv with a prefix, .mzid, .pin, .pepXML,
    .dat without the Mascot F-number pattern) are NOT used. Executables/libraries are dropped.
H8  Documentation: description/title length, keyword count, protocol lengths, README, SDRF,
    experimental-design tables (author-written vs tool-written), search parameter/log files;
    trends and descriptive association with the results form (Habits/Dataset-form.csv.gz).

Rates are reported with group-bootstrap 95% intervals (hashed pi_group; datasets without a group
id are their own cluster; seed 20261001). Run from work/stack-and-reuse/.
"""
import json, re, sys
from collections import Counter
import numpy as np, pandas as pd
from common import R, O, SEED, features, cluster_boot
from software_registry import TOOLS
sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H5-H6-H8-strategies-software-docs'; OUT.mkdir(parents=True, exist_ok=True)
AC, HS = O / 'Archive-contents', O / 'Header-sniff'
YEARS = list(range(2013, 2026))
PERIOD_CUTS, PERIOD_LABS = [2012, 2016, 2019, 2022, 2025], ['2013-16', '2017-19', '2020-22', '2023-25']
STRIP = r'\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$'
BINARY = r'\.(exe|dll|jar|so|class|pyc|pdb|dylib|msi)$'
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
SUMMARY = {}


def sq(s):
    return "'" + s.replace("'", "''") + "'"


def boot(values, groups, reps=1000):
    """Rate (%) with group-bootstrap 95% CI (common.cluster_boot, seed fixed)."""
    v = np.asarray(values, float)
    if len(v) == 0:
        return np.nan, np.nan, np.nan
    lo, hi = cluster_boot(v, np.asarray(groups), reps=reps, seed=SEED)
    return round(100 * v.mean(), 1), round(100 * lo, 1), round(100 * hi, 1)


def rate_row(mask, d, label, **extra):
    x = d[mask] if mask is not None else d
    return dict(extra, subset=label, n=len(x), **{k: v for k, v in zip(['pct', 'lo', 'hi'], boot(x.value, x.grp))})


# =============================================================================================
# 0. Base table
# =============================================================================================
F = features()
F['grp'] = F.pi_group.fillna(F.accession)
con = connect()
con.execute("set threads=1; set memory_limit='3GB'; set preserve_insertion_order=false")
T = con.execute("""select s.accession, lower(coalesce(s.title,'')) title, lower(coalesce(s.description,'')) description,
  lower(coalesce(s.keywords,'')) keywords, lower(coalesce(s.processing,'')) processing, lower(coalesce(s.quant,'')) quant,
  lower(coalesce(s.experiments,'')) experiments, lower(coalesce(s.software,'')) software,
  lower(coalesce(b.sampleProcessingProtocol,'')) sample_protocol, lower(coalesce(b.dataProcessingProtocol,'')) data_protocol,
  lower(coalesce(array_to_string(list_transform(b.identifiedPTMStrings, x -> x.name), ' | '), '')) ptm_struct,
  coalesce(p.keyword_count, 0) n_keywords, coalesce(s.description,'') description_raw
  from study s left join pride_bronze_projects b on b.accession = s.accession
  left join pride_projects p on p.project_accession = s.accession""").df()
assert len(T) == 36249 and T.accession.is_unique
D = F[['accession', 'release_year', 'grp', 'pi_group', 'dia_expanded', 'tmt_expanded', 'lfq_expanded', 'silac_expanded',
       'family', 'vendor_single', 'standards', 'has_software_annotation', 'submission_type', 'group_seq', 'group_first_year', 'submission_date']
      + [c for c in F.columns if c.startswith('sw_')]].merge(T, on='accession', validate='one_to_one')
form = pd.read_csv(O / 'Habits/Dataset-form.csv.gz')[['accession', 'form', 'recognisable']]
D = D.merge(form, on='accession', how='left', validate='one_to_one')
D['period'] = pd.cut(D.release_year, PERIOD_CUTS, labels=PERIOD_LABS).astype(str)
D['design_text'] = D.title + ' ' + D.description + ' ' + D.keywords
D['protocol_text'] = D.sample_protocol + ' ' + D.data_protocol
D['all_text'] = D.design_text + ' ' + D.protocol_text + ' ' + D.quant + ' ' + D.experiments + ' ' + D.software

# =============================================================================================
# 1. File-name patterns (applied to lower-cased path, one compression suffix removed)
# =============================================================================================
# H6: conservative tool signatures. Spot-checked against the most frequent matching names
# (written to H6-signature-hits.csv). Deliberately excluded as not tool-specific: bare
# peptides.txt/summary.txt/parameters.txt, psm/peptide/protein/ion.tsv, *_report.tsv, .mzid, .pin,
# .pepXML (written by Comet, MSFragger, X!Tandem...), .dat without an F-number, .blib, '.peaks*'.
SIG = {
    'MaxQuant': r'(^|/)(proteingroups|evidence|msms|allpeptides|modificationspecificpeptides|msmsscans)\.txt$|(^|/)mqpar[^/]*\.xml$|\.apl$'
                r'|(^|/)[^/]*\([^/)]*\) ?sites\.txt$|(^|/)combined/(txt|proc|andromeda)/',
    'Proteome Discoverer': r'\.(msf|pdresult|pdstudy|pdresultview|msfview|pdprocessingwf|pdconsensuswf)$',
    'Mascot': r'(^|/)f\d{5,6}\.dat$|mascot',
    'PEAKS': r'(^|/)(db search psm|de novo peptides|protein-peptides|peaks[^/]*)\.(csv|txt|tsv|xlsx)$|(^|[/_ .-])peaks[ _-]?(studio|xpro|online|db)\b',
    'Scaffold': r'\.(sf3|sfd)$|scaffold',
    'ProteinPilot': r'\.group$|proteinpilot|protein ?pilot',
    'Progenesis': r'progenesis',
    'Byonic': r'\.byrslt$|byonic',
    'PLGS': r'final_peptide|final_fragment|(^|[/_ .-])plgs',
    'DIA-NN': r'(^|/)report(-first-pass)?\.(tsv|parquet|stats\.tsv|log\.txt|manifest\.txt|protein_description\.tsv)$'
              r'|\.(pr|pg|gg|unique_genes)_matrix\.(tsv|parquet)$|\.speclib$|unimod_\d+_sites_\d+\.tsv$|diann|dia-nn|dia_nn',
    'Spectronaut': r'\.(sne|kit|bgsfasta)$|spectronaut',
    'FragPipe': r'(^|/)combined_(protein|peptide|ion|modified_peptide|site_[^/]*)\.tsv$|fragpipe|msfragger|fragger\.params$|\.fp-manifest$'
                r'|(^|/)(abundance|ratio)_(protein|gene|peptide|single-site|multi-site)_[a-z]+\.tsv$|(^|/)ptm-shepherd',
    'Skyline': r'\.(sky|skyd|skyl|skyr)$',
    'OpenSWATH': r'\.(osw|oswpq|pqp|sqmass)$|pyprophet|openswath',
    'EncyclopeDIA': r'\.(elib|dlib)$|encyclopedia',
    'MS-GF+': r'msgf',
    'MetaMorpheus': r'\.psmtsv$|metamorpheus|(^|/)allquantified',
    'Comet': r'(^|[/_.-])comet([/_.-]|$)',
    'X!Tandem': r'\.t\.xml$|x!?tandem',
    'PeptideShaker': r'\.(cpsx|psdb)$|peptideshaker|searchgui',
    'OpenMS': r'\.(idxml|consensusxml|featurexml)$',
    'TPP': r'(^|/)interact[^/]*\.(pep|prot|ipro)\.xml$|\.prot\.xml$|\.protxml$|iprophet',
    'Percolator': r'\.pout$|percolator',
    'Sage': r'results\.sage\.(tsv|parquet)$|\.sage\.',
    'pFind': r'pfind|(^|[/_.-])plink([/_.-]|$)',
}
assert set(SIG) <= set(TOOLS), set(SIG) - set(TOOLS)
DESIGN_TOOL = (r'(^|/)experimentaldesign(template)?[^/]*\.(txt|tsv|xlsx|csv)$|(^|/)experimental design template[^/]*$|(^|/)filelist_[^/]*\.txt$'
               r'|(^|/)experiment_annotation\.tsv$|\.fp-manifest$|(^|/)report(-first-pass)?\.manifest\.txt$')
OTHER_PAT = {
    # H5
    'site_table': r'(^|/)[^/]*sites\.(txt|tsv|csv|xlsx|xls)$|(^|/)combined_site_[^/]*\.tsv$|unimod_\d+_sites_\d+\.tsv$|(^|/)(abundance|ratio)_(single|multi)-site[^/]*\.tsv$',
    'mq_core_txt': r'(^|/)(proteingroups|evidence)\.txt$',
    'mq_oxidation_sites': r'(^|/)oxidation \(m\) ?sites\.txt$',
    'mq_proc_lfq': r'(^|/)proc/label-free[ _]quantification',
    'mq_proc_ibaq': r'(^|/)proc/ibaq',
    'mq_proc_any': r'(^|/)proc/[^/]*\.finished\.txt$',
    'tmt_integrator': r'(^|/)(abundance|ratio)_(protein|gene|peptide|single-site|multi-site)_[a-z]+\.tsv$',
    'fragpipe_diann_folder': r'(^|/)diann-output/',
    # H8 documentation
    'doc_readme': (r'(^|/)[^/]*read[ _-]?me[^/]*$', r'licen[cs]e|third-party'),
    # bounded so that peptide-sequence file names containing 'sdrf' are not matched; quantms-derived files excluded
    'doc_sdrf': (r'(^|/)([^/]*[_. -])?sdrf([_. -][^/]*)?\.(tsv|txt|csv|xlsx|xls)$', r'openms|msstats'),
    # tool-written sample-to-run mappings: MaxQuant experimentalDesign(Template), FragPipe manifest/annotation/filelist, DIA-NN manifest
    'doc_design_tool': DESIGN_TOOL,
    'doc_design_author': (r'(^|/)[^/]*(experiment(al)?[ _-]?design|study[ _-]?design|sample[ _-]?(info|information|list|key|description|annotation|annotations|sheet|table|metadata|map|mapping|overview)|samplesheet|file[ _-]?(description|list|key|annotation)|metadata|(^|[_ -])design)[^/]*\.(txt|tsv|csv|xlsx|xls|docx|doc|pdf)$', DESIGN_TOOL + '|metadatasection'),
    'doc_params': r'(^|/)mqpar[^/]*\.xml$|\.workflow$|fragger\.params$|comet[^/]*\.params$|(^|/)report(-first-pass)?\.log\.txt$|\.pd(processing|consensus)wf$|(^|/)searchgui_input|open-pfind\.cfg$|(^|/)txt/parameters\.txt$',
    'doc_document': r'\.(pdf|docx?|pptx?|md|rtf)$',
}
ALL_PAT = {**{'tool:' + k: v for k, v in SIG.items()}, **OTHER_PAT}
inc = lambda v: v[0] if isinstance(v, tuple) else v
union = '|'.join(f'(?:{inc(v)})' for v in ALL_PAT.values())

con.execute(f"""create temp table top as
  select project_accession accession, lower(regexp_replace(file_name, {sq(STRIP)}, '', 'i')) p, coalesce(file_category_value, '') cat
  from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
  and project_accession in (select accession from study) and not regexp_matches(lower(file_name), {sq(BINARY)})""")
con.execute(f"""create temp table zm as
  select a.accession, lower(regexp_replace(m.member, {sq(STRIP)}, '', 'i')) p
  from read_parquet('{AC / 'Members.parquet'}') m
  join (select url, any_value(accession) accession from read_csv('{AC / 'Archives.csv.gz'}', sample_size=-1) group by url) a using (url)
  where not ends_with(m.member, '/') and not regexp_matches(m.member, '(^|/)__MACOSX/|(^|/)\\._')
    and not regexp_matches(lower(m.member), {sq(BINARY)})
    and regexp_matches(lower(regexp_replace(m.member, {sq(STRIP)}, '', 'i')), {sq(union)})""")
cols = ', '.join(f'bool_or(regexp_matches(p, {sq(r[0])}) and not regexp_matches(p, {sq(r[1])})) "{k}"' if isinstance(r, tuple)
                 else f'bool_or(regexp_matches(p, {sq(r)})) "{k}"' for k, r in ALL_PAT.items())
FT = con.execute(f'select accession, {cols} from top group by 1').df().set_index('accession')
FZ = con.execute(f'select accession, {cols} from zm group by 1').df().set_index('accession')
FT = FT.reindex(D.accession).fillna(False).astype(bool); FZ = FZ.reindex(D.accession).fillna(False).astype(bool)
A = pd.read_csv(AC / 'Archives.csv.gz')
listed = set(A[A.status.isin(['ok', 'partial'])].accession)
D['zip_listed'] = D.accession.isin(listed)
# Run-name evidence (weak): top-level raw files whose names mention TMT/iTRAQ or phospho enrichment
RAWX = r'\.(raw|wiff|wiff2|d|tdf|baf|lcd|mzml|mzxml|mgf)$'
runs = con.execute(f"""select accession, bool_or(regexp_matches(p, '(tmt|itraq)')) run_isobaric,
  bool_or(regexp_matches(p, 'phospho|tio2|imac|fe-?nta|ptyr|psty')) run_phospho, count(*) n_runs
  from top where cat = 'RAW' or regexp_matches(p, {sq(RAWX)}) group by 1""").df().set_index('accession').reindex(D.accession)
D['has_runs'] = runs.n_runs.notna().to_numpy()
D['run_isobaric'] = runs.run_isobaric.fillna(False).astype(bool).to_numpy()
D['run_phospho'] = runs.run_phospho.fillna(False).astype(bool).to_numpy()

# signature spot-check sheet: most frequent matching names per tool and source
hits = []
for k, r in SIG.items():
    for src in ['top', 'zm']:
        h = con.execute(f"""select regexp_extract(p, '[^/]*/?[^/]*$') as name, count(distinct accession) as datasets from {src}
                            where regexp_matches(p, {sq(r)}) group by 1 order by datasets desc, name limit 15""").df()
        h.insert(0, 'source', 'top-level' if src == 'top' else 'inside zip'); h.insert(0, 'tool', k); hits.append(h)
pd.concat(hits).to_csv(OUT / 'H6-signature-hits.csv', index=False)
# PEAKS pitfall check: how many datasets would a naive '.peaks' rule have credited?
naive = con.execute(f"""select count(distinct a.accession) from read_parquet('{AC / 'Members.parquet'}') m
  join (select url, any_value(accession) accession from read_csv('{AC / 'Archives.csv.gz'}', sample_size=-1) group by url) a using (url)
  where regexp_matches(lower(m.member), '\\.peaks(_t2)?$')""").fetchone()[0]
naive_p0 = con.execute(f"""select count(distinct a.accession) from read_parquet('{AC / 'Members.parquet'}') m
  join (select url, any_value(accession) accession from read_csv('{AC / 'Archives.csv.gz'}', sample_size=-1) group by url) a using (url)
  where regexp_matches(lower(m.member), '(^|/)p0/[^/]*\\.peaks(_t2)?$')""").fetchone()[0]
SUMMARY['peaks_pitfall'] = {'datasets_with_.peaks_members_in_zips': int(naive), 'of_which_in_p0_folders': int(naive_p0),
                            'note': "'.peaks' members are MaxQuant intermediates; not counted as PEAKS evidence"}

# =============================================================================================
# 2. Table headers
# =============================================================================================
HEADER_SIG = {  # conservative subset of header_sniff.SIGNATURES, mapped to registry tool names
    'DIA-NN': [['protein.group', 'protein.ids', 'genes'], ['protein.group', 'protein.names'], ['precursor.id', 'stripped.sequence'],
               ['file.name', 'protein.group', 'precursor.id']],
    'Spectronaut': [['pg.proteingroups'], ['pg.proteinaccessions'], ['eg.precursorid'], ['r.filename', 'pg.genes']],
    'MaxQuant': [['majority protein ids'], ['razor + unique peptides'], ['leading razor protein', 'raw file'],
                 ['protein ids', 'peptide counts (all)'], ['positions within proteins', 'localization prob']],
    'Proteome Discoverer': [['master', 'accession'], ['annotated sequence', '# protein groups'], ['confidence', 'annotated sequence'],
                            ['identifying node']],
    'FragPipe': [['protein id', 'protein probability'], ['spectrum file', 'modified peptide'], ['protein id', 'combined spectral count']],
    'PEAKS': [['-10lgp']],
    'Scaffold': [['exclusive unique peptide count'], ['protein identification probability']],
    'Mascot': [['prot_acc', 'pep_seq'], ['prot_hit_num']],
    'PLGS': [['protein.key', 'protein.accession'], ['peptide.seq', 'protein.accession']],
}
LAB = r'(11[3-9]|121|12[6-9][nclh]?|13[0-5][nclh]?)'   # 127L/127H (isobarQuant) = 127N/127C
CH_PD = re.compile(r'^abundances?(?: \([a-z ]+\))?: (?:f\d+: )?' + LAB + r'(?![0-9.])')
CH_GEN = re.compile(r'(?:tmt|itraq|channel|reporter|ion|signal_sum|sn)[^a-z0-9]{0,3}(?:\d{1,2}plex|pro|plex)?[^a-z0-9]{0,3}' + LAB + r'(?![0-9.])')
CH_EXACT = re.compile(r'^(?:tmt[a-z0-9]*[ _-])?' + LAB + r'$')
MQ_REP = re.compile(r'^reporter intensity (?:corrected )?(\d{1,2})(?![0-9.])')
RUNCOL = re.compile(r'\.(raw|d|wiff|wiff2|mzml|dia)$')
PTMCOL = re.compile(r'probabilities$|site ids$|site positions$|score diffs$|^number of |phosphors|ptmrs|site probabilit|ptm\.')


def plex_from_labels(L):
    it = {l for l in L if l in ('113', '114', '115', '116', '117', '118', '119', '121')}
    tmt = L - it
    if len(it) >= 3 and len(tmt) <= 1:
        return 'iTRAQ 8-plex' if it & {'113', '118', '119', '121'} else 'iTRAQ 4-plex'
    if tmt & {'134c', '135n'}:
        return 'TMTpro 18-plex'
    if {l for l in tmt if l[:3] in ('132', '133', '134')}:
        return 'TMTpro 16-plex'
    if '131c' in tmt:
        return 'TMT 11-plex'
    if any(l[-1] in 'nc' for l in tmt):
        return 'TMT 10-plex'
    if len(tmt) >= 6:
        return 'TMT 6-plex'
    return 'isobaric, unclear'


MODS = [  # (plain modification name, regex) first match wins; applied to site-table names and PTM columns
    ('search defaults/artefacts (oxidation, deamidation, alkylation)',
     r'oxidation|deamid|carbamidomethyl|\bnem\b|^nem|carbamyl|pyro|methylthio|unimod_35_|unimod_4_|unimod_7_|_15\.99|_57\.02|_0\.98|ammonia|dehydrat'),
    ('phosphorylation', r'phospho|unimod_21_|79\.96|sty_79'),
    ('SUMO remnant', r'qqtgg|sumo'),
    ('ubiquitin remnant (GlyGly)', r'glygly|gly-?gly|digly|\bgg ?\(k\)|unimod_121_|114\.04|ubiquit'),
    ('lysine acetylation', r'acetyl[^a-z]*k([^a-z]|sites|$)|^ac ?\(k\)|unimod_1_|k_42\.01'),
    ('glycosylation', r'hexnac|hex ?\(|glcnac|galnac|glyco|n_1[0-9]{3}\.|unimod_43_'),
    ('other lysine acylations', r'succ|croton|crot ?\(|crot\b|lactyl|^la ?\(k\)|2-?ohib|hydroxyisobut|malon|propion|butyr|glutar|benzoyl'),
    ('methylation', r'methyl|unimod_3[467]_'),
    ('ADP-ribosylation', r'adp-?ribos|^adpr'),
    ('citrullination', r'citrul'),
]


def mod_of(s):
    for name, rx in MODS:
        if re.search(rx, s):
            return name
    return 'other'


def parse_header(cols):
    s = [c.lower().strip() for c in cols]
    S = set(s)
    tools = [t for t, alts in HEADER_SIG.items() if any(all(x in S for x in alt) for alt in alts)]
    labels = set()
    for c in s:
        for rx in (CH_PD, CH_GEN, CH_EXACT):
            m = rx.search(c)
            if m:
                labels.add(m.group(1).replace('l', 'n').replace('h', 'c')); break
    mq = {int(m.group(1)) for c in s if (m := MQ_REP.search(c))}
    lfq_cols = [c for c in s if c.startswith('lfq intensity ')]
    runs = sum(bool(RUNCOL.search(c)) for c in s)
    spn_runs = sum(bool(re.match(r'^\[\d+\] .*\.(pg|pep|eg)\.', c)) for c in s)
    ptm = {mod_of(c) for c in s if PTMCOL.search(c)}
    ptm |= {'phosphorylation' for c in s if c.startswith(('phosphors', 'ptmrs')) and 'phospho' in c}
    ptm.discard('search defaults/artefacts (oxidation, deamidation, alkylation)')
    ptm.discard('other')
    pd_ptmrs = any(c.startswith(('ptmrs', 'phosphors')) for c in s)
    return dict(
        h_tools=tools,
        h_iso_labels=len(labels) >= 4, h_iso_n=len(labels) if len(labels) >= 4 else 0,
        h_iso_kit=plex_from_labels(labels) if len(labels) >= 4 else None,
        h_mq_reporter_n=len(mq) if len(mq) >= 2 else 0,
        h_lfq=bool(lfq_cols) or any('maxlfq' in c for c in s), h_lfq_n=len({re.sub(r'___\d+$', '', c) for c in lfq_cols}),
        h_ibaq=any(c.startswith('ibaq') for c in s),
        h_silac=any(re.match(r'^ratio [hml]/[hml]', c) for c in s) or any(re.search(r'abundance ratio.*\((heavy|light|medium)\)', c) for c in s)
                or any(re.match(r'^abundances?(?: \([a-z ]+\))?: (f\d+: )?(heavy|light|medium)\b', c) for c in s),
        h_spc=any('spectral count' in c for c in s),
        h_runs=max(runs, spn_runs, len(lfq_cols)),
        h_ptm=sorted(ptm), h_ptmrs=pd_ptmrs)


H = pd.DataFrame([json.loads(l) for l in open(HS / 'Headers.jsonl')]).drop_duplicates('url', keep='last')
H = H[H.url.isin(pd.read_csv(HS / 'Targets.csv.gz', usecols=['url']).url)]
H = H[H.status.isin(['ok', 'partial']) & (H['columns'].str.len() > 0)].reset_index(drop=True)
P = pd.DataFrame([parse_header(c) for c in H['columns']]); P.insert(0, 'accession', H.accession); P.insert(1, 'file_name', H.file_name)
SUMMARY['headers'] = {'tables_read': len(P), 'datasets_with_a_read_table': int(P.accession.nunique())}
P['h_iso'] = P.h_iso_labels | (P.h_mq_reporter_n >= 2)
P['h_iso_channels'] = np.where(P.h_iso_labels, P.h_iso_n, P.h_mq_reporter_n)
HD = P.groupby('accession').agg(h_iso=('h_iso', 'any'), h_iso_channels=('h_iso_channels', 'max'), h_lfq=('h_lfq', 'any'),
                                h_ibaq=('h_ibaq', 'any'), h_silac=('h_silac', 'any'), h_spc=('h_spc', 'any'), h_runs=('h_runs', 'max'),
                                h_ptmrs=('h_ptmrs', 'any'),
                                h_iso_kit=('h_iso_kit', lambda x: Counter(x.dropna()).most_common(1)[0][0] if x.notna().any() else None),
                                h_ptm=('h_ptm', lambda x: sorted(set(sum(x, [])))), h_tools=('h_tools', lambda x: sorted(set(sum(x, [])))))
D = D.merge(HD, left_on='accession', right_index=True, how='left')
D['header_read'] = D.accession.isin(HD.index)
for c in ['h_iso', 'h_lfq', 'h_ibaq', 'h_silac', 'h_spc', 'h_ptmrs']:
    D[c] = D[c].fillna(False).astype(bool)
for c in ['h_ptm', 'h_tools']:
    D[c] = D[c].apply(lambda v: v if isinstance(v, list) else [])
# spot-check sheet: example parsed tables per signal
spot = []
for col in ['h_iso_labels', 'h_mq_reporter_n', 'h_lfq', 'h_silac', 'h_ibaq']:
    x = P[P[col].astype(bool)].sample(min(8, int(P[col].astype(bool).sum())), random_state=SEED)
    for i, r in x.iterrows():
        spot.append({'signal': col, 'accession': r.accession, 'file_name': r.file_name, 'kit': r.h_iso_kit, 'channels': r.h_iso_channels,
                     'matched_columns': ' | '.join(c for c in H.loc[i, 'columns'] if re.search(
                         {'h_iso_labels': '(11[3-9]|121|12[6-9]|13[0-5])', 'h_mq_reporter_n': '(?i)reporter intensity', 'h_lfq': '(?i)lfq|maxlfq',
                          'h_silac': '(?i)ratio|heavy|light', 'h_ibaq': '(?i)ibaq'}[col], c))[:600]})
pd.DataFrame(spot).to_csv(OUT / 'H5-header-parse-spotcheck.csv', index=False)

# =============================================================================================
# 3. Site tables (H5)
# =============================================================================================
site_rx = OTHER_PAT['site_table']
ST = pd.concat([con.execute(f"""select distinct accession, '{lab}' as source, regexp_extract(p, '[^/]*$') as name from {src}
                                where regexp_matches(p, {sq(site_rx)})""").df() for src, lab in [('top', 'top-level'), ('zm', 'inside zip')]])
ST['modification'] = ST.name.map(mod_of)
ST.loc[ST.name.str.contains(r'^(abundance|ratio)_(single|multi)-site'), 'modification'] = 'site quantification, modification not named'
ST.groupby(['name', 'modification']).agg(datasets=('accession', 'nunique')).reset_index().sort_values('datasets', ascending=False) \
    .to_csv(OUT / 'H5-site-table-names.csv', index=False)
site_mods = ST.groupby('accession').modification.agg(lambda x: sorted(set(x)))
D['site_mods'] = D.accession.map(site_mods).apply(lambda v: v if isinstance(v, list) else [])
D['any_site_table'] = D.accession.isin(ST.accession)
D['mq_sites_visible'] = FT.mq_oxidation_sites.to_numpy() | FZ.mq_oxidation_sites.to_numpy()   # MaxQuant writes Oxidation (M)Sites.txt by default
D['mq_core'] = FT.mq_core_txt.to_numpy() | FZ.mq_core_txt.to_numpy()

# metadata PTM text
PTM_TEXT = {
    'phosphorylation': (r'phospho', r'phospho', r'tio2|titanium dioxide|\bimac\b|fe-?nta|fe\(iii\)|fe3\+|ti4\+|ti\(iv\)|ti-?imac|zr-?imac|phospho-?(peptide )?enrich|anti-?phospho|phosphotyrosine|p-?tyr|ptmscan'),
    'ubiquitin remnant (GlyGly)': (r'ubiquitinom|ubiquitylom|ubiquitin remnant|ubiquitination site|ubiquitylation site|digly|di-gly|gly-gly|glygly|k-?[εe]-?gg|\bkgg\b', r'glygly|gly-gly|ubiquitin', r'k-?[εe]-?gg|digly|di-gly|ubiquitin remnant|ptmscan ubiquit'),
    'lysine acetylation': (r'acetylom|acetyl-?lysine|lysine acetyl|\bkac\b|acetylation', r'n6-acetyl|acetyl-?l?-?lysine|acetyllysine', r'acetyl-?lysine antibod|anti-?acetyl|acetyl-?lysine enrich|ptmscan acetyl'),
    'glycosylation': (r'glycoproteom|glycopeptide|glycosylat|glycan|glcnac|galnac|n-glyco|o-glyco|sialyl', r'glyco|hexnac|hex\(|glcnac', r'hilic|lectin|glycopeptide enrich|zic-hilic|hydrazide|pngase'),
}
for k, (dx, sx, ex) in PTM_TEXT.items():
    D[f'design_{k}'] = D.design_text.str.contains(dx, regex=True)
    D[f'struct_{k}'] = D.ptm_struct.str.contains(sx, regex=True)
    D[f'enrich_{k}'] = D.sample_protocol.str.contains(ex, regex=True)
    D[f'site_{k}'] = D.site_mods.apply(lambda v, k=k: k in v)
    D[f'hdr_{k}'] = D.h_ptm.apply(lambda v, k=k: k in v)
D['hdr_phosphorylation'] |= D.h_ptmrs   # PD ptmRS/phosphoRS localisation columns (mostly phospho)

# =============================================================================================
# 4. H5 analyses
# =============================================================================================
# 4a. site tables by modification
rows = []
mq_vis = D[D.mq_sites_visible]
for mod in [m for m, _ in MODS] + ['other', 'site quantification, modification not named']:
    has = D.site_mods.apply(lambda v: mod in v)
    rows.append({'modification': mod, 'datasets_with_site_table': int(has.sum()),
                 'top_level': int(D.accession.isin(ST[(ST.modification == mod) & (ST.source == 'top-level')].accession).sum()),
                 'inside_zip': int(D.accession.isin(ST[(ST.modification == mod) & (ST.source == 'inside zip')].accession).sum()),
                 'pct_of_datasets_with_MaxQuant_site_tables_visible': round(100 * mq_vis.site_mods.apply(lambda v: mod in v).mean(), 1)})
pd.DataFrame(rows).to_csv(OUT / 'H5-site-tables-by-modification.csv', index=False)
SUMMARY['H5_site_tables'] = {'datasets_with_any_site_table': int(D.any_site_table.sum()),
                             'datasets_with_MaxQuant_Oxidation_M_Sites_visible': int(D.mq_sites_visible.sum()),
                             'datasets_with_MaxQuant_proteinGroups_or_evidence_visible': int(D.mq_core.sum())}

# 4b. agreement between site tables and metadata, within MaxQuant datasets whose site tables are visible
agree = []
for k in PTM_TEXT:
    x = D[D.mq_sites_visible].copy()
    site = x[f'site_{k}']
    for meta in ['design', 'struct', 'enrich']:
        m = x[f'{meta}_{k}']
        agree.append({'modification': k, 'population': 'MaxQuant output with site tables visible', 'denominator': len(x),
                      'metadata_source': {'design': 'title/description/keywords', 'struct': 'PRIDE identified-PTM field',
                                          'enrich': 'enrichment named in sample protocol'}[meta],
                      'site_table_and_metadata': int((site & m).sum()), 'site_table_only': int((site & ~m).sum()),
                      'metadata_only': int((~site & m).sum()), 'neither': int((~site & ~m).sum()),
                      'pct_metadata_given_site_table': round(100 * m[site].mean(), 1) if site.any() else np.nan,
                      'pct_site_table_given_metadata': round(100 * site[m].mean(), 1) if m.any() else np.nan})
    # like-for-like stricter: site table AND an enrichment term vs design text
pd.DataFrame(agree).to_csv(OUT / 'H5-ptm-site-vs-metadata.csv', index=False)
# bootstrap CIs for the key phospho comparisons
x = D[D.mq_sites_visible].assign(value=lambda d: d['design_phosphorylation'].astype(float))
ci_rows = [rate_row(x.site_phosphorylation, x, 'design text mentions phospho | Phospho site table present'),
           rate_row(~x.site_phosphorylation, x, 'design text mentions phospho | no Phospho site table')]
x = x.assign(value=lambda d: d['enrich_phosphorylation'].astype(float))
ci_rows += [rate_row(x.site_phosphorylation, x, 'enrichment named in sample protocol | Phospho site table present'),
            rate_row(~x.site_phosphorylation, x, 'enrichment named in sample protocol | no Phospho site table')]
x = x.assign(value=lambda d: d['site_phosphorylation'].astype(float))
ci_rows += [rate_row(x.design_phosphorylation & x.enrich_phosphorylation, x, 'Phospho site table | design text AND enrichment named'),
            rate_row(~x.design_phosphorylation & ~x.enrich_phosphorylation & ~x.struct_phosphorylation, x,
                     'Phospho site table | no phospho in design text, structured PTMs or enrichment')]
pd.DataFrame(ci_rows).to_csv(OUT / 'H5-phospho-agreement-CIs.csv', index=False)

# 4c. strategy by year: metadata (all datasets) and file evidence (with its own denominators)
D['itraq_text'] = D.all_text.str.contains(r'itraq', regex=True)
D['isobaric_meta'] = D.tmt_expanded | D.itraq_text
D['silac_meta'] = D.silac_expanded | D.design_text.str.contains('silac')
D['any_site_phospho_file'] = D.site_phosphorylation | D.hdr_phosphorylation
yr = []
for y, x in D.groupby('release_year'):
    mq = x[x.mq_sites_visible]; hd = x[x.header_read]; rn = x[x.has_runs]
    yr.append({'release_year': y, 'datasets': len(x),
               'meta_tmt_expanded_pct': 100 * x.tmt_expanded.mean(), 'meta_itraq_text_pct': 100 * x.itraq_text.mean(),
               'meta_silac_expanded_pct': 100 * x.silac_expanded.mean(), 'meta_lfq_expanded_pct': 100 * x.lfq_expanded.mean(),
               'meta_dia_expanded_pct': 100 * x.dia_expanded.mean(),
               **{f'design_{k}_pct': 100 * x[f'design_{k}'].mean() for k in PTM_TEXT},
               'struct_phospho_pct': 100 * x.struct_phosphorylation.mean(),
               'mq_sites_visible_n': len(mq), **{f'mq_site_{k}_pct': 100 * mq[f'site_{k}'].mean() if len(mq) else np.nan for k in PTM_TEXT},
               'header_read_n': len(hd), 'hdr_isobaric_pct': 100 * hd.h_iso.mean() if len(hd) else np.nan,
               'hdr_lfq_pct': 100 * hd.h_lfq.mean() if len(hd) else np.nan, 'hdr_silac_pct': 100 * hd.h_silac.mean() if len(hd) else np.nan,
               'hdr_ibaq_pct': 100 * hd.h_ibaq.mean() if len(hd) else np.nan,
               'with_raw_runs_n': len(rn), 'run_names_isobaric_pct': 100 * rn.run_isobaric.mean() if len(rn) else np.nan,
               'run_names_phospho_pct': 100 * rn.run_phospho.mean() if len(rn) else np.nan})
pd.DataFrame(yr).round(1).to_csv(OUT / 'H5-strategy-by-year.csv', index=False)

# 4d. isobaric: header/run-name evidence vs metadata; plex from text vs headers
q = []
for lab, pop, ev in [('header read', D.header_read, 'h_iso'), ('raw runs listed', D.has_runs, 'run_isobaric')]:
    x = D[pop]
    e, m = x[ev], x.isobaric_meta
    q.append({'evidence': 'isobaric reporter columns' if ev == 'h_iso' else 'TMT/iTRAQ in raw run names', 'population': lab,
              'denominator': len(x), 'evidence_and_metadata': int((e & m).sum()), 'evidence_only': int((e & ~m).sum()),
              'metadata_only': int((~e & m).sum()), 'neither': int((~e & ~m).sum())})
for lab, ev, meta in [('SILAC ratio/heavy-light columns', 'h_silac', 'silac_meta'), ('LFQ intensity / MaxLFQ columns', 'h_lfq', 'lfq_expanded'),
                      ('iBAQ columns', 'h_ibaq', 'lfq_expanded')]:
    x = D[D.header_read]; e, m = x[ev], x[meta]
    q.append({'evidence': lab, 'population': 'header read', 'metadata_flag': meta, 'denominator': len(x),
              'evidence_and_metadata': int((e & m).sum()), 'evidence_only': int((e & ~m).sum()),
              'metadata_only': int((~e & m).sum()), 'neither': int((~e & ~m).sum())})
D['lfq_wide'] = D.all_text.str.contains(r'label.free|\blfq\b|maxlfq|ibaq|intensity.based absolute', regex=True)
for lab, ev in [('LFQ intensity / MaxLFQ columns', 'h_lfq'), ('iBAQ columns', 'h_ibaq')]:
    x = D[D.header_read]; e, m = x[ev], x.lfq_wide
    q.append({'evidence': lab, 'population': 'header read', 'metadata_flag': 'label-free/LFQ/MaxLFQ/iBAQ anywhere in text (wider)',
              'denominator': len(x), 'evidence_and_metadata': int((e & m).sum()), 'evidence_only': int((e & ~m).sum()),
              'metadata_only': int((~e & m).sum()), 'neither': int((~e & ~m).sum())})
# MaxQuant processing steps vs lfq flag (datasets with a MaxQuant proc folder)
proc = FT.mq_proc_any.to_numpy() | FZ.mq_proc_any.to_numpy()
D['mq_proc'] = proc; D['mq_proc_lfq'] = FT.mq_proc_lfq.to_numpy() | FZ.mq_proc_lfq.to_numpy()
D['mq_proc_ibaq'] = FT.mq_proc_ibaq.to_numpy() | FZ.mq_proc_ibaq.to_numpy()
x = D[D.mq_proc]; e, m = x.mq_proc_lfq, x.lfq_expanded
q.append({'evidence': 'MaxQuant ran label-free quantification (proc/ step files)', 'population': 'MaxQuant proc folder deposited',
          'metadata_flag': 'lfq_expanded', 'denominator': len(x), 'evidence_and_metadata': int((e & m).sum()),
          'evidence_only': int((e & ~m).sum()), 'metadata_only': int((~e & m).sum()), 'neither': int((~e & ~m).sum())})
QA = pd.DataFrame(q)
QA['pct_metadata_given_evidence'] = (100 * QA.evidence_and_metadata / (QA.evidence_and_metadata + QA.evidence_only)).round(1)
QA['pct_evidence_given_metadata'] = (100 * QA.evidence_and_metadata / (QA.evidence_and_metadata + QA.metadata_only)).round(1)
QA.to_csv(OUT / 'H5-quant-evidence-vs-metadata.csv', index=False)

# text plex
TPLEX = [(r'tmt ?pro[^a-z0-9]{0,3}(18|eighteen)|18[ -]?plex|tmt18|tmt-18|tmt 18', 18), (r'tmt ?pro[^a-z0-9]{0,3}16|16[ -]?plex|tmt16|tmt-16|tmt 16', 16),
         (r'tmt ?pro[^a-z0-9]{0,3}(32|35)|3[25][ -]?plex', 35), (r'11[ -]?plex|tmt11|tmt-11|tmt 11\b', 11), (r'10[ -]?plex|tmt10|tmt-10|tmt 10\b', 10),
         (r'6[ -]?plex|tmt6|tmt-6|tmt 6\b|sixplex', 6), (r'itraq[^.]{0,25}8[ -]?plex|8[ -]?plex[^.]{0,25}itraq|itraq ?8', 8),
         (r'itraq[^.]{0,25}4[ -]?plex|4[ -]?plex[^.]{0,25}itraq|itraq ?4', 4), (r'(duplex|2[ -]?plex)[^.]{0,25}tmt|tmt ?(duplex|2[ -]?plex)', 2)]


def text_plex(t):
    """Highest TMT plex stated (the Unimod label name 'TMT6plex' is also used for 10/11-plex reagents, so 6 is
    kept only when nothing higher is stated); iTRAQ plex if only iTRAQ is stated."""
    found = sorted({n for rx, n in TPLEX if re.search(rx, t)})
    tmt = [n for n in found if n in (2, 6, 10, 11, 16, 18, 35)]
    itq = [n for n in found if n in (4, 8)]
    if tmt and itq and re.search('itraq', t) and re.search(r'tmt|tandem mass tag', t):
        return 'TMT and iTRAQ both stated', found
    if tmt:
        return f'TMT {max(tmt)}', found
    if itq:
        return f'iTRAQ {max(itq)}', found
    if re.search(r'tmt ?pro', t):
        return 'TMTpro, plex not stated', found
    return None, found


tp = D.all_text.apply(text_plex)
D['text_plex'] = tp.str[0]; D['text_plex_set'] = tp.str[1]
iso = D[D.isobaric_meta]
pt = pd.crosstab(iso.release_year, iso.text_plex.fillna('plex not stated'))
pt['isobaric_datasets'] = iso.groupby('release_year').size()
pt.to_csv(OUT / 'H5-isobaric-plex-text-by-year.csv')
hp = D[D.h_iso]
hk = pd.crosstab(hp.release_year, hp.h_iso_channels.astype(int).astype(str) + ' channels')
hk['datasets_with_reporter_columns'] = hp.groupby('release_year').size()
hk.to_csv(OUT / 'H5-isobaric-channels-header-by-year.csv')
both = D[D.h_iso & D.text_plex_set.str.len().gt(0)].copy()
# header plex: reagent kit inferred from channel labels (PD/other tables) when available, else MaxQuant channel count
KITN = {'iTRAQ 4-plex': 'iTRAQ 4', 'iTRAQ 8-plex': 'iTRAQ 8', 'TMT 6-plex': 'TMT 6', 'TMT 10-plex': 'TMT 10', 'TMT 11-plex': 'TMT 11',
        'TMTpro 16-plex': 'TMT 16', 'TMTpro 18-plex': 'TMT 18'}
both['header_plex'] = [KITN.get(k) if isinstance(k, str) and k in KITN else f'{int(c)} channels (MaxQuant or unclear)'
                       for k, c in zip(both.h_iso_kit, both.h_iso_channels)]
num = lambda v: int(re.search(r'(\d+)', v).group(1)) if isinstance(v, str) and re.search(r'\d+', v) else np.nan
both['hp'] = both.header_plex.map(num); both['tp'] = both.text_plex.map(num)
both['match'] = both.hp == both.tp
both['header_lt_text'] = both.hp < both.tp
both['header_gt_text'] = both.hp > both.tp
pd.crosstab(both.header_plex, both.text_plex, margins=True).to_csv(OUT / 'H5-isobaric-plex-text-vs-header.csv')
both[['accession', 'release_year', 'h_iso_kit', 'h_iso_channels', 'header_plex', 'text_plex', 'match']].to_csv(OUT / 'H5-isobaric-plex-text-vs-header-datasets.csv', index=False)
SUMMARY['H5_isobaric'] = {
    'datasets_isobaric_metadata': int(D.isobaric_meta.sum()), 'of_which_plex_stated_in_text': int(iso.text_plex.notna().sum()),
    'datasets_with_reporter_columns_in_headers': int(D.h_iso.sum()), 'header_and_text_plex_both': len(both),
    'header_plex_equals_text_plex_pct': round(100 * both.match.mean(), 1), 'header_below_text_pct': round(100 * both.header_lt_text.mean(), 1),
    'header_above_text_pct': round(100 * both.header_gt_text.mean(), 1),
    'kit_from_labels_only': {'n': int(both.h_iso_kit.notna().sum()), 'match_pct': round(100 * both[both.h_iso_kit.notna()].match.mean(), 1)},
    'maxquant_channel_count_only': {'n': int(both.h_iso_kit.isna().sum()), 'match_pct': round(100 * both[both.h_iso_kit.isna()].match.mean(), 1)}}
# by-period plex shares among datasets with a stated plex (text)
st = iso[iso.text_plex.notna() & (iso.text_plex != 'TMT and iTRAQ both stated')]
pp = pd.crosstab(st.period, st.text_plex, normalize='index').mul(100).round(1); pp['n'] = st.groupby('period').size()
pp.to_csv(OUT / 'H5-isobaric-plex-text-by-period.csv')
# samples per quantitative matrix (header tables): median runs/sample columns by year and by kind
sm = P[P.h_runs >= 2].merge(D[['accession', 'release_year', 'dia_expanded']], on='accession')
sm['kind'] = np.where(sm.h_lfq_n >= 2, 'MaxQuant LFQ matrix', np.where(sm.dia_expanded, 'per-run matrix, DIA dataset', 'per-run matrix, other'))
smy = sm.groupby(['accession', 'release_year', 'kind']).h_runs.max().reset_index()
smy.groupby(['kind', pd.cut(smy.release_year, PERIOD_CUTS, labels=PERIOD_LABS)], observed=True).h_runs \
    .describe(percentiles=[.25, .5, .75]).round(1).to_csv(OUT / 'H5-samples-per-matrix-by-period.csv')

# =============================================================================================
# 5. H6 declared vs evidenced software
# =============================================================================================
# Registry text regexes used as-is, except two that match non-software phrases in this cohort (checked by reading contexts):
#   EncyclopeDIA: 'encyclopedia' is almost always 'Kyoto Encyclopedia of Genes and Genomes' (KEGG) or the Cancer Cell Line Encyclopedia
#   Sage: '\bsage\b' mostly matches 'Sage-N Research' (vendor of SORCERER/SEQUEST)
# Sensitivity: spelling variants seen in protocols of evidenced-but-unnamed datasets ('max quant', 'protein discoverer', 'pd1.3')
VARIANT = {'MaxQuant': r'max ?quant|andromeda', 'Proteome Discoverer': r'prot(eome|ein)s? ?discover|\bpd ?v?\d\.\d',
           'FragPipe': r'ms ?fragger|frag ?pipe', 'Mascot': r'mascot', 'DIA-NN': r'dia[-_ ]?nn', 'Spectronaut': r'spectronaut'}
TEXT_OVERRIDE = {'EncyclopeDIA': r'(?<!kyoto )(?<!line )encyclopedia(?! of)', 'Sage': r'\bsage\b(?! group)(?!-n)'}
nm_proc = {t: D.processing.str.contains(TOOLS[t]['text'], regex=True) for t in ['MaxQuant', 'Proteome Discoverer', 'Mascot', 'DIA-NN', 'Spectronaut', 'FragPipe']}
nm_sw = {t: D.software.str.contains(TOOLS[t]['text'], regex=True) for t in nm_proc}
SUMMARY['H6_structured_field_vs_protocol'] = {t: {'protocol': int(nm_proc[t].sum()), 'structured_field': int(nm_sw[t].sum()),
                                                  'protocol_not_structured': int((nm_proc[t] & ~nm_sw[t]).sum()),
                                                  'structured_not_protocol': int((nm_sw[t] & ~nm_proc[t]).sum())} for t in nm_proc}
SUMMARY['H6_text_regex_overrides'] = {t: {'registry_matches': int(D.software.str.cat(D.processing, sep=' ').str.contains(TOOLS[t]['text'], regex=True).sum()),
                                          'override_matches': int(D.software.str.cat(D.processing, sep=' ').str.contains(r, regex=True).sum())}
                                      for t, r in TEXT_OVERRIDE.items()}
ev_rows, yr_rows = [], []
for t in SIG:
    trx = TEXT_OVERRIDE.get(t, TOOLS[t]['text'])
    nm = D.software.str.cat(D.processing, sep=' ').str.contains(trx, regex=True)   # as common.software_flags
    wide = D.all_text.str.contains(trx, regex=True)
    ev_top, ev_zip = FT['tool:' + t].to_numpy(), FZ['tool:' + t].to_numpy()
    ev_hdr = D.h_tools.apply(lambda v: t in v).to_numpy()
    ev = ev_top | ev_zip | ev_hdr
    D[f'named_{t}'] = nm.to_numpy(); D[f'ev_{t}'] = ev; D[f'wide_{t}'] = wide.to_numpy()
for t in SIG:
    D[f'ev_other_than_{t}'] = D[[f'ev_{u}' for u in SIG if u != t]].any(axis=1)
D['ev_any_tool'] = D[[f'ev_{t}' for t in SIG]].any(axis=1)
D['n_tools_evidenced'] = D[[f'ev_{t}' for t in SIG]].sum(axis=1)
D['n_tools_named'] = D[[f'named_{t}' for t in SIG]].sum(axis=1)
for t in SIG:
    ev, nm, wide = D[f'ev_{t}'], D[f'named_{t}'], D[f'wide_{t}']
    evnn = D[ev & ~nm]
    insteadnamed = Counter(sum([[u for u in SIG if r[f'named_{u}']] for _, r in evnn.iterrows()], []))
    nne = D[nm & ~ev]
    other_ev = Counter(sum([[u for u in SIG if r[f'ev_{u}']] for _, r in nne.iterrows()], []))
    ev_rows.append({'tool': t, 'named': int(nm.sum()), 'evidenced': int(ev.sum()),
                    'evidenced_top_level': int(FT['tool:' + t].sum()), 'evidenced_in_zip': int(FZ['tool:' + t].sum()),
                    'evidenced_by_header': int(D.h_tools.apply(lambda v: t in v).sum()),
                    'both': int((ev & nm).sum()), 'evidenced_not_named': int((ev & ~nm).sum()),
                    'pct_evidenced_not_named': round(100 * (ev & ~nm).sum() / ev.sum(), 1) if ev.sum() else np.nan,
                    'pct_evidenced_not_named_anywhere_in_text': round(100 * (ev & ~wide).sum() / ev.sum(), 1) if ev.sum() else np.nan,
                    'pct_evidenced_not_named_spelling_variants_allowed': round(100 * (ev & ~nm & ~D.all_text.str.contains(VARIANT[t], regex=True)).sum() / ev.sum(), 1)
                    if t in VARIANT and ev.sum() else np.nan,
                    'evidenced_not_named_but_no_tool_named': int((ev & ~nm & (D.n_tools_named == 0)).sum()),
                    'evidenced_not_named_most_often_named_instead': '; '.join(f'{k} {v}' for k, v in insteadnamed.most_common(3)),
                    'named_not_evidenced': int((nm & ~ev).sum()),
                    'pct_named_not_evidenced': round(100 * (nm & ~ev).sum() / nm.sum(), 1) if nm.sum() else np.nan,
                    'named_not_evidenced_but_other_tool_evidenced': int((nm & ~ev & D[f'ev_other_than_{t}']).sum()),
                    'named_not_evidenced_other_tools_evidenced': '; '.join(f'{k} {v}' for k, v in other_ev.most_common(3)),
                    'export': TOOLS[t]['export'],
                    'note': 'pepXML/protXML are TPP formats but are also written by Proteome Discoverer exports and FragPipe (Philosopher); '
                            'file evidence here means "TPP-format files", not TPP software' if t == 'TPP' else
                            ('text regex overridden (see Summary.json)' if t in TEXT_OVERRIDE else '')})
EV = pd.DataFrame(ev_rows).sort_values('evidenced', ascending=False)
EV.to_csv(OUT / 'H6-declared-vs-evidenced-by-tool.csv', index=False)
MAIN6 = ['MaxQuant', 'Proteome Discoverer', 'Mascot', 'DIA-NN', 'Spectronaut', 'FragPipe', 'PEAKS', 'Skyline']
for t in MAIN6:
    x = D[D[f'ev_{t}']].assign(value=lambda d: (~d[f'named_{t}']).astype(float))
    for y, g in x.groupby('release_year'):
        yr_rows.append({'tool': t, 'release_year': y, 'evidenced': len(g), 'pct_evidenced_not_named': round(100 * g.value.mean(), 1)})
    for p, g in x.groupby('period'):
        yr_rows.append(dict(rate_row(None, g, 'period'), tool=t, release_year=p, evidenced=len(g)))
    xn = D[D[f'named_{t}']].assign(value=lambda d: d[f'ev_{t}'].astype(float))
    for p, g in xn.groupby('period'):
        r = rate_row(None, g, 'period, evidenced given named'); r.update(tool=t, release_year=p, named=len(g)); yr_rows.append(r)
    for p_, g in x.groupby('period'):
        for ann, gg in g.groupby('has_software_annotation'):
            r = rate_row(None, gg, f'period, structured software field {"filled" if ann else "empty"}'); r.update(tool=t, release_year=p_, evidenced=len(gg))
            yr_rows.append(r)
pd.DataFrame(yr_rows).to_csv(OUT / 'H6-evidenced-not-named-by-year.csv', index=False)
sa = D.groupby('release_year').agg(datasets=('accession', 'size'), structured_software_field_filled_pct=('has_software_annotation', 'mean'),
                                   no_tool_named_pct=('family', lambda f: (f == 'No tool named').mean()))
sa[['structured_software_field_filled_pct', 'no_tool_named_pct']] *= 100
sa.round(1).to_csv(OUT / 'H6-software-annotation-by-year.csv')
# DIA-NN evidenced but not named: is FragPipe (which runs DIA-NN internally) named or evidenced?
dn = D[D['ev_DIA-NN'] & ~D['named_DIA-NN']]
SUMMARY['H6_diann_not_named'] = {'n': len(dn), 'FragPipe_named_pct': round(100 * dn.named_FragPipe.mean(), 1),
                                 'FragPipe_evidenced_pct': round(100 * dn.ev_FragPipe.mean(), 1),
                                 'fragpipe_diann_output_folder_pct': round(100 * (FT.fragpipe_diann_folder | FZ.fragpipe_diann_folder)[dn.accession].mean(), 1),
                                 'Spectronaut_named_pct': round(100 * dn.named_Spectronaut.mean(), 1),
                                 'no_tool_named_pct': round(100 * (dn.n_tools_named == 0).mean(), 1)}
# pipelines: tool sets evidenced in files
ev_sets = D[D.ev_any_tool].apply(lambda r: ' + '.join(t if t != 'TPP' else 'TPP-format files' for t in SIG if r[f'ev_{t}']), axis=1)
combos = ev_sets[D.loc[ev_sets.index, 'n_tools_evidenced'] >= 2].value_counts()
combos.rename('datasets').to_frame().to_csv(OUT / 'H6-multi-tool-evidence-combinations.csv')
pairs = Counter()
for s in ev_sets[D.loc[ev_sets.index, 'n_tools_evidenced'] >= 2]:
    ts = s.split(' + ')
    for i in range(len(ts)):
        for j in range(i + 1, len(ts)):
            pairs[(ts[i], ts[j])] += 1
pd.DataFrame([{'tool_a': a, 'tool_b': b, 'datasets': n} for (a, b), n in pairs.most_common()]).to_csv(OUT / 'H6-tool-pairs-evidenced.csv', index=False)
D['fp_diann_folder'] = (FT.fragpipe_diann_folder | FZ.fragpipe_diann_folder).to_numpy()
mt = []
for y, x in D[D.ev_any_tool].groupby('release_year'):
    mt.append({'release_year': y, 'datasets_with_evidenced_tool': len(x), 'pct_two_or_more_tools': round(100 * (x.n_tools_evidenced >= 2).mean(), 1),
               'pct_FragPipe_and_DIA-NN': round(100 * (x.ev_FragPipe & x['ev_DIA-NN']).mean(), 2),
               'n_FragPipe_diann-output_folder': int(x.fp_diann_folder.sum()),
               'pct_two_or_more_tools_excl_TPP_Percolator': round(100 * (x[[f'ev_{t}' for t in SIG if t not in ('TPP', 'Percolator')]].sum(axis=1) >= 2).mean(), 1)})
pd.DataFrame(mt).to_csv(OUT / 'H6-tools-per-dataset-by-year.csv', index=False)
# 'No tool named' datasets: what is evidenced
nt = D[D.family == 'No tool named']
SUMMARY['H6_no_tool_named'] = {'datasets': len(nt), 'with_evidenced_tool': int(nt.ev_any_tool.sum()),
                               'pct': round(100 * nt.ev_any_tool.mean(), 1),
                               'top_evidenced': {t: int(nt[f'ev_{t}'].sum()) for t in SIG if nt[f'ev_{t}'].sum() >= 20}}

# =============================================================================================
# 6. H8 documentation
# =============================================================================================
D['desc_chars'] = D.description_raw.str.len(); D['desc_words'] = D.description_raw.str.split().str.len()
D['title_chars'] = D.title.str.len(); D['sample_protocol_chars'] = D.sample_protocol.str.len(); D['data_protocol_chars'] = D.data_protocol.str.len()
for k in ['doc_readme', 'doc_sdrf', 'doc_design_tool', 'doc_design_author', 'doc_params', 'doc_document']:
    D[k + '_top'] = FT[k].to_numpy(); D[k + '_zip'] = FZ[k].to_numpy(); D[k] = D[k + '_top'] | D[k + '_zip']
D['doc_document'] = D.doc_document_top   # inside zips, documents are mostly tool-written (MaxQuant tables.pdf in 5,369 datasets)
# every zipped 'sdrf.tsv' found sits next to FragPipe outputs (fragpipe.workflow/manifest/combined_protein.tsv): FragPipe-written
D['doc_sdrf'] = D.doc_sdrf_top
SUMMARY['H8_sdrf_inside_zips'] = {'datasets': int(D.doc_sdrf_zip.sum()), 'zip_only': int((D.doc_sdrf_zip & ~D.doc_sdrf_top).sum()),
                                  'note': 'zipped sdrf.tsv: 47 of these datasets have it inside FragPipe output folders (FragPipe-written); SDRF counted top-level only'}
D['sdrf_frozen_layer'] = F.set_index('accession').fmt_sdrf.reindex(D.accession).fillna(False).astype(bool).to_numpy()
dup = D.description_raw.str.strip().str.lower()
D['desc_duplicated'] = dup.map(dup.value_counts()).gt(1) & dup.str.len().gt(0)
dupg = D[D.desc_duplicated].groupby(dup[D.desc_duplicated]).agg(n=('accession', 'size'), groups=('grp', 'nunique'))
DOCS = ['doc_readme', 'doc_sdrf', 'doc_design_author', 'doc_design_tool', 'doc_params', 'doc_document']
by = D.groupby('release_year').agg(datasets=('accession', 'size'), desc_chars_median=('desc_chars', 'median'),
                                   desc_chars_p25=('desc_chars', lambda x: x.quantile(.25)), desc_chars_p75=('desc_chars', lambda x: x.quantile(.75)),
                                   keywords_median=('n_keywords', 'median'), sample_protocol_chars_median=('sample_protocol_chars', 'median'),
                                   data_protocol_chars_median=('data_protocol_chars', 'median'), title_chars_median=('title_chars', 'median'),
                                   desc_duplicated_pct=('desc_duplicated', 'mean'), **{k + '_pct': (k, 'mean') for k in DOCS},
                                   doc_readme_top_pct=('doc_readme_top', 'mean'), doc_sdrf_top_pct=('doc_sdrf_top', 'mean'),
                                   zip_listed_pct=('zip_listed', 'mean'))
for c in by.columns:
    if c.endswith('_pct'):
        by[c] = 100 * by[c]
by.round(1).to_csv(OUT / 'H8-documentation-by-year.csv')
# minimum-length floors (PRIDE submission requires minimum text lengths)
SUMMARY['H8_lengths'] = {'desc_chars_min': int(D.desc_chars.min()), 'desc_chars_lt_100_pct': round(100 * (D.desc_chars < 100).mean(), 1),
                         'sample_protocol_lt_100_pct': round(100 * (D.sample_protocol_chars < 100).mean(), 1),
                         'data_protocol_lt_100_pct': round(100 * (D.data_protocol_chars < 100).mean(), 1),
                         'descriptions_shared_with_another_dataset_pct': round(100 * D.desc_duplicated.mean(), 1),
                         'shared_description_texts': len(dupg), 'shared_texts_spanning_2plus_groups': int((dupg.groups > 1).sum()),
                         'sdrf_bounded_rule_vs_frozen_layer_top': {'bounded_rule': int(D.doc_sdrf_top.sum()), 'frozen_fmt_sdrf': int(D.sdrf_frozen_layer.sum())}}
# relation to results form: documentation by form (2019-2025), and form by documentation, year-standardized
R19 = D[D.release_year.between(2019, 2025) & D.form.notna()].copy()
R19['desc_q'] = R19.groupby('release_year').desc_chars.transform(lambda x: pd.qcut(x.rank(method='first'), 4, labels=['Q1 shortest', 'Q2', 'Q3', 'Q4 longest'])).astype(str)
bf = R19.groupby('form').agg(datasets=('accession', 'size'), desc_chars_median=('desc_chars', 'median'), keywords_median=('n_keywords', 'median'),
                             sample_protocol_chars_median=('sample_protocol_chars', 'median'), data_protocol_chars_median=('data_protocol_chars', 'median'),
                             **{k + '_pct': (k, 'mean') for k in DOCS})
for c in bf.columns:
    if c.endswith('_pct'):
        bf[c] = (100 * bf[c]).round(1)
bf.to_csv(OUT / 'H8-documentation-by-results-form-2019-2025.csv')
w = R19.release_year.value_counts(normalize=True)


def std_rate(x, outcome, mask):
    """Year-standardized rate (2019-2025 year weights) of `outcome` among rows where mask; group bootstrap."""
    xs = x[mask]
    rates = xs.groupby('release_year')[outcome].mean()
    est = 100 * (rates * w).sum() / w[rates.index].sum()
    gi = xs.groupby('grp').indices; gs = np.array(list(gi)); rng = np.random.default_rng(SEED)
    yv, ov = xs.release_year.to_numpy(), xs[outcome].to_numpy(float)
    bs = []
    for _ in range(500):
        idx = np.concatenate([gi[g] for g in gs[rng.integers(0, len(gs), len(gs))]])
        r = pd.Series(ov[idx]).groupby(yv[idx]).mean()
        bs.append(100 * (r * w).sum() / w[r.index].sum())
    return round(est, 1), round(np.percentile(bs, 2.5), 1), round(np.percentile(bs, 97.5), 1), len(xs)


rr = []
R19['standards_f'] = R19.standards.astype(float); R19['recognisable_f'] = R19.recognisable.astype(float)
for outcome in ['recognisable_f', 'standards_f']:
    for k in ['doc_readme', 'doc_sdrf', 'doc_design_author', 'doc_params']:
        for val in [True, False]:
            e, lo, hi, n = std_rate(R19, outcome, R19[k] == val)
            rr.append({'outcome': outcome[:-2], 'documentation': k, 'present': val, 'n': n, 'year_std_pct': e, 'lo': lo, 'hi': hi,
                       'note': 'partly definitional: parameter files (mqpar.xml, fragpipe.workflow) are themselves tool output'
                       if k == 'doc_params' and outcome == 'recognisable_f' else ''})
    for qv in ['Q1 shortest', 'Q2', 'Q3', 'Q4 longest']:
        e, lo, hi, n = std_rate(R19, outcome, R19.desc_q == qv)
        rr.append({'outcome': outcome[:-2], 'documentation': 'description length quartile (within year)', 'present': qv, 'n': n,
                   'year_std_pct': e, 'lo': lo, 'hi': hi})
pd.DataFrame(rr).to_csv(OUT / 'H8-results-by-documentation-2019-2025.csv', index=False)
# within-group persistence of README / SDRF / author design table (consecutive datasets of the same group)
G = D[D.pi_group.notna()].sort_values(['pi_group', 'submission_date', 'accession']).copy()
pr = []
for k in ['doc_readme', 'doc_sdrf', 'doc_design_author', 'doc_params']:
    G['prev'] = G.groupby('pi_group')[k].shift()
    x = G[G.prev.notna()]
    base = x.groupby('release_year')[k].mean()
    for pv in [True, False]:
        s = x[x.prev == pv]
        exp = s.release_year.map(base).mean()
        pr.append({'documentation': k, 'previous_dataset_had_it': pv, 'pairs': len(s), 'groups': s.pi_group.nunique(),
                   'pct_now': round(100 * s[k].mean(), 1), 'pct_expected_same_years': round(100 * exp, 1)})
pd.DataFrame(pr).to_csv(OUT / 'H8-documentation-persistence-within-groups.csv', index=False)
# turnover: first-year vs returning groups, by period
D['who'] = np.where(D.pi_group.isna(), 'no group id', np.where(D.group_first_year == D.release_year, 'first-year group', 'returning group'))
tv = D.groupby(['period', 'who']).agg(datasets=('accession', 'size'), desc_chars_median=('desc_chars', 'median'),
                                      sample_protocol_chars_median=('sample_protocol_chars', 'median'),
                                      **{k + '_pct': (k, 'mean') for k in ['doc_readme', 'doc_sdrf', 'doc_design_author']})
for c in tv.columns:
    if c.endswith('_pct'):
        tv[c] = (100 * tv[c]).round(2)
tv.to_csv(OUT / 'H8-documentation-first-year-vs-returning.csv')
# SDRF/README by tool family (2023-2025) to separate tool ecosystems from time
fam = D[D.release_year >= 2023].groupby('family').agg(datasets=('accession', 'size'), sdrf_pct=('doc_sdrf', 'mean'), readme_pct=('doc_readme', 'mean'),
                                                       design_author_pct=('doc_design_author', 'mean'), params_pct=('doc_params', 'mean'),
                                                       desc_chars_median=('desc_chars', 'median'))
for c in fam.columns:
    if c.endswith('_pct'):
        fam[c] = (100 * fam[c]).round(1)
fam.to_csv(OUT / 'H8-documentation-by-software-family-2023-2025.csv')

# period rates with group-bootstrap CIs for the main trends (H5 and H8)
pr_rows = []
TREND = [('isobaric metadata (tmt_expanded or iTRAQ)', 'isobaric_meta', None), ('iTRAQ named in text', 'itraq_text', None),
         ('silac_expanded', 'silac_expanded', None), ('phospho in title/description/keywords', 'design_phosphorylation', None),
         ('Phospho site table | MaxQuant site tables visible', 'site_phosphorylation', 'mq_sites_visible'),
         ('GlyGly site table | MaxQuant site tables visible', 'site_ubiquitin remnant (GlyGly)', 'mq_sites_visible'),
         ('SDRF (top-level)', 'doc_sdrf', None), ('README (top-level or zip)', 'doc_readme', None),
         ('author-written design/sample table', 'doc_design_author', None), ('search parameter/log file', 'doc_params', None),
         ('structured software field filled', 'has_software_annotation', None)]
for lab, col, pop in TREND:
    for per, x in D.groupby('period'):
        x = x[x[pop]] if pop else x
        pr_rows.append(rate_row(None, x.assign(value=x[col].astype(float)), lab, period=per, population=pop or 'all datasets'))
pd.DataFrame(pr_rows).to_csv(OUT / 'Trends-by-period-CIs.csv', index=False)
# composition check for the phospho site-table trend: software family and DIA among phospho-titled datasets
ph = D[D.design_phosphorylation]
comp = pd.crosstab(ph.period, ph.family, normalize='index').mul(100).round(1)
comp['DIA_pct'] = ph.groupby('period').dia_expanded.mean().mul(100).round(1); comp['datasets'] = ph.groupby('period').size()
comp.to_csv(OUT / 'H5-phospho-titled-datasets-composition-by-period.csv')

# =============================================================================================
# 7. denominators and summary
# =============================================================================================
den = {'datasets': len(D), 'with_top_level_files': int(D.has_runs.sum() + 0), 'datasets_with_raw_runs_listed': int(D.has_runs.sum()),
       'datasets_with_listed_results_zip': int(D.zip_listed.sum()), 'datasets_with_a_read_table_header': int(D.header_read.sum()),
       'datasets_with_MaxQuant_site_tables_visible': int(D.mq_sites_visible.sum()), 'datasets_with_MaxQuant_proc_folder': int(D.mq_proc.sum()),
       'zip_members_scanned_after_prefilter': int(con.execute('select count(*) from zm').fetchone()[0])}
den['with_top_level_files'] = int(con.execute('select count(distinct accession) from top').fetchone()[0])
pd.Series(den).to_csv(OUT / 'Denominators.csv', header=['n'])
SUMMARY['denominators'] = den
SUMMARY['H6_overall'] = {'datasets_with_any_evidenced_tool': int(D.ev_any_tool.sum()),
                         'pct_two_or_more_tools_among_evidenced': round(100 * (D[D.ev_any_tool].n_tools_evidenced >= 2).mean(), 1)}
keep = ['accession', 'release_year', 'pi_group', 'dia_expanded', 'tmt_expanded', 'isobaric_meta', 'text_plex', 'h_iso', 'h_iso_channels', 'h_iso_kit',
        'h_lfq', 'h_ibaq', 'h_silac', 'header_read', 'mq_sites_visible', 'mq_proc', 'mq_proc_lfq', 'any_site_table', 'run_isobaric', 'run_phospho', 'zip_listed']
keep += [f'{a}_{k}' for k in PTM_TEXT for a in ['site', 'hdr', 'design', 'struct', 'enrich']]
keep += [f'{a}_{t}' for t in SIG for a in ['named', 'ev']] + ['n_tools_evidenced', 'n_tools_named']
keep += ['desc_chars', 'n_keywords', 'sample_protocol_chars', 'data_protocol_chars', 'desc_duplicated'] + DOCS + ['form']
X = D[keep].copy(); X['site_mods'] = D.site_mods.str.join('; '); X['h_ptm'] = D.h_ptm.str.join('; ')
X.to_csv(OUT / 'Dataset-strategy-software-docs.csv.gz', index=False)
(OUT / 'Summary.json').write_text(json.dumps(SUMMARY, indent=2, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))
print(json.dumps(SUMMARY, indent=1, default=str))
